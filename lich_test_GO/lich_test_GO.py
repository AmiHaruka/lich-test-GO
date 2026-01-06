'''
Author:  Kanbe 
Date: 2025-12-02 15:15:07
LastEditors:  Kanbe 
LastEditTime: 2025-12-02 16:24:09
FilePath: /Water/lich_test/multi_frame_lich.py
Description: 
'''
#!/usr/bin/env python3
"""
Multi-frame driver for the LICH-TEST ice classifier.

Features
--------
- Reads Top/Trajectory (prmtop + md.nc) via MDAnalysis.
- Lets you choose which frames to analyse (start/stop/step).
- Lets you define which atom names count as water oxygens (default: OW/mW).
- Reuses the original LICH-TEST pipeline: coord, coordlayer, box (xdim, ydim, zdim),
  watoms, latoms, line2in are prepared per frame to keep parity with the single-frame script.
- Writes optional per-frame labelled xyz and a counts CSV.
"""

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np
import MDAnalysis as mda
from find_neigh_dir import find_neigh_dir
from ice_labelling_func import ice_labelling_func
from neighlistcell import neighbours
from temp_matching_func import temp_matching_func


LABEL_ORDER: List[Tuple[str, int, str]] = [
    ("ncubic", 1, "C"),
    ("nhex", 2, "H"),
    ("nmix", 3, "IM"),
    ("nci", 4, "IC"),
    ("nhi", 5, "IH"),
    ("nch", 6, "CH"),
    ("nint", 7, "In"),
]


def wrap_positions(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Apply minimum image convention to wrap coords into the primary box."""
    return coords - box * np.round(coords / box)


def build_selection(oxygen_names: Iterable[str]) -> str:
    return " or ".join([f"name {nm}" for nm in oxygen_names])


def count_labels(labels: np.ndarray) -> dict:
    counts = {name: int(np.count_nonzero(labels == code)) for name, code, _ in LABEL_ORDER}
    counts["nliq"] = int(labels.size - sum(counts.values()))
    return counts


def classify_frame(coord: np.ndarray, box_lengths: np.ndarray, min_score: float) -> Tuple[np.ndarray, dict]:
    xdim, ydim, zdim = box_lengths
    neigh_list, N, neigh_num = neighbours(coord, xdim, ydim, zdim)
    neigh_dir = find_neigh_dir(coord, neigh_list, xdim, ydim, zdim, N)
    stg_score, ecl_score = temp_matching_func(neigh_dir, neigh_list, neigh_num, N, min_score)
    labels = ice_labelling_func(stg_score, ecl_score, neigh_list, neigh_num, N)
    counts = count_labels(labels)
    return labels, counts


def write_label_xyz(path: Path, coord: np.ndarray, labels: np.ndarray, box_lengths: np.ndarray, line2in: str) -> None:
    """Write per-atom labels in the same style as script.py's `_ice` output."""
    xdim, ydim, zdim = box_lengths
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fw:
        fw.write(f"{coord.shape[0]} \n")
        fw.write(line2in)
        for i in range(coord.shape[0]):
            code = "L"
            for _, lbl_code, lbl_char in LABEL_ORDER:
                if labels[i] == lbl_code:
                    code = lbl_char
                    break
            fw.write(f"{code} {coord[i,0]} {coord[i,1]} {coord[i,2]} \n")


def write_counts_csv(path: Path, rows: List[dict]) -> None:
    headers = ["frame"] + [name for name, _, _ in LABEL_ORDER] + ["nliq"]
    with path.open("w") as fw:
        fw.write(",".join(headers) + "\n")
        for row in rows:
            fw.write(",".join(str(row[h]) for h in headers) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Multi-frame LICH-TEST-GO runner")
    parser.add_argument("--top", required=True, help="Topology file (e.g., system.prmtop)")
    parser.add_argument("--traj", required=True, help="Trajectory file (e.g., md.nc; md.dcd)")
    parser.add_argument("--start", type=int, default=None, help="Start frame")
    parser.add_argument("--stop", type=int, default=None, help="Stop frame")
    parser.add_argument("--step", type=int, default=50, help="Frame stride")
    parser.add_argument(
        "--oxygen-names",
        nargs="+",
        default=["OW", "mW","O"],
        help="Atom names to treat as water oxygens",
    )
    parser.add_argument("--min-score", type=float, default=0.5, help="Minimum similarity score")
    parser.add_argument(
        "--write-label-prefix",
        default=None,
        help="If set, write labelled xyz per frame using this prefix (prefix_frameXXXX_ice.xyz)",
    )
    parser.add_argument("--counts-csv", default="./frameCount.csv", help="Path to write counts CSV")
    args = parser.parse_args()

    sel_str = build_selection(args.oxygen_names)
    u = mda.Universe(args.top, args.traj)
    sel = u.select_atoms(sel_str)
    if sel.n_atoms == 0:
        sys.exit(f"No atoms found for selection: {sel_str}")

    results = []
    label_prefix = Path(args.write_label_prefix) if args.write_label_prefix else None

    print(f"Analysing frames with selection '{sel_str}', total oxygens: {sel.n_atoms}")
    for ts in u.trajectory[args.start:args.stop:args.step]:
        frame_idx = ts.frame
        box_lengths = np.array(ts.dimensions[:3], dtype=float)
        coord = sel.positions.copy().astype(float)
        coord = wrap_positions(coord, box_lengths)

        line2in = f"frame {frame_idx} cell {box_lengths[0]} {box_lengths[1]} {box_lengths[2]}\n"
        watoms = coord.shape[0]
        latoms = u.atoms.n_atoms - watoms
        coordlayer = np.empty((latoms, 3))  # placeholder to keep parity with original signature
        _ = (coordlayer, watoms, latoms)  # quiet linters; available if needed

        labels, counts = classify_frame(coord, box_lengths, args.min_score)
        counts_row = {"frame": frame_idx, **counts}
        results.append(counts_row)

        print(
            f"Frame {frame_idx}: "
            f"cubic={counts['ncubic']} hex={counts['nhex']} mixed={counts['nmix']} "
            f"ci={counts['nci']} hi={counts['nhi']} ch={counts['nch']} int={counts['nint']} "
            f"liq={counts['nliq']}"
        )

        if label_prefix:
            out_path = label_prefix.parent / f"{label_prefix.name}_frame{frame_idx:06d}_ice.xyz"
            write_label_xyz(out_path, coord, labels, box_lengths, line2in)

    if args.counts_csv:
        counts_path = Path(args.counts_csv)
        counts_path.parent.mkdir(parents=True, exist_ok=True)
        write_counts_csv(counts_path, results)
        print(f"Counts written to {args.counts_csv}")


if __name__ == "__main__":
    main()
