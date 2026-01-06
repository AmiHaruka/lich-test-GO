#!/usr/bin/env python3
"""
Cluster growth analysis for ice-like molecules using LICH-TEST labels.

Per frame it:
- runs LICH-TEST labeling
- selects ice-like labels (configurable, default: cubic=1, hex=2, HI=5; optional INT=7)
- builds an O–O neighbor graph with cutoff (PBC-aware via FastNS)
- finds connected components (clusters) and reports:
  * total ice count
  * number of clusters
  * largest (Smax) and second-largest (S2) cluster sizes
  * composition of the largest cluster (fractions of hex/cubic/hi/int)

Outputs a CSV summarizing these metrics per frame.
"""

import argparse
import sys
from pathlib import Path
from typing import Iterable, List

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.nsgrid import FastNS
except ImportError as exc:  # pragma: no cover - dependency check
    sys.exit(
        "MDAnalysis is required (pip install MDAnalysis).\n"
        f"Import error: {exc}"
    )

from find_neigh_dir import find_neigh_dir
from ice_labelling_func import ice_labelling_func
from neighlistcell import neighbours
from temp_matching_func import temp_matching_func


class DSU:
    """Disjoint set union for connected components."""

    def __init__(self, n: int):
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, x: int) -> int:
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]


def wrap_positions(coords: np.ndarray, box_lengths: np.ndarray) -> np.ndarray:
    return coords - box_lengths * np.round(coords / box_lengths)


def build_selection(oxygen_names: Iterable[str]) -> str:
    return " or ".join([f"name {nm}" for nm in oxygen_names])


def classify_frame(coord: np.ndarray, box_lengths: np.ndarray, min_score: float) -> np.ndarray:
    xdim, ydim, zdim = box_lengths
    neigh_list, N, neigh_num = neighbours(coord, xdim, ydim, zdim)
    neigh_dir = find_neigh_dir(coord, neigh_list, xdim, ydim, zdim, N)
    stg_score, ecl_score = temp_matching_func(neigh_dir, neigh_list, neigh_num, N, min_score)
    labels = ice_labelling_func(stg_score, ecl_score, neigh_list, neigh_num, N)
    return labels


def write_csv(path: Path, rows: List[dict], headers: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fw:
        fw.write(",".join(headers) + "\n")
        for row in rows:
            fw.write(",".join(str(row[h]) for h in headers) + "\n")


def cluster_stats(ice_coords: np.ndarray, ice_labels: np.ndarray, cutoff: float, box_lengths: np.ndarray) -> dict:
    """Return cluster metrics for a set of ice-like atoms."""
    n = ice_coords.shape[0]
    if n == 0:
        return {
            "n_ice": 0,
            "n_clusters": 0,
            "Smax": 0,
            "S2": 0,
            "n_hex_max": 0,
            "n_cubic_max": 0,
            "n_hi_max": 0,
            "n_int_max": 0,
        }

    # neighbor search with cutoff (PBC-aware)
    ns = FastNS(cutoff, ice_coords, box_lengths)
    pairs = ns.self_search().get_pairs()

    dsu = DSU(n)
    for i, j in pairs:
        dsu.union(int(i), int(j))

    # component sizes
    roots = np.array([dsu.find(i) for i in range(n)], dtype=int)
    sizes = np.bincount(roots)
    if sizes.size == 0:
        sizes = np.array([n], dtype=int)
    order = np.argsort(sizes)[::-1]
    smax = int(sizes[order[0]]) if sizes.size > 0 else 0
    s2 = int(sizes[order[1]]) if sizes.size > 1 else 0
    n_clusters = int(np.count_nonzero(sizes))

    # composition of largest cluster
    if sizes.size == 0:
        n_hex_max = n_cubic_max = n_hi_max = n_int_max = 0
    else:
        largest_root = order[0]
        mask_max = roots == largest_root
        n_hex_max = int(np.sum(ice_labels[mask_max] == 2))
        n_cubic_max = int(np.sum(ice_labels[mask_max] == 1))
        n_hi_max = int(np.sum(ice_labels[mask_max] == 5))
        n_int_max = int(np.sum(ice_labels[mask_max] == 7))

    return {
        "n_ice": n,
        "n_clusters": n_clusters,
        "Smax": smax,
        "S2": s2,
        "n_hex_max": n_hex_max,
        "n_cubic_max": n_cubic_max,
        "n_hi_max": n_hi_max,
        "n_int_max": n_int_max,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Cluster growth analysis from LICH-TEST labels")
    parser.add_argument("--top", required=True, help="Topology file (prmtop)")
    parser.add_argument("--traj", required=True, help="Trajectory file (e.g., md.nc)")
    parser.add_argument("--start", type=int, default=None, help="Start frame (0-based, inclusive)")
    parser.add_argument("--stop", type=int, default=None, help="Stop frame (exclusive)")
    parser.add_argument("--step", type=int, default=1, help="Frame stride")
    parser.add_argument("--oxygen-names", nargs="+", default=["OW", "mW"], help="Atom names for water oxygens")
    parser.add_argument("--min-score", type=float, default=0.5, help="LICH-TEST similarity cutoff")
    parser.add_argument("--cutoff", type=float, default=3.5, help="O–O cutoff for cluster connectivity (Å)")
    parser.add_argument("--ice-labels", nargs="+", type=int, default=[1, 2, 5], help="Labels to treat as ice-like (e.g., 1 cubic, 2 hex, 5 HI)")
    parser.add_argument("--include-int", action="store_true", help="Also include label 7 (interfacial) in ice-like set")
    parser.add_argument("--dt-ps", type=float, default=None, help="Time per frame (ps). If omitted, uses ts.dt if available else 1.0")
    parser.add_argument("--out-csv", default="cluster_growth.csv", help="Output CSV path")
    args = parser.parse_args()

    ice_labels_set = set(args.ice_labels)
    if args.include_int:
        ice_labels_set.add(7)

    sel_str = build_selection(args.oxygen_names)
    u = mda.Universe(args.top, args.traj)
    sel = u.select_atoms(sel_str)
    if sel.n_atoms == 0:
        sys.exit(f"No atoms found for selection: {sel_str}")

    rows = []
    headers = [
        "frame",
        "time_ps",
        "n_ice",
        "n_clusters",
        "Smax",
        "S2",
        "n_hex_total",
        "n_cubic_total",
        "n_hi_total",
        "n_int_total",
        "frac_hex_max",
        "frac_cubic_max",
        "frac_hi_max",
        "frac_int_max",
    ]

    times_ps = []
    print(f"Analysing frames {args.start}:{args.stop}:{args.step}, selection '{sel_str}'")
    for ts in u.trajectory[args.start:args.stop:args.step]:
        box_lengths = np.array(ts.dimensions[:3], dtype=float)
        coord = wrap_positions(sel.positions.copy().astype(float), box_lengths)
        labels = classify_frame(coord, box_lengths, args.min_score)
        times_ps.append(ts.time if ts.dt is not None else np.nan)

        ice_mask = np.isin(labels, list(ice_labels_set))
        ice_coords = coord[ice_mask]
        ice_labels = labels[ice_mask]

        stats = cluster_stats(ice_coords, ice_labels, args.cutoff, box_lengths)
        n_hex_total = int(np.sum(labels == 2))
        n_cubic_total = int(np.sum(labels == 1))
        n_hi_total = int(np.sum(labels == 5))
        n_int_total = int(np.sum(labels == 7))

        denom = stats["n_hex_max"] + stats["n_cubic_max"] + stats["n_hi_max"] + stats["n_int_max"]
        if denom > 0:
            frac_hex = stats["n_hex_max"] / denom
            frac_cubic = stats["n_cubic_max"] / denom
            frac_hi = stats["n_hi_max"] / denom
            frac_int = stats["n_int_max"] / denom
        else:
            frac_hex = frac_cubic = frac_hi = frac_int = 0.0

        rows.append({
            "frame": ts.frame,
            "time_ps": times_ps[-1],
            "n_ice": stats["n_ice"],
            "n_clusters": stats["n_clusters"],
            "Smax": stats["Smax"],
            "S2": stats["S2"],
            "n_hex_total": n_hex_total,
            "n_cubic_total": n_cubic_total,
            "n_hi_total": n_hi_total,
            "n_int_total": n_int_total,
            "frac_hex_max": frac_hex,
            "frac_cubic_max": frac_cubic,
            "frac_hi_max": frac_hi,
            "frac_int_max": frac_int,
        })

    # time handling
    times_ps = np.array(times_ps)
    if args.dt_ps is not None:
        dt_ps = args.dt_ps
    else:
        if np.all(np.isfinite(times_ps)) and len(times_ps) > 1:
            dt_ps = float(np.nanmean(np.diff(times_ps)))
        else:
            dt_ps = 1.0
            # fill times by frame index if missing
            for i, row in enumerate(rows):
                row["time_ps"] = i * dt_ps

    write_csv(Path(args.out_csv), rows, headers)
    print(f"Saved cluster metrics to {args.out_csv}")


if __name__ == "__main__":
    main()
