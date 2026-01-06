## LICH-TEST-GO

**LICH-TEST-GO** is a functionally extended, “out-of-the-box” workflow built on top of the official **LICH-TEST** project for classifying ice-related local structures in atomistic simulations. The original LICH-TEST (Liquid water, Interfacial, Cubic and Hexagonal Ice Classification through Eclipsed and Staggered Conformation Template Matching) is designed to classify different ice types in simulations via template matching. ([GitHub][1])

This GO version focuses on practical usability for trajectory processing:

* Command-line, batch-friendly execution for **topology + trajectory** inputs.
* Frame subsampling via **frame stride** (e.g., `--step 50`) to accelerate long-trajectory analysis.
* A standardized **per-frame count table** output (`.csv`) for downstream analysis.
* A simple plotting step to turn the `.csv` into a `.png` time series.

---
<p align="center">
  <img src="/lich_test_GO/frameCount.png?raw=true" alt="LICHTEST" width="600">
</p>

## 🚀 How to Use 🚀

### 1) Requirements

You need a working Python environment with:

```bash
mamba install mdanalysis -c conda-forge
```


MDAnalysis supports DCD trajectories via dedicated readers (commonly produced by CHARMM/NAMD/LAMMPS). ([docs.mdanalysis.org][3])

### 2) Run classification (produces `.csv`)

Typical usage (example: Amber `prmtop` + `dcd` trajectory; **frame stride = 5**):

```bash
python lich_test_GO.py \
  --top YOUR_TOPFILE \
  --traj YOUR_TRAJFILE \
  --step 50 \
  --counts-csv ./frameCount.csv
```

Conceptually, this follows the MDAnalysis model `Universe(topology, trajectory)` where the topology provides atom identities/connectivity and the trajectory provides coordinates (and, when present, box dimensions). ([userguide.mdanalysis.org][2])

### 3) Plot your results (produces `.png`)

If you use the companion plotting script:

```bash
# Run in a working directory that contains frame_counts.csv
python lich_test_GO_draw.py
```

This step reads the per-frame counts table and creates a PNG plot (useful for quickly visualizing ice fraction evolution).

---

## Outputs

### 1) CSV output (`.csv`)

The main artifact is a per-frame count table with the following header:

```text
frame,ncubic,nhex,nmix,nci,nhi,nch,nint,nliq
```

Example rows:

```csv
frame,ncubic,nhex,nmix,nci,nhi,nch,nint,nliq
0,0,0,0,0,0,0,0,256
50,0,0,0,0,0,0,0,256
100,0,0,0,0,0,1,3,252
150,0,0,0,0,0,0,0,256
200,0,0,0,0,0,0,0,256
250,0,0,0,0,0,2,6,248
300,0,0,0,0,0,2,6,248
```

**Column meanings (counts per analyzed frame):**

* `frame`: the trajectory frame index that was analyzed (note: with `--step`, frames are subsampled)
* `ncubic`: cubic ice–like (Ice Ic–like) molecules
* `nhex`: hexagonal ice–like (Ice Ih–like) molecules
* `nmix`: mixed / stacking-disordered–like ice
* `nci`: cubic-interfacial–like
* `nhi`: hex-interfacial–like
* `nch`: clathrate-hydrate–like
* `nint`: other interfacial/intermediate structures
* `nliq`: liquid-like (unclassified into the above solid/interfacial categories)

### 2) Console output (runtime logs)

During execution you will typically see:

**(A) Per-frame summary lines**, e.g.

```text
Frame 0: cubic=0 hex=0 mixed=0 ci=0 hi=0 ch=0 int=0 liq=256
```

These lines provide an immediate sanity check that the classification is evolving plausibly across frames.

**(B) Neighbour-list construction / progress logs**, e.g.

```text
Neighbourlist generation:
Generating a cell list with 64 cells, x:4, y:4, z:4
...
Progress: 0 %
...
Progress: 100 %
Avg. # of neighbours: 4.0
Type and size of neigh_list: <class 'numpy.ndarray'> (256, 4)
Calculating neighbour directions
Time spent on UtV_func, score_func: 0.0166 0.0164
```

These logs are intended for transparency and debugging (e.g., verifying neighbor counts, monitoring progress on large systems, and profiling time spent in key kernels).

---

## Acknowledgements

* Official LICH-TEST repository (original method and implementation). ([GitHub][1])
* LICH-TEST method description and scope (template matching; liquid/Ic/Ih/clathrate/interfacial recognition; DOI: 10.1021/acs.jpcb.1c01926). ([ACS JPCB][3])
* MDAnalysis documentation for loading `Universe(topology, trajectory)` and DCD trajectory support. ([userguide.mdanalysis.org][2])

[1]: https://github.com/opakarin/lich-test "GitHub - opakarin/lich-test: lich-test is an analysis tool for ice ..."
[2]: https://userguide.mdanalysis.org/stable/universe.html "Universe — MDAnalysis User Guide documentation"
[3]: https://pubs.acs.org/doi/10.1021/acs.jpcb.1c01926 "Liquid Water and Interfacial, Cubic, and Hexagonal Ice ... - Aalto"
