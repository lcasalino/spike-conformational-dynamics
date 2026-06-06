# Spike conformational dynamics analyses

*Developed and written by Lorenzo Casalino (UC San Diego).*

VMD/Tcl scripts used to quantify large-scale conformational motions of the
SARS-CoV-2 spike ectodomain. The same set of observables is computed for two
datasets:

* **RAV** — the *Respiratory Aerosol Virion* all-atom simulation, which
  contains many spikes; analysed spike-by-spike.
* **Single-spike** — the isolated-spike simulations (open / closed / mutant),
  each with several replicas; analysed replica-by-replica.

The two datasets differ only in how the chains A/B/C are named in the topology;
the geometric definition of every observable is identical and is written once
(in `rav_selections.tcl`) and applied to both.

## Observables

| Observable | Description |
|------------|-------------|
| Ankle / Hip / Knee tilting angle | Tilting angle of three points along the spike stem (angle at an origin COM between two arm COMs). |
| NTD triangle area | Area of the triangle formed by the COMs of the three N-terminal domains (residues 13–291). |
| RBD–central-helix (CH) distance | Distance from each receptor-binding domain to the central helix (RBD opening). |

## Files

Shared library and configuration:

| File | Role |
|------|------|
| `rav_common.tcl`     | Geometry primitives, output formatters, per-trajectory helpers, and the two dataset drivers (`rav_run_spikes`, `rav_run_replicas`). |
| `rav_selections.tcl` | The atom selections and per-frame measurement for each observable (the only analysis-specific code). |
| `rav_paths.tcl`      | **Edit this** — input/output paths and dataset parameters, in one place. |

Runner scripts (one per observable per dataset):

| RAV | Single-spike |
|-----|--------------|
| `ankle_tilting_rav.tcl`        | `ankle_tilting_singlespike.tcl` |
| `hip_tilting_rav.tcl`          | `hip_tilting_singlespike.tcl` |
| `knee_tilting_rav.tcl`         | `knee_tilting_singlespike.tcl` |
| `ntd_area_rav.tcl`             | `ntd_area_singlespike.tcl` |
| `rbd_ch_distance_rav.tcl`      | `rbd_ch_distance_singlespike.tcl` |

Utilities (not required to reproduce the analysis):

| File | Role |
|------|------|
| `plotting_for_RAV_paper.ipynb` | Notebook that reads the per-trajectory `.txt` files and produces the paper figures (KDE distributions; RBD distance in 4 time blocks). Reads only the value columns, recomputing the time axis from the frame index. |
| `verify_refactor.tcl` | Self-test (plain `tclsh`, no VMD): checks the selection strings and the driver logic. Run `tclsh verify_refactor.tcl`. |

## Requirements

* **VMD 1.9.x** (developed and used with VMD 1.9.4, embedded Tcl 8.5).
* No external Tcl packages. The only non-core commands used,
  `tcl::mathfunc::min` / `tcl::mathfunc::max`, are Tcl 8.5+ built-ins; all
  other commands are provided by VMD.

## Input data

**RAV** — one pair of files per spike (`NN` = `00`..`29`):

```
<RAV_INPPATH>/<RAV_NAME><NN>.alignedBackbone.psf
<RAV_INPPATH>/<RAV_NAME><NN>.alignedBackbone.dcd
```

Each spike has six protein segments `A1<NN> A2<NN> B1<NN> B2<NN> C1<NN> C2<NN>`
(two per protomer). 2530 saved frames at **0.2 ns/frame** = **506 ns/spike**.
Spike index **6 does not exist in the RAV system**, so it is skipped (`RAV_SKIP`
in `rav_paths.tcl`); the analyses therefore cover 29 spikes (`00`–`29` minus
`06`).

**Single-spike** — one topology per system plus one dcd per replica:

```
<SS_INPPATH>/<system>.psf
<SS_INPPATH>/<system>_<rep>.dcd
```

The single spike's six segments are `AS1 AS2 BS1 BS2 CS1 CS2`. Systems and
replica counts: open (6), closed (3), mutant (6). Frame spacing **0.1 ns/frame**.

The trajectories were pre-aligned on the spike backbone, so the scripts perform
no alignment.

## Configuration

Edit `rav_paths.tcl` only. It defines the input/output directories, the
trajectory basenames, the spike/replica counts, and the frame spacing (`dt`)
for each dataset.

## Running

```bash
vmd -dispdev text -e ankle_tilting_rav.tcl
vmd -dispdev text -e ankle_tilting_singlespike.tcl
# ... and likewise for hip / knee / ntd_area / rbd_ch_distance
```

Keep `rav_common.tcl`, `rav_selections.tcl` and `rav_paths.tcl` in the same
directory as the runner scripts (each runner sources them by path).

## Regenerating all data

By default the analyses write to a self-contained `results/` folder next to the
scripts (set in `rav_paths.tcl`), so re-running never overwrites the original
output folders.

To regenerate everything with one command:

```bash
cd source_code
VMD=/home/lcasalino/Software/vmd-1.9.4a57/bin/vmd ./run_all.sh        # all (RAV + single-spike)
./run_all.sh rav           # only the RAV analyses
./run_all.sh singlespike   # only the single-spike analyses
```

`run_all.sh` runs each of the ten analyses with `vmd -dispdev text -e`, writing
data to `results/<ANALYSIS>/` and per-script logs to `results/logs/`. Set the
`VMD` variable to your VMD binary (the default points at VMD 1.9.4a57). **The
RAV runs read ~2 GB per spike across 29 spikes — expect several hours and heavy
I/O.** Any single analysis can also be run on its own, e.g.
`vmd -dispdev text -e ankle_tilting_rav.tcl`.

To then make the figures from the regenerated data, the notebook's `path_out`
points at the `results` folder; run it in Jupyter (or
`jupyter nbconvert --to notebook --execute plotting_for_RAV_paper.ipynb` if
`nbconvert` is installed).

The per-trajectory time-series outputs are provided alongside the code as
`results_timeseries.zip` (the `.txt` files only — not the large PDB snapshots).
Unzip it in this folder to recreate `results/`, then run the notebook to
reproduce the figures **without** re-running VMD:

```bash
unzip results_timeseries.zip       # creates results/<ANALYSIS>/*.txt
```

### Memory and parallelism

VMD loads the whole trajectory into RAM, so peak memory per job is roughly the
trajectory size plus ~0.5 GB of VMD overhead. The single-spike jobs need much
more memory than the RAV jobs because their replicas have up to ~10 000 frames
(vs 2530 for RAV):

| Job type | Loaded trajectory | Peak RAM per job |
|----------|-------------------|------------------|
| RAV (each of 5)          | 2530 frames, ~2.0 GB                                   | **~2.5 GB** |
| single-spike (each of 5) | one replica at a time; largest ~8.1 GB (open/mutant), ~4–5 GB (closed) | **~8.5 GB** |

Total peak RAM is `(concurrent jobs) × (per-job peak)`. Worst-case totals:

| Command | Concurrent | Peak RAM |
|---------|-----------|----------|
| `./run_all.sh rav -j 5`         | 5 × 2.5 GB | ~13 GB |
| `./run_all.sh singlespike -j 2` | 2 × 8.5 GB | ~17 GB |
| `./run_all.sh singlespike -j 3` | 3 × 8.5 GB | ~26 GB |
| `./run_all.sh singlespike -j 5` | 5 × 8.5 GB | ~43 GB |
| `./run_all.sh all -j 10`        | 5 RAV + 5 SS | ~55 GB |

Notes:
* The five single-spike scripts iterate the systems in the same order, so if
  launched together they tend to load the same large replica simultaneously —
  plan for the worst-case total above, and expect heavy `results/` filesystem
  I/O.
* Per-job RAM is fixed by the trajectory size; the only knob is `-j`. Check the
  node first with `free -g` and leave ~20% headroom.
* Rough rule of thumb: RAV set is comfortable on a 32 GB node with `-j 5`; the
  single-spike set wants `-j 2–3` on a typical node and `-j 5` only on 64 GB+;
  `all -j 10` needs a 64 GB+ node.
* On a large-memory node, memory stops being the limit (the whole set peaks at
  ~55 GB): set `-j` by the available CPU core count and filesystem bandwidth
  instead, and run everything at once with `./run_all.sh all -j 10`. Each job is
  effectively single-threaded, so there is no benefit to `-j` beyond the number
  of analyses (10).
* Splitting across machines (e.g. the RAV set on one node, the single-spike set
  on another, all writing to the shared `results/`) keeps each node's peak
  modest and halves the wall-clock.

## Output

Results are written under `<outpath>/<ANALYSIS>/`:

| File | Contents |
|------|----------|
| `spike<i>.<ANALYSIS>.txt` / `rep<r>.<ANALYSIS>.<system>.txt` | Per-trajectory time series. Tilting: `time  angle  (angle − angle@frame0)`. NTD: `time \t area`. RBD: `time \t dRBD_A dRBD_B dRBD_C`. |
| `spike<i>.max.frame<f>.pdb` / `spike<i>.min...` (RAV); `rep<r>.max.frame<f>.<system>.pdb` / `rep<r>.min...` (single-spike) | Conformation at each trajectory's own extremum (one per spike/replica). |
| `absolute.max.pdb` / `absolute.min.pdb` (RAV); `absolute.max.<system>.pdb` / `absolute.min.<system>.pdb` (single-spike) | The global extremum conformation — a single file, overwritten as new running records are found. The spike/replica and frame are recorded in the `MAX/MIN OVERALL` line of the min/max log. |
| `maxmin.spikes.txt` / `maxmin.reps.txt` | Per-trajectory and overall min/max log. |
| `concatenated_spikes_*.txt` / `concatenated_reps_*.txt` | All per-trajectory series concatenated. |

Exported PDBs contain the spike's six segments plus any glycan residue with at
least one atom within 2 Å.

## Plotting (notebook)

`plotting_for_RAV_paper.ipynb` builds the paper figures from the per-trajectory
`.txt` files. It reads only the value columns and recomputes the time axis from
the frame index, so it is unaffected by the time column in the input files.

**Python dependencies:** `numpy`, `matplotlib`, `seaborn`, `scipy`
(e.g. `pip install numpy matplotlib seaborn scipy`). Developed with Python 3 /
Jupyter.

**Fonts:** the figures use the *Libre Franklin* sans-serif font. If it is not
installed, matplotlib falls back to a default font silently (the plots are
still correct, only the typography differs). Install Libre Franklin and clear
the matplotlib cache (`rm -rf ~/.cache/matplotlib`) to reproduce the exact
appearance.

**Time window for the RAV vs single-spike comparison:** for the tilting-angle
and NTD-area distributions, only the **last 1760 frames (≈352 ns)** of each RAV
spike are used. This approximately matches the timescale sampled by the
individual single-spike simulations (which range from ~400 to ~1000 ns per
replica), so the two distributions are compared over comparable amounts of
simulation time. The RBD–CH distance figures instead split each full RAV spike
trajectory into 4 equal time blocks.

## Figure-to-script mapping

Each paper figure (saved by the notebook into `path_out`) and the scripts whose
output it is built from:

| Figure file | Notebook section | RAV data (script → output folder) | Single-spike data (script → output folder) |
|-------------|------------------|-----------------------------------|--------------------------------------------|
| `tilting_HIP_ss_vs_RAVlast352_wide_v3.png`   | Tilting / NTD | `hip_tilting_rav.tcl`   → `HIP_TILTING_ANGLE`   | `hip_tilting_singlespike.tcl`   → `HIP_TILTING_ANGLE_SINGLESPIKE`   |
| `tilting_KNEE_ss_vs_RAVlast352_wide_v3.png`  | Tilting / NTD | `knee_tilting_rav.tcl`  → `KNEE_TILTING_ANGLE`  | `knee_tilting_singlespike.tcl`  → `KNEE_TILTING_ANGLE_SINGLESPIKE`  |
| `tilting_ANKLE_ss_vs_RAVlast352_wide_v3.png` | Tilting / NTD | `ankle_tilting_rav.tcl` → `ANKLE_TILTING_ANGLE` | `ankle_tilting_singlespike.tcl` → `ANKLE_TILTING_ANGLE_SINGLESPIKE` |
| `NTD_area_ss_vs_RAVlast352_wide_v3.png`      | Tilting / NTD | `ntd_area_rav.tcl`      → `NTD_AREA_RAV`         | `ntd_area_singlespike.tcl`      → `NTD_AREA_SINGLESPIKE`            |
| `RBD_CH_RAV_4timeblocks_chainA_viridis_lines.png` | RBD–CH distance | `rbd_ch_distance_rav.tcl` → `RBD_CH_DIST_RAV` (col 1, chain A) | not used |
| `RBD_CH_RAV_4timeblocks_chainB_plasmaR_lines.png` | RBD–CH distance | `rbd_ch_distance_rav.tcl` → `RBD_CH_DIST_RAV` (col 2, chain B) | not used |
| `RBD_CH_RAV_4timeblocks_chainC_plasmaR_lines.png` | RBD–CH distance | `rbd_ch_distance_rav.tcl` → `RBD_CH_DIST_RAV` (col 3, chain C) | not used |

The four tilting/NTD figures compare the single-spike distribution (full
trajectories) against the RAV distribution (last ≈352 ns of each spike). The
three RBD figures use only the RAV data, split into 4 time blocks, one figure
per protomer (chains A/B/C = output columns 1/2/3).

## Note on the angle convention

`rav_angle` returns the angle at a vertex using the **absolute value** of the
dot product, i.e. an unsigned acute angle in [0°, 90°]. This is the convention
used to generate the published data and is preserved here unchanged (including
the original 57.2958 ≈ 180/π factor), so results are reproduced exactly.
