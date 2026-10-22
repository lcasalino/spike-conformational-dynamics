# Sampling and temporal robustness

*Developed and written by Lorenzo Casalino (UC San Diego).*

Python analyses quantifying how well the spike conformational observables are
sampled within the RAV trajectory, and how the 29 spikes of the RAV compare
with the single-spike simulations. These produce the two supplementary figures
and the numbers quoted in the "Sampling and temporal robustness" section of the
Supporting Information.

The scripts read only the per-trajectory `.txt` time series written by the
VMD/Tcl analyses in the parent folder (`../results/<ANALYSIS>/`); no
trajectories are re-read.

## What is computed

| Quantity | Description |
|----------|-------------|
| Integrated autocorrelation time | Per trajectory, Sokal automatic windowing. RAV values are window-limited and reported as lower bounds. |
| Block-averaged uncertainty | Flyvbjerg–Petersen block averaging, compared against the SEM across spikes and the naive per-frame SEM. |
| Drift | Least-squares slope fitted per spike over the analysis window, tested against zero by one-sample *t*-test. |
| Block-wise stationarity | Distributions over four sequential 88-ns blocks; Jensen–Shannon distance between the first and last block against a spike-resampling noise floor. |
| Variance decomposition | Between-spike vs within-spike variance and the intraclass correlation, corrected for the sampling noise in each per-spike mean. |
| Fluctuation amplitude | Per-trajectory standard deviation over the matched 154–506 ns interval, RAV vs single-spike (Mann–Whitney *U*). |
| Mean differences | RAV vs single-spike, with 95% CIs from a hierarchical bootstrap over spikes and replicas. |
| Departure from the initial structure | First-frame value vs window mean, per spike, and their correlation across spikes. |

## Analysis window

The first 154 ns of each RAV spike (770 frames) contain the initial
conformational relaxation and are excluded; all quantities use the final 1760
frames (352 ns), the same window as Fig. 5 of the paper. This gives
29 × 352 ns = 10.21 µs, matched to the 10.24 µs of the single-spike set.

## Branched single-spike replicas

Replicas 4–6 of the open and mutant systems were started from frames of earlier
replicas rather than from the common starting structure. Their branch times
(taken from `extractframes.tcl` in the respective simulation directories) are
hard-coded in `SS_BRANCH` at the top of `sampling_robustness.py` and are used
whenever a common time axis is required.

## Files

| File | Role |
|------|------|
| `sampling_robustness.py` | Main analysis; writes the summary, the two tables and both figures. |
| `per_molecule_table.py`  | Per-trajectory table of first-frame value, last-frame value, mean and SD. |
| `si_table_per_spike.py`  | Supplementary table of per-spike summary statistics (CSV, Markdown, LaTeX, HTML and .docx). |
| `si_fig_per_spike_hinges.py` | Supplementary figure: per-spike hinge-angle distributions, one panel per spike. |

## Requirements

Python 3 with NumPy, SciPy and Matplotlib. Developed with Python 3.6,
NumPy 1.16, SciPy 1.2 and Matplotlib 2.2; no features newer than these are
used. The figures use the *Libre Franklin* font if installed and fall back to
the default sans-serif otherwise.

## Running

```bash
python3 sampling_robustness.py            # reads ../results
python3 sampling_robustness.py <results>  # or an explicit results directory
python3 per_molecule_table.py
python3 si_table_per_spike.py
python3 si_fig_per_spike_hinges.py
```

`sampling_robustness.py` takes a few minutes; the hierarchical bootstrap is the
slow part. `per_molecule_table.py` is immediate.

## Output

| File | Contents |
|------|----------|
| `sampling_robustness_summary.txt` | All quantities, one block per observable, human-readable. |
| `summary_table.csv`      | The same numbers, one row per observable. |
| `per_spike_table.csv`    | Per-spike mean, SD and correlation time for every observable. |
| `per_molecule_values.csv`| Per-trajectory first-frame value, last-frame value, mean and SD, for the 29 RAV spikes and the 15 single-spike replicas. |
| `SI_table_per_spike.{csv,md,tex,html,docx}` | Supplementary table: per-spike value at t = 0 and mean ± SD over the analysis window for each observable, plus the chain-A RBD state. |
| `SI_per_spike_hinges.{png,pdf}` | Supplementary figure: hip, knee and ankle distributions for each of the 29 spikes. |
| `SI_spike_sampling.{png,pdf}` | Supplementary figure: sampling and temporal stability within the RAV trajectory. |
| `SI_spike_rav_vs_single.{png,pdf}` | Supplementary figure: comparison of the RAV spikes with the single-spike simulations. |

Random seeds are fixed, so repeated runs reproduce the published numbers
exactly.
