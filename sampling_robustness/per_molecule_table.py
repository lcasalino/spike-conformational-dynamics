#!/usr/bin/env python3
# ---------------------------------------------------------------------------
#  Per-molecule summary table: for each RAV spike and each single-spike replica,
#  value at the first frame, value at the last frame, mean and SD of the hip,
#  knee and ankle tilt angles and the NTD triangle area.
#
#  RAV spikes: mean and SD over the analysis window (last 352 ns, as in Fig. 5);
#  first/last frame of the full 506-ns trajectory.
#  Single-spike replicas: mean and SD over the whole replica (as pooled in Fig. 5);
#  first/last frame of the replica. Branched replicas start from a frame of their
#  parent (branch time on the parent's clock in start_ns).
#
#  Output: per_molecule_values.csv
#  Usage:  python3 per_molecule_table.py
# ---------------------------------------------------------------------------
import os
import csv
import numpy as np

from sampling_robustness import (OBS, RAV_DT, SS_DT, SS_SYSTEMS, SS_BRANCH, SPIKES,
                                 N_LAST, load_rav, load_ss, HERE)

NAMES = ("hip", "knee", "ankle", "ntd")
PARENT = {("spike_open_prot_glyc_amarolab", 4): 1, ("spike_open_prot_glyc_amarolab", 5): 3,
          ("spike_open_prot_glyc_amarolab", 6): 2, ("spike_mutant_prot_glyc_amarolab", 4): 2,
          ("spike_mutant_prot_glyc_amarolab", 5): 2, ("spike_mutant_prot_glyc_amarolab", 6): 3}


def main():
    obs = [o for o in OBS if o[0] in NAMES]
    rav = {o[0]: load_rav(o[3], o[5]) for o in obs}
    ss = {o[0]: load_ss(o[4], o[5]) for o in obs}

    header = ["system", "molecule", "branched_from", "start_ns", "length_ns", "mean_sd_window"]
    for n in NAMES:
        header += [f"{n}_frame0", f"{n}_last", f"{n}_mean", f"{n}_sd"]
    rows = []
    for k, i in enumerate(SPIKES):
        row = ["RAV", f"spike {i}", "", 0.0, rav["hip"].shape[1] * RAV_DT, "last 352 ns"]
        for n in NAMES:
            x = rav[n][k]
            w = x[-N_LAST:]
            row += [x[0], x[-1], w.mean(), w.std(ddof=1)]
        rows.append(row)
    k = 0
    for name, nrep in SS_SYSTEMS:
        label = name.split("_")[1]                            # open / closed / mutant
        for r in range(1, nrep + 1):
            key = (name, r)
            par = f"rep {PARENT[key]}" if key in PARENT else ""
            row = [f"single spike ({label})", f"rep {r}", par, SS_BRANCH.get(key, 0.0),
                   len(ss["hip"][k]) * SS_DT, "whole replica"]
            for n in NAMES:
                x = ss[n][k]
                row += [x[0], x[-1], x.mean(), x.std(ddof=1)]
            rows.append(row)
            k += 1

    out = os.path.join(HERE, "per_molecule_values.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in rows:
            w.writerow([f"{v:.2f}" if isinstance(v, float) else v for v in row])
    print(f"wrote {out} ({len(rows)} molecules)")


if __name__ == "__main__":
    main()
