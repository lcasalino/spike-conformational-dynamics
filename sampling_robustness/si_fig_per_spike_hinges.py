#!/usr/bin/env python3
# ---------------------------------------------------------------------------
#  Supplementary figure: per-spike distributions of the three hinge angles.
#
#  One panel per RAV spike (29 panels), each showing kernel density estimates
#  of the hip, knee and ankle tilt angles over the analysis window (final 1760
#  frames, ~352 ns). Dashed vertical lines mark the modal value of each
#  distribution; solid vertical lines mark the value at the first production
#  frame of that spike. All panels share identical axis ranges.
#
#  Colours are those used for the same observables in Fig. 5 and the other
#  supplementary figures (hip red, knee yellow, ankle orange).
#
#  Usage:  python3 si_fig_per_spike_hinges.py
# ---------------------------------------------------------------------------
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy import stats

from sampling_robustness import OBS, N_LAST, SPIKES, load_rav, HERE

HINGES = ("hip", "knee", "ankle")
LABEL = {"hip": "Hip", "knee": "Knee", "ankle": "Ankle"}
XLIM = (0, 55)
NCOL = 6


def main():
    obs = {o[0]: o for o in OBS}
    col = {h: obs[h][6] for h in HINGES}
    full = {h: load_rav(obs[h][3], obs[h][5]) for h in HINGES}
    win = {h: v[:, -N_LAST:] for h, v in full.items()}

    plt.rcParams.update({"font.sans-serif": "Libre Franklin", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    nrow = int(np.ceil(len(SPIKES) / NCOL))
    fig, axes = plt.subplots(nrow, NCOL, figsize=(3.0 * NCOL, 2.5 * nrow),
                             sharex=True, sharey=True)
    axes = axes.ravel()
    g = np.linspace(*XLIM, 600)

    ymax = 0
    for k in range(len(SPIKES)):
        ax = axes[k]
        t0 = []
        for h in HINGES:
            y = stats.gaussian_kde(win[h][k])(g)
            ymax = max(ymax, y.max())
            ax.fill_between(g, y, color=col[h], alpha=0.35, lw=0)
            ax.plot(g, y, color=col[h], lw=1.4)
            ax.vlines(g[y.argmax()], 0, y.max(), color=col[h], linestyles="dashed", lw=1.2)
            t0.append((full[h][k, 0], col[h]))
        for x, c in t0:                      # value at t = 0: stub at the baseline
            ax.vlines(x, 0, 0.16 * ymax, color=c, lw=2.0, zorder=6)
        ax.text(0.97, 0.93, f"spike {SPIKES[k]}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8.5, fontweight="bold")
        ax.set_xlim(XLIM)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))

    for ax in axes[len(SPIKES):]:
        ax.axis("off")
    for k in range(len(SPIKES)):
        ax = axes[k]
        ax.set_ylim(0, 1.12 * ymax)
        ax.set_xlabel("Tilt angle (°)")
        ax.set_ylabel("Density")
        ax.tick_params(labelbottom=True, labelleft=True, labelsize=8)
        ax.tick_params(which="minor", length=2)

    handles = [plt.Line2D([], [], color=col[h], lw=3, label=LABEL[h]) for h in HINGES]
    handles.append(plt.Line2D([], [], color="0.4", ls="--", lw=1.2, label="modal value"))
    handles.append(plt.Line2D([], [], color="0.4", lw=2.0, label="value at t = 0"))
    fig.legend(handles=handles, frameon=False, fontsize=14, ncol=5,
               loc="lower center", bbox_to_anchor=(0.5, 0.0))

    fig.tight_layout(rect=(0, 0.045, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"SI_per_spike_hinges.{ext}"), dpi=300)
    print("wrote SI_per_spike_hinges.{png,pdf}")


if __name__ == "__main__":
    main()
