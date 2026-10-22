# Developed and written by Lorenzo Casalino (UC San Diego).
# =============================================================================
#  sampling_robustness.py
#
#  Sampling and temporal-robustness analysis of the spike conformational
#  observables (hinge tilting angles: ankle / hip / knee; NTD "breathing"
#  triangle area; RBD-CH distance) for the revision SI.
#
#  Reads the per-trajectory .txt series produced by the VMD scripts in
#  ../ (results/<ANALYSIS>/) and reports, for each observable:
#    1. integrated autocorrelation time (Sokal automatic window) per spike /
#       replica and the resulting number of effectively independent samples;
#    2. variance decomposition: between-spike (molecule-to-molecule) vs
#       within-spike (temporal) variability, intraclass correlation (ICC);
#    3. uncertainty of the ensemble mean, with spikes as independent units,
#       cross-checked by Flyvbjerg-Petersen block averaging of each spike;
#    4. sequential time blocks: block means +- SEM across spikes, a paired
#       per-spike drift test, and Jensen-Shannon distances between block
#       distributions against a spike-bootstrap noise floor;
#    5. RAV vs single-spike comparison of mean and KDE mode, with
#       hierarchical bootstrap CIs (resampling spikes / replicas).
#
#  Usage:  python3 sampling_robustness.py [results_dir]
#  Output: sampling_robustness_summary.txt, *.csv and SI figure PNG/PDF,
#          written next to this script.
# =============================================================================
import os
import sys

import numpy as np
from scipy import stats
from scipy.spatial.distance import jensenshannon

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "results")

RAV_DT = 0.2                      # ns/frame
SS_DT = 0.1                       # ns/frame
N_LAST = 1760                     # last 352 ns of each RAV spike (as in Fig. 5)
SPIKES = [i for i in range(30) if i != 6]
SS_SYSTEMS = [("spike_open_prot_glyc_amarolab", 6),
              ("spike_closed_prot_glyc_amarolab", 3),
              ("spike_mutant_prot_glyc_amarolab", 6)]
# replicas 4-6 of the open and mutant systems were branched from replicas 1-3.
# Start time (ns) on the parent's clock = 0-based frame index x 0.1 ns, from
# SPIKE_CLEAVED/OPEN/4-MD/extractframes.tcl and SPIKE_MUTATED_CLEAVED/4-MD/extractframes.tcl
# (frames counted from the start of pr.4, as in the analysed time series):
#   open   rep4 <- rep1 frame 2716 (MIN), rep5 <- rep3 frame 715 (MIDDLE), rep6 <- rep2 frame 3118 (MAX)
#   mutant rep4 <- rep2 frame 2523 (MIN), rep5 <- rep2 frame 859 (MIDDLE), rep6 <- rep3 frame 2088 (MAX;
#   the restart file is named "r2" but the script loads REPLICA_3)
SS_BRANCH = {("spike_open_prot_glyc_amarolab", 4): 271.6,
             ("spike_open_prot_glyc_amarolab", 5): 71.5,
             ("spike_open_prot_glyc_amarolab", 6): 311.8,
             ("spike_mutant_prot_glyc_amarolab", 4): 252.3,
             ("spike_mutant_prot_glyc_amarolab", 5): 85.9,
             ("spike_mutant_prot_glyc_amarolab", 6): 208.8}
SS_OFFSET = [SS_BRANCH.get((name, r), 0.0) for name, n in SS_SYSTEMS for r in range(1, n + 1)]
N_BOOT = 2000
RNG = np.random.RandomState(20261016)

# name, label, unit, RAV folder, SS folder, column, colour (as in the paper figures)
OBS = [
    ("hip",   "Hip tilt angle",      "°",  "HIP_TILTING_ANGLE",   "HIP_TILTING_ANGLE_SINGLESPIKE",   1, (0.85975, 0.360588, 0.406917)),
    ("knee",  "Knee tilt angle",     "°",  "KNEE_TILTING_ANGLE",  "KNEE_TILTING_ANGLE_SINGLESPIKE",  1, (0.91, 0.774, 0.0)),
    ("ankle", "Ankle tilt angle",    "°",  "ANKLE_TILTING_ANGLE", "ANKLE_TILTING_ANGLE_SINGLESPIKE", 1, (0.977856, 0.602051, 0.241387)),
    ("ntd",   "NTD triangle area",   "Å²", "NTD_AREA_RAV",        "NTD_AREA_SINGLESPIKE",            1, (0.828, 0.657, 1.0)),
    ("rbdA",  "RBD–CH distance (A)", "Å",  "RBD_CH_DIST_RAV",     "RBD_CH_DIST_SINGLESPIKE",         1, (0.4, 0.4, 0.4)),
    ("rbdB",  "RBD–CH distance (B)", "Å",  "RBD_CH_DIST_RAV",     "RBD_CH_DIST_SINGLESPIKE",         2, (0.4, 0.4, 0.4)),
    ("rbdC",  "RBD–CH distance (C)", "Å",  "RBD_CH_DIST_RAV",     "RBD_CH_DIST_SINGLESPIKE",         3, (0.4, 0.4, 0.4)),
]


# ---------------------------------------------------------------- loading ---
def load_rav(folder, col):
    return np.array([np.loadtxt(f"{RES}/{folder}/spike{i}.{folder}.txt", usecols=col)
                     for i in SPIKES])                       # (29, 2530)


def load_ss(folder, col):
    reps = []
    for name, n in SS_SYSTEMS:
        for r in range(1, n + 1):
            reps.append(np.loadtxt(f"{RES}/{folder}/rep{r}.{folder}.{name}.txt", usecols=col))
    return reps                                               # list of 15 arrays


# ------------------------------------------------------------- statistics ---
def tau_int(x, dt, c=5.0):
    """Integrated autocorrelation time (ns), Sokal automatic windowing."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    if np.allclose(x, 0):
        return dt * 0.5
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf /= acf[0]
    tau = 0.5
    for m in range(1, n):
        tau += acf[m]
        if m >= c * tau:
            break
    return max(tau, 0.5) * dt


def fp_block_sem(x):
    """Flyvbjerg-Petersen block-averaging SEM: max over the doubling sequence
    (a conservative plateau estimate)."""
    x = np.asarray(x, float)
    sems = []
    while len(x) >= 16:
        n = len(x)
        sems.append(np.std(x, ddof=1) / np.sqrt(n))
        x = 0.5 * (x[: n // 2 * 2 : 2] + x[1: n // 2 * 2 : 2])
    return max(sems)


def kde_mode(x, grid):
    k = stats.gaussian_kde(x)
    return grid[np.argmax(k(grid))]


def hist_js(a, b, edges):
    pa, _ = np.histogram(a, edges)
    pb, _ = np.histogram(b, edges)
    return jensenshannon(pa + 1e-12, pb + 1e-12, base=2)


# ---------------------------------------------------------------- analysis ---
def analyse(o):
    name, label, unit, fr, fs, col, _ = o
    full = load_rav(fr, col)
    win = full[:, -N_LAST:]
    ss = load_ss(fs, col)
    out = {"name": name, "label": label, "unit": unit, "full": full, "win": win, "ss": ss}

    # 1. correlation times / effective samples
    tau = np.array([tau_int(x, RAV_DT) for x in win])
    tau_ss = np.array([tau_int(x, SS_DT) for x in ss])
    out["tau"], out["tau_ss"] = tau, tau_ss
    out["neff"] = np.sum(N_LAST * RAV_DT / (2 * tau))
    out["neff_ss"] = np.sum([len(x) * SS_DT / (2 * t) for x, t in zip(ss, tau_ss)])

    # 2. variance decomposition (one-way random effects, spikes as groups)
    means = win.mean(1)
    within_var = win.var(1, ddof=1).mean()
    between_var = means.var(ddof=1)
    out["means"], out["sds"] = means, win.std(1, ddof=1)
    out["sd_between"], out["sd_within"] = np.sqrt(between_var), np.sqrt(within_var)
    # ICC corrected for the finite-sampling noise in each spike mean
    noise = np.mean([fp_block_sem(x) ** 2 for x in win])
    out["icc"] = max(between_var - noise, 0) / (max(between_var - noise, 0) + within_var)
    # persistence of spike-specific offsets: correlation of per-spike means
    # between window halves and between the first and last 88-ns blocks
    out["r_halves"] = stats.pearsonr(win[:, :N_LAST // 2].mean(1), win[:, N_LAST // 2:].mean(1))[0]
    out["r_ends"] = stats.pearsonr(win[:, :N_LAST // 4].mean(1), win[:, -N_LAST // 4:].mean(1))[0]

    # 2b. per-molecule fluctuation amplitude over the same time interval:
    # SD within each RAV spike's analysed window (154-506 ns) vs SD within the same
    # 154-506 ns of each single-spike replica on the common clock from the initial
    # structure (branched replicas offset by SS_BRANCH; their segment can be shorter).
    # Check: SD over the last 352 ns of each replica (same length, different stage).
    t_lo, t_hi = (full.shape[1] - N_LAST) * RAV_DT, full.shape[1] * RAV_DT
    seg = []
    for x, off in zip(ss, SS_OFFSET):
        t = off + np.arange(len(x)) * SS_DT
        seg.append(x[(t >= t_lo) & (t < t_hi)])
    n_ss = int(round(N_LAST * RAV_DT / SS_DT))
    out["fsd_rav"] = win.std(1, ddof=1)
    out["fsd_ss"] = np.array([y.std(ddof=1) for y in seg])
    out["fsd_ss_len"] = np.array([len(y) * SS_DT for y in seg])
    out["fsd_p"] = stats.mannwhitneyu(out["fsd_rav"], out["fsd_ss"], alternative="two-sided").pvalue
    out["fsd_ss_last"] = np.array([x[-n_ss:].std(ddof=1) for x in ss])
    out["fsd_p_last"] = stats.mannwhitneyu(out["fsd_rav"], out["fsd_ss_last"], alternative="two-sided").pvalue

    # 3. uncertainty of the ensemble mean
    out["grand"] = win.mean()
    out["sem_spikes"] = means.std(ddof=1) / np.sqrt(len(means))
    out["sem_fp"] = np.sqrt(np.sum([fp_block_sem(x) ** 2 for x in win])) / len(win)
    out["sem_naive"] = win.std(ddof=1) / np.sqrt(win.size)

    # 4a. sequential blocks over the analysed window (4 x 88 ns) and the full run (5 x ~101 ns)
    for key, arr, nb in (("blk_win", win, 4), ("blk_full", full, 5)):
        L = arr.shape[1] // nb
        bm = np.array([arr[:, b * L:(b + 1) * L].mean(1) for b in range(nb)])   # (nb, 29)
        out[key] = {"L": L, "mean": bm.mean(1),
                    "sem": bm.std(1, ddof=1) / np.sqrt(bm.shape[1]), "per_spike": bm}
    # 4b. paired drift test within the analysed window: per-spike linear slope
    t = np.arange(N_LAST) * RAV_DT
    slopes = np.array([np.polyfit(t, x, 1)[0] for x in win]) * 100      # per 100 ns
    out["slope"], out["slope_sem"] = slopes.mean(), slopes.std(ddof=1) / np.sqrt(len(slopes))
    out["slope_p"] = stats.ttest_1samp(slopes, 0).pvalue
    out["drift_352"] = slopes.mean() * 3.52
    # 4c. JS distance between the first and last window block vs spike-bootstrap noise floor
    L = out["blk_win"]["L"]
    lo, hi = np.percentile(win, [0.5, 99.5])
    edges = np.linspace(lo, hi, 41)
    b0, b3 = win[:, :L], win[:, 3 * L:4 * L]
    out["js_first_last"] = hist_js(b0.ravel(), b3.ravel(), edges)
    # excluded relaxation phase (0-154 ns) vs last window block, for reference (not a matched block)
    excl = np.clip(full[:, :full.shape[1] - N_LAST].ravel(), edges[0], edges[-1])
    out["js_excl_last"] = hist_js(excl, b3.ravel(), edges)
    nf = []
    for _ in range(500):
        i1 = RNG.randint(0, len(win), len(win))
        i2 = RNG.randint(0, len(win), len(win))
        blk = RNG.randint(0, 4)
        nf.append(hist_js(win[i1, blk * L:(blk + 1) * L].ravel(),
                          win[i2, blk * L:(blk + 1) * L].ravel(), edges))
    out["js_floor"] = np.percentile(nf, 95)

    # 5. RAV vs single-spike, hierarchical bootstrap over spikes / replicas
    grid = np.linspace(min(win.min(), min(x.min() for x in ss)),
                       max(win.max(), max(x.max() for x in ss)), 600)
    ss_all = np.concatenate(ss)
    out["mode_rav"], out["mode_ss"] = kde_mode(win.ravel(), grid), kde_mode(ss_all, grid)
    out["mean_ss"] = ss_all.mean()
    bm_r, bm_s, bmo_r, bmo_s = [], [], [], []
    for _ in range(N_BOOT):
        r = win[RNG.randint(0, len(win), len(win))].ravel()
        s = np.concatenate([ss[j] for j in RNG.randint(0, len(ss), len(ss))])
        bm_r.append(r.mean()); bm_s.append(s.mean())
        if _ < 300:   # KDE is slower; 300 resamples are enough for a 95% CI
            bmo_r.append(kde_mode(r[::4], grid)); bmo_s.append(kde_mode(s[::8], grid))
    d = np.array(bm_r) - np.array(bm_s)
    out["dmean"], out["dmean_ci"] = out["grand"] - out["mean_ss"], np.percentile(d, [2.5, 97.5])
    out["mode_rav_ci"], out["mode_ss_ci"] = np.percentile(bmo_r, [2.5, 97.5]), np.percentile(bmo_s, [2.5, 97.5])
    out["ss_rep_means"] = np.array([x.mean() for x in ss])
    return out


# ----------------------------------------------------------------- report ---
def report(R):
    lines = []
    p = lines.append
    p("Sampling and temporal robustness of spike conformational observables")
    p(f"RAV: {len(SPIKES)} spikes, last {N_LAST} frames ({N_LAST*RAV_DT:.0f} ns) each "
      f"= {len(SPIKES)*N_LAST*RAV_DT/1000:.2f} us aggregate")
    nss = sum(len(x) for x in R[0]["ss"])
    p(f"Single-spike: 15 replicas (open 6, closed 3, mutant 6), {nss} frames "
      f"= {nss*SS_DT/1000:.2f} us aggregate")
    p("All uncertainties are SEM (or 95% CI) with spikes / replicas as the independent unit.\n")
    for o in R:
        u = o["unit"]
        p(f"=== {o['label']} ({u}) ===")
        p(f"  tau_int RAV   : median {np.median(o['tau']):.1f} ns "
          f"(IQR {np.percentile(o['tau'],25):.1f}-{np.percentile(o['tau'],75):.1f}; "
          f"range {o['tau'].min():.1f}-{o['tau'].max():.1f})  -> N_eff ~ {o['neff']:.0f}")
        p(f"  tau_int SS    : median {np.median(o['tau_ss']):.1f} ns "
          f"(range {o['tau_ss'].min():.1f}-{o['tau_ss'].max():.1f})  -> N_eff ~ {o['neff_ss']:.0f}")
        p(f"  RAV mean      : {o['grand']:.2f} ± {o['sem_spikes']:.2f} (SEM over spikes); "
          f"block-averaging SEM {o['sem_fp']:.2f}; naive per-frame SEM {o['sem_naive']:.3f}")
        p(f"  variability   : between-spike SD {o['sd_between']:.2f}, within-spike SD "
          f"{o['sd_within']:.2f}, ICC {o['icc']:.2f}")
        p(f"  fluctuation   : mean per-trajectory SD, RAV {o['fsd_rav'].mean():.2f} ± "
          f"{o['fsd_rav'].std(ddof=1)/np.sqrt(len(o['fsd_rav'])):.2f} (n={len(o['fsd_rav'])}) vs SS "
          f"{o['fsd_ss'].mean():.2f} ± {o['fsd_ss'].std(ddof=1)/np.sqrt(len(o['fsd_ss'])):.2f} "
          f"(n={len(o['fsd_ss'])}); ratio {o['fsd_rav'].mean()/o['fsd_ss'].mean():.2f}, "
          f"Mann-Whitney p = {o['fsd_p']:.3g} [SS over 154-506 ns, segments "
          f"{o['fsd_ss_len'].min():.0f}-{o['fsd_ss_len'].max():.0f} ns]")
        p(f"                  check, SS last 352 ns: {o['fsd_ss_last'].mean():.2f}, "
          f"p = {o['fsd_p_last']:.3g}")
        p(f"  persistence   : r(per-spike mean, 1st vs 2nd half) = {o['r_halves']:.2f}; "
          f"r(first vs last 88 ns) = {o['r_ends']:.2f}")
        p(f"  per-spike means range {o['means'].min():.2f}-{o['means'].max():.2f}; "
          f"SS replica means range {o['ss_rep_means'].min():.2f}-{o['ss_rep_means'].max():.2f}")
        bw, bf = o["blk_win"], o["blk_full"]
        p("  blocks (last 352 ns, 4 x {:.0f} ns): ".format(bw['L']*RAV_DT) +
          ", ".join(f"{m:.2f}±{s:.2f}" for m, s in zip(bw['mean'], bw['sem'])))
        p("  blocks (full 506 ns, 5 x {:.0f} ns): ".format(bf['L']*RAV_DT) +
          ", ".join(f"{m:.2f}±{s:.2f}" for m, s in zip(bf['mean'], bf['sem'])))
        p(f"  drift in window: {o['slope']:+.3f} ± {o['slope_sem']:.3f} per 100 ns "
          f"(p = {o['slope_p']:.3g}); total over 352 ns {o['drift_352']:+.2f} "
          f"= {abs(o['drift_352'])/o['sd_within']:.2f} within-spike SD")
        p(f"  JS distance first vs last 88-ns block: {o['js_first_last']:.3f} "
          f"(spike-resampling 95% noise floor {o['js_floor']:.3f})")
        p(f"  RAV vs SS     : mean {o['grand']:.2f} vs {o['mean_ss']:.2f}, "
          f"diff {o['dmean']:+.2f} [95% CI {o['dmean_ci'][0]:+.2f}, {o['dmean_ci'][1]:+.2f}]")
        p(f"                  KDE mode {o['mode_rav']:.2f} [{o['mode_rav_ci'][0]:.2f}, {o['mode_rav_ci'][1]:.2f}] "
          f"vs {o['mode_ss']:.2f} [{o['mode_ss_ci'][0]:.2f}, {o['mode_ss_ci'][1]:.2f}]\n")
    txt = "\n".join(lines)
    open(os.path.join(HERE, "sampling_robustness_summary.txt"), "w").write(txt + "\n")
    print(txt)

    # machine-readable tables
    with open(os.path.join(HERE, "summary_table.csv"), "w") as f:
        f.write("observable,unit,tau_med_ns,tau_q1,tau_q3,neff_rav,neff_ss,mean_rav,sem_spikes,"
                "sem_blockavg,sd_between,sd_within,icc,drift_per100ns,drift_sem,drift_p,"
                "js_first_last,js_floor,mean_ss,dmean,dmean_lo,dmean_hi,mode_rav,mode_rav_lo,"
                "mode_rav_hi,mode_ss,mode_ss_lo,mode_ss_hi,fluct_sd_rav,fluct_sd_ss,fluct_p\n")
        for o in R:
            f.write(",".join(str(v) for v in [
                o["name"], o["unit"], *np.round([np.median(o["tau"]), *np.percentile(o["tau"], [25, 75]),
                o["neff"], o["neff_ss"], o["grand"], o["sem_spikes"], o["sem_fp"], o["sd_between"],
                o["sd_within"], o["icc"], o["slope"], o["slope_sem"]], 4), f"{o['slope_p']:.3g}",
                *np.round([o["js_first_last"], o["js_floor"], o["mean_ss"], o["dmean"], *o["dmean_ci"],
                o["mode_rav"], *o["mode_rav_ci"], o["mode_ss"], *o["mode_ss_ci"],
                o["fsd_rav"].mean(), o["fsd_ss"].mean()], 4), f"{o['fsd_p']:.3g}"]) + "\n")
    with open(os.path.join(HERE, "per_spike_table.csv"), "w") as f:
        f.write("spike," + ",".join(f"{o['name']}_mean,{o['name']}_sd,{o['name']}_tau_ns" for o in R) + "\n")
        for k, s in enumerate(SPIKES):
            f.write(f"{s}," + ",".join(f"{o['means'][k]:.3f},{o['sds'][k]:.3f},{o['tau'][k]:.2f}"
                                       for o in R) + "\n")


# ----------------------------------------------------------------- figure ---
def figure(R):
    rows = [o for o in R if o["name"] in ("hip", "knee", "ankle", "ntd")]
    plt.rcParams.update({"font.sans-serif": "Libre Franklin", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    # figure 1: sampling of the RAV trajectory (time course, correlation time, block stability)
    # figure 2: RAV vs single spike (spike-to-spike variability, fluctuation, evolution from t = 0)
    fig1, ax1 = plt.subplots(len(rows), 3, figsize=(12, 3.0 * len(rows)),
                             gridspec_kw={"width_ratios": [1.25, 1, 1.25]})
    fig2, ax2 = plt.subplots(len(rows), 3, figsize=(12, 3.0 * len(rows)),
                             gridspec_kw={"width_ratios": [1, 0.75, 1.25]})
    axes = np.concatenate([ax1, ax2], axis=1)
    blues = plt.get_cmap("plasma_r")(np.linspace(0.15, 1.0, 4))
    t0 = (2530 - N_LAST) * RAV_DT
    for r, o in enumerate(rows):
        c, u = OBS[[x[0] for x in OBS].index(o["name"])][6], o["unit"]
        full = o["full"]
        # (a) time course: 10-ns block means, median and IQR across spikes
        nb = 10
        L = int(nb / RAV_DT)
        k = full.shape[1] // L
        bm = full[:, :k * L].reshape(len(full), k, L).mean(2)
        tc = (np.arange(k) + 0.5) * nb
        ax = axes[r, 0]
        for x in bm:
            ax.plot(tc, x, color=c, lw=0.5, alpha=0.35)
        ax.fill_between(tc, *np.percentile(bm, [25, 75], axis=0), color=c, alpha=0.25, lw=0)
        ax.plot(tc, bm.mean(0), color="k", lw=2, label="mean of 29 spikes")
        # sequential ~101-ns block means over the full run, ± SEM across spikes
        bf = o["blk_full"]
        bc = (np.arange(len(bf["mean"])) + 0.5) * bf["L"] * RAV_DT
        ax.errorbar(bc, bf["mean"], yerr=bf["sem"], fmt="o", ms=6, mfc="white", mec="k",
                    mew=1.5, ecolor="k", elinewidth=1.5, capsize=3, zorder=5,
                    label=f"{bf['L'] * RAV_DT:.0f}-ns block mean ± SEM")
        ax.axvspan(0, t0, color="0.85", lw=0, zorder=0)
        ax.text(t0 / 2, 0.97, "excluded", transform=ax.get_xaxis_transform(),
                ha="center", va="top", color="0.4", fontsize=8)
        ax.set_xlim(0, 506)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel(f"{o['label']} ({u})")
        if r == 0:
            ax.legend(frameon=False, loc="lower right", fontsize=8)
        # (b) autocorrelation-time distribution
        ax = axes[r, 1]
        tb = np.arange(0, 125, 5)
        ax.hist(o["tau"], bins=tb, color=c, alpha=0.75, edgecolor="white", lw=1,
                label="RAV spikes")
        ax.hist(o["tau_ss"], bins=tb, histtype="step", color="k", lw=1.2,
                label="single-spike simulations")
        ax.set_xlim(0, 120)
        ax.axvline(np.median(o["tau"]), color="k", ls="--", lw=1)
        ax.set_xlabel(r"$\tau_{\mathrm{int}}$ (ns)")
        ax.set_ylabel("Count")
        ax.text(0.97, 0.95, f"median {np.median(o['tau']):.0f} ns\n"
                f"$N_{{\\mathrm{{eff}}}}$ ≈ {o['neff']:.0f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8)
        if r == 0:
            ax.legend(frameon=False, loc="center right", fontsize=8)
        # (c) sequential 88-ns blocks of the analysed window
        ax = axes[r, 2]
        win = o["win"]
        Lw = o["blk_win"]["L"]
        excl = full[:, :full.shape[1] - N_LAST]
        lo, hi = np.percentile(np.r_[win.ravel(), excl.ravel()], [0.2, 99.8])
        g = np.linspace(lo, hi, 400)
        ax.plot(g, stats.gaussian_kde(excl.ravel())(g), color="0.45", lw=1.5, ls="--",
                label=f"0–{t0:.0f} ns (excluded)")
        for b in range(4):
            x = win[:, b * Lw:(b + 1) * Lw].ravel()
            a0 = t0 + b * Lw * RAV_DT
            ax.plot(g, stats.gaussian_kde(x)(g), color=blues[b], lw=2,
                    label=f"{a0:.0f}–{a0 + Lw * RAV_DT:.0f} ns")
        ax.set_ylim(0, 1.3 * ax.get_ylim()[1])                 # headroom for the legend
        ax.set_xlabel(f"{o['label']} ({u})")
        ax.set_ylabel("Density")
        ax.legend(frameon=False, fontsize=7.5, loc="upper right")
        # (d) molecule-to-molecule variability: per-spike mean ± SD (sorted) and value at frame 0
        ax = axes[r, 3]
        order = np.argsort(o["means"])
        y = np.arange(len(order))
        ax.errorbar(o["means"][order], y, xerr=o["sds"][order], fmt="o", ms=3.5, color=c,
                    ecolor=c, elinewidth=1, capsize=0, label="RAV spike mean ± SD (row colour)")
        ax.plot(full[order, 0], y, "|", ms=6, mew=1.5, color="k", label="value at first frame (t = 0)")
        ax.axvline(o["grand"], color="k", lw=1, ls="--", label="RAV ensemble mean")
        ax.set_yticks([])
        ax.set_ylabel("Spikes (sorted by mean)")
        ax.set_xlabel(f"{o['label']} ({u})")
        ax.text(0.97, 0.03, f"ICC = {o['icc']:.2f}\nr(halves) = {o['r_halves']:.2f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
        ax.spines["left"].set_visible(False)
        if r == len(rows) - 1:
            ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=1,
                      bbox_to_anchor=(0.5, -0.22))
        # (e) per-molecule fluctuation amplitude over the same 154-506 ns interval
        ax = axes[r, 4]
        for xpos, vals, kw in ((0, o["fsd_rav"], dict(color=c, alpha=0.8, mec="white", mew=0.5)),
                               (1, o["fsd_ss"], dict(mfc="none", mec="k"))):
            jit = RNG.uniform(-0.15, 0.15, len(vals))
            ax.plot(xpos + jit, vals, "o" if xpos == 0 else "s", ms=4, ls="none",
                    label=("RAV spike (row colour)" if xpos == 0 else "single-spike simulation"), **kw)
            sem = vals.std(ddof=1) / np.sqrt(len(vals))
            ax.errorbar(xpos + 0.32, vals.mean(), yerr=sem, fmt="_", ms=12, mew=2,
                        color="k", elinewidth=1.5, capsize=0, label="mean ± SEM" if xpos == 0 else None)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["RAV\nspikes", "single-spike\nsimulations"])
        ax.set_xlim(-0.5, 1.6)
        ax.set_ylim(0, 1.3 * max(o["fsd_rav"].max(), o["fsd_ss"].max()))
        ax.set_ylabel(f"SD over 154–506 ns ({u})")
        ptxt = f"p = {o['fsd_p']:.2g}" if o["fsd_p"] >= 0.001 else "p < 0.001"
        ax.text(0.5, 0.97, f"ratio {o['fsd_rav'].mean() / o['fsd_ss'].mean():.2f}, {ptxt}",
                transform=ax.transAxes, ha="center", va="top", fontsize=8)
        if r == len(rows) - 1:
            ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=1,
                      bbox_to_anchor=(0.5, -0.22))
        # (f) evolution from the initial structure (10-ns block means), RAV spikes vs SS replicas
        ax = axes[r, 5]
        ax.plot(tc, bm.mean(0), color=c, lw=2, zorder=5, label="RAV, mean of 29 spikes (row colour)")
        ax.fill_between(tc, bm.mean(0) - bm.std(0, ddof=1) / np.sqrt(len(bm)),
                        bm.mean(0) + bm.std(0, ddof=1) / np.sqrt(len(bm)), color=c, alpha=0.35, lw=0, zorder=5)
        # single-spike replicas on a common clock from the initial structure: reps 1-3 start at
        # t = 0, branched reps 4-6 enter at their branch time; mean ± SEM over the replicas
        # running in each 10-ns block (dotted where < 9 replicas, cut where < 3)
        Ls = int(nb / SS_DT)
        sb = [x[:len(x) // Ls * Ls].reshape(-1, Ls).mean(1) for x in o["ss"]]
        sh = [int(round(off / nb)) for off in SS_OFFSET]
        kmax = max(len(x) + k0 for x, k0 in zip(sb, sh))
        sm_t, se_t, n_t = [], [], []
        for j in range(kmax):
            v = np.array([x[j - k0] for x, k0 in zip(sb, sh) if 0 <= j - k0 < len(x)])
            n_t.append(len(v))
            sm_t.append(v.mean())
            se_t.append(v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0)
        sm_t, se_t, n_t = map(np.array, (sm_t, se_t, n_t))
        keep = n_t >= 3
        sm_t, se_t, n_t = sm_t[keep], se_t[keep], n_t[keep]
        ts = (np.arange(kmax)[keep] + 0.5) * nb
        many = n_t >= 9
        ax.fill_between(ts, sm_t - se_t, sm_t + se_t, color="0.5", alpha=0.3, lw=0)
        ax.plot(np.where(many, ts, np.nan), sm_t, color="k", lw=2,
                label="single-spike simulations, mean of 9–15")
        few = ~many | np.r_[~many[1:], False]             # overlap one block so the line is unbroken
        ax.plot(np.where(few, ts, np.nan), sm_t, color="k", lw=1.5, ls=":",
                label="single-spike simulations, mean of 3–8")
        ax.axvspan(t0, 506, color=c, alpha=0.08, lw=0, zorder=0)
        ax.text((t0 + 506) / 2, 0.97, "RAV analysis\nwindow", transform=ax.get_xaxis_transform(),
                ha="center", va="top", color="0.4", fontsize=8)
        ax.set_xlim(0, 1000)
        ax.set_xlabel("Time from initial structure (ns)")
        ax.set_ylabel(f"{o['label']} ({u})")
        if r == len(rows) - 1:
            ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=1,
                      bbox_to_anchor=(0.5, -0.22))
    titles = ["Time course (29 RAV spikes)", "Autocorrelation time",
              "Stability across time blocks", "Spike-to-spike variability",
              "Fluctuation within each trajectory", "Evolution from the initial structure"]
    for j, (t, title) in enumerate(zip("abcabc", titles)):
        axes[0, j].set_title(f"{t}   {title}", loc="left", fontweight="bold", fontsize=11, pad=10)
    fig1.tight_layout()
    fig2.tight_layout()
    fig2.subplots_adjust(bottom=0.11)
    for ext in ("png", "pdf"):
        fig1.savefig(os.path.join(HERE, f"SI_spike_sampling.{ext}"), dpi=300)
        fig2.savefig(os.path.join(HERE, f"SI_spike_rav_vs_single.{ext}"), dpi=300)


if __name__ == "__main__":
    R = [analyse(o) for o in OBS]
    report(R)
    figure(R)
