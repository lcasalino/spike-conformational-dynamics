#!/usr/bin/env bash
# Developed and written by Lorenzo Casalino (UC San Diego).
# =============================================================================
#  run_all.sh — regenerate analysis data with VMD.
#
#  Output goes to ./results/<ANALYSIS>/ (set in rav_paths.tcl), so the original
#  output folders are left untouched.
#
#  Usage:
#     ./run_all.sh [all|rav|singlespike] [-j N]
#       all|rav|singlespike   which set of analyses to run (default: all)
#       -j N                  run up to N analyses concurrently (default: 1 = serial)
#
#  The 10 analyses are independent and each writes to its own results/<ANALYSIS>/
#  subdir, so they are safe to run in parallel and/or split across machines that
#  share the results/ filesystem. NOTE: each RAV job streams ~2 GB per spike
#  across 29 spikes, so a high -j on one node multiplies I/O and RAM — keep -j
#  modest for the RAV set. Single-spike jobs are much lighter.
#
#  To split across machines instead, just run individual scripts, e.g.:
#     machineA$ vmd -dispdev text -e ankle_tilting_rav.tcl
#     machineB$ vmd -dispdev text -e hip_tilting_rav.tcl
#
#  VMD location: override by exporting VMD, e.g.
#     VMD=/home/lcasalino/Software/vmd-1.9.4a57/bin/vmd ./run_all.sh
# =============================================================================
set -euo pipefail

cd "$(dirname "$(readlink -f "$0")")"

# ---- locate VMD: honour $VMD, else 'vmd' on PATH (e.g. after `module load vmd`),
#      else the local install. First candidate that resolves wins. ----
resolve_vmd() {
    local v="$1"
    [[ -n "$v" && -x "$v" ]] && { printf '%s\n' "$v"; return 0; }
    [[ -n "$v" ]] && command -v "$v" >/dev/null 2>&1 && { command -v "$v"; return 0; }
    return 1
}
VMD_FOUND=""
for cand in "${VMD:-}" vmd /home/lcasalino/Software/vmd-1.9.4a57/bin/vmd; do
    if VMD_FOUND=$(resolve_vmd "$cand"); then break; fi
done
if [[ -z "$VMD_FOUND" ]]; then
    echo "ERROR: VMD not found. Tried \$VMD='${VMD:-}', 'vmd' on PATH, and the local install." >&2
    echo "       Load a module (e.g. 'module load vmd') or set VMD=/path/to/vmd and retry." >&2
    exit 1
fi
VMD="$VMD_FOUND"

# ---- parse args: one positional (which) + optional -j N ----
WHICH="all"
JOBS=1
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j) JOBS="${2:-}"; shift 2 ;;
        -j*) JOBS="${1#-j}"; shift ;;
        all|rav|singlespike) WHICH="$1"; shift ;;
        *) echo "Usage: $0 [all|rav|singlespike] [-j N]" >&2; exit 2 ;;
    esac
done
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: -j must be a positive integer" >&2; exit 2; }

RAV_SCRIPTS=(ankle_tilting_rav.tcl hip_tilting_rav.tcl knee_tilting_rav.tcl ntd_area_rav.tcl rbd_ch_distance_rav.tcl)
SS_SCRIPTS=(ankle_tilting_singlespike.tcl hip_tilting_singlespike.tcl knee_tilting_singlespike.tcl ntd_area_singlespike.tcl rbd_ch_distance_singlespike.tcl)
case "$WHICH" in
    all)         SCRIPTS=("${RAV_SCRIPTS[@]}" "${SS_SCRIPTS[@]}") ;;
    rav)         SCRIPTS=("${RAV_SCRIPTS[@]}") ;;
    singlespike) SCRIPTS=("${SS_SCRIPTS[@]}") ;;
esac

mkdir -p results/logs
echo "VMD: $VMD"
echo "Running ${#SCRIPTS[@]} analysis script(s), up to $JOBS at a time; output -> ./results/  logs -> ./results/logs/"

run_one() {
    local s="$1" log="results/logs/${1%.tcl}.log"
    echo ">>> start $s   (log: $log)"
    if "$VMD" -dispdev text -e "$s" >"$log" 2>&1; then
        echo "<<< OK    $s"
    else
        echo "<<< FAIL  $s — see $log" >&2
        return 1
    fi
}

fail=0
if [[ "$JOBS" -eq 1 ]]; then
    for s in "${SCRIPTS[@]}"; do run_one "$s" || fail=1; done
else
    # bounded-concurrency pool
    for s in "${SCRIPTS[@]}"; do
        run_one "$s" &
        while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n || fail=1; done
    done
    while (( $(jobs -rp | wc -l) > 0 )); do wait -n || fail=1; done
fi

if [[ "$fail" -eq 0 ]]; then
    echo "All requested analyses finished. Data is in ./results/"
else
    echo "One or more analyses FAILED — check results/logs/" >&2
    exit 1
fi
