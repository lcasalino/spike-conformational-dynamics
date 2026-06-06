# Developed and written by Lorenzo Casalino (UC San Diego).
# =============================================================================
#  rav_paths.tcl
#
#  >>> EDIT THESE PATHS for your environment <<<
#
#  Central configuration of input/output locations and dataset parameters,
#  sourced by every runner script so paths are defined in one place.
# =============================================================================

# ---- where results are written (a sub-folder per analysis is created here) --
# Default: a self-contained "results/" folder next to these scripts, so that
# re-running regenerates all data inside source_code/ without touching the
# original output folders. Set this to an absolute path to write elsewhere.
set _SRCDIR     [file dirname [file normalize [info script]]]
set RAV_OUTPATH [file join $_SRCDIR results]

# ---- RAV dataset (many spikes, one trajectory per spike) --------------------
# Files: $RAV_INPPATH/${RAV_NAME}<NN>.alignedBackbone.{psf,dcd}, NN = 00..29
set RAV_INPPATH "/net/gpfs-amarolab/nwauer/GORDONBELL/RAV_v2/Analysis/Spike_Align/Trajs"
set RAV_NAME    "RAV_v2_all_spike-glycan"
set RAV_NSPIKE  29          ;# highest spike index, inclusive
set RAV_SKIP    {6}         ;# spike 6 does not exist in the RAV system, so it is skipped
set RAV_DT      0.2         ;# ns per saved frame (2530 frames = 506 ns/spike)

# ---- Single-spike dataset (isolated spike: open / closed / mutant) ----------
# Files (per system): <inppath>/<name>.psf and <inppath>/<name>_<rep>.dcd
set SS_OUTPATH  [file join $_SRCDIR results]
set SS_NAMES    {spike_open_prot_glyc_amarolab spike_closed_prot_glyc_amarolab spike_mutant_prot_glyc_amarolab}
set SS_INPPATHS {
    /net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/OPEN
    /net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/CLOSED
    /net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/MUTANT
}
set SS_NREPS    {6 3 6}     ;# replica counts, parallel to SS_NAMES (open closed mutant)
set SS_DT       0.1         ;# ns per saved frame
