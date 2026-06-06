# Developed and written by Lorenzo Casalino (UC San Diego).
# RBD -> central-helix distance, RAV dataset (per spike).
#   vmd -dispdev text -e rbd_ch_distance_rav.tcl
set _dir [file dirname [file normalize [info script]]]
source [file join $_dir rav_common.tcl]
source [file join $_dir rav_selections.tcl]
source [file join $_dir rav_paths.tcl]

rav_run_spikes [list \
    name     $RAV_NAME \
    analysis "RBD_CH_DIST_RAV" \
    inppath  $RAV_INPPATH \
    outpath  $RAV_OUTPATH \
    nspike   $RAV_NSPIKE \
    skip     $RAV_SKIP \
    dt       $RAV_DT \
    catfile  "concatenated_spikes_distance.txt" \
    make     rbd_make \
    measure  rbd_measure \
    emit     rav_emit_columns]

quit
