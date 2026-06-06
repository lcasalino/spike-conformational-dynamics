# Developed and written by Lorenzo Casalino (UC San Diego).
# Knee stem tilting angle, RAV dataset (per spike).
#   vmd -dispdev text -e knee_tilting_rav.tcl
set _dir [file dirname [file normalize [info script]]]
source [file join $_dir rav_common.tcl]
source [file join $_dir rav_selections.tcl]
source [file join $_dir rav_paths.tcl]

rav_run_spikes [list \
    name     $RAV_NAME \
    analysis "KNEE_TILTING_ANGLE" \
    inppath  $RAV_INPPATH \
    outpath  $RAV_OUTPATH \
    nspike   $RAV_NSPIKE \
    skip     $RAV_SKIP \
    dt       $RAV_DT \
    catfile  "concatenated_spikes_tilting.txt" \
    make     knee_make \
    measure  knee_measure \
    emit     rav_emit_tilting]

quit
