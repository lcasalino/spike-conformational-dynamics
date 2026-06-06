# Developed and written by Lorenzo Casalino (UC San Diego).
# Ankle stem tilting angle, single-spike dataset (open / closed / mutant).
#   vmd -dispdev text -e ankle_tilting_singlespike.tcl
set _dir [file dirname [file normalize [info script]]]
source [file join $_dir rav_common.tcl]
source [file join $_dir rav_selections.tcl]
source [file join $_dir rav_paths.tcl]

rav_run_replicas [list \
    names    $SS_NAMES \
    inppaths $SS_INPPATHS \
    nreps    $SS_NREPS \
    outpath  $SS_OUTPATH \
    analysis "ANKLE_TILTING_ANGLE_SINGLESPIKE" \
    dt       $SS_DT \
    catfile  "concatenated_reps_tilting.txt" \
    make     ankle_make \
    measure  ankle_measure \
    emit     rav_emit_tilting]

quit
