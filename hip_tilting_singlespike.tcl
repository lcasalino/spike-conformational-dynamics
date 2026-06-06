# Developed and written by Lorenzo Casalino (UC San Diego).
# Hip stem tilting angle, single-spike dataset (open / closed / mutant).
#   vmd -dispdev text -e hip_tilting_singlespike.tcl
set _dir [file dirname [file normalize [info script]]]
source [file join $_dir rav_common.tcl]
source [file join $_dir rav_selections.tcl]
source [file join $_dir rav_paths.tcl]

rav_run_replicas [list \
    names    $SS_NAMES \
    inppaths $SS_INPPATHS \
    nreps    $SS_NREPS \
    outpath  $SS_OUTPATH \
    analysis "HIP_TILTING_ANGLE_SINGLESPIKE" \
    dt       $SS_DT \
    catfile  "concatenated_reps_tilting.txt" \
    make     hip_make \
    measure  hip_measure \
    emit     rav_emit_tilting]

quit
