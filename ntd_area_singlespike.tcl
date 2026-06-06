# Developed and written by Lorenzo Casalino (UC San Diego).
# NTD triangle area, single-spike dataset (open / closed / mutant).
#   vmd -dispdev text -e ntd_area_singlespike.tcl
set _dir [file dirname [file normalize [info script]]]
source [file join $_dir rav_common.tcl]
source [file join $_dir rav_selections.tcl]
source [file join $_dir rav_paths.tcl]

rav_run_replicas [list \
    names    $SS_NAMES \
    inppaths $SS_INPPATHS \
    nreps    $SS_NREPS \
    outpath  $SS_OUTPATH \
    analysis "NTD_AREA_SINGLESPIKE" \
    dt       $SS_DT \
    catfile  "concatenated_reps_area.txt" \
    make     ntd_make \
    measure  ntd_measure \
    emit     rav_emit_columns]

quit
