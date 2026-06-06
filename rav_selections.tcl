# Developed and written by Lorenzo Casalino (UC San Diego).
# =============================================================================
#  rav_selections.tcl
#
#  Atom selections and per-frame measurements for each observable. These are
#  the ONLY analysis-specific definitions, and they are dataset-agnostic: the
#  driver supplies the segment names, so the identical geometric definition is
#  applied to both the RAV and the single-spike data (which differ only in how
#  the chains A/B/C are named).
#
#  Each observable provides:
#    <obs>_make    {mol segnames chains} -> list of selection-text strings
#    <obs>_measure {sels frame}          -> list of values (element 0 drives extrema)
#
#  NOTE: _make returns selection *text*, not atomselect handles. VMD deletes
#  atomselect handles when the proc that created them returns, so the handles
#  must be created by the caller (rav_measure_series) whose call frame stays
#  alive for the whole measurement loop. _measure receives the live handles.
#
#  'segnames' is the full six-segment selection for the spike; 'chains' is a
#  three-element list giving the segment name(s) of protomers A, B and C.
#
#  Sourced by the runner scripts together with rav_common.tcl.
# =============================================================================


# -----------------------------------------------------------------------------
#  Stem tilting angles (ankle / hip / knee)
#
#  Three CA selections per spike: an origin vertex and the two arms whose
#  angle (at the origin) is reported. Defined over the full segment set, so
#  'chains' is unused here.
# -----------------------------------------------------------------------------

proc _tilt_make {segnames origin_resid arm1_resid arm2_resid} {
    return [list \
        "protein and segname $segnames and resid $origin_resid and name CA" \
        "protein and segname $segnames and resid $arm1_resid and name CA" \
        "protein and segname $segnames and resid $arm2_resid and name CA"]
}

proc _tilt_measure {sels frame} {
    lassign $sels origin arm1 arm2
    $origin frame $frame
    $arm1   frame $frame
    $arm2   frame $frame
    set com_origin [measure center $origin weight mass]
    set com_arm1   [measure center $arm1   weight mass]
    set com_arm2   [measure center $arm2   weight mass]
    return [format "%.2f" [rav_angle $com_arm1 $com_origin $com_arm2]]
}

proc ankle_make {mol segnames chains} { return [_tilt_make $segnames 1213 "1169 to 1206" "1212 to 1239"] }
proc ankle_measure {sels frame}       { return [_tilt_measure $sels $frame] }

proc hip_make {mol segnames chains}   { return [_tilt_make $segnames "1136 to 1140" "1141 to 1161" "816 to 1135"] }
proc hip_measure {sels frame}         { return [_tilt_measure $sels $frame] }

proc knee_make {mol segnames chains}  { return [_tilt_make $segnames "1161 to 1168" "1169 to 1206" "1141 to 1161"] }
proc knee_measure {sels frame}        { return [_tilt_measure $sels $frame] }


# -----------------------------------------------------------------------------
#  NTD triangle area
#
#  One NTD CA selection per protomer (chains A/B/C, residues 13-291); the area
#  of the triangle formed by the three COMs is reported.
# -----------------------------------------------------------------------------

proc ntd_make {mol segnames chains} {
    lassign $chains a b c
    return [list \
        "name CA and segname $a and (resid 13 to 291)" \
        "name CA and segname $b and (resid 13 to 291)" \
        "name CA and segname $c and (resid 13 to 291)"]
}

proc ntd_measure {sels frame} {
    lassign $sels a b c
    $a frame $frame
    $b frame $frame
    $c frame $frame
    set com_a [measure center $a weight mass]
    set com_b [measure center $b weight mass]
    set com_c [measure center $c weight mass]
    return [format "%.2f" [rav_triangle_area $com_a $com_b $com_c]]
}


# -----------------------------------------------------------------------------
#  RBD -> central-helix (CH) distance
#
#  Per protomer (chains A/B/C) the COM of a set of RBD residues; its distance
#  to the COM of the shared central-helix selection is reported. Three values
#  per frame (A, B, C); chain A drives the extrema.
# -----------------------------------------------------------------------------

proc rbd_make {mol segnames chains} {
    lassign $chains a b c
    set rbd "((resid 375 to 380) or (resid 394 to 404) or (resid 431 to 438) or (resid 508 to 517))"
    set ch  "((resid 747 to 784) or (resid 946 to 967) or (resid 986 to 1034))"
    return [list \
        "name CA and segname $a and $rbd" \
        "name CA and segname $b and $rbd" \
        "name CA and segname $c and $rbd" \
        "name CA and segname $segnames and $ch"]
}

proc rbd_measure {sels frame} {
    lassign $sels a b c ch
    foreach s $sels { $s frame $frame }
    set com_ch [measure center $ch weight mass]
    return [list \
        [format "%.2f" [vecdist [measure center $a weight mass] $com_ch]] \
        [format "%.2f" [vecdist [measure center $b weight mass] $com_ch]] \
        [format "%.2f" [vecdist [measure center $c weight mass] $com_ch]]]
}
