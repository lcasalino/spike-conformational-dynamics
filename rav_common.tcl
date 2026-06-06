# =============================================================================
#  rav_common.tcl
#
#  Shared library for the geometric analyses of the SARS-CoV-2 spike
#  ectodomain. The same observables are computed for two datasets:
#
#    * RAV  : the Respiratory Aerosol Virion all-atom simulation, which
#             contains many spikes; analysed spike-by-spike (rav_run_spikes).
#    * SS   : the isolated single-spike simulations (open / closed / mutant),
#             each with several replicas; analysed replica-by-replica
#             (rav_run_replicas).
#
#  This file provides:
#    * geometry primitives        -> rav_angle, rav_triangle_area
#    * output line formatters      -> rav_emit_tilting, rav_emit_columns
#    * per-trajectory helpers      -> rav_measure_series, rav_extrema,
#                                     rav_write_series, rav_write_pdb
#    * two dataset drivers         -> rav_run_spikes, rav_run_replicas
#
#  The science-specific code (atom selections and the measured quantity) lives
#  in rav_selections.tcl as make/measure callbacks that are dataset-agnostic:
#  the driver passes the segment names so that the SAME definition is applied
#  to both datasets.
#
#  Requirements: VMD 1.9.x (used with VMD 1.9.4, embedded Tcl 8.5). No external
#  Tcl packages (tcl::mathfunc::min/max are Tcl 8.5+ built-ins; vecsub/vecdot/
#  veclength/vecdist/measure/atomselect/mol are provided by VMD).
#
#  Developed and written by Lorenzo Casalino (UC San Diego).
# =============================================================================


# -----------------------------------------------------------------------------
#  Geometry primitives
# -----------------------------------------------------------------------------

# rav_angle T1 O T2 -> angle (deg) at vertex O between the points T1 and T2.
#
# NOTE (preserved from the original analysis code): the ABSOLUTE VALUE of the
# dot product is used, so the result is the unsigned acute angle in [0, 90]
# degrees. This is the convention that produced the published data; the
# 57.2958 (= 180/pi) factor is kept at the original precision so results are
# reproduced exactly.
proc rav_angle {T1 O T2} {
    set v1 [vecsub $O $T1]
    set v2 [vecsub $O $T2]

    set v1mag [expr {double([veclength $v1])}]
    set v2mag [expr {double([veclength $v2])}]

    set dotprod [expr {abs(double([vecdot $v1 $v2]))}]

    # Degenerate guard kept from the original code (effectively never fires).
    if {$v1mag == $v2mag} {
        return [expr {double(57.2958 * acos(1))}]
    }
    return [expr {double(57.2958 * acos($dotprod / ($v1mag * $v2mag)))}]
}

# rav_triangle_area com1 com2 com3 -> area of the triangle (Heron's formula).
proc rav_triangle_area {com1 com2 com3} {
    set a [veclength [vecsub $com2 $com1]]
    set b [veclength [vecsub $com3 $com2]]
    set c [veclength [vecsub $com1 $com3]]
    set s [expr {($a + $b + $c) / 2.0}]
    return [expr {sqrt($s * ($s - $a) * ($s - $b) * ($s - $c))}]
}


# -----------------------------------------------------------------------------
#  Output line formatters
#  Each receives: channel, time, list of values for the frame, values@frame0.
# -----------------------------------------------------------------------------

# Tilting: "<time>  <angle>  <angle - angle(frame0)>"
proc rav_emit_tilting {fh t values values0} {
    set a  [lindex $values 0]
    set a0 [lindex $values0 0]
    set d  [expr {$a - $a0}]
    puts $fh "[format {%.2f} $t]  [format {%12.5f} $a]  [format {%12.5f} $d]"
}

# Area / distance: "<time> \t <v1> <v2> ..."
proc rav_emit_columns {fh t values values0} {
    puts $fh "[format {%.2f} $t] \t [join $values { }]"
}


# -----------------------------------------------------------------------------
#  Per-trajectory helpers (shared by both drivers)
# -----------------------------------------------------------------------------

# Build selections once, then measure every frame. Returns a list of rows,
# each row being the value list returned by the measure callback.
#
# The make callback returns selection-TEXT strings; the atomselect handles are
# created here so they live in this call frame for the whole loop (VMD deletes
# atomselect handles when the proc that created them returns).
proc rav_measure_series {mol make measure segnames chains numframes} {
    set sels {}
    foreach text [$make $mol $segnames $chains] {
        lappend sels [atomselect $mol $text]
    }
    set rows {}
    for {set f 0} {$f < $numframes} {incr f} {
        lappend rows [$measure $sels $f]
    }
    foreach s $sels { catch {$s delete} }
    return $rows
}

# Extrema of column 0 across all rows. Returns {max fmax min fmin}.
proc rav_extrema {rows} {
    set col0 {}
    foreach r $rows { lappend col0 [lindex $r 0] }
    set mx [tcl::mathfunc::max {*}$col0]
    set mn [tcl::mathfunc::min {*}$col0]
    return [list $mx [lsearch -exact $col0 $mx] $mn [lsearch -exact $col0 $mn]]
}

# Write the per-trajectory time series.
proc rav_write_series {path rows dt stride emit} {
    set fh [open $path w]
    set v0 [lindex $rows 0]
    set n  [llength $rows]
    for {set f 0} {$f < $n} {incr f} {
        set t [expr {double($f * $dt * $stride)}]
        $emit $fh $t [lindex $rows $f] $v0
    }
    close $fh
}

# Write a PDB of the requested segments plus glycans within 2 A, at one frame.
proc rav_write_pdb {mol segnames frame path} {
    set sel [atomselect $mol \
        "segname $segnames or (glycan and same segname as within 2 of (segname $segnames))" \
        frame $frame]
    $sel writepdb $path
    $sel delete
}


# -----------------------------------------------------------------------------
#  Driver 1 of 2 : RAV (many spikes, one trajectory per spike)
# -----------------------------------------------------------------------------

# rav_run_spikes cfg
#
# Required cfg keys:
#   name      trajectory basename (files <inppath>/<name><NN>.alignedBackbone.{psf,dcd})
#   analysis  output sub-folder / label
#   inppath   directory with the per-spike psf/dcd files
#   outpath   directory under which <analysis>/ is created
#   catfile   name of the concatenated output file
#   make      selection callback {make mol segnames chains} -> sels
#   measure   measurement callback {measure sels frame} -> value list (col 0 = extrema)
#   emit      output formatter (rav_emit_tilting | rav_emit_columns)
# Optional keys (defaults): nspike 29, skip {6}, stride 1, dt 0.2
proc rav_run_spikes {cfg} {
    array set C {nspike 29 skip 6 stride 1 dt 0.2}
    array set C $cfg

    set outdir "$C(outpath)/$C(analysis)"
    file mkdir $outdir
    set maxmin "$outdir/maxmin.spikes.txt"
    if {[file exists $maxmin]} { file delete $maxmin }
    puts "INFO: \[$C(analysis)\] RAV spikes -> $outdir"

    set abs_max ""; set abs_min ""; set first 1

    for {set spike 0} {$spike <= $C(nspike)} {incr spike} {
        if {[lsearch -exact $C(skip) $spike] >= 0} { puts "INFO: skipping spike $spike"; continue }

        set id [format "%02d" $spike]
        set segnames "A1$id A2$id B1$id B2$id C1$id C2$id"
        set chains [list "A1$id A2$id" "B1$id B2$id" "C1$id C2$id"]
        puts "\n#### SPIKE $spike : segnames $segnames ####"

        cd $C(inppath)
        mol new     $C(name)$id.alignedBackbone.psf type psf \
            first 0 last -1 step 1          filebonds 1 autobonds 1 waitfor all
        mol addfile $C(name)$id.alignedBackbone.dcd type dcd \
            first 0 last -1 step $C(stride) filebonds 1 autobonds 1 waitfor all
        set mol top
        set numframes [molinfo $mol get numframes]
        puts "INFO: $numframes frames"

        set rows [rav_measure_series $mol $C(make) $C(measure) $segnames $chains $numframes]
        lassign [rav_extrema $rows] mx fmx mn fmn
        puts "INFO: spike $spike  MAX $mx (frame $fmx)  MIN $mn (frame $fmn)"

        rav_write_pdb $mol $segnames $fmx "$outdir/spike$spike.max.frame$fmx.pdb"
        rav_write_pdb $mol $segnames $fmn "$outdir/spike$spike.min.frame$fmn.pdb"

        if {$first} { set abs_max $mx; set abs_min $mn; set first 0 }
        if {$mx >= $abs_max} { set abs_max $mx; set abs_max_at "spike $spike frame $fmx"; rav_write_pdb $mol $segnames $fmx "$outdir/absolute.max.pdb" }
        if {$mn <= $abs_min} { set abs_min $mn; set abs_min_at "spike $spike frame $fmn"; rav_write_pdb $mol $segnames $fmn "$outdir/absolute.min.pdb" }

        set fh2 [open $maxmin a+]
        puts $fh2 "# SPIKE $spike"
        puts $fh2 "MAX: $mx frame $fmx "
        puts $fh2 "MIN: $mn frame $fmn \n"
        close $fh2

        rav_write_series "$outdir/spike$spike.$C(analysis).txt" $rows $C(dt) $C(stride) $C(emit)

        mol delete $mol
        catch {exec chmod g+rwX -R $outdir}
        puts "INFO: done with spike $spike"
    }

    set fh2 [open $maxmin a+]
    puts $fh2 "MAX OVERALL: $abs_max ($abs_max_at)  MIN OVERALL: $abs_min ($abs_min_at)"
    close $fh2

    rav_concat $outdir "spike*.$C(analysis).txt" "$outdir/$C(catfile)"
}


# -----------------------------------------------------------------------------
#  Driver 2 of 2 : single spike (several systems, several replicas each)
# -----------------------------------------------------------------------------

# rav_run_replicas cfg
#
# Required cfg keys:
#   names     list of system basenames (each has <inppath>/<name>.psf and
#             <inppath>/<name>_<rep>.dcd)
#   inppaths  list of input directories, parallel to names
#   nreps     list of replica counts, parallel to names
#   outpath   directory under which <analysis>/ is created
#   analysis  output sub-folder / label
#   catfile   name of the concatenated output file
#   make / measure / emit   as for rav_run_spikes
# Optional keys (defaults): stride 1, dt 0.1
#
# The single spike's six segments are AS1 AS2 BS1 BS2 CS1 CS2; for the
# per-chain observables (NTD area, RBD distance) the original analysis used
# only the first sub-segment of each chain (AS1 / BS1 / CS1), which is
# reproduced here. Absolute extrema and the OVERALL line are per system.
proc rav_run_replicas {cfg} {
    array set C {stride 1 dt 0.1}
    array set C $cfg

    set segnames "AS1 AS2 BS1 BS2 CS1 CS2"
    set chains   [list "AS1" "BS1" "CS1"]

    set outdir "$C(outpath)/$C(analysis)"
    file mkdir $outdir
    set maxmin "$outdir/maxmin.reps.txt"
    if {[file exists $maxmin]} { file delete $maxmin }
    puts "INFO: \[$C(analysis)\] single-spike replicas -> $outdir"

    foreach name $C(names) inppath $C(inppaths) nrep $C(nreps) {
        puts "\n######## SYSTEM $name ($nrep replicas) ########"
        set abs_max ""; set abs_min ""; set first 1

        for {set rep 1} {$rep <= $nrep} {incr rep} {
            puts "\n#### $name REPLICA $rep : segnames $segnames ####"

            cd $inppath
            mol new     $name.psf type psf \
                first 0 last -1 step 1          filebonds 1 autobonds 1 waitfor all
            mol addfile ${name}_$rep.dcd type dcd \
                first 0 last -1 step $C(stride) filebonds 1 autobonds 1 waitfor all
            set mol top
            set numframes [molinfo $mol get numframes]
            puts "INFO: $numframes frames"

            set rows [rav_measure_series $mol $C(make) $C(measure) $segnames $chains $numframes]
            lassign [rav_extrema $rows] mx fmx mn fmn
            puts "INFO: $name rep $rep  MAX $mx (frame $fmx)  MIN $mn (frame $fmn)"

            rav_write_pdb $mol $segnames $fmx "$outdir/rep$rep.max.frame$fmx.$name.pdb"
            rav_write_pdb $mol $segnames $fmn "$outdir/rep$rep.min.frame$fmn.$name.pdb"

            if {$first} { set abs_max $mx; set abs_min $mn; set first 0 }
            if {$mx >= $abs_max} { set abs_max $mx; set abs_max_at "rep $rep frame $fmx"; rav_write_pdb $mol $segnames $fmx "$outdir/absolute.max.$name.pdb" }
            if {$mn <= $abs_min} { set abs_min $mn; set abs_min_at "rep $rep frame $fmn"; rav_write_pdb $mol $segnames $fmn "$outdir/absolute.min.$name.pdb" }

            set fh2 [open $maxmin a+]
            puts $fh2 "# SYSTEM $name REPLICA $rep"
            puts $fh2 "MAX: $mx frame $fmx "
            puts $fh2 "MIN: $mn frame $fmn \n"
            close $fh2

            rav_write_series "$outdir/rep$rep.$C(analysis).$name.txt" $rows $C(dt) $C(stride) $C(emit)

            mol delete $mol
            catch {exec chmod g+rwX -R $outdir}
            puts "INFO: done with $name rep $rep"
        }

        set fh2 [open $maxmin a+]
        puts $fh2 "MAX OVERALL ($name): $abs_max ($abs_max_at)  MIN OVERALL ($name): $abs_min ($abs_min_at)"
        close $fh2
    }

    rav_concat $outdir "rep*.$C(analysis).*.txt" "$outdir/$C(catfile)"
}


# Concatenate matching per-trajectory files (natural order) into one file.
proc rav_concat {outdir pattern catpath} {
    if {[file exists $catpath]} { file delete $catpath }
    set out [open $catpath w]
    foreach f [lsort -dictionary [glob -nocomplain "$outdir/$pattern"]] {
        set in [open $f r]; fcopy $in $out; close $in
    }
    close $out
    puts "INFO: wrote $catpath"
}
