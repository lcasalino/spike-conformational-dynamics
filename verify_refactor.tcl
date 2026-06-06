# Developed and written by Lorenzo Casalino (UC San Diego).
# Verification harness (plain tclsh, no VMD). Run: tclsh verify_refactor.tcl
# Stubs the VMD commands, then exercises the refactored library end-to-end.

set DIR [file dirname [file normalize [info script]]]

# ---- pure-Tcl VMD vector ops ----
proc vecsub {a b} { set r {}; foreach x $a y $b { lappend r [expr {$x-$y}] }; return $r }
proc vecdot {a b} { set s 0.0; foreach x $a y $b { set s [expr {$s+$x*$y}] }; return $s }
proc veclength {a} { set s 0.0; foreach x $a { set s [expr {$s+$x*$x}] }; return [expr {sqrt($s)}] }
proc vecdist {a b} { return [veclength [vecsub $a $b]] }

# ---- selection-recording atomselect (for Test B) and object stub (for Test C) ----
set ::SELTEXTS {}
set ::SELCNT 0
proc atomselect {mol text args} {
    lappend ::SELTEXTS $text
    set name "::sel[incr ::SELCNT]"
    proc $name {sub args} {
        switch -- $sub {
            writepdb { set p [lindex $args 0]; set f [open $p w]; puts $f "REMARK stub"; close $f }
            frame    { }
            delete   { }
            default  { }
        }
    }
    return $name
}

set FAIL 0
proc check {label cond} {
    global FAIL
    if {$cond} { puts "  ok   $label" } else { puts "  FAIL $label"; incr FAIL }
}
proc norm {s} { regsub -all {\s+} [string trim $s] " " s; return $s }

# =====================================================================
# Load the library (Test A: syntax)
# =====================================================================
puts "== Test A: load library =="
source [file join $DIR rav_common.tcl]
source [file join $DIR rav_selections.tcl]
puts "  ok   rav_common.tcl + rav_selections.tcl loaded"

# =====================================================================
# Test B: selection strings, token-identical to originals (RAV and SS)
# =====================================================================
puts "== Test B: selection strings =="
proc collect {makeproc segnames chains} {
    # _make now returns selection-text strings directly
    set out {}; foreach t [$makeproc mol $segnames $chains] { lappend out [norm $t] }
    return $out
}
# RAV spike 0 naming and SS naming
set RAV_SEG "A100 A200 B100 B200 C100 C200"
set RAV_CH  [list "A100 A200" "B100 B200" "C100 C200"]
set SS_SEG  "AS1 AS2 BS1 BS2 CS1 CS2"
set SS_CH   [list "AS1" "BS1" "CS1"]

# expected (normalized) original strings, parameterized by segnames/chain list
proc exp_tilt {seg o a1 a2} {
    return [list \
      "protein and segname $seg and resid $o and name CA" \
      "protein and segname $seg and resid $a1 and name CA" \
      "protein and segname $seg and resid $a2 and name CA"]
}
proc exp_ntd {ch} {
    lassign $ch a b c
    return [list \
      "name CA and segname $a and (resid 13 to 291)" \
      "name CA and segname $b and (resid 13 to 291)" \
      "name CA and segname $c and (resid 13 to 291)"]
}
proc exp_rbd {seg ch} {
    lassign $ch a b c
    set r "((resid 375 to 380) or (resid 394 to 404) or (resid 431 to 438) or (resid 508 to 517))"
    set h "((resid 747 to 784) or (resid 946 to 967) or (resid 986 to 1034))"
    return [list \
      "name CA and segname $a and $r" \
      "name CA and segname $b and $r" \
      "name CA and segname $c and $r" \
      "name CA and segname $seg and $h"]
}
foreach {obj o a1 a2} {ankle 1213 {1169 to 1206} {1212 to 1239} hip {1136 to 1140} {1141 to 1161} {816 to 1135} knee {1161 to 1168} {1169 to 1206} {1141 to 1161}} {
    check "RAV ${obj}" [string equal [collect ${obj}_make $RAV_SEG $RAV_CH] [exp_tilt $RAV_SEG $o $a1 $a2]]
    check "SS  ${obj}" [string equal [collect ${obj}_make $SS_SEG  $SS_CH ] [exp_tilt $SS_SEG  $o $a1 $a2]]
}
check "RAV ntd" [string equal [collect ntd_make $RAV_SEG $RAV_CH] [exp_ntd $RAV_CH]]
check "SS  ntd" [string equal [collect ntd_make $SS_SEG  $SS_CH ] [exp_ntd $SS_CH ]]
check "RAV rbd" [string equal [collect rbd_make $RAV_SEG $RAV_CH] [exp_rbd $RAV_SEG $RAV_CH]]
check "SS  rbd" [string equal [collect rbd_make $SS_SEG  $SS_CH ] [exp_rbd $SS_SEG  $SS_CH ]]

# =====================================================================
# Test C: mock end-to-end for both drivers
# =====================================================================
puts "== Test C: driver end-to-end (mock VMD) =="
set ::FRAMES 5
proc mol {sub args} { if {$sub eq "new"} { return 0 } ; return 0 }
proc molinfo {m get what} { return $::FRAMES }
proc measure {args} { return {0 0 0} }
# deterministic per-frame values: max at frame 2, min at frame 1
set ::VALS {2.00 0.00 4.00 1.00 3.00}
proc test_make {mol segnames chains} { return [list "stub"] }
proc test_measure {sels frame} { return [format "%.2f" [lindex $::VALS $frame]] }

set TMP [file join /tmp rav_verify_[pid]]
file mkdir [file join $TMP in]
file mkdir [file join $TMP in_open] [file join $TMP in_closed] [file join $TMP in_mut]
set OUT [file join $TMP out]

# expected per-trajectory time series for dt and the VALS above
proc expected_series {dt} {
    set v0 [lindex $::VALS 0]
    set s ""
    for {set f 0} {$f < 5} {incr f} {
        set t [expr {double($f*$dt)}]
        set a [lindex $::VALS $f]
        append s "[format {%.2f} $t]  [format {%12.5f} $a]  [format {%12.5f} [expr {$a-$v0}]]\n"
    }
    return $s
}

# ---- spikes driver ----
rav_run_spikes [list name traj analysis A_SPK inppath [file join $TMP in] outpath $OUT \
    nspike 3 skip {2} stride 1 dt 0.2 catfile cat_spk.txt \
    make test_make measure test_measure emit rav_emit_tilting]
set sd [file join $OUT A_SPK]
check "spike txt exists"        [file exists [file join $sd spike0.A_SPK.txt]]
check "spike skipped 2"         [expr {![file exists [file join $sd spike2.A_SPK.txt]]}]
check "spike max pdb frame2"    [file exists [file join $sd spike0.max.frame2.pdb]]
check "spike min pdb frame1"    [file exists [file join $sd spike0.min.frame1.pdb]]
check "spike min NOT frame2 (bug fixed)" [expr {![file exists [file join $sd spike0.min.frame2.pdb]]}]
check "spike absolute.max.pdb fixed name"  [file exists [file join $sd absolute.max.pdb]]
check "spike absolute.min.pdb fixed name"  [file exists [file join $sd absolute.min.pdb]]
check "spike NO frame-named absolute pdbs" [expr {[llength [glob -nocomplain [file join $sd absolute.*.frame*.pdb]]] == 0}]
check "spike series content"    [string equal [read [set h [open [file join $sd spike0.A_SPK.txt]]]][close $h] [expected_series 0.2]]
check "spike concat lines=15"   [expr {[regexp -all "\n" [read [set h [open [file join $sd cat_spk.txt]]]][close $h]] == 15}]
check "spike maxmin OVERALL"    [expr {[string match "*MAX OVERALL: 4.00*MIN OVERALL: 0.00*" [read [set h [open [file join $sd maxmin.spikes.txt]]]][close $h]]}]

# ---- replicas driver ----
rav_run_replicas [list names {sysO sysC} inppaths [list [file join $TMP in_open] [file join $TMP in_closed]] \
    nreps {2 1} outpath $OUT analysis A_REP dt 0.1 catfile cat_rep.txt \
    make test_make measure test_measure emit rav_emit_tilting]
set rd [file join $OUT A_REP]
check "rep txt naming incl system" [file exists [file join $rd rep1.A_REP.sysO.txt]]
check "rep2 only for sysO"         [expr {[file exists [file join $rd rep2.A_REP.sysO.txt]] && ![file exists [file join $rd rep2.A_REP.sysC.txt]]}]
check "rep max pdb frame2"         [file exists [file join $rd rep1.max.frame2.sysO.pdb]]
check "rep min pdb frame1"         [file exists [file join $rd rep1.min.frame1.sysO.pdb]]
check "rep absolute.max.<sys>.pdb" [file exists [file join $rd absolute.max.sysO.pdb]]
check "rep NO frame-named absolute pdbs" [expr {[llength [glob -nocomplain [file join $rd absolute.*.frame*.pdb]]] == 0}]
check "rep series content (dt0.1)" [string equal [read [set h [open [file join $rd rep1.A_REP.sysO.txt]]]][close $h] [expected_series 0.1]]
check "rep concat lines=15"        [expr {[regexp -all "\n" [read [set h [open [file join $rd cat_rep.txt]]]][close $h]] == 15}]
check "rep maxmin per-system"      [expr {[string match "*MAX OVERALL (sysO):*MAX OVERALL (sysC):*" [read [set h [open [file join $rd maxmin.reps.txt]]]][close $h]]}]

file delete -force $TMP
puts "================================"
puts [expr {$FAIL==0 ? "ALL CHECKS PASSED" : "$FAIL CHECK(S) FAILED"}]
