# MODEL_STAGE boundaries and per-stage KIND of recorder ladruno, classic Tcl.
#
# Usage:  OpenSees ladruno_stage_kind_model.tcl <out.ladruno> <case>
# Prints "STAGE_KIND_MODEL OK <case>" as its last line on success (the
# interpreter exits 0 even when a deck fails, so the caller keys on the marker).
#
# One 2-D cantilever (elastic beam, tip mass) and `recorder ladruno` with NO
# -kind unless the case says so. Every analysis records every step (-T dt 0).
# The cases and the MODEL_STAGE list (KIND, recorded rows) expected of them:
#
#   dc_static_transient  static 4 | loadConst, wipeAnalysis, domainChange | transient 5
#                        -> [static 4, transient 5]   (user barrier + kind flip)
#   flip_static_transient  the same with NO domainChange
#                        -> [static 4, transient 5]   (kind flip alone)
#   dc_static_static     static 4 | loadConst, new load pattern, domainChange | static 3
#                        -> [static 4, static 3]      (user barrier alone)
#   repeat_static        static 4 | wipeAnalysis, a new Static analysis | static 3
#                        -> [static 7]                (a repeated kind is no boundary)
#   internal_restamp     static 2 | SP pattern | static 2 | eleLoad pattern | static 2
#                        -> [static 6]                (internal re-stamps are folded)
#   dc_first             domainChange before the first analyze | static 4
#                        -> [static 4]                (no stage for the first record)
#   transient_only       transient 5 -> [transient 5]
#   record_before_analysis  `record` with no analysis | transient 5
#                        -> [transient 6]             (provisional KIND rewritten)
#   override_kind        -kind static; static 4 | wipeAnalysis | transient 5 (no dc)
#                        -> [static 4, static 5]      (the flip still splits; -kind
#                                                      only sets the label)
#   override_kind_dc     -kind transient; static 4 | new pattern, domainChange | static 3
#                        -> [transient 4, transient 3] (the barrier splits too)
#   eigen_snapshot       eigen 1; record | gravity, analysis Static, analyze 4
#                        -> [static 5]   (the transient analysis classic `eigen`
#                                         builds is not a running analysis)
#   eigen_snapshot_transient  eigen 1; record | analysis Transient, analyze 5
#                        -> [transient 6]
#   flip_then_topology   static 4 | wipeAnalysis | transient 5 | add a node + element | transient 2
#                        -> [static 4, transient 5, transient 2], ids strictly increasing
#                           (a flip that did not move the stamp must not collide with
#                            the stamp of the later topology change)
#   freefield_shape      a free-field deck shape sampled on -T dt 0.02: gravity static 4 |
#                        eigen with the static analysis alive | loadConst, wipeAnalysis |
#                        a new node + element, domainChange, rayleigh, UniformExcitation |
#                        transient 10 x 0.01
#                        -> [static 4, transient 5]  (5 = the -T dt 0.02 samples of 10 steps)

if {$argc < 2} { puts "usage: ladruno_stage_kind_model.tcl <out.ladruno> <case>"; exit 1 }
set out [lindex $argv 0]
set case [lindex $argv 1]
if {[file exists $out]} { file delete $out }

wipe
model basic -ndm 2 -ndf 3
node 1 0.0 0.0
node 2 0.0 3.0
fix 1 1 1 1
mass 2 10.0 10.0 0.0
geomTransf Linear 1
element elasticBeamColumn 1 1 2 0.1 2.0e8 1.0e-3 1

set kind_opt {}
if {$case == "override_kind"} { set kind_opt {-kind static} }
if {$case == "override_kind_dc"} { set kind_opt {-kind transient} }
# the free-field case samples on a -T dt grid; every other case records every step.
set t_opt {-T dt 0.0}
if {$case == "freefield_shape"} { set t_opt {-T dt 0.02} }
recorder ladruno $out -N displacement {*}$kind_opt {*}$t_opt

proc static_chain {} {
    constraints Plain
    numberer RCM
    system BandGen
    test NormDispIncr 1.0e-10 10
    algorithm Newton
    integrator LoadControl 0.25
    analysis Static
}
proc transient_chain {} {
    constraints Plain
    numberer RCM
    system BandGen
    test NormDispIncr 1.0e-10 10
    algorithm Newton
    integrator Newmark 0.5 0.25
    analysis Transient
}
proc gravity {} {
    timeSeries Linear 1
    pattern Plain 1 1 { load 2 0.0 -100.0 0.0 }
}
proc lateral {} {
    timeSeries Sine 2 0.0 10.0 0.5
    pattern Plain 2 2 { load 2 100.0 0.0 0.0 }
}

switch -- $case {
    dc_static_transient {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        wipeAnalysis
        domainChange
        lateral; transient_chain; analyze 5 0.01
    }
    flip_static_transient {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        wipeAnalysis
        lateral; transient_chain; analyze 5 0.01
    }
    dc_static_static {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        timeSeries Linear 3
        pattern Plain 3 3 { load 2 50.0 0.0 0.0 }
        domainChange
        analyze 3
    }
    repeat_static {
        gravity; static_chain; analyze 4
        wipeAnalysis
        static_chain; analyze 3
    }
    internal_restamp {
        gravity; static_chain; analyze 2
        # an SP pattern moves the domain-change stamp (Domain::addLoadPattern)
        pattern Plain 5 1 { sp 2 3 0.0 }
        analyze 2
        # an element load (eleLoad) on the beam
        pattern Plain 6 1 { eleLoad -ele 1 -type -beamUniform -1.0 }
        analyze 2
    }
    dc_first {
        domainChange
        gravity; static_chain; analyze 4
    }
    transient_only {
        lateral; transient_chain; analyze 5 0.01
    }
    record_before_analysis {
        record
        lateral; transient_chain; analyze 5 0.01
    }
    override_kind {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        wipeAnalysis
        lateral; transient_chain; analyze 5 0.01
    }
    override_kind_dc {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        timeSeries Linear 3
        pattern Plain 3 3 { load 2 50.0 0.0 0.0 }
        domainChange
        analyze 3
    }
    eigen_snapshot {
        eigen 1
        record
        gravity; static_chain; analyze 4
    }
    eigen_snapshot_transient {
        eigen 1
        record
        lateral; transient_chain; analyze 5 0.01
    }
    flip_then_topology {
        gravity; static_chain; analyze 4
        loadConst -time 0.0
        wipeAnalysis
        lateral; transient_chain; analyze 5 0.01
        node 3 1.0 3.0
        fix 3 1 1 1
        element elasticBeamColumn 2 2 3 0.1 2.0e8 1.0e-3 1
        analyze 2 0.01
    }
    freefield_shape {
        # a free-field deck, reduced. The recorder was created before both
        # stages and samples on -T dt 0.02 (set above). Gravity runs Static;
        # `eigen` feeds rayleigh while the static analysis is alive (no
        # eigen-only mark); the stage close is loadConst + wipeAnalysis; the
        # earthquake stage adds a node and an element (stamp AND fingerprint
        # move), fires the domainChange barrier, binds rayleigh, drives a
        # UniformExcitation and runs Transient.
        gravity; static_chain; analyze 4
        set lambda [eigen -fullGenLapack 1]
        loadConst -time 0.0
        wipeAnalysis
        node 3 1.0 3.0
        mass 3 5.0 5.0 0.0
        element elasticBeamColumn 2 2 3 0.1 2.0e8 1.0e-3 1
        domainChange
        set w [expr {sqrt([lindex $lambda 0])}]
        rayleigh [expr {0.05 * $w}] 0.0 [expr {0.05 / $w}] 0.0
        timeSeries Path 4 -dt 0.01 -values {0.0 0.5 1.0 0.5 0.0 -0.5 -1.0 -0.5 0.0 0.5 1.0 0.5}
        pattern UniformExcitation 4 1 -accel 4
        transient_chain; analyze 10 0.01
    }
    default {
        puts "unknown case $case"
        exit 1
    }
}
remove recorders
puts "STAGE_KIND_MODEL OK $case"
