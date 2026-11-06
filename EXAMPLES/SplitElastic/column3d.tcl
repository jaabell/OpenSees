# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, Jose Antonio Abell. All rights reserved.
# A four-brick paired elastic column. Run serially, or with 2/4 OpenSeesMP ranks,
# in an external scratch directory. Compare column3d-*.txt between executions.
wipe
set rank [getPID]
set parts [getNP]
if {$parts != 1 && $parts != 2 && $parts != 4} {error "use 1, 2 or 4 ranks"}
set first [expr {$rank*4/$parts}]
set last [expr {($rank+1)*4/$parts}]
for {set i $first} {$i <= $last} {incr i} {
    if {$i == 0} {model BasicBuilder -ndm 3 -ndf 3} else {model BasicBuilder -ndm 3 -ndf 6}
    set x [expr {$i*.25}]
    set j 0
    foreach {y z} {0 0 1 0 1 1 0 1} {
        set n [expr {1+4*$i+$j}]
        node $n $x $y $z
        if {$i == 0} {
            fix $n 0 1 1
        } else {
            fix $n [expr {$i == 4}] 1 1 0 1 1
        }
        incr j
    }
}
# The element reads the actual node DOFs, including after switching builders.
model BasicBuilder -ndm 3 -ndf 3
for {set i $first} {$i < $last} {incr i} {
    set a [expr {1+4*$i}]
    set b [expr {$a+4}]
    set f0 [expr {1.-$i/4.}]
    set f1 [expr {1.-($i+1)/4.}]
    element SplitElastic [expr {$i+1}] $a $b [expr {$b+1}] [expr {$a+1}] \
        [expr {$a+3}] [expr {$b+3}] [expr {$b+2}] [expr {$a+2}] \
        200.0 0.2 2.0 -coupling 10.0 -couplingFactors $f0 $f1 $f1 $f0 $f0 $f1 $f1 $f0
}
rayleigh 0.1 0.0001 0 0
timeSeries Trig 1 0 0.1 0.1
pattern Plain 1 1 {
    if {$rank == 0} {
        foreach n {1 2 3 4} {load $n 0.25 0 0}
    }
}
for {set i [expr {$first+1}]} {$i <= $last} {incr i} {
    set n [expr {1+4*$i}]
    recorder Node -file column3d-$i.txt -precision 16 -time -dT .01 -node $n -dof 1 4 disp
}
constraints Transformation
if {$parts == 1} {
    numberer Plain
    system Diagonal
} else {
    numberer ParallelPlain
    system MPIDiagonal
}
algorithm Linear
integrator ExplicitDifference
analysis Transient
if {[analyze 1000 .001] != 0} {error "3D analysis failed"}
remove recorders
wipe
puts "COLUMN3D_OK $rank"
exit
