"""DOF numbering with MP constraints, serial and in OpenSeesMP.

The serial cases use openseespy with constraints Plain, where the numberer
gives every MP-constrained DOF (marked -4) the equation number of its
retained DOF.

The MPI cases run OpenSeesMP through mpiexec and are skipped unless the
environment variable OPENSEESMP_EXE points to an OpenSeesMP executable and
mpiexec is found (on PATH, or through the environment variable MPIEXEC).
"""

import os
import shutil
import signal
import subprocess
from pathlib import Path

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


# ---------------------------------------------------------------------------
# serial: equalDOF-heavy model under constraints Plain
# ---------------------------------------------------------------------------

NUM_MASTERS = 40
NUM_SLAVES = 400
OVERLAP_NODE = 3001


def _equaldof_model(numberer):
    """Masters on springs to ground, many slaves tied to them by equalDOF,
    and one node tied to two masters with an overlapping DOF."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ele = 1
    for i in range(1, NUM_MASTERS + 1):
        ops.node(i, float(i), 0.0)
        ops.node(1000 + i, float(i), -1.0)
        ops.fix(1000 + i, 1, 1)
        ops.element("zeroLength", ele, 1000 + i, i, "-mat", 1, 1, "-dir", 1, 2)
        ele += 1
    for k in range(NUM_SLAVES):
        tag = 2001 + k
        master = 1 + k % NUM_MASTERS
        ops.node(tag, float(master), 1.0)
        ops.equalDOF(master, tag, 1, 2)
        if k > 0:
            ops.element("zeroLength", ele, tag - 1, tag, "-mat", 1, 1, "-dir", 1, 2)
            ele += 1
    # DOF 1 is first tied to master 1, then again to master 2; the DOF ends
    # up with the equation of the constraint defined last
    ops.node(OVERLAP_NODE, 0.0, 2.0)
    ops.equalDOF(1, OVERLAP_NODE, 1, 2)
    ops.equalDOF(2, OVERLAP_NODE, 1)
    ops.element("zeroLength", ele, 3, OVERLAP_NODE, "-mat", 1, 1, "-dir", 1, 2)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(NUM_MASTERS, 1.0, 0.5)
    ops.constraints("Plain")
    ops.numberer(numberer)
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


@pytest.mark.parametrize("numberer", ["Plain", "RCM"])
def test_equaldof_numbering(numberer):
    try:
        _equaldof_model(numberer)
        dofs = {tag: list(ops.nodeDOFs(tag)) for tag in ops.getNodeTags()}

        # free DOFs of the unconstrained nodes: a bijection onto 0..neq-1
        free = [eq for tag, ids in dofs.items()
                if tag <= NUM_MASTERS for eq in ids]
        assert sorted(free) == list(range(len(free)))
        assert ops.systemSize() == len(free)

        # every slave DOF carries the equation of its retained DOF
        for k in range(NUM_SLAVES):
            master = 1 + k % NUM_MASTERS
            assert dofs[2001 + k] == dofs[master]

        # several MP_Constraints on one node are applied in order
        assert dofs[OVERLAP_NODE] == [dofs[2][0], dofs[1][1]]
    finally:
        ops.wipe()


# ---------------------------------------------------------------------------
# MPI: ParallelNumberer in OpenSeesMP
# ---------------------------------------------------------------------------

MP_EXE = os.environ.get("OPENSEESMP_EXE")
MPIEXEC = os.environ.get("MPIEXEC") or shutil.which("mpiexec")
needs_mp = pytest.mark.skipif(
    not MP_EXE or not os.path.isfile(MP_EXE) or not MPIEXEC,
    reason="set OPENSEESMP_EXE (and have mpiexec on PATH or in MPIEXEC)")

TIMEOUT = 120


def _run_mp(workdir, deck, nproc, env_extra):
    """Run deck on nproc processes. Returns (completed, output). A run that
    does not finish within TIMEOUT is killed with its children, since a
    process stopping in the middle of a collective step can leave the other
    processes waiting."""
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "deck.tcl").write_text(deck, encoding="ascii")
    env = dict(os.environ)
    env.update(env_extra)
    popen_kw = {}
    if os.name != "nt":
        popen_kw["start_new_session"] = True
    proc = subprocess.Popen([MPIEXEC, "-n", str(nproc), MP_EXE, "deck.tcl"],
                            cwd=workdir, env=env, stdin=subprocess.DEVNULL,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, **popen_kw)
    try:
        out, _ = proc.communicate(timeout=TIMEOUT)
        return True, out
    except subprocess.TimeoutExpired:
        if os.name == "nt":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                           stdin=subprocess.DEVNULL, capture_output=True)
        else:
            os.killpg(proc.pid, signal.SIGKILL)
        out, _ = proc.communicate()
        return False, out


# A truss chain along x split over the processes: process p owns elements
# p*EPR+1 .. (p+1)*EPR, the end nodes of each segment are shared. Node 1 is
# fixed on process 0 and the last process adds a grounded spring at the free
# end, so that process holds a constraint of its own. Exact end displacement
# 1/(k_chain + k_spring).
LAGRANGE_DECK = r"""
set pid [getPID]
set np  [getNP]
set EPR 5
model BasicBuilder -ndm 1 -ndf 1
uniaxialMaterial Elastic 1 1000.0
set n0 [expr $pid*$EPR]
set n1 [expr ($pid+1)*$EPR]
for {set i $n0} {$i <= $n1} {incr i} {
    node [expr $i+1] [expr double($i)]
}
if {$pid == 0} { fix 1 1 }
for {set e [expr $n0+1]} {$e <= $n1} {incr e} {
    element truss $e $e [expr $e+1] 1.0 1
}
if {$pid == $np-1} {
    node 1000 [expr double($n1+1)]
    fix 1000 1
    element truss 1000 [expr $n1+1] 1000 1.0 1
}
pattern Plain 1 Linear {
    if {$pid == $np-1} { load [expr $n1+1] 1.0 }
}
constraints $::env(TEST_HANDLER)
numberer ParallelRCM
system Mumps
test NormDispIncr 1e-10 10
algorithm Linear
integrator LoadControl 1.0
analysis Static
set ok [analyze 1]
if {$pid == $np-1} {
    puts "RESULT $ok [format %.15g [nodeDisp [expr $n1+1] 1]]"
}
"""


def _result(out):
    for line in out.splitlines():
        if line.startswith("RESULT"):
            _, ok, disp = line.split()
            return int(ok), float(disp)
    return None


@needs_mp
def test_lagrange_multipliers_not_fused(tmp_path):
    """Under constraints Lagrange every constraint has a node-less DOF_Group.
    The merge used to fuse them across processes and the analysis returned a
    wrong displacement without any warning. It must either give the exact
    answer or stop with an error naming the cause."""
    exact = 1.0 / (1000.0 / 10 + 1000.0)  # 10 elements in series + spring

    done, out = _run_mp(tmp_path / "transformation", LAGRANGE_DECK, 2,
                        {"TEST_HANDLER": "Transformation"})
    assert done, out[-2000:]
    ok, disp = _result(out)
    assert ok == 0 and disp == pytest.approx(exact, rel=1e-10), out[-2000:]

    done, out = _run_mp(tmp_path / "lagrange", LAGRANGE_DECK, 2,
                        {"TEST_HANDLER": "Lagrange"})
    res = _result(out)
    if res is not None and res[0] == 0:
        assert res[1] == pytest.approx(exact, rel=1e-10), out[-2000:]
    else:
        assert "has no node" in out, out[-2000:]


# A 2-D truss lattice of 3 rows split by columns over the processes, fixed
# at the left edge except for the vertical DOF of node 1, whose DOF_Group has
# tag 0. The last process adds a node tied by equalDOF to the top
# right node (a -4 DOF under constraints Plain). system Mumps keeps the
# global equation numbers the numberer assigns, which every process writes
# for its own nodes.
LATTICE_DECK = r"""
set pid [getPID]
set np  [getNP]
set C 3
model BasicBuilder -ndm 2 -ndf 2
uniaxialMaterial Elastic 1 1000.0
set c0 [expr $pid*$C]
set c1 [expr ($pid+1)*$C]
proc nt {c r} { return [expr $c*3 + $r + 1] }
for {set c $c0} {$c <= $c1} {incr c} {
    for {set r 0} {$r < 3} {incr r} {
        node [nt $c $r] [expr double($c)] [expr double($r)]
        if {$c == 0 && $r == 0} { fix [nt $c $r] 1 0 }
        if {$c == 0 && $r > 0} { fix [nt $c $r] 1 1 }
    }
}
set e [expr $pid*100 + 1]
for {set c $c0} {$c <= $c1} {incr c} {
    for {set r 0} {$r < 2} {incr r} {
        if {$c > $c0 || $pid == 0} {
            element truss $e [nt $c $r] [nt $c [expr $r+1]] 1.0 1; incr e
        }
    }
}
for {set c $c0} {$c < $c1} {incr c} {
    for {set r 0} {$r < 3} {incr r} {
        element truss $e [nt $c $r] [nt [expr $c+1] $r] 1.0 1; incr e
        if {$r < 2} {
            element truss $e [nt $c $r] [nt [expr $c+1] [expr $r+1]] 1.0 1; incr e
        }
    }
}
set top [nt $c1 2]
if {$pid == $np-1} {
    node 9001 [expr double($c1)+0.5] 2.0
    equalDOF $top 9001 1 2
    element truss $e 9001 [nt $c1 0] 1.0 1
}
pattern Plain 1 Linear {
    if {$pid == $np-1} { load $top 1.0 -0.5 }
}
constraints Plain
numberer $::env(TEST_NUMBERER)
system Mumps
test NormDispIncr 1e-10 10
algorithm Linear
integrator LoadControl 1.0
analysis Static
set ok [analyze 1]
set fd [open dofs_$pid.txt w]
foreach n [lsort -integer [getNodeTags]] {
    puts $fd "$n [nodeDOFs $n]"
}
close $fd
puts "RESULT $ok"
"""

# Equation numbers written by the code before the linear-time merge,
# {nodeTag: [eq, eq]}. The merge must reproduce them exactly. With
# ParallelPlain the vertex with tag 0 (node 1) is numbered last.
EXPECTED = {
    ('ParallelRCM', 2): {
        1: [-1, 36], 2: [-1, -1], 3: [-1, -1], 4: [34, 35], 5: [32, 33],
        6: [30, 31], 7: [28, 29], 8: [26, 27], 9: [24, 25], 10: [22, 23],
        11: [20, 21], 12: [18, 19], 13: [16, 17], 14: [14, 15], 15: [12, 13],
        16: [10, 11], 17: [8, 9], 18: [6, 7], 19: [4, 5], 20: [2, 3],
        21: [0, 1], 9001: [0, 1],
    },
    ('ParallelRCM', 3): {
        1: [-1, 54], 2: [-1, -1], 3: [-1, -1], 4: [52, 53], 5: [50, 51],
        6: [48, 49], 7: [46, 47], 8: [44, 45], 9: [42, 43], 10: [40, 41],
        11: [38, 39], 12: [36, 37], 13: [34, 35], 14: [32, 33], 15: [30, 31],
        16: [28, 29], 17: [26, 27], 18: [24, 25], 19: [22, 23], 20: [20, 21],
        21: [18, 19], 22: [16, 17], 23: [14, 15], 24: [12, 13], 25: [10, 11],
        26: [8, 9], 27: [6, 7], 28: [4, 5], 29: [2, 3], 30: [0, 1],
        9001: [0, 1],
    },
    ('ParallelPlain', 2): {
        1: [-1, 36], 2: [-1, -1], 3: [-1, -1], 4: [24, 25], 5: [26, 27],
        6: [28, 29], 7: [30, 31], 8: [32, 33], 9: [34, 35], 10: [0, 1],
        11: [2, 3], 12: [4, 5], 13: [6, 7], 14: [8, 9], 15: [10, 11],
        16: [12, 13], 17: [14, 15], 18: [16, 17], 19: [18, 19], 20: [20, 21],
        21: [22, 23], 9001: [22, 23],
    },
    ('ParallelPlain', 3): {
        1: [-1, 54], 2: [-1, -1], 3: [-1, -1], 4: [42, 43], 5: [44, 45],
        6: [46, 47], 7: [48, 49], 8: [50, 51], 9: [52, 53], 10: [0, 1],
        11: [2, 3], 12: [4, 5], 13: [6, 7], 14: [8, 9], 15: [10, 11],
        16: [12, 13], 17: [14, 15], 18: [16, 17], 19: [18, 19], 20: [20, 21],
        21: [22, 23], 22: [24, 25], 23: [26, 27], 24: [28, 29], 25: [30, 31],
        26: [32, 33], 27: [34, 35], 28: [36, 37], 29: [38, 39], 30: [40, 41],
        9001: [40, 41],
    },
}


def _read_dofs(path):
    dofs = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if parts:
            dofs[int(parts[0])] = [int(x) for x in parts[1:]]
    return dofs


@needs_mp
@pytest.mark.parametrize("nproc", [2, 3])
@pytest.mark.parametrize("numberer", ["ParallelRCM", "ParallelPlain"])
def test_parallel_numbering_unchanged(tmp_path, numberer, nproc):
    done, out = _run_mp(tmp_path, LATTICE_DECK, nproc,
                        {"TEST_NUMBERER": numberer})
    assert done, out[-2000:]
    ranks = {}
    for p in range(nproc):
        path = tmp_path / ("dofs_%d.txt" % p)
        assert path.exists(), out[-2000:]
        ranks[p] = _read_dofs(path)

    # shared nodes carry the same numbers on every process, and the numbers
    # of all unconstrained nodes form a bijection onto 0..neq-1
    merged = {}
    for dofs in ranks.values():
        for tag, ids in dofs.items():
            if tag in merged:
                assert merged[tag] == ids, (tag, merged[tag], ids)
            merged[tag] = ids
    free = sorted(eq for tag, ids in merged.items() if tag != 9001
                  for eq in ids if eq >= 0)
    assert free == list(range(len(free)))
    top = 3 * (3 * nproc) + 3
    assert merged[9001] == merged[top]

    assert merged == EXPECTED[(numberer, nproc)]
