"""LoadControl -tangentPredictor: first-iteration forcing -K*du_p for
non-homogeneous sp constraints under constraints Transformation.

Model: a column of 20 stdBrick elements (10 x 10 x 100, nu = 0), base fixed,
lateral dofs fixed, top face driven axially by sp to delta = 0.15. With
J2Plasticity (fy = 379.5, H = 2000) the converged stress is 300 < fy, so
nothing yields physically and the exact midpoint displacement is delta / 2.
The standard predictor drives only the top layer in the first iteration,
which yields spuriously and costs iterations and cutbacks. The tangent
predictor must remove that cost without changing the converged answer.

The midpoint is probed, not the driven face: the driven face reports the
prescribed value whether or not the interior moved.
"""

import os
import subprocess
from pathlib import Path

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


L, AX, E, FY, H, DELTA, N = 100.0, 10.0, 200000.0, 379.5, 2000.0, 0.15, 20
MID_NODE = 4 * (N // 2) + 1

STOCK = ()
TANGENT = ("-tangentPredictor",)


def _build(mat, drive="disp", handler=("Transformation",), max_iter=40):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    h = L / N
    for k in range(N + 1):
        z, b = k * h, 4 * k
        ops.node(b + 1, 0.0, 0.0, z)
        ops.node(b + 2, AX, 0.0, z)
        ops.node(b + 3, AX, AX, z)
        ops.node(b + 4, 0.0, AX, z)
    if mat == "j2":
        ops.nDMaterial("J2Plasticity", 1, E / 3.0, E / 2.0, FY, FY, 0.0, H)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, 0.0)
    for k in range(N):
        b = 4 * k
        ops.element("stdBrick", k + 1, *[b + j for j in range(1, 9)], 1)
    for k in range(N + 1):
        for j in range(1, 5):
            if k == 0:
                ops.fix(4 * k + j, 1, 1, 1)
            else:
                ops.fix(4 * k + j, 1, 1, 0)
    top = 4 * N
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(1, 5):
        if drive == "disp":
            ops.sp(top + j, 3, DELTA)
        else:
            ops.load(top + j, 0.0, 0.0, E * AX * AX * DELTA / L / 4.0)
    ops.constraints(*handler)
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-10, max_iter, 0)
    ops.algorithm("KrylovNewton")
    ops.analysis("Static")


def _fixed_march(mat, option, drive="disp", handler=("Transformation",),
                 steps=10):
    """Equal steps, integrator re-issued every step. Returns (iters, u_mid)."""
    _build(mat, drive, handler, max_iter=40)
    iters = 0
    for _ in range(steps):
        ops.integrator("LoadControl", 1.0 / steps, *option)
        assert ops.analyze(1) == 0
        iters += ops.testIter()
    return iters, ops.nodeDisp(MID_NODE, 3)


def _adaptive_march(mat, option, max_iter=5):
    """Halve on failure, double after two successes, cap 0.2.

    Returns (increments, cutbacks, iterations, u_mid).
    """
    _build(mat, max_iter=max_iter)
    lam, dl, inc, cuts, iters, good = 0.0, 0.1, 0, 0, 0, 0
    while lam < 1.0 - 1.0e-12:
        dl = min(dl, 1.0 - lam)
        ops.integrator("LoadControl", dl, *option)
        ok = ops.analyze(1)
        iters += ops.testIter()
        if ok == 0:
            lam += dl
            inc += 1
            good += 1
            if good >= 2:
                dl, good = min(2.0 * dl, 0.2), 0
        else:
            cuts += 1
            good = 0
            dl /= 2.0
            assert dl >= 1.0e-4, "march stalled"
    return inc, cuts, iters, ops.nodeDisp(MID_NODE, 3)


def teardown_function(_):
    ops.wipe()


def test_option_absent_is_unchanged():
    """Without the option LoadControl reproduces the reference path exactly.

    Reference values recorded with the integrator before the option existed.
    """
    inc, cuts, iters, u_mid = _adaptive_march("j2", STOCK)
    assert (inc, cuts, iters) == (43, 23, 224)
    assert u_mid.hex() == ADAPTIVE_U_MID_HEX

    iters, u_mid = _fixed_march("j2", STOCK)
    assert iters == 60
    assert u_mid.hex() == FIXED_U_MID_HEX


ADAPTIVE_U_MID_HEX = "0x1.3333333333331p-4"
FIXED_U_MID_HEX = "0x1.3333333333333p-4"


def test_tangent_predictor_removes_cutbacks():
    """Adaptive march, maxIter 5: the option must match the elastic control."""
    stock = _adaptive_march("j2", STOCK)
    elastic = _adaptive_march("elastic", STOCK)
    tangent = _adaptive_march("j2", TANGENT)

    assert stock[1] > 0, "premise: the standard predictor must cut back here"
    assert tangent[1] == 0
    assert tangent[:3] == elastic[:3]
    assert tangent[2] < stock[2]
    assert tangent[3] == pytest.approx(0.5 * DELTA, rel=1e-9)
    assert tangent[3] == pytest.approx(stock[3], rel=1e-9)


def test_tangent_predictor_fixed_march():
    """Ten equal steps: iterations drop to the elastic control, answer kept."""
    it_stock, u_stock = _fixed_march("j2", STOCK)
    it_elastic, _ = _fixed_march("elastic", STOCK)
    it_tangent, u_tangent = _fixed_march("j2", TANGENT)

    assert it_tangent == it_elastic
    assert it_tangent < it_stock
    assert u_tangent == pytest.approx(0.5 * DELTA, rel=1e-9)
    assert u_tangent == pytest.approx(u_stock, rel=1e-9)


def test_tangent_predictor_after_the_triple():
    """The option is also accepted after numIter/minLambda/maxLambda."""
    _build("j2", max_iter=40)
    iters = 0
    for _ in range(10):
        ops.integrator("LoadControl", 0.1, 1, 0.1, 0.1, "-tangentPredictor")
        assert ops.analyze(1) == 0
        iters += ops.testIter()
    it_elastic, _ = _fixed_march("elastic", STOCK)
    assert iters == it_elastic
    assert ops.nodeDisp(MID_NODE, 3) == pytest.approx(0.5 * DELTA, rel=1e-9)


@pytest.mark.parametrize("drive, handler", [("load", "Transformation"),
                                            ("disp", "Penalty")])
def test_no_contributing_element_falls_back(capfd, drive, handler):
    """No element supplies K*du_p: the SPs must be enforced at once.

    Without the fallback the first iteration of each step would see no
    prescribed motion, dU would vanish and NormDispIncr would accept an
    unmoved step. Covers a load-driven model (no non-homogeneous sp) and a
    constraint handler other than Transformation.
    """
    if handler == "Penalty":
        handler_args = ("Penalty", 1.0e12, 1.0e12)
    else:
        handler_args = (handler,)

    capfd.readouterr()
    it_tangent, u_tangent = _fixed_march("j2", TANGENT, drive, handler_args)
    err = capfd.readouterr().err
    it_stock, u_stock = _fixed_march("j2", STOCK, drive, handler_args)

    assert "-tangentPredictor: no element supplied" in err
    assert it_tangent == it_stock
    assert u_tangent == u_stock
    assert u_tangent == pytest.approx(0.5 * DELTA, rel=1e-6)


TCL_SCRIPT = """
set L 100.0; set AX 10.0; set E 200000.0; set N 20
model basic -ndm 3 -ndf 3
set h [expr $L/$N]
for {set k 0} {$k <= $N} {incr k} {
    set z [expr $k*$h]; set b [expr 4*$k]
    node [expr $b+1] 0.0 0.0 $z
    node [expr $b+2] $AX 0.0 $z
    node [expr $b+3] $AX $AX $z
    node [expr $b+4] 0.0 $AX $z
}
nDMaterial J2Plasticity 1 [expr $E/3.0] [expr $E/2.0] 379.5 379.5 0.0 2000.0
for {set k 0} {$k < $N} {incr k} {
    set b [expr 4*$k]
    element stdBrick [expr $k+1] [expr $b+1] [expr $b+2] [expr $b+3] [expr $b+4] \\
        [expr $b+5] [expr $b+6] [expr $b+7] [expr $b+8] 1
}
for {set k 0} {$k <= $N} {incr k} {
    for {set j 1} {$j <= 4} {incr j} {
        if {$k == 0} { fix [expr 4*$k+$j] 1 1 1 } else { fix [expr 4*$k+$j] 1 1 0 }
    }
}
timeSeries Linear 1
pattern Plain 1 1 {
    for {set j 1} {$j <= 4} {incr j} { sp [expr 4*$N+$j] 3 0.15 }
}
constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-10 40 0
algorithm KrylovNewton
analysis Static
set iters 0
for {set s 0} {$s < 10} {incr s} {
    integrator LoadControl 0.1 %s
    if {[analyze 1] != 0} { puts "RESULT failed"; exit }
    set iters [expr $iters + [testIter]]
}
puts "RESULT $iters [nodeDisp %d 3]"
"""


def _tcl_executable():
    here = Path(ops.__file__).resolve().parent
    for name in ("OpenSees.exe", "OpenSees"):
        exe = here / name
        if exe.is_file():
            return exe
    return None


@pytest.mark.skipif(_tcl_executable() is None,
                    reason="no OpenSees Tcl executable next to the module")
def test_tcl_parser(tmp_path):
    exe = _tcl_executable()

    def run(flag):
        script = tmp_path / "chain.tcl"
        script.write_text(TCL_SCRIPT % (flag, MID_NODE))
        res = subprocess.run([str(exe), str(script)], capture_output=True,
                             text=True, cwd=tmp_path, timeout=300,
                             env=dict(os.environ))
        line = [x for x in (res.stdout + res.stderr).splitlines()
                if x.startswith("RESULT")]
        assert line, res.stdout + res.stderr
        fields = line[-1].split()
        assert fields[1] != "failed"
        return int(fields[1]), float(fields[2])

    it_stock, u_stock = run("")
    it_tangent, u_tangent = run("-tangentPredictor")
    assert it_tangent < it_stock
    assert u_tangent == pytest.approx(u_stock, rel=1e-9)
    assert u_tangent == pytest.approx(0.075, rel=1e-9)
