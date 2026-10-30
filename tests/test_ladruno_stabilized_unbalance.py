"""LadrunoStabilizedUnbalance convergence test.

Model: the shallow von Mises truss of test_ladruno_arclength.py (two
corotTruss bars, E A = 1e4, half span 1.0, rise 0.1, apex reference load
1.0, limit load about 3.80).

Checks:
  - with LoadControl (no stabilization) the test behaves as NormUnbalance:
    three steps of 0.05 converge to the same apex displacement (1e-9);
  - under LadrunoArcLength -stabilize (load increment 0.1, c = 0.01) the
    stock NormUnbalance accepts the first step at 1e-8, because it norms
    the residual with the viscous force subtracted; LadrunoStabilizedUnbalance
    rejects the same step at 1e-8, because the static unbalance stays of the
    size of the viscous force, and accepts it at 0.1;
  - the full argument form (tol maxIter printFlag normType) is accepted
    by Python and, when an OpenSees Tcl executable is found, by classic Tcl.
"""
import os
import subprocess
from pathlib import Path

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


LSU = "LadrunoStabilizedUnbalance"
APEX = 2


def _build_truss():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, -1.0, 0.0)
    ops.node(APEX, 0.0, 0.1)
    ops.node(3, 1.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(APEX, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)
    ops.element("corotTruss", 1, 1, APEX, 1.0, 1)
    ops.element("corotTruss", 2, APEX, 3, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(APEX, 0.0, -1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")


def _loadcontrol_tip(test_cmd):
    _build_truss()
    ops.test(*test_cmd)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")
    ok = ops.analyze(3)
    return ok, ops.nodeDisp(APEX, 2)


def test_fallback_equals_normunbalance():
    ok_ref, uy_ref = _loadcontrol_tip(("NormUnbalance", 1.0e-8, 50))
    ok_lsu, uy_lsu = _loadcontrol_tip((LSU, 1.0e-8, 50))
    assert ok_ref == 0
    assert ok_lsu == 0
    assert abs(uy_ref - uy_lsu) < 1.0e-9


def _stabilized_step(test_cmd, cVisc=0.01):
    _build_truss()
    ops.test(*test_cmd)
    ops.algorithm("Newton")
    ops.integrator("LadrunoArcLength", 0.1, 1.0, "-stabilize", 2.0e-4,
                   "-cVisc", cVisc)
    ops.analysis("Static")
    return ops.analyze(1)


def test_norms_the_residual_without_the_viscous_force():
    assert _stabilized_step(("NormUnbalance", 1.0e-8, 100)) == 0
    assert _stabilized_step((LSU, 1.0e-8, 100)) != 0
    assert _stabilized_step((LSU, 0.1, 100)) == 0


def test_full_argument_form():
    _build_truss()
    ops.test(LSU, 1.0e-8, 50, 0, 2)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


TCL_SCRIPT = """
model basic -ndm 2 -ndf 2
node 1 -1.0 0.0
node 2  0.0 0.1
node 3  1.0 0.0
fix 1 1 1
fix 3 1 1
fix 2 1 0
uniaxialMaterial Elastic 1 1.0e4
element corotTruss 1 1 2 1.0 1
element corotTruss 2 2 3 1.0 1
timeSeries Linear 1
pattern Plain 1 1 { load 2 0.0 -1.0 }
constraints Plain
numberer Plain
system BandGeneral
test LadrunoStabilizedUnbalance 1.0e-8 50 0 2
algorithm Newton
integrator LoadControl 0.05
analysis Static
set ok [analyze 3]
puts "RESULT $ok [nodeDisp 2 2]"
"""


def _tcl_executable():
    env = os.environ.get("OPENSEES_EXE")
    if env and Path(env).is_file():
        return Path(env)
    here = Path(ops.__file__).resolve().parent
    for name in ("OpenSees.exe", "OpenSees"):
        exe = here / name
        if exe.is_file():
            return exe
    return None


@pytest.mark.skipif(_tcl_executable() is None,
                    reason="no OpenSees Tcl executable (set OPENSEES_EXE)")
def test_classic_tcl_parser(tmp_path):
    _ok, uy_ref = _loadcontrol_tip(("NormUnbalance", 1.0e-8, 50))
    script = tmp_path / "stabilized_unbalance.tcl"
    script.write_text(TCL_SCRIPT)
    res = subprocess.run([str(_tcl_executable()), str(script)],
                         capture_output=True, text=True, cwd=tmp_path,
                         timeout=300, stdin=subprocess.DEVNULL,
                         env=dict(os.environ))
    out = res.stdout + res.stderr
    line = [x for x in out.splitlines() if x.startswith("RESULT")]
    assert line, out
    fields = line[-1].split()
    assert int(fields[1]) == 0, out
    assert float(fields[2]) == pytest.approx(uy_ref, rel=1e-9)
