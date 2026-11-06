"""LadrunoArcLength static integrator.

Model: a shallow von Mises truss, two corotTruss bars (E A = 1e4, half span
1.0, rise 0.1) with the apex loaded downward by a reference load of 1.0. The
snap-through limit load is about 3.80, and the far stable branch lies below
uy = -2 rise.

Checks:
  - without options LadrunoArcLength reproduces stock ArcLength step for step
    (load factor and apex displacement to 1e-7 relative);
  - -adapt (Ramm iteration-count rule, Jd = 8) grows the radius from 0.004
    and reaches load factor 3.0 in fewer steps than the fixed radius;
  - with Newton limited to 25 iterations per step, stock LoadControl
    (increment 0.1) stops at the limit point (the uncontrolled jump to the
    far branch takes about 76 iterations), while -stabilize (viscous
    regularization, c = 0.03, same increment and Newton settings) passes the
    limit load with at most 10 iterations per step and ends on the far
    branch at load factor 6.0;
  - after a failed step the converged per-step state (deltaLambdaStep,
    sign, ||deltaUstep||, load factor, apex displacement) is restored, and a
    script loop of `ladrunoArcLength reduceStep 0.5` + analyze recovers;
  - the documented argument forms are accepted;
  - the classic Tcl interpreter answers `ladrunoArcLength` (run only when
    an OpenSees Tcl executable is found).
"""
import math
import os
import subprocess
from pathlib import Path

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


LAL = "LadrunoArcLength"
H = 0.10
SPAN = 1.0
E = 1.0e4
A = 1.0
APEX = 2
LIMIT = 3.80


def _build_truss():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, -SPAN, 0.0)
    ops.node(APEX, 0.0, H)
    ops.node(3, SPAN, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(APEX, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.element("corotTruss", 1, 1, APEX, A, 1)
    ops.element("corotTruss", 2, APEX, 3, A, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(APEX, 0.0, -1.0)


def _setup_solver(tol=1.0e-10, iters=50):
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")


def _trace(integrator_cmd, nsteps):
    _build_truss()
    _setup_solver()
    ops.integrator(*integrator_cmd)
    ops.analysis("Static")
    path = []
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        path.append((ops.getTime(), ops.nodeDisp(APEX, 2)))
    return path


def test_identity_with_arclength():
    arc, alpha, nsteps = 0.01, 0.1, 25
    ref = _trace(("ArcLength", arc, alpha), nsteps)
    lad = _trace((LAL, arc, alpha), nsteps)
    assert len(ref) > 10
    assert len(lad) == len(ref)
    for i, ((lr, ur), (ll, ul)) in enumerate(zip(ref, lad)):
        assert abs(ll - lr) <= 1e-9 + 1e-7 * abs(lr), (i, ll, lr)
        assert abs(ul - ur) <= 1e-9 + 1e-7 * abs(ur), (i, ul, ur)


def _steps_to_reach(path, target):
    for i, (lam, _uy) in enumerate(path):
        if lam >= target:
            return i + 1
    return None


def test_adapt_reaches_target_in_fewer_steps():
    arc, alpha, target = 0.004, 0.1, 3.0
    fixed = _trace((LAL, arc, alpha), 400)
    adapt = _trace((LAL, arc, alpha, "-adapt", 8, 1.0e-3, 5.0e-2), 400)
    n_fixed = _steps_to_reach(fixed, target)
    n_adapt = _steps_to_reach(adapt, target)
    assert n_fixed is not None
    assert n_adapt is not None
    assert n_adapt < n_fixed, (n_adapt, n_fixed)


# Newton settings shared by the control and the stabilized run. With an
# elastic truss an unlimited Newton iteration at the first load step above
# the limit can still land on the far branch (measured: 76 iterations at
# lambda = 3.9 with LoadControl 0.1), so the control limits the iterations.
STAB_ITERS = 25


def _load_steps(integrator_cmd, nsteps=60):
    _build_truss()
    _setup_solver(1.0e-8, STAB_ITERS)
    ops.integrator(*integrator_cmd)
    ops.analysis("Static")
    lam, uy, maxit = 0.0, 0.0, 0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        lam = ops.getTime()
        uy = ops.nodeDisp(APEX, 2)
        maxit = max(maxit, ops.testIter())
    return lam, uy, maxit


def test_loadcontrol_stops_at_limit_point():
    # the reference failure that -stabilize removes
    # last converged step lambda = 3.8 on the near branch; the step to 3.9
    # fails
    lam, uy, _it = _load_steps(("LoadControl", 0.1))
    assert lam < LIMIT + 0.05, lam
    assert uy > -H, uy


def test_stabilize_clears_limit_point():
    lam, uy, maxit = _load_steps(
        (LAL, 0.1, 1.0, "-stabilize", 2.0e-4, "-cVisc", 0.03))
    assert lam == pytest.approx(6.0), lam
    assert uy < -2.0 * H, uy
    assert maxit <= 10, maxit


def test_revert_restores_committed_state_and_retry_recovers():
    arc, alpha = 0.05, 0.1
    _build_truss()
    _setup_solver()
    ops.integrator(LAL, arc, alpha)
    ops.analysis("Static")

    # a fixed radius this large fails near the limit point
    committed = None
    failed = False
    for _ in range(2000):
        if ops.analyze(1) != 0:
            failed = True
            break
        committed = dict(
            lam=ops.getTime(),
            uy=ops.nodeDisp(APEX, 2),
            dlds=ops.ladrunoArcLength("deltaLambdaStep"),
            sign=ops.ladrunoArcLength("sign"),
            nrm=ops.ladrunoArcLength("deltaUstepNorm"),
            al=ops.ladrunoArcLength("arcLength"),
        )
    assert failed
    assert committed is not None
    assert committed["lam"] > 3.0
    assert abs(committed["al"] - arc) <= 1e-12
    assert committed["nrm"] > 0.0

    assert ops.getTime() == pytest.approx(committed["lam"])
    assert ops.nodeDisp(APEX, 2) == pytest.approx(committed["uy"])
    assert ops.ladrunoArcLength("deltaLambdaStep") == pytest.approx(committed["dlds"])
    assert ops.ladrunoArcLength("sign") == committed["sign"]
    assert ops.ladrunoArcLength("deltaUstepNorm") == pytest.approx(committed["nrm"])

    recovered = False
    for _ in range(14):
        ops.ladrunoArcLength("reduceStep", 0.5)
        if ops.analyze(1) == 0:
            recovered = True
            break
    assert recovered
    assert ops.ladrunoArcLength("arcLength") < committed["al"]
    moved = (ops.getTime() != pytest.approx(committed["lam"])) or (
        ops.nodeDisp(APEX, 2) != pytest.approx(committed["uy"]))
    assert moved


def test_command_mutators_and_guards():
    _build_truss()
    _setup_solver()
    ops.integrator(LAL, 0.02, 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    al = ops.ladrunoArcLength()
    assert al == pytest.approx(0.02)
    assert ops.ladrunoArcLength("reduceStep", 0.5) == pytest.approx(0.01)
    assert ops.ladrunoArcLength("increaseStep", 3.0) == pytest.approx(0.03)
    assert ops.ladrunoArcLength("setArcLength", 0.015) == pytest.approx(0.015)
    with pytest.raises(Exception):
        ops.ladrunoArcLength("notASubcommand")
    ops.integrator("LoadControl", 0.1)
    with pytest.raises(Exception):
        ops.ladrunoArcLength("arcLength")


@pytest.mark.parametrize("cmd", [
    (LAL, 0.01, 1.0),
    (LAL, 0.01, 0.1, "-adapt", 5, 1.0e-3, 1.0e-1),
    (LAL, 0.01, 0.1, "-adapt", 5, 1.0e-3, 1.0e-1, "-p", 0.5),
    (LAL, 0.1, 1.0, "-stabilize", 2.0e-4),
    (LAL, 0.1, 1.0, "-stabilize", 2.0e-4, "-cVisc", 0.5),
    (LAL, 0.1, 1.0, "-stabilize", 2.0e-4, "-mass", "gershgorin"),
    (LAL, 0.1, 1.0, "-stabilize", "-adaptStab", "-mass", "unity"),
])
def test_parser_forms(cmd):
    _build_truss()
    _setup_solver()
    ops.integrator(*cmd)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, cmd
    assert math.isfinite(ops.nodeDisp(APEX, 2))


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
test NormDispIncr 1.0e-10 50 0
algorithm Newton
integrator LadrunoArcLength 0.02 1
analysis Static
analyze 1
set al [ladrunoArcLength arcLength]
set bare [ladrunoArcLength]
set half [ladrunoArcLength reduceStep 0.5]
set bad [catch {ladrunoArcLength notASubcommand} msg]
puts "RESULT $al $bare $half $bad"
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
def test_classic_tcl_command(tmp_path):
    script = tmp_path / "arclength.tcl"
    script.write_text(TCL_SCRIPT)
    res = subprocess.run([str(_tcl_executable()), str(script)],
                         capture_output=True, text=True, cwd=tmp_path,
                         timeout=300, stdin=subprocess.DEVNULL,
                         env=dict(os.environ))
    out = res.stdout + res.stderr
    assert "invalid command name" not in out, out
    line = [x for x in out.splitlines() if x.startswith("RESULT")]
    assert line, out
    al, bare, half, bad = (float(v) for v in line[-1].split()[1:5])
    assert al == pytest.approx(0.02)
    assert bare == pytest.approx(al)
    assert half == pytest.approx(0.5 * al)
    assert bad == 1.0
