"""LadrunoIndirectControl static integrator.

Indirect displacement control prescribes the increment of a weighted sum of
DOFs, zeta = c . U = sum_k coef_k U(node_k, dof_k), per step.

Checks:
  - with one entry and coef = 1 the path on the von Mises truss (two
    corotTruss bars, E A = 1e4, half span 1.0, rise 0.1, apex load -1.0,
    through the snap-through limit point) equals DisplacementControl step for
    step (1e-8);
  - the constraint holds exactly, coef * dU = Delta-zeta (1e-9), for
    coef = 1, 2 and 0.5;
  - a relative two-DOF quantity u(3,1) - u(2,1) (the elongation of the second
    bar of a two-bar chain, a crack-mouth-opening type control) increases by
    the prescribed amount each step;
  - the -iter option is accepted with and without dmin dmax, and before or
    after the -dof entries;
  - a missing -dof entry or a nonexistent control node is rejected.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


LIC = "LadrunoIndirectControl"
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


def _setup_solver():
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
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


def test_single_dof_equals_displacement_control():
    incr, nsteps = -0.004, 40
    ref = _trace(("DisplacementControl", APEX, 2, incr), nsteps)
    lic = _trace((LIC, incr, "-dof", APEX, 2, 1.0), nsteps)
    assert len(ref) > 20
    assert len(lic) == len(ref)
    for k, ((lr, ur), (ll, ul)) in enumerate(zip(ref, lic)):
        assert abs(lr - ll) < 1.0e-8, (k, lr, ll)
        assert abs(ur - ul) < 1.0e-8, (k, ur, ul)
    # the path goes through the limit point: the load factor falls again
    lams = [p[0] for p in ref]
    assert max(lams) > lams[-1]


@pytest.mark.parametrize("coef", [1.0, 2.0, 0.5])
def test_constraint_is_exact(coef):
    incr = -0.004
    _build_truss()
    _setup_solver()
    ops.integrator(LIC, incr, "-dof", APEX, 2, coef)
    ops.analysis("Static")
    prev = ops.nodeDisp(APEX, 2)
    nconv = 0
    for _ in range(15):
        if ops.analyze(1) != 0:
            break
        now = ops.nodeDisp(APEX, 2)
        assert abs(coef * (now - prev) - incr) < 1.0e-9
        prev = now
        nconv += 1
    assert nconv >= 10


def _build_chain():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 2.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(3, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.element("truss", 2, 2, 3, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0, 0.0)
    _setup_solver()


def test_relative_two_dof_control():
    _build_chain()
    incr = 0.002
    ops.integrator(LIC, incr, "-dof", 3, 1, 1.0, "-dof", 2, 1, -1.0)
    ops.analysis("Static")
    prev = ops.nodeDisp(3, 1) - ops.nodeDisp(2, 1)
    nconv = 0
    for _ in range(10):
        if ops.analyze(1) != 0:
            break
        now = ops.nodeDisp(3, 1) - ops.nodeDisp(2, 1)
        assert abs((now - prev) - incr) < 1.0e-9
        prev = now
        nconv += 1
    assert nconv >= 8
    # elongation 0.002 per step of a bar with E A = 1e4, length 1: load 20
    assert ops.getTime() == pytest.approx(20.0 * nconv, rel=1e-9)


@pytest.mark.parametrize("args", [
    (-0.004, "-dof", APEX, 2, 1.0, "-iter", 4, -0.02, -0.001),
    (-0.004, "-dof", APEX, 2, 1.0, "-iter", 4),
    (-0.004, "-iter", 4, "-dof", APEX, 2, 1.0),
    (-0.004, "-iter", 4, -0.02, -0.001, "-dof", APEX, 2, 1.0),
])
def test_parser_forms(args):
    _build_truss()
    _setup_solver()
    ops.integrator(LIC, *args)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(APEX, 2) < 0.0


@pytest.mark.parametrize("args", [
    (-0.004,),
    (-0.004, "-dof", 99, 2, 1.0),
    (-0.004, "-dof", APEX, 3, 1.0),
])
def test_parser_rejects_bad_input(args):
    _build_truss()
    with pytest.raises(Exception):
        ops.integrator(LIC, *args)
