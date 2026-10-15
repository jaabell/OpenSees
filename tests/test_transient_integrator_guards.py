"""Transient integrator guards: zero time step and the Hall tangent.

Model: one-DOF oscillator, truss of stiffness k = E*A/L = 100 between a
fixed node and a free node of mass m = 1, under a constant load P = 1.

Zero time step: TRBDF2, TRBDF3, Houbolt and BackwardEuler form their
integration constants as 2/dt and 4/dt^2. With dt = 0.0 the constants are
inf, the predicted acceleration -4/dt * 0 is NaN, and with algorithm Linear
analyze(1, 0.0) returns 0 and commits a NaN state. Each integrator must
instead fail the step (negative return code), leave the committed state
finite (zero), and accept a following step with dt > 0.

Hall tangent: BackwardEuler, Newmark1 and Collocation must assemble the
element stiffness for algorithm Newton -hall, so that a linear model gives
the same displacement history as with the default current tangent.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


K = 100.0
M = 1.0
P = 1.0


def _build_sdof(integrator, algorithm):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.mass(2, M)
    ops.uniaxialMaterial("Elastic", 1, K)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, P)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm(*algorithm)
    ops.integrator(*integrator)
    ops.analysis("Transient")


def _state_is_finite(tag):
    values = (ops.nodeDisp(tag, 1), ops.nodeVel(tag, 1), ops.nodeAccel(tag, 1))
    return all(math.isfinite(v) for v in values)


@pytest.mark.parametrize("integrator", [
    ("TRBDF2",),
    ("TRBDF3",),
    ("Houbolt",),
    ("BackwardEuler",),
])
def test_zero_time_step_is_rejected(integrator):
    try:
        _build_sdof(integrator, ("Linear",))
        assert ops.analyze(1, 0.0) < 0
        assert _state_is_finite(2)
        assert ops.nodeDisp(2, 1) == 0.0
        assert ops.analyze(1, 0.01) == 0
        assert _state_is_finite(2)
        assert ops.nodeDisp(2, 1) > 0.0
    finally:
        ops.wipe()


def _history(integrator, algorithm, num_steps=10, dt=0.5):
    _build_sdof(integrator, algorithm)
    disp = []
    for _ in range(num_steps):
        assert ops.analyze(1, dt) == 0
        disp.append(ops.nodeDisp(2, 1))
    return disp


@pytest.mark.parametrize("integrator", [
    ("BackwardEuler",),
    ("Newmark1", 0.5, 0.25),
    ("Collocation", 1.0, 0.25, 0.5),
])
def test_hall_tangent_matches_current_tangent(integrator):
    """Newton -hall with iFactor + cFactor = 1 on a linear model.

    For an elastic truss Ki = Kt = k, so the Hall tangent
    c1*(0.1*Ki + 0.9*Kt) equals the current tangent and both runs must
    give the same displacement history. Before the fix the three
    integrators add no element term for HALL_TANGENT, so the tangent is
    only the nodal c3*m = 16 (c3 = 4/dt^2 = 1/(beta*dt^2) = 16 at
    dt = 0.5) against the true 116; the Newton correction is amplified by
    k/(c3*m) = 6.25 per iteration and the first step fails to converge.
    """
    try:
        reference = _history(integrator, ("Newton",))
        hall = _history(integrator, ("Newton", "-hall", 0.1, 0.9))
        assert reference[-1] != 0.0
        assert hall == pytest.approx(reference, rel=1.0e-10, abs=1.0e-14)
    finally:
        ops.wipe()
