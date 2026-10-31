"""constraints LadrunoProjection: equalDOF ties under CentralDifferenceLadruno.

The handler keeps every constrained DOF in the equation set with its own
diagonal mass, and CentralDifferenceLadruno projects the starting and every
solved acceleration onto the constraint manifold with
a_proj = L (L^T M L)^-1 L^T M a_raw.

Checks:
- two masses (2 and 3) tied by equalDOF under a constant force of 10: the
  tie holds to 1e-12, the acceleration is F/(m1 + m2) and the linear
  momentum equals the impulse F t to 1e-12 (relative);
- a tied spring-mass oscillator reproduces constraints Transformation (exact
  for equalDOF, C = I) to 1e-9;
- the tie does not move the stability limit: stable at 0.99 dt_cr, divergent
  at 1.02 dt_cr, as the untied oscillator with the combined mass;
- refusals, each making analyze() fail: a DOF constrained twice, a constraint
  chain, a constrained DOF that is also fixed, a massless constrained DOF
  (refused at handle(), whatever the order of the commands), a system other
  than Diagonal, an integrator that does not apply the projection;
- initial conditions off the constraint manifold make the analysis fail, and
  -projectICs sets them on the manifold;
- removing the last constraint mid-run rebuilds with no projector and runs.
"""

import math

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"


def _analysis(handler=(PROJ,), system=("Diagonal",), integrator=(CDL,)):
    ops.constraints(*handler)
    ops.numberer("Plain")
    ops.system(*system)
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(*integrator)
    ops.analysis("Transient")


def test_tie_combined_mass_and_momentum():
    m1, m2, F = 2.0, 3.0, 10.0
    a_exact = F / (m1 + m2)
    dt, nsteps = 1.0e-3, 400

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.equalDOF(1, 2, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, F)
    _analysis()

    t = 0.0
    tie_err = mom_err = a_err = 0.0
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        t += dt
        u1, u2 = ops.nodeDisp(1, 1), ops.nodeDisp(2, 1)
        v1, v2 = ops.nodeVel(1, 1), ops.nodeVel(2, 1)
        tie_err = max(tie_err, abs(u1 - u2), abs(v1 - v2))
        mom_err = max(mom_err, abs(m1 * v1 + m2 * v2 - F * t))
        a_err = max(a_err, abs(ops.nodeAccel(1, 1) - a_exact))

    assert tie_err < 1e-12
    assert a_err < 1e-10
    assert mom_err < 1e-12 * (1.0 + F * t)


def _run_tied_shm(handler_args, m1, m2, k, d0, dt, nsteps):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
    ops.equalDOF(1, 2, 1)
    # both tied DOFs start at d0 (the constraint was built at u = 0)
    ops.setNodeDisp(1, 1, d0, "-commit")
    ops.setNodeDisp(2, 1, d0, "-commit")
    _analysis(handler=handler_args)
    u = [ops.nodeDisp(1, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        u.append(ops.nodeDisp(1, 1))
    return np.array(u)


def test_matches_transformation_reference():
    m1, m2, k, d0 = 2.0, 3.0, 100.0, 1.0
    dt = 0.2 * (2.0 / math.sqrt(k / (m1 + m2)))
    nsteps = 300

    u_proj = _run_tied_shm((PROJ,), m1, m2, k, d0, dt, nsteps)
    u_ref = _run_tied_shm(("Transformation",), m1, m2, k, d0, dt, nsteps)

    assert len(u_proj) == nsteps + 1 and len(u_ref) == nsteps + 1
    scale = np.max(np.abs(u_ref))
    assert np.max(np.abs(u_proj - u_ref)) < 1e-9 * (1.0 + scale)
    assert scale > 0.5 * d0


def _peak_single(M, k, d0, dt, nsteps):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)
    ops.node(1, 0.0); ops.mass(1, M)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
    ops.setNodeDisp(1, 1, d0, "-commit")
    _analysis(handler=("Plain",))
    peak = abs(d0)
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return float("inf")
        peak = max(peak, abs(ops.nodeDisp(1, 1)))
        if not math.isfinite(peak) or peak > 1e6 * abs(d0):
            return float("inf")
    return peak


def _peak_tied(m1, m2, k, d0, dt, nsteps):
    u = _run_tied_shm((PROJ,), m1, m2, k, d0, dt, nsteps)
    if len(u) < nsteps + 1:
        return float("inf")
    peak = float(np.max(np.abs(u)))
    if not math.isfinite(peak) or peak > 1e6 * abs(d0):
        return float("inf")
    return peak


def test_tie_does_not_change_stability_limit():
    m1, m2, k, d0 = 2.0, 3.0, 100.0, 1.0
    dt_cr = 2.0 / math.sqrt(k / (m1 + m2))
    nsteps = 600

    a_stable = _peak_single(m1 + m2, k, d0, 0.99 * dt_cr, nsteps)
    a_unstable = _peak_single(m1 + m2, k, d0, 1.02 * dt_cr, nsteps)
    b_stable = _peak_tied(m1, m2, k, d0, 0.99 * dt_cr, nsteps)
    b_unstable = _peak_tied(m1, m2, k, d0, 1.02 * dt_cr, nsteps)

    assert math.isfinite(a_stable) and a_stable < 10.0 * d0
    assert not math.isfinite(a_unstable)
    assert math.isfinite(b_stable) and b_stable < 10.0 * d0
    assert not math.isfinite(b_unstable)


def _require_handler():
    # the refusal must come from the handler, not from a build that lacks it
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.constraints("LadrunoProjection")
    ops.wipe()


def _runs_ok(build):
    _require_handler()
    try:
        build()
        return ops.analyze(1, 1e-4) == 0
    except Exception:
        return False


def _two_masses(m2=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, 1.0)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.equalDOF(1, 2, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 1.0)


def test_double_constraint_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        for n in (1, 2, 3):
            ops.node(n, float(n)); ops.mass(n, 1.0)
        ops.equalDOF(1, 3, 1)
        ops.equalDOF(2, 3, 1)
        _analysis()
    assert not _runs_ok(build)


def test_chain_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        for n in (1, 2, 3):
            ops.node(n, float(n)); ops.mass(n, 1.0)
        ops.equalDOF(1, 2, 1)
        ops.equalDOF(2, 3, 1)
        _analysis()
    assert not _runs_ok(build)


def test_fixed_constrained_dof_refused():
    def build():
        _two_masses()
        ops.fix(2, 1)
        _analysis()
    assert not _runs_ok(build)


def test_massless_constrained_dof_refused():
    def build():
        _two_masses(m2=0.0)
        _analysis()
    assert not _runs_ok(build)


def test_massless_constrained_dof_refused_handler_first():
    # the handler is defined before the constraint and the mass
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        _analysis()
        ops.node(1, 0.0); ops.mass(1, 1.0)
        ops.node(2, 1.0)
        ops.equalDOF(1, 2, 1)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, 1.0)
    assert not _runs_ok(build)


def test_non_diagonal_system_refused():
    def build():
        _two_masses()
        _analysis(system=("FullGeneral",))
    assert not _runs_ok(build)


def test_integrator_without_projection_refused():
    def build():
        _two_masses()
        _analysis(integrator=("Newmark", 0.5, 0.25))
    assert not _runs_ok(build)


def test_initial_conditions_off_manifold():
    def build(handler):
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, 2.0)
        ops.node(2, 1.0); ops.mass(2, 3.0)
        ops.equalDOF(1, 2, 1)
        ops.setNodeDisp(1, 1, 1.0, "-commit")
        ops.setNodeDisp(2, 1, 2.0, "-commit")   # the tie requires u1 == u2
        _analysis(handler=handler)

    build((PROJ,))
    assert ops.analyze(1, 1e-3) != 0

    build((PROJ, "-projectICs"))
    assert ops.analyze(1, 1e-3) == 0
    assert abs(ops.nodeDisp(1, 1) - ops.nodeDisp(2, 1)) < 1e-9


def test_rebuild_without_constraints():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, 2.0)
    ops.node(2, 1.0); ops.mass(2, 3.0)
    ops.equalDOF(1, 2, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 10.0)
    _analysis()
    assert ops.analyze(5, 1e-3) == 0
    ops.remove("mp", 2)
    assert ops.analyze(5, 1e-3) == 0

