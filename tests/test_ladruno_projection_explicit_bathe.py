"""constraints LadrunoProjection under ExplicitBathe.

ExplicitBathe reads diag(M) for the projector on the first step after a
change, projects the committed acceleration there, and projects both
sub-step accelerations (t + p dt and t + dt) before they form velocities
and displacements. Without a projector the scheme is unchanged
(tests/test_explicit_bathe_default_path.py checks it bit for bit).

Model: masses 2 and 3 tied by equalDOF, spring k = 100 to ground, both
starting at u = 1, dt = 0.2 (2/omega), p = 0.54.

Checks:
- the tie holds to 1e-12 in displacement, velocity and acceleration;
- the response matches constraints Transformation (exact for equalDOF) to
  1e-9, with and without local non-viscous damping (-lnvd 0);
- with -lnvd 0.8 the tie still holds to 1e-12 and the motion decays;
- refusal of a massless constrained DOF.
"""

import math

import numpy as np

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


PROJ = "LadrunoProjection"
_M1, _M2, _K, _D0 = 2.0, 3.0, 100.0, 1.0
_DT = 0.2 * (2.0 / math.sqrt(_K / (_M1 + _M2)))
_N = 300


def _run_tied_shm(integrator_args, handler_args, m2=_M2):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)
    ops.node(1, 0.0); ops.mass(1, _M1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
    ops.equalDOF(1, 2, 1)
    ops.setNodeDisp(1, 1, _D0, "-commit")
    ops.setNodeDisp(2, 1, _D0, "-commit")
    ops.constraints(*handler_args)
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    u1 = [ops.nodeDisp(1, 1)]
    tie = 0.0
    for _ in range(_N):
        if ops.analyze(1, _DT) != 0:
            break
        u1.append(ops.nodeDisp(1, 1))
        tie = max(tie, abs(ops.nodeDisp(1, 1) - ops.nodeDisp(2, 1)),
                  abs(ops.nodeVel(1, 1) - ops.nodeVel(2, 1)),
                  abs(ops.nodeAccel(1, 1) - ops.nodeAccel(2, 1)))
    return np.array(u1), tie


def test_tie_exact():
    u, tie = _run_tied_shm(("ExplicitBathe", 0.54), (PROJ,))
    assert len(u) == _N + 1
    assert tie < 1e-12, tie
    assert np.max(np.abs(u)) > 0.5 * _D0


def _match_transformation(integ):
    u_proj, _ = _run_tied_shm(integ, (PROJ,))
    u_ref, _ = _run_tied_shm(integ, ("Transformation",))
    assert len(u_proj) == _N + 1 and len(u_ref) == _N + 1
    scale = float(np.max(np.abs(u_ref)))
    assert float(np.max(np.abs(u_proj - u_ref))) < 1e-9 * (1.0 + scale)
    assert scale > 0.5 * _D0


def test_matches_transformation():
    _match_transformation(("ExplicitBathe", 0.54))


def test_matches_transformation_lnvd_zero():
    _match_transformation(("ExplicitBathe", 0.54, "-lnvd", 0.0))


def test_tie_exact_with_local_damping():
    u, tie = _run_tied_shm(("ExplicitBathe", 0.54, "-lnvd", 0.8), (PROJ,))
    assert len(u) == _N + 1
    assert tie < 1e-12, tie
    assert abs(u[-1]) < abs(u[0]) + 1e-9


def test_massless_constrained_dof_refused():
    # the refusal must come from the handler, not from a build that lacks it
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.constraints(PROJ)
    ops.wipe()
    try:
        u, _ = _run_tied_shm(("ExplicitBathe", 0.54), (PROJ,), m2=0.0)
        ran = len(u) == _N + 1
    except Exception:
        ran = False
    assert not ran
