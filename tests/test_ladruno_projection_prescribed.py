"""constraints LadrunoProjection: prescribed motion (non-homogeneous SP and
imposedMotion) under CentralDifferenceLadruno.

A prescribed DOF is removed from the equation set and its displacement is
imposed each step by the handler's applyLoad(), as
TransformationDOF_Group::enforceSPs() does; ImposedMotionSP sets the
velocity and acceleration itself. A constrained DOF tied only to prescribed
DOFs is removed as well and its displacement, velocity and acceleration are
set from them each step.

Checks (references: constraints Transformation with system FullGeneral):
- base excitation by imposedMotion on a support: response of the mass and
  of the support match the reference to 1e-6 (relative RMS);
- a constant prescribed settlement: response matches the reference to 1e-6;
- imposedMotion and a constant SP on an equalDOF retained DOF: the
  constrained DOF follows it to 1e-10, the downstream DOF matches the
  reference to 1e-9;
- rigidLink -bar with only the retained ux prescribed: the constrained ux
  follows it to 1e-10 while the uy tie (both free) is held to 1e-10;
- refusals: a prescribed constrained DOF, a constrained DOF tied to both a
  free and a prescribed DOF, a model with no free equation.
"""

import numpy as np

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"
_M, _K = 1.0, 100.0


def _rms(a):
    return float(np.sqrt(np.mean(a * a)))


def _rel(a, b):
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    d = _rms(a) + _rms(b)
    return 0.0 if d == 0.0 else _rms(a - b) / (0.5 * d)


def _analysis(handler, system):
    ops.constraints(*handler)
    ops.numberer("Plain")
    ops.system(*system)
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")


def _ground_motion_on(node):
    ops.timeSeries("Trig", 1, 0.0, 5.0, 0.4, "-factor", 2.0)
    ops.pattern("MultipleSupport", 1)
    ops.groundMotion(1, "Plain", "-accel", 1)
    ops.imposedMotion(node, 1, 1)


def _record(nodes, nsteps, dt):
    hist = [[ops.nodeDisp(n, d)] for (n, d) in nodes]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None
        for h, (n, d) in zip(hist, nodes):
            h.append(ops.nodeDisp(n, d))
    return [np.array(h) for h in hist]


def _run_base_excitation(handler, system):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0); ops.mass(2, _M)
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    _ground_motion_on(1)
    _analysis(handler, system)
    return _record([(2, 1), (1, 1)], 1000, 0.005)


def test_base_excitation_matches_reference():
    proj = _run_base_excitation((PROJ,), ("Diagonal",))
    ref = _run_base_excitation(("Transformation",), ("FullGeneral",))
    assert proj is not None and ref is not None
    assert _rms(proj[0]) > 1e-3
    assert _rms(proj[1]) > 1e-3
    assert _rel(proj[0], ref[0]) < 1e-6, _rel(proj[0], ref[0])
    assert _rel(proj[1], ref[1]) < 1e-6, _rel(proj[1], ref[1])


def _run_settlement(handler, system):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0); ops.mass(2, _M)
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, 1, 0.5, "-const")
    _analysis(handler, system)
    return _record([(2, 1)], 800, 0.005)


def test_constant_settlement_matches_reference():
    proj = _run_settlement((PROJ,), ("Diagonal",))
    ref = _run_settlement(("Transformation",), ("FullGeneral",))
    assert proj is not None and ref is not None
    assert _rms(proj[0]) > 1e-3
    assert _rel(proj[0], ref[0]) < 1e-6, _rel(proj[0], ref[0])


def _equaldof_on_prescribed(handler, system, motion):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)                     # prescribed retained DOF
    ops.node(2, 0.0); ops.mass(2, _M)    # tied to node 1
    ops.node(3, 0.0); ops.mass(3, _M)    # free DOF, spring to node 2
    ops.equalDOF(1, 2, 1)
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.element("zeroLength", 1, 2, 3, "-mat", 1, "-dir", 1)
    if motion == "ground":
        _ground_motion_on(1)
        nsteps = 1000
    else:
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.sp(1, 1, 0.4, "-const")
        nsteps = 800
    _analysis(handler, system)
    return _record([(1, 1), (2, 1), (3, 1)], nsteps, 0.005)


def test_imposed_motion_on_retained_dof_drives_constrained_dof():
    proj = _equaldof_on_prescribed((PROJ,), ("Diagonal",), "ground")
    ref = _equaldof_on_prescribed(("Transformation",), ("FullGeneral",), "ground")
    assert proj is not None and ref is not None
    p1, p2, p3 = proj
    assert _rms(p1) > 1e-3
    assert np.max(np.abs(p2 - p1)) < 1e-10
    assert _rms(p3) > 1e-3
    assert _rel(p3, ref[2]) < 1e-9, _rel(p3, ref[2])


def test_constant_sp_on_retained_dof():
    proj = _equaldof_on_prescribed((PROJ,), ("Diagonal",), "constant")
    ref = _equaldof_on_prescribed(("Transformation",), ("FullGeneral",), "constant")
    assert proj is not None and ref is not None
    p1, p2, p3 = proj
    assert abs(p1[-1] - 0.4) < 1e-10
    assert np.max(np.abs(p2 - p1)) < 1e-10
    assert _rel(p3, ref[2]) < 1e-9, _rel(p3, ref[2])


def test_rigidlink_bar_with_prescribed_ux():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.mass(1, _M, _M)   # ux prescribed, uy free
    ops.node(2, 1.0, 0.0); ops.mass(2, _M, _M)   # tied to node 1 by rigidLink -bar
    ops.node(3, 1.0, 0.0); ops.mass(3, _M, _M)
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.node(4, 0.0, 0.0); ops.fix(4, 1, 1)
    ops.element("zeroLength", 1, 4, 1, "-mat", 1, "-dir", 2)
    ops.element("zeroLength", 2, 2, 3, "-mat", 1, "-dir", 2)
    ops.node(5, 1.0, 0.0); ops.fix(5, 1, 1)
    ops.element("zeroLength", 3, 5, 3, "-mat", 1, "-dir", 2)
    ops.rigidLink("bar", 1, 2)
    _ground_motion_on(1)                         # ux of node 1 only
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(1, 0.0, 5.0)
    _analysis((PROJ,), ("Diagonal",))
    out = _record([(2, 1), (2, 2), (1, 1), (1, 2), (3, 2)], 800, 0.005)
    assert out is not None
    sx, sy, mx, my, u3y = out
    assert _rms(mx) > 1e-3
    assert _rms(my) > 1e-3
    assert _rms(u3y) > 1e-4
    assert np.max(np.abs(sx - mx)) < 1e-10
    assert np.max(np.abs(sy - my)) < 1e-10


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
        return ops.analyze(1, 1e-3) == 0
    except Exception:
        return False


def test_prescribed_constrained_dof_refused():
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, _M)
        ops.node(2, 1.0); ops.mass(2, _M)
        ops.equalDOF(1, 2, 1)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.sp(2, 1, 0.3, "-const")
        _analysis((PROJ,), ("Diagonal",))
    assert not _runs_ok(build)


def test_mixed_free_and_prescribed_retained_refused():
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, _M)
        ops.node(2, 0.0); ops.mass(2, _M)
        ops.node(3, 0.0); ops.mass(3, _M)
        # u3 = 0.5 u1 + 0.5 u2 with u1 prescribed and u2 free
        ops.equationConstraint(3, 1, 1.0, 1, 1, -0.5, 2, 1, -0.5)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.sp(1, 1, 0.3, "-const")
        _analysis((PROJ,), ("Diagonal",))
    assert not _runs_ok(build)


def test_no_free_equation_refused():
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, _M)
        ops.node(2, 0.0); ops.mass(2, _M)
        ops.equalDOF(1, 2, 1)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.sp(1, 1, 0.3, "-const")
        _analysis((PROJ,), ("Diagonal",))
    assert not _runs_ok(build)
