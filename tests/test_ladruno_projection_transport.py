"""constraints LadrunoProjection: general constraint matrices (rigidLink -beam,
rigidDiaphragm, equationConstraint) under CentralDifferenceLadruno.

With constraints Transformation and system Diagonal, the off-diagonal terms of
the condensed mass T^T M T (the translation-rotation coupling of an offset
mass) are dropped at assembly. The projection keeps every DOF with its own
diagonal mass and enforces the constraint on the acceleration, so the
diagonal recipe gives the coupled answer.

Checks:
- an offset mass tied by rigidLink -beam and driven by a moment matches
  constraints Transformation with system FullGeneral to 1e-6 (relative RMS),
  stays within 0.1 of implicit Newmark, holds the tie u_s,y = u_m,y + d
  theta_m to 1e-12, and matches the closed-form modal response of the
  condensed 2-DOF system to 2e-2;
- the same model in free vibration (smaller moment, longer record) matches
  the FullGeneral reference to 1e-6;
- a 3D rigidDiaphragm with four corner masses under lateral load and torque
  matches the reference to 1e-6 and stays rigid to 1e-10;
- an equationConstraint u_2 = u_1 is enforced to 1e-12 with the combined
  acceleration F/(m1 + m2);
- a tied undamped oscillator conserves energy to 1e-3 over 2000 steps;
- refusals: a tied rotation with zero rotational mass, a diaphragm-tied DOF
  that is also fixed, a constraint to a node that is not in the domain.
"""

import math

import numpy as np

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"


def _rms(a):
    return float(np.sqrt(np.mean(a * a)))


def _rel(a, b):
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    d = _rms(a) + _rms(b)
    return 0.0 if d == 0.0 else _rms(a - b) / (0.5 * d)


def _analysis(handler, system, integrator):
    ops.constraints(*handler)
    ops.numberer("Plain")
    ops.system(*system)
    ops.test("NormDispIncr", 1e-10, 50)
    if integrator == "cdl":
        ops.algorithm("Linear")
        ops.integrator(CDL)
    else:
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


# 2D master with an offset slave mass tied by rigidLink -beam
_M0, _IRZ, _MS, _MSRZ, _D = 1.0, 1.0, 1.0, 0.2, 2.0
_KX = _KY = _KR = 100.0


def _run_rigidlink_beam(handler, system, integrator, moment=5.0, nsteps=400, dt=0.01):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(10, 0.0, 0.0); ops.fix(10, 1, 1, 1)
    ops.node(1, 0.0, 0.0); ops.mass(1, _M0, _M0, _IRZ)
    ops.node(2, _D, 0.0); ops.mass(2, _MS, _MS, _MSRZ)   # the tied rotation needs mass
    for i, k in ((1, _KX), (2, _KY), (3, _KR)):
        ops.uniaxialMaterial("Elastic", i, k)
    ops.element("zeroLength", 1, 10, 1, "-mat", 1, 2, 3, "-dir", 1, 2, 3)
    ops.rigidLink("beam", 1, 2)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, moment)
    _analysis(handler, system, integrator)
    uy = [ops.nodeDisp(1, 2)]
    gap = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None, None
        uy.append(ops.nodeDisp(1, 2))
        gap = max(gap, abs(ops.nodeDisp(2, 2) - (ops.nodeDisp(1, 2) + _D * ops.nodeDisp(1, 3))))
    return np.array(uy), gap


def _analytic_uy(t, moment=5.0):
    """Modal response of the condensed (uy, rz) master system, assembled by
    hand: M = [[M0 + MS, MS d], [MS d, IRZ + MS d^2 + MSRZ]], K = diag(KY, KR),
    constant moment applied from rest."""
    t = np.asarray(t, dtype=float)
    M22 = np.array([[_M0 + _MS, _MS * _D],
                    [_MS * _D, _IRZ + _MS * _D * _D + _MSRZ]])
    K22 = np.array([[_KY, 0.0], [0.0, _KR]])
    u_s = np.linalg.solve(K22, np.array([0.0, moment]))
    w2, phi = np.linalg.eig(np.linalg.solve(M22, K22))
    w = np.sqrt(w2.real)
    phi = phi.real
    for i in range(phi.shape[1]):
        phi[:, i] /= math.sqrt(phi[:, i] @ M22 @ phi[:, i])
    uy = np.full_like(t, u_s[0])
    for i in range(len(w)):
        uy += phi[0, i] * (-(phi[:, i] @ M22 @ u_s)) * np.cos(w[i] * t)
    return uy


def test_rigidlink_beam_matches_dense_reference():
    proj, gap = _run_rigidlink_beam((PROJ,), ("Diagonal",), "cdl")
    ref, _ = _run_rigidlink_beam(("Transformation",), ("FullGeneral",), "cdl")
    imp, _ = _run_rigidlink_beam(("Transformation",), ("FullGeneral",), "imp")
    assert proj is not None and ref is not None and imp is not None
    assert _rms(proj) > 1e-3
    assert _rel(proj, ref) < 1e-6, _rel(proj, ref)
    assert _rel(proj, imp) < 0.1, _rel(proj, imp)
    assert gap < 1e-12, gap


def test_rigidlink_beam_matches_closed_form():
    dt, nsteps = 0.01, 400
    proj, _ = _run_rigidlink_beam((PROJ,), ("Diagonal",), "cdl", moment=5.0, nsteps=nsteps, dt=dt)
    assert proj is not None
    anchor = _analytic_uy(np.arange(nsteps + 1) * dt, moment=5.0)
    assert _rel(proj, anchor) < 2e-2, _rel(proj, anchor)


def test_rigidlink_beam_free_vibration():
    proj, _ = _run_rigidlink_beam((PROJ,), ("Diagonal",), "cdl", moment=2.0, nsteps=600)
    ref, _ = _run_rigidlink_beam(("Transformation",), ("FullGeneral",), "cdl",
                                 moment=2.0, nsteps=600)
    assert proj is not None and ref is not None
    assert _rel(proj, ref) < 1e-6, _rel(proj, ref)


def _run_diaphragm(handler, system, nsteps=300, dt=0.005):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(100, 0.0, 0.0, 0.0)
    ops.mass(100, 4.0, 4.0, 1e-9, 1e-9, 1e-9, 8.0)
    ops.fix(100, 0, 0, 1, 1, 1, 0)
    gnd = 10
    for tag, x, y in ((1, 1.0, 1.0), (2, -1.0, 1.0), (3, -1.0, -1.0), (4, 1.0, -1.0)):
        ops.node(tag, x, y, 0.0)
        # the diaphragm ties rz of the corners: it needs a (small) rotational mass
        ops.mass(tag, 1.0, 1.0, 1e-9, 1e-9, 1e-9, 0.01)
        ops.fix(tag, 0, 0, 1, 1, 1, 0)
        ops.node(gnd, x, y, -1.0); ops.fix(gnd, 1, 1, 1, 1, 1, 1)
        ops.uniaxialMaterial("Elastic", tag, 50.0)
        ops.element("zeroLength", tag, gnd, tag, "-mat", tag, tag, "-dir", 1, 2)
        gnd += 1
    ops.rigidDiaphragm(3, 100, 1, 2, 3, 4)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(100, 3.0, 0.0, 0.0, 0.0, 0.0, 2.0)
    _analysis(handler, system, "cdl")
    ux = [ops.nodeDisp(100, 1)]
    gap = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None, None
        ux.append(ops.nodeDisp(100, 1))
        # corner 1 (x = 1, y = 1): ux = ux_m - rz_m (y1 - y_m)
        gap = max(gap, abs(ops.nodeDisp(1, 1) - (ops.nodeDisp(100, 1) - ops.nodeDisp(100, 6))))
    return np.array(ux), gap


def test_rigid_diaphragm_matches_reference_and_stays_rigid():
    proj, gap = _run_diaphragm((PROJ,), ("Diagonal",))
    ref, _ = _run_diaphragm(("Transformation",), ("FullGeneral",))
    assert proj is not None and ref is not None
    assert _rms(proj) > 1e-3
    assert _rel(proj, ref) < 1e-6, _rel(proj, ref)
    assert gap < 1e-10, gap


def test_equation_constraint_enforced():
    m1, m2, F = 2.0, 3.0, 10.0
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.equationConstraint(2, 1, 1.0, 1, 1, -1.0)   # u_2 - u_1 = 0
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, F)
    ops.constraints(PROJ)
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL); ops.analysis("Transient")

    tie_err = 0.0
    for _ in range(200):
        assert ops.analyze(1, 1e-3) == 0
        tie_err = max(tie_err, abs(ops.nodeDisp(1, 1) - ops.nodeDisp(2, 1)))
    assert tie_err < 1e-12
    assert abs(ops.nodeAccel(1, 1) - F / (m1 + m2)) < 1e-9


def test_energy_conserved():
    m1, m2, k, d0 = 2.0, 3.0, 100.0, 1.0
    dt = 0.02 * (2.0 / math.sqrt(k / (m1 + m2)))
    nsteps = 2000

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
    ops.equalDOF(1, 2, 1)
    ops.setNodeDisp(1, 1, d0, "-commit")
    ops.setNodeDisp(2, 1, d0, "-commit")
    ops.constraints(PROJ)
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL); ops.analysis("Transient")

    E0 = 0.5 * k * d0 * d0
    energies = []
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        u = ops.nodeDisp(1, 1)
        v1, v2 = ops.nodeVel(1, 1), ops.nodeVel(2, 1)
        energies.append(0.5 * m1 * v1 * v1 + 0.5 * m2 * v2 * v2 + 0.5 * k * u * u)
    energies = np.array(energies)
    assert len(energies) == nsteps
    assert np.max(np.abs(energies - E0)) / E0 < 1e-3


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


def test_massless_rotational_tie_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(10, 0.0, 0.0); ops.fix(10, 1, 1, 1)
        ops.node(1, 0.0, 0.0); ops.mass(1, 1.0, 1.0, 1.0)
        ops.node(2, 2.0, 0.0); ops.mass(2, 1.0, 1.0, 0.0)    # tied rz without mass
        for i in (1, 2, 3):
            ops.uniaxialMaterial("Elastic", i, 100.0)
        ops.element("zeroLength", 1, 10, 1, "-mat", 1, 2, 3, "-dir", 1, 2, 3)
        ops.rigidLink("beam", 1, 2)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0., 0., 5.)
        _analysis((PROJ,), ("Diagonal",), "cdl")
    assert not _runs_ok(build)


def test_fixed_diaphragm_dof_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(100, 0., 0., 0.); ops.mass(100, 4., 4., 1e-9, 1e-9, 1e-9, 8.)
        ops.fix(100, 0, 0, 1, 1, 1, 0)
        ops.node(1, 1., 1., 0.); ops.mass(1, 1., 1., 1e-9, 1e-9, 1e-9, 0.01)
        ops.fix(1, 0, 0, 1, 1, 1, 1)              # rz fixed but tied by the diaphragm
        ops.node(10, 1., 1., -1.); ops.fix(10, 1, 1, 1, 1, 1, 1)
        ops.uniaxialMaterial("Elastic", 1, 50.)
        ops.element("zeroLength", 1, 10, 1, "-mat", 1, 1, "-dir", 1, 2)
        ops.rigidDiaphragm(3, 100, 1)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.load(100, 3., 0., 0., 0., 0., 2.)
        _analysis((PROJ,), ("Diagonal",), "cdl")
    assert not _runs_ok(build)


def test_constraint_to_missing_node_refused():
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, 1.0)
        ops.node(2, 1.0); ops.mass(2, 1.0)
        ops.equalDOF(1, 2, 1)
        ops.equalDOF(1, 99, 1)      # node 99 does not exist
        ops.constraints(PROJ)
        ops.numberer("Plain"); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")
    # either equalDOF or the handler refuses; the analysis must not run
    assert not _runs_ok(build)
