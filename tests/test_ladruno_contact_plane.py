"""Frictionless penalty contact of slave nodes against a rigid analytical
plane (contactPlane), 3D.

Checks:

- Static penetration against the plane x = 0 equals P/kn to 1e-6; against an
  inclined plane n = (0.6, 0.8, 0) with the slave free in x only it equals
  P/(kn*nx^2) to 1e-6.
- Explicit impact (NewmarkExplicit): restitution 0.9 < e <= 1.05 and peak
  penetration v*sqrt(m/kn) = 0.002 to 5 percent along an axis; on the
  inclined plane the slave rebounds with e > 0.9.
- A model whose only finite element is the contact runs twice in one process
  and tears down cleanly.
- A contact tag already used by contact is refused by contactPlane.
- contactPlane -soft is refused.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


KN = 1.0e6


def _anchor_truss(slave, X):
    """Truss of negligible axial stiffness from the slave to a fixed node,
    so the model has a regular element next to the contact."""
    a = 9990 + slave
    ops.node(a, X[0], X[1], X[2] + 1.0)
    ops.fix(a, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 9999, 1.0)
    ops.element("Truss", a, slave, a, 1.0e-9, 9999)


def _static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-12, 30, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _explicit():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("NewmarkExplicit", 0.5)
    ops.algorithm("Linear")
    ops.analysis("Transient")


def _static_press(nx, ny, P=1.0e3):
    """Slave 1e-8 behind a plane through the origin with normal (nx, ny, 0),
    free in x only, loaded -P in x. Returns the x penetration."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, -1.0e-8, 0.0, 0.0)
    ops.fix(1, 0, 1, 1)
    _anchor_truss(1, [-1.0e-8, 0.0, 0.0])
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, nx, ny, 0.0, 0.0, 0.0, 0.0, KN)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, -P, 0.0, 0.0)
    _static()
    assert ops.analyze(1) == 0
    return -(-1.0e-8 + ops.nodeDisp(1, 1))


def test_static_penetration_axis_aligned():
    P = 1.0e3
    pen = _static_press(1.0, 0.0, P)
    assert pen > 0 and abs(pen - P / KN) / (P / KN) < 1e-6, pen


def test_static_penetration_inclined():
    P, nx = 1.0e3, 0.6
    pen = _static_press(nx, 0.8, P)
    expect = P / (KN * nx * nx)
    assert pen > 0 and abs(pen - expect) / expect < 1e-6, (pen, expect)


def _impact(normal, v=2.0, m=1.0, gap0=0.01, nsteps=4000):
    """Mass gap0 outside the plane thrown along -n. Returns
    (ok, restitution, max normal penetration)."""
    nl = math.sqrt(sum(c * c for c in normal))
    n = [c / nl for c in normal]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    X0 = [gap0 * c for c in n]
    ops.node(1, *X0)
    ops.mass(1, m, m, m)
    _anchor_truss(1, X0)
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, n[0], n[1], n[2], 0.0, 0.0, 0.0, KN)
    for d in range(3):
        ops.setNodeVel(1, d + 1, -v * n[d], "-commit")
    _explicit()
    dt = math.sqrt(m / KN)
    max_pen, v_out = 0.0, 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return False, 0.0, max_pen
        gap = sum(n[d] * (X0[d] + ops.nodeDisp(1, d + 1)) for d in range(3))
        vn = sum(n[d] * ops.nodeVel(1, d + 1) for d in range(3))
        if gap < 0:
            max_pen = max(max_pen, -gap)
        if gap > 0.5 * gap0 and vn > 0:
            v_out = vn
    return True, abs(v_out / v), max_pen


def test_explicit_impact_axis_aligned():
    ok, e, max_pen = _impact([1, 0, 0])
    assert ok and 0.9 < e <= 1.05, e
    assert abs(max_pen - 0.002) / 0.002 < 0.05, max_pen


def test_explicit_impact_inclined():
    ok, e, max_pen = _impact([0.6, 0.8, 0.0])
    assert ok and 0.9 < e <= 1.05 and max_pen < 0.02, (e, max_pen)


def test_contact_only_model_runs_twice_and_tears_down():
    """No structural element: the contact is the only finite element."""
    dt = math.sqrt(1.0 / KN)
    for _ in range(2):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.node(1, 0.01, 0.0, 0.0)
        ops.fix(1, 0, 1, 1)
        ops.mass(1, 1.0, 1.0, 1.0)
        ops.contactSurface(1, "-slave", 1)
        ops.contactPlane(1, 1, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, KN)
        ops.setNodeVel(1, 1, -2.0, "-commit")
        _explicit()
        for _ in range(2000):
            assert ops.analyze(1, dt) == 0
        assert 0.01 + ops.nodeDisp(1, 1) > 0.0     # bounced back
    ops.wipe()


def _two_surfaces():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.node(5, 0.0, 0.0, 1.0e-4)
    ops.mass(5, 1.0, 1.0, 1.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)


def test_duplicate_tag_refused_across_contact_and_plane():
    _two_surfaces()
    ops.contact(7, 1, 2, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    with pytest.raises(ops.OpenSeesError):
        ops.contactPlane(7, 2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN)
    ops.contactPlane(8, 2, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, KN)


def test_soft_option_refused():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 1.0e-4)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(1, "-slave", 1)
    try:
        ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN, "-soft", 0.1)
    except Exception:
        return
    _explicit()
    assert ops.analyze(1, 1.0e-4) != 0
