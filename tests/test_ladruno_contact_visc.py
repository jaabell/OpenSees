"""Viscous normal contact stabilization (-visc muc): p_visc = muc * gap rate,
active only while in contact, on contactPlane and on node-to-segment contact.

Checks (Kelvin-Voigt contact of a unit mass, zeta = muc / (2 sqrt(kn m))):

- Explicit impact on a plane: e > 0.95 with muc = 0, e decreasing in muc,
  e < 0.7 at muc = 500; e matches exp(-zeta pi / sqrt(1 - zeta^2)) to 0.06
  for muc = 200 and 500.
- Static press with muc = 5e3 is bitwise identical to muc = 0 and equals
  -P/kn to 1e-9.
- A static solve with a stale committed velocity is bitwise identical with
  and without -visc, and the disablement is reported.
- Implicit Newmark impact converges with -visc and rebounds slower than
  without.
- Node-to-segment impact on a fixed facet: e > 0.9 with muc = 0, e < 0.7 with
  muc = 500.
"""
import math

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


KN = 1.0e6
M = 1.0


def _zeta(muc):
    return muc / (2.0 * math.sqrt(KN * M))


def _e_kelvin_voigt(muc):
    z = _zeta(muc)
    return math.exp(-z * math.pi / math.sqrt(1.0 - z * z)) if z < 1.0 else 0.0


def _explicit():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("NewmarkExplicit", 0.5)
    ops.algorithm("Linear")
    ops.analysis("Transient")


def _plane_impact(muc, v=2.0, gap0=0.01):
    """Unit mass gap0 above the plane z = 0, thrown down at v, explicit.
    Returns (ok, restitution)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, gap0)
    ops.mass(1, M, M, M)
    ops.fix(1, 1, 1, 0)
    ops.node(901, 0.0, 0.0, gap0 + 1.0)
    ops.fix(901, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 900, 1.0e-6)
    ops.element("truss", 901, 901, 1, 1.0, 900)
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN, "-visc", muc)
    ops.setNodeVel(1, 3, -v, "-commit")
    T = 2.0 * math.pi * math.sqrt(M / KN)
    dt = T / 60.0
    _explicit()
    v_out = 0.0
    for _ in range(int(40.0 * T / dt) + 200):
        if ops.analyze(1, dt) != 0:
            return False, 0.0
        z = gap0 + ops.nodeDisp(1, 3)
        vz = ops.nodeVel(1, 3)
        if z > 0.5 * gap0 and vz > 0.0:
            v_out = vz
    return True, abs(v_out / v)


def test_explicit_plane_impact_dissipation_monotone():
    ok0, e0 = _plane_impact(0.0)
    ok1, e1 = _plane_impact(200.0)
    ok2, e2 = _plane_impact(500.0)
    assert ok0 and ok1 and ok2
    assert e0 > 0.95, e0
    assert e1 < e0 and e2 < e1, (e0, e1, e2)
    assert e2 < 0.7, e2


def test_explicit_plane_restitution_kelvin_voigt():
    for muc in (200.0, 500.0):
        ok, e = _plane_impact(muc)
        assert ok
        assert abs(e - _e_kelvin_voigt(muc)) < 0.06, (muc, e, _e_kelvin_voigt(muc))


def _static_plane_press(muc, stale_vz=0.0):
    """Node on the plane z = 0, x and y fixed, a 0.1 z spring to ground,
    pressed by -1e3. Returns the converged z displacement."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 0)
    ops.node(90, 0.0, 0.0, 0.0)
    ops.fix(90, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 90, 0.1)
    ops.element("zeroLength", 90, 90, 1, "-mat", 90, "-dir", 3)
    if stale_vz != 0.0:
        ops.setNodeVel(1, 3, stale_vz, "-commit")
    ops.contactSurface(1, "-slave", 1)
    args = [1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN]
    if muc > 0.0:
        args += ["-visc", muc]
    ops.contactPlane(*args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return ops.nodeDisp(1, 3)


def test_static_unchanged_by_visc():
    z_off = _static_plane_press(0.0)
    z_on = _static_plane_press(5.0e3)
    assert z_off == z_on, (z_off, z_on)
    assert abs(z_off - (-1.0e3 / KN)) < 1.0e-9, z_off


def test_static_with_stale_velocity_ignores_visc(capfd):
    z_off = _static_plane_press(0.0, stale_vz=-5.0)
    z_on = _static_plane_press(50.0, stale_vz=-5.0)
    err = capfd.readouterr().err
    assert z_on == z_off, (z_on, z_off)
    assert z_off > -0.1, z_off             # contact engaged, not the spring alone
    assert "-visc" in err and "static integrator" in err, err[-300:]


def test_implicit_newmark_converges_and_dissipates():
    def impact(muc):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        gap0 = 0.01
        ops.node(1, 0.0, 0.0, gap0)
        ops.mass(1, M, M, M)
        ops.fix(1, 1, 1, 0)
        ops.contactSurface(1, "-slave", 1)
        ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN, "-visc", muc)
        ops.setNodeVel(1, 3, -2.0, "-commit")
        T = 2.0 * math.pi * math.sqrt(M / KN)
        dt = T / 40.0
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.test("NormDispIncr", 1.0e-10, 50, 0)
        ops.algorithm("Newton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        v_out = 0.0
        for _ in range(int(60.0 * T / dt)):
            if ops.analyze(1, dt) != 0:
                return None
            z = gap0 + ops.nodeDisp(1, 3)
            vz = ops.nodeVel(1, 3)
            if z > 0.5 * gap0 and vz > 0.0:
                v_out = vz
        return abs(v_out / 2.0)

    e_off, e_on = impact(0.0), impact(500.0)
    assert e_off is not None and e_on is not None
    assert e_on < e_off, (e_on, e_off)


def test_segment_impact_dissipates():
    def impact(muc):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
            ops.node(t, float(x), float(y), 0.0)
            ops.fix(t, 1, 1, 1)
        gap0 = 0.01
        ops.node(5, 0.0, 0.0, gap0)
        ops.mass(5, M, M, M)
        ops.fix(5, 1, 1, 0)
        ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
        ops.contactSurface(2, "-slave", 5)
        ops.contact(1, 1, 2, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0, "-visc", muc)
        ops.setNodeVel(5, 3, -2.0, "-commit")
        T = 2.0 * math.pi * math.sqrt(M / KN)
        dt = T / 60.0
        _explicit()
        v_out = 0.0
        for _ in range(int(40.0 * T / dt) + 200):
            if ops.analyze(1, dt) != 0:
                return None
            z = gap0 + ops.nodeDisp(5, 3)
            vz = ops.nodeVel(5, 3)
            if z > 0.5 * gap0 and vz > 0.0:
                v_out = vz
        return abs(v_out / 2.0)

    e0, e1 = impact(0.0), impact(500.0)
    assert e0 is not None and e1 is not None
    assert e0 > 0.9, e0
    assert e1 < 0.7, e1
