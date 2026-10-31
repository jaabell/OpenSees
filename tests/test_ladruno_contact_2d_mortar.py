"""2D mortar contact: penalty, thickness, augmentation, mesh tying and friction.

A slave segment (or a chain of slave segments) seeded 1e-8 below a fixed master line of length L,
-outward 0 1 (the interface is flush, so the orientation is given).

Checks:

- Penalty: a matched pair pressed by P settles at the uniform penetration P/(epsN L h) (1e-10) for
  h = 1 and with -thickness 2 (half the penetration); a non-matching interface of three slave
  segments on two master segments, loaded by the consistent nodal loads of a uniform pressure,
  settles at the uniform penetration P/(epsN L) (1e-8): the mortar lane passes the patch test.
- Augmentation: the held-load loop (ladrunoBeginAugment, analyze at a zero load increment,
  ladrunoEndAugment) drives ladrunoMortarPenetration below 1e-9 for epsN = 1e4 and 1e6.
- Mesh tie (on a master longer than the slave, so the overlap stays the whole slave): loads
  (Qx, Qy) on a tied slave segment give the displacements Q/(epsTie L) in both
  directions (1e-10).
- Friction (same longer master): with springs k in x on both slave nodes, a load above the cone slips to
  (Q - mu P)/(2 k) (1e-9) and a load below it sticks at a small slip.

Every test fails on the base with an unknown command.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


EPS = 1.0e6
SEED = 1.0e-8
P = 100.0
L = 2.0


def _analyze_augmented(maxAug=40, augTol=1.0e-9):
    ops.integrator("LoadControl", 0.0)
    try:
        ops.ladrunoBeginAugment()
        pen = ops.ladrunoMortarPenetration()
        for k in range(1, int(maxAug) + 1):
            rc = ops.analyze(1)
            if rc != 0:
                return rc, k, pen
            pen = ops.ladrunoMortarPenetration()
            if pen < augTol:
                return 0, k, pen
        return 1, int(maxAug), pen
    finally:
        ops.ladrunoEndAugment()


def _static(system="FullGeneral"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _pair(extra=(), fix_x=True, loads=(0.0, -P / 2.0), master_half=1.0):
    """Master 1-2 on y = 0 (fixed, from -master_half to master_half), slave segment 5-6
    (from -1 to 1) seeded below it."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, x in ((1, -master_half), (2, master_half)):
        ops.node(t, x, 0.0)
        ops.fix(t, 1, 1)
    for t, x in ((5, -1.0), (6, 1.0)):
        ops.node(t, x, -SEED)
        if fix_x:
            ops.fix(t, 1, 0)
    ops.contactSurface(10, "-master", 2, 1, 2)
    ops.contactSurface(20, "-slave-segments", 2, 5, 6)
    ops.contact(1, 10, 20, "-mortar", "-epsN", EPS, "-outward", 0.0, 1.0, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, *loads)
    _static()


def test_penalty_penetration():
    _pair()
    assert ops.analyze(1) == 0
    g = P / (EPS * L)
    for t in (5, 6):
        assert abs(ops.nodeDisp(t, 2) - SEED + g) < 1.0e-10 * g


def test_thickness_scales_the_penalty():
    _pair(extra=("-thickness", 2.0))
    assert ops.analyze(1) == 0
    g = P / (EPS * L * 2.0)
    for t in (5, 6):
        assert abs(ops.nodeDisp(t, 2) - SEED + g) < 1.0e-10 * g


def test_nonmatching_interface_transmits_the_load():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    mx = [-1.0, 0.0, 1.0]
    for i, x in enumerate(mx):
        ops.node(1 + i, x, 0.0)
        ops.fix(1 + i, 1, 1)
    sx = [-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0]
    for i, x in enumerate(sx):
        ops.node(5 + i, x, -SEED)
        ops.fix(5 + i, 1, 0)
    ops.contactSurface(10, "-master", 2, 1, 2, 2, 3)
    ops.contactSurface(20, "-slave-segments", 2, 5, 6, 6, 7, 7, 8)
    ops.contact(1, 10, 20, "-mortar", "-epsN", EPS, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # consistent nodal loads of a uniform pressure P / L
    le = L / 3.0
    for t, w in ((5, le / 2.0), (6, le), (7, le), (8, le / 2.0)):
        ops.load(t, 0.0, -P * w / L)
    _static()
    assert ops.analyze(1) == 0
    g = P / (EPS * L)
    for t in (5, 6, 7, 8):
        assert abs(ops.nodeDisp(t, 2) - SEED + g) < 1.0e-8 * g


@pytest.mark.parametrize("epsN", [1.0e4, 1.0e6])
def test_augmentation_drives_the_penetration_to_zero(epsN):
    global EPS
    saved, EPS = EPS, epsN
    try:
        _pair()
    finally:
        EPS = saved
    assert ops.analyze(1) == 0
    assert ops.ladrunoMortarPenetration() > 1.0e-9
    status, _, pen = _analyze_augmented()
    assert status == 0 and pen < 1.0e-9, (status, pen)


def test_tie_binds_both_components():
    qx, qy = 3.0, -5.0
    # a longer master keeps the displaced slave fully covered, so the overlap
    # (and with it D) does not change with the tangential displacement
    _pair(extra=("-tie",), fix_x=False, loads=(qx / 2.0, qy / 2.0), master_half=2.0)
    assert ops.analyze(1) == 0
    for t in (5, 6):
        assert abs(ops.nodeDisp(t, 1) - qx / (EPS * L)) < 1.0e-10 * abs(qx / (EPS * L))
        assert abs(ops.nodeDisp(t, 2) - qy / (EPS * L)) < 1.0e-10 * abs(qy / (EPS * L))


MU = 0.3
EPS_T = 1.0e4
K_X = 50.0


def _friction(Q):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, x in ((1, -2.0), (2, 2.0)):          # covers the slave after it slides
        ops.node(t, x, 0.0)
        ops.fix(t, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, K_X)
    for t, x in ((5, -1.0), (6, 1.0)):
        ops.node(t, x, -SEED)
        ops.node(t + 10, x, -SEED)
        ops.fix(t + 10, 1, 1)
        ops.element("zeroLength", t, t + 10, t, "-mat", 1, "-dir", 1)
    ops.contactSurface(10, "-master", 2, 1, 2)
    ops.contactSurface(20, "-slave-segments", 2, 5, 6)
    ops.contact(1, 10, 20, "-mortar", "-epsN", EPS, "-mu", MU, "-epsT", EPS_T,
                "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, Q / 2.0, -P / 2.0)
    _static()
    assert ops.analyze(1) == 0
    return [ops.nodeDisp(t, 1) for t in (5, 6)]


def test_friction_slip():
    Q = 50.0                                  # > mu P = 30
    for u in _friction(Q):
        assert abs(u - (Q - MU * P) / (2.0 * K_X)) < 1.0e-9


def test_friction_stick():
    Q = 10.0                                  # < mu P = 30
    for u in _friction(Q):
        assert 0.0 < u < Q / (EPS_T * L) * 1.01
