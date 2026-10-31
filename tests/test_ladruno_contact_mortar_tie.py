"""Mortar mesh tying (contact ... -mortar -tie): a permanent bond across a non-matching interface.

Checks:

- Patch test: a stdBrick column split at z = 1 into one brick (master face) and two bricks
  (slave face, non-matching, shared edge nodes) reproduces the monolithic bar u = P z / (E A)
  to 1e-6 under tension (penalty epsTie = 1e8), and compression mirrors tension (1e-9).
- Static Newton converges on the symmetric tie tangent and the tie adds no equation
  (16 free DOFs).
- Listing the slave facets in either order gives the same solution (1e-12).
- A held-load augmentation loop drives ladrunoMortarTieResidual below 1e-11 for
  epsTie = 1e6 and 1e8.
- Proximity guard: a slave facet coincident with one master facet and parallel to a second
  master facet one facet size away (the far side of a box) is tied to the near facet only:
  the bond stiffness is epsTie A (1e-10), not twice that.
- Refusals: -tie without -mortar.

Every test fails on the base with an unknown command.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E, NU = 2.0e4, 0.0
PTOT = 20.0


def _analyze_augmented(maxAug=20, augTol=1.0e-9, query=None):
    """Held-load Uzawa loop (see tests/test_ladruno_contact_mortar.py)."""
    if query is None:
        query = ops.ladrunoMortarPenetration
    ops.integrator("LoadControl", 0.0)
    try:
        ops.ladrunoBeginAugment()
        pen = query()
        for k in range(1, int(maxAug) + 1):
            rc = ops.analyze(1)
            if rc != 0:
                return rc, k, pen
            pen = query()
            if pen < augTol:
                return 0, k, pen
        return 1, int(maxAug), pen
    finally:
        ops.ladrunoEndAugment()


def _static_setup():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _split_column(epsTie, swap_facets=False, tension=True):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    bot = [(1, 0, 0, 0), (2, 1, 0, 0), (3, 1, 1, 0), (4, 0, 1, 0),
           (5, 0, 0, 1), (6, 1, 0, 1), (7, 1, 1, 1), (8, 0, 1, 1)]
    for t, x, y, z in bot:
        ops.node(t, float(x), float(y), float(z))
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    top = [(9, 0, 0, 1), (10, 0.5, 0, 1), (11, 1, 0, 1), (12, 0, 1, 1), (13, 0.5, 1, 1),
           (14, 1, 1, 1), (15, 0, 0, 2), (16, 0.5, 0, 2), (17, 1, 0, 2), (18, 0, 1, 2),
           (19, 0.5, 1, 2), (20, 1, 1, 2)]
    for t, x, y, z in top:
        ops.node(t, float(x), float(y), float(z))
    ops.element("stdBrick", 2, 9, 10, 13, 12, 15, 16, 19, 18, 1)
    ops.element("stdBrick", 3, 10, 11, 14, 13, 16, 17, 20, 19, 1)
    for t in range(5, 21):
        ops.fix(t, 1, 1, 0)                      # uniaxial
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    facetA, facetB = [9, 10, 13, 12], [10, 11, 14, 13]
    slave = (facetB + facetA) if swap_facets else (facetA + facetB)
    ops.contactSurface(2, "-slave-segments", 4, *slave)
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", epsTie, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    sgn = 1.0 if tension else -1.0
    for t in (15, 17, 18, 20):
        ops.load(t, 0.0, 0.0, sgn * PTOT / 8.0)
    for t in (16, 19):
        ops.load(t, 0.0, 0.0, sgn * PTOT / 4.0)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    if ops.analyze(1) != 0:
        return None
    return dict(sysN=ops.systemSize(), u_slave=ops.nodeDisp(9, 3),
                u_master=ops.nodeDisp(5, 3), u_tip=ops.nodeDisp(15, 3))


def test_patch_test_split_matches_monolithic():
    r = _split_column(1.0e8)
    assert r is not None
    assert abs(r["u_slave"] - r["u_master"]) < 1.0e-6
    assert abs(r["u_slave"] - PTOT / E) < 1.0e-6
    assert abs(r["u_tip"] - 2.0 * PTOT / E) < 1.0e-6


def test_tie_adds_no_equation():
    r = _split_column(1.0e7)
    assert r is not None and r["sysN"] == 16


def test_facet_order_independent():
    fwd = _split_column(1.0e7, swap_facets=False)
    rev = _split_column(1.0e7, swap_facets=True)
    assert fwd is not None and rev is not None
    assert abs(fwd["u_tip"] - rev["u_tip"]) < 1.0e-12
    assert abs(fwd["u_slave"] - rev["u_slave"]) < 1.0e-12


def test_compression_mirrors_tension():
    rt = _split_column(1.0e8, tension=True)
    rc = _split_column(1.0e8, tension=False)
    assert rt is not None and rc is not None
    assert abs(rc["u_tip"] + rt["u_tip"]) < 1.0e-9


def _unit_square(tags, z):
    for t, x, y in zip(tags, (0.0, 1.0, 1.0, 0.0), (0.0, 0.0, 1.0, 1.0)):
        ops.node(t, x, y, z)


@pytest.mark.parametrize("epsTie", [1.0e6, 1.0e8])
def test_held_load_augmentation_closes_the_tie(epsTie):
    augTol = 1.0e-11
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    _unit_square([5, 6, 7, 8], 1.0e-4)          # a gap that a contact would open
    for t in (5, 6, 7, 8):
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", epsTie, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, 400.0 / 4.0)       # tension: only the tie holds it
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    pen = ops.ladrunoMortarTieResidual()
    status, _, res = _analyze_augmented(maxAug=80, augTol=augTol,
                                        query=ops.ladrunoMortarTieResidual)
    assert status == 0 and res < augTol
    assert pen > res


def _box_tie(extra=()):
    """Master: two parallel fixed facets, z = 0 and z = 1 (bottom and top of a box).
    Slave: one facet coincident with the top facet, free in z, pulled up by P = 1.
    Returns the slave z-displacement."""
    epsTie = 1.0e4
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)
    _unit_square([5, 6, 7, 8], 1.0)
    for t in range(1, 9):
        ops.fix(t, 1, 1, 1)
    _unit_square([9, 10, 11, 12], 1.0)
    for t in (9, 10, 11, 12):
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave-segments", 4, 9, 10, 11, 12)
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", epsTie, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (9, 10, 11, 12):
        ops.load(t, 0.0, 0.0, 0.25)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    return ops.nodeDisp(9, 3), 1.0 / epsTie


def test_tie_is_not_bonded_to_the_far_side_of_a_closed_surface():
    u, u_one_bond = _box_tie()
    assert u == pytest.approx(u_one_bond, rel=1e-10)


def test_tie_requires_mortar():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(1 + i, float(x), float(y), 0.0)
        ops.node(5 + i, float(x), float(y), 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    with pytest.raises(ops.OpenSeesError):
        ops.contact(2, 1, 2, "-tie", "-epsTie", 1.0e5)
