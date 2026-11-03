"""Edge-to-edge contact (-edgeedge on a mortar contact): penalty, friction, augmented Lagrangian.

Two crossed bars: a slave facet in the x-z plane whose bottom edge runs along x, over a fixed
master facet in the y-z plane whose top edge runs along y. Only the two crossing edges touch: no
slave node projects onto the master face and the perpendicular facets have no clip overlap, so
neither the node-to-segment nor the mortar lane transmits a force. The common normal of the edges
is +z.

Checks:

- Penalty: the slave bottom edge, pressed down by P on two soft springs, settles at the penetration
  P/(epsN + 2k) (1e-7), symmetric (1e-9), with static Newton converging.
- The same model with node-to-segment contact of the slave nodes transmits no force: the edge sinks
  to the spring depth P/(2k) (1e-6).
- Oblique facets (normals about 60 degrees apart, the mortar clip not empty): the mortar lane stands
  down, so the penetration is the edge-to-edge one P/(epsN + 2k) (1e-6), not a superposition.
- Bars far apart: a Newmark run with -edgeedge is bitwise equal to the run without it.
- Friction: an explicit raking edge (NewmarkExplicit) is decelerated by Coulomb friction and coasts
  at its initial velocity without it; a sub-cone static load sticks with a small slip; the symmetric
  tangent gives the same answer on ProfileSPD and FullGeneral (1e-10); -edgeMu 0 is bitwise equal
  to no friction options.
- Augmented Lagrangian (-edgeAlm): a held-load augmentation drives ladrunoEdgePenetration below
  1e-9 for epsN = 1e5 and 1e7, while the penalty penetration of the first solve depends on epsN
  (ratio > 10); a reversed load lifts the bar off and the query returns 0; the first solve with
  -edgeAlm is bitwise equal to the one without.
- Refusals: -edgeedge without -mortar, edge options without -edgeedge, -edgeedge with -tie,
  -edgeTauMax without -edgeMu or -edgeCohesion, and -edgeSoft.

Every test fails on the base with an unknown command.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


EPS = 1.0e6      # edge penalty epsN
KSPRING = 500.0  # soft spring on each slave bottom node
PTOT = 400.0     # total downward load on the slave bottom edge
ZC = -1.0e-4     # slave bottom edge seeded slightly below the master top edge


def _analyze_augmented(maxAug=20, augTol=1.0e-9, query=None):
    """Held-load Uzawa loop: re-solve and commit at a zero load increment until the
    measure drops below augTol. Returns (status, augmentations, measure)."""
    if query is None:
        query = ops.ladrunoEdgePenetration
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


def _build_bars(xoff=0.0, zc=ZC):
    """Master (fixed) facet in the y-z plane, top edge 1-2 along y at z = 0. Slave facet in the
    x-z plane (shifted by xoff), bottom edge 5-6 along x at z = zc."""
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    for t, x, y, z in [(5, -1.0 + xoff, 0.0, zc), (6, 1.0 + xoff, 0.0, zc),
                       (7, 1.0 + xoff, 0.0, zc + 0.3), (8, -1.0 + xoff, 0.0, zc + 0.3)]:
        ops.node(t, x, y, z)


def _z_springs(zg=ZC, k=KSPRING):
    """z-springs on the slave bottom nodes 5, 6: the edge pair shares one gap, so the
    differential mode of the two nodes needs a stiffness of its own."""
    ops.uniaxialMaterial("Elastic", 1, k)
    for gnd, sn, x in [(305, 5, -1.0), (306, 6, 1.0)]:
        ops.node(gnd, x, 0.0, zg)
        ops.fix(gnd, 1, 1, 1)
        ops.element("zeroLength", 900 + sn, gnd, sn, "-mat", 1, "-dir", 3)


def _surfaces():
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)


def _load_down():
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 0.0, -PTOT / 2.0)


def _static_setup(system="FullGeneral", tol=1.0e-11, maxit=50):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", tol, maxit, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _pressed_bars():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)
    _z_springs()


# ----------------------------------------------------------------- penalty
def test_crossed_bars_restrained():
    _pressed_bars()
    _surfaces()
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge",
                "-edgeBand", 0.15, "-outward", 0.0, 0.0, 1.0)
    _load_down()
    _static_setup()
    assert ops.analyze(1) == 0
    delta = PTOT / (EPS + 2.0 * KSPRING)
    z5 = ZC + ops.nodeDisp(5, 3)
    z6 = ZC + ops.nodeDisp(6, 3)
    assert abs(z5 - z6) < 1.0e-9
    assert abs(z5 + delta) < 1.0e-7
    assert abs(z5) < 1.0e-3


def test_node_to_segment_passes_through():
    _pressed_bars()
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(3, "-slave", 5, 6)
    ops.contact(1, 1, 3, EPS, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    _load_down()
    _static_setup()
    assert ops.analyze(1) == 0
    z5 = ZC + ops.nodeDisp(5, 3)
    assert abs(z5 - (ZC - PTOT / (2.0 * KSPRING))) < 1.0e-6
    assert abs(z5) > 100.0 * PTOT / (EPS + 2.0 * KSPRING)


def test_oblique_pair_is_not_counted_twice():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # master top edge along y at z = 0, facet tilted 45 degrees; slave bottom edge along x at
    # z = ZC, facet tilted 45 degrees the other way: normals about 60 degrees apart
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.3, 1.0, -0.3), (4, 0.3, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    for t, x, y, z in [(5, -1.0, 0.0, ZC), (6, 1.0, 0.0, ZC),
                       (7, 1.0, 0.3, ZC + 0.3), (8, -1.0, 0.3, ZC + 0.3)]:
        ops.node(t, x, y, z)
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)
    _z_springs()
    _surfaces()
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge",
                "-edgeBand", 0.15, "-outward", 0.0, 0.0, 1.0)
    _load_down()
    _static_setup()
    assert ops.analyze(1) == 0
    z5 = ZC + ops.nodeDisp(5, 3)
    assert abs(z5 + PTOT / (EPS + 2.0 * KSPRING)) < 1.0e-6


def test_edgeedge_far_apart_is_bitwise_inert():
    def run(edgeedge):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _build_bars(xoff=10.0)
        for t in (5, 6, 7, 8):
            ops.fix(t, 1, 1, 0)
            ops.mass(t, 1.0, 1.0, 1.0)
            ops.setNodeVel(t, 3, 0.2, "-commit")
        _surfaces()
        args = [1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0]
        if edgeedge:
            args += ["-edgeedge", "-edgeBand", 0.15]
        ops.contact(*args)
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25, 0)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        out = []
        for _ in range(8):
            assert ops.analyze(1, 1.0e-3) == 0
            out.append(tuple(ops.nodeDisp(t, 3) for t in (5, 6, 7, 8)))
        return out
    assert run(True) == run(False)


# ----------------------------------------------------------------- friction
EPS_F = 1.0e5    # normal penalty of the friction cases: N = EPS_F |ZC| = 10
KT = 1.0e4


def _friction_contact(mu, friction=True, cohesion=0.0, extra=()):
    _surfaces()
    args = [1, 1, 2, "-mortar", "-epsN", EPS_F, "-edgeedge", "-edgeBand", 0.15,
            "-outward", 0.0, 0.0, 1.0]
    if friction:
        if mu > 0.0:
            args += ["-edgeMu", mu]
        if cohesion > 0.0:
            args += ["-edgeCohesion", cohesion]
        args += ["-edgeKt", KT]
    args += list(extra)
    ops.contact(*args)


def _rake(mu, friction=None, nsteps=60, dt=1.0e-3, v0=0.5, extra=()):
    """Bottom edge held at z = ZC (constant normal force), free in y with unit masses and an
    initial +y velocity; returns the y-velocity history of node 5."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6):
        ops.fix(t, 1, 0, 1)
        ops.mass(t, 1.0, 1.0, 1.0)
        ops.setNodeVel(t, 2, v0, "-commit")
    _friction_contact(mu, friction=(mu > 0.0) if friction is None else friction, extra=extra)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("NewmarkExplicit", 0.5)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    vy = [ops.nodeVel(5, 2)]
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        vy.append(ops.nodeVel(5, 2))
    return vy


def test_explicit_rake_friction_opposes_motion():
    v0 = 0.5
    vy_fric = _rake(0.5, v0=v0)
    vy_free = _rake(0.0, v0=v0)
    assert abs(vy_free[-1] - v0) < 1.0e-6
    assert vy_fric[-1] < 0.8 * v0
    assert all(vy_fric[i + 1] <= vy_fric[i] + 1.0e-9 for i in range(len(vy_fric) - 1))
    assert vy_fric[-1] >= 0.0


def _stick(system):
    """Static sub-cone tangential load: cap = mu N + c = 0.3 * 10 + 20 = 23 > 10, so the pair
    sticks. y-springs carry the differential mode of the two bottom nodes."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6):
        ops.fix(t, 1, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    for gnd, sn, x in [(405, 5, -1.0), (406, 6, 1.0)]:
        ops.node(gnd, x, 0.0, ZC)
        ops.fix(gnd, 1, 1, 1)
        ops.element("zeroLength", 800 + sn, gnd, sn, "-mat", 1, "-dir", 2)
    _friction_contact(0.3, cohesion=20.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 5.0, 0.0)
    _static_setup(system=system, tol=1.0e-10, maxit=40)
    assert ops.analyze(1) == 0
    return ops.nodeDisp(5, 2)


def test_implicit_stick_converges():
    y = _stick("FullGeneral")
    assert 0.0 < y < 0.05


def test_symmetric_friction_tangent_is_solver_safe():
    assert abs(_stick("ProfileSPD") - _stick("FullGeneral")) < 1.0e-10


def test_zero_mu_is_bitwise_frictionless():
    assert _rake(0.0, friction=False) == _rake(0.0, friction=False, extra=("-edgeMu", 0.0))


# --------------------------------------------------------- augmented Lagrangian
ZSPRING = 0.01   # spring anchors above the master: removing the load opens the contact


def _alm_model(epsN, alm):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)
    _z_springs(zg=ZSPRING)
    _surfaces()
    args = [1, 1, 2, "-mortar", "-epsN", epsN, "-edgeedge", "-edgeBand", 0.15,
            "-outward", 0.0, 0.0, 1.0]
    if alm:
        args += ["-edgeAlm", "-edgeAugTol", 1.0e-9]
    ops.contact(*args)
    _load_down()
    _static_setup()


def test_augmentation_is_epsN_independent():
    pen_penalty, pen_alm = {}, {}
    for epsN in (1.0e5, 1.0e7):
        _alm_model(epsN, alm=True)
        assert ops.analyze(1) == 0
        pen_penalty[epsN] = ops.ladrunoEdgePenetration()   # lambdaN = 0 on the first solve
        status, _, pen = _analyze_augmented(maxAug=40)
        assert status == 0, (epsN, status, pen)
        pen_alm[epsN] = pen
    assert pen_penalty[1.0e5] / pen_penalty[1.0e7] > 10.0
    assert max(pen_alm.values()) < 1.0e-9


def test_release_after_augmentation():
    _alm_model(1.0e6, alm=True)
    assert ops.analyze(1) == 0
    status, _, pen = _analyze_augmented(maxAug=40)
    assert status == 0 and pen < 1.0e-9
    ops.pattern("Plain", 2, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 0.0, PTOT)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    _analyze_augmented(maxAug=20)
    assert ZC + ops.nodeDisp(5, 3) > 0.0
    assert ops.ladrunoEdgePenetration() < 1.0e-9


def test_first_solve_with_edgealm_is_bitwise_penalty():
    def first(alm):
        _alm_model(1.0e6, alm)
        assert ops.analyze(1) == 0
        return tuple(ops.nodeDisp(t, 3) for t in (5, 6))
    assert first(True) == first(False)


def test_edge_penetration_is_zero_without_edgealm():
    _alm_model(1.0e6, alm=False)
    assert ops.analyze(1) == 0
    assert ops.ladrunoEdgePenetration() == 0.0


# ----------------------------------------------------------------- refusals
@pytest.mark.parametrize("args", [
    (1, 1, 2, EPS, "-edgeedge"),                                   # no -mortar
    (1, 1, 2, "-mortar", "-epsN", EPS, "-edgeKn", EPS),            # edge option, no -edgeedge
    (1, 1, 2, "-mortar", "-epsN", EPS, "-tie", "-edgeedge"),       # tie
    (1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge", "-edgeTauMax", 5.0),
    (1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge", "-edgeSoft", 0.1),
])
def test_refusals(args):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    _surfaces()
    if not isinstance(args[3], str):
        ops.contactSurface(3, "-slave", 5, 6)
        args = (args[0], args[1], 3) + tuple(args[3:])
    with pytest.raises(ops.OpenSeesError):
        ops.contact(*args)
