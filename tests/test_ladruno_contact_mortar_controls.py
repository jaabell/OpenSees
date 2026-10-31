"""Controls of mortar contact: -augment, -maxGap, -gapOffset, -adjust; consistent friction
tangent at shared multi-pair nodes.

Checks:

- -augment never: a linear penalty tie (split stdBrick column) gives the same tip
  displacement for 1, 2 and 5 load steps (1e-12) and converges to the monolithic bar as
  1/epsTie; the default (commit) drifts with the step count exactly as the series-spring
  Uzawa recursion (2.2e-3, 2.1e-3, 2.04e-3 to 1e-10).
- -gapOffset -delta on two stacked blocks gives the interference pressure
  delta / (2L/E + 1/epsN) (1e-10); -augment request keeps it on a physical step and
  reaches delta E / (2L) inside the held-load bracket (1e-8); commit augments at every
  physical step; never stays on the penalty value inside the bracket.
- Without a shift a coincident interface carries no load; a cohesive bond with an
  interference prestress under -augment never is step-count independent (1e-10).
- A single tie around a closed cylinder equals the 4-sector split (1e-6) with the default
  proximity condition and with -maxGap 0.1; a -maxGap that admits the far side of the ring
  is more than 1 percent stiffer. A slave facet between two parallel master facets is
  bonded twice when -maxGap admits the far one.
- -adjust on a faceted tube in a faceted ring (spurious gaps and penetrations as meshed)
  starts with exactly zero displacement and penetration; the tolerance branch adjusts
  only nodes with |gbar_ref| <= tol; the adjusted interface carries a lateral load and a
  shrink fit moves every tube node radially inward with no net lateral motion.
- Shared multi-pair ridge nodes of a creased stiff roof with mu = 0.2: with -consistanttan
  the contact force equals the analytic rigid descent (2e-4) in at most 6 iterations per
  step, with no lateral force, and the converged state equals that of the symmetric tangent
  when that one converges (1e-6).
- Refusals: bad values, the options without -mortar, a gap shift on a tie.

Every test fails on the base with an unknown command.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


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


def _static(tol=1.0e-12, maxit=60, system="FullGeneral", test="NormDispIncr"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.test(test, tol, maxit, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


# ---------------------------------------------------------------------- -augment
E_COL, P_COL = 2.0e4, 20.0          # sigma = 20, u_tip(monolithic) = 2*P/E = 2e-3


def _split_column(epsTie, nsteps, extra=()):
    """The split column: bottom brick [0,1]^2x[0,1] (MASTER face z=1) tied to a top block of
    TWO bricks split at x=0.5 (SLAVE face, non-matching). Uniaxial (x,y fixed, nu=0). End load P in
    `nsteps` equal LoadControl increments. Returns the tip z-displacement."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E_COL, 0.0)
    bot = [(1, 0, 0, 0), (2, 1, 0, 0), (3, 1, 1, 0), (4, 0, 1, 0),
           (5, 0, 0, 1), (6, 1, 0, 1), (7, 1, 1, 1), (8, 0, 1, 1)]
    top = [(9, 0, 0, 1), (10, .5, 0, 1), (11, 1, 0, 1), (12, 0, 1, 1), (13, .5, 1, 1), (14, 1, 1, 1),
           (15, 0, 0, 2), (16, .5, 0, 2), (17, 1, 0, 2), (18, 0, 1, 2), (19, .5, 1, 2), (20, 1, 1, 2)]
    for t, x, y, z in bot + top:
        ops.node(t, float(x), float(y), float(z))
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    for t in range(5, 21):
        ops.fix(t, 1, 1, 0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 2, 9, 10, 13, 12, 15, 16, 19, 18, 1)
    ops.element("stdBrick", 3, 10, 11, 14, 13, 16, 17, 20, 19, 1)
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave-segments", 4, 9, 10, 13, 12, 10, 11, 14, 13)
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", epsTie, "-outward", 0.0, 0.0, 1.0, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (15, 17, 18, 20):
        ops.load(t, 0.0, 0.0, P_COL / 8.0)
    for t in (16, 19):
        ops.load(t, 0.0, 0.0, P_COL / 4.0)
    _static()
    ops.integrator("LoadControl", 1.0 / nsteps)
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
    return ops.nodeDisp(15, 3)


def test_augment_never_linear_tie_is_step_count_independent():
    """`-augment never`: a linear penalty tie is a LINEAR problem, so 1, 2 and 5 load steps give the
    same answer (to Newton round-off). The default augments at every commit, so the same model drifts
    with the step count."""
    u = [_split_column(1.0e5, n, ("-augment", "never")) for n in (1, 2, 5)]
    for ui in u[1:]:
        assert abs(ui - u[0]) <= 1.0e-12 * abs(u[0]), f"never: step-count dependent {u}"
    d = [_split_column(1.0e5, n) for n in (1, 2, 5)]
    assert d[0] == pytest.approx(u[0], rel=1e-14)       # 1 step: nothing to augment before it
    # the column is exactly 1-D (nu = 0), so the drift of the default is the series-spring Uzawa
    # recursion: 2.2e-3 / 2.1e-3 / 2.04e-3 for 1 / 2 / 5 steps.
    assert d == pytest.approx([2.2e-3, 2.1e-3, 2.04e-3], rel=1e-10), d
    assert u == pytest.approx([2.2e-3] * 3, rel=1e-12), u


def test_augment_never_converges_to_the_tie_limit():
    """Pure penalty: the tip error against the monolithic bar (2P/E) is O(1/epsTie) -- one decade
    of epsTie removes one decade of error -- and vanishes in the limit."""
    exact = 2.0 * P_COL / E_COL
    err = [abs(_split_column(eps, 1, ("-augment", "never")) - exact) for eps in (1e5, 1e6, 1e7, 1e8)]
    for a, b in zip(err, err[1:]):
        assert 8.0 < a / b < 12.0, f"penalty error not O(1/epsTie): {err}"
    assert err[-1] < 2.0e-4 * exact          # measured 1.0e-4 relative at epsTie = 1e8


# ------------------------------------------------------------------------------------ blocks
E_BLK, DELTA = 2.0e4, 1.0e-3        # two unit blocks, prescribed interference delta


def _blocks(extra, epsN=1.0e6, top_free_x=False, shear=0.0):
    """Two stacked unit bricks, bottom z=0 clamped, top z=2 held in z. Uniaxial in z (x,y fixed,
    nu=0) unless top_free_x (then the top block's x DOFs are free and `shear` is applied in x at the
    top face -- the interface carries it by cohesion). Frictionless mortar contact across the
    MATCHING, coincident z=1 interface, gap shifted by the options in `extra`."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E_BLK, 0.0)
    sq = [(0, 0), (1, 0), (1, 1), (0, 1)]
    for k, z in enumerate((0.0, 1.0)):
        for i, (x, y) in enumerate(sq):
            ops.node(1 + 4 * k + i, float(x), float(y), z)
    for k, z in enumerate((1.0, 2.0)):
        for i, (x, y) in enumerate(sq):
            ops.node(9 + 4 * k + i, float(x), float(y), z)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6, 7, 8):
        ops.fix(t, 1, 1, 0)
    for t in (9, 10, 11, 12):
        ops.fix(t, 0 if top_free_x else 1, 1, 0)
    for t in (13, 14, 15, 16):
        ops.fix(t, 0 if top_free_x else 1, 1, 1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 2, 9, 10, 11, 12, 13, 14, 15, 16, 1)
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave-segments", 4, 9, 10, 11, 12)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (13, 14, 15, 16):
        ops.load(t, shear / 4.0, 0.0, 0.0)
    _static()


def _base_pressure():
    ops.reactions()
    return sum(ops.nodeReaction(t, 3) for t in (1, 2, 3, 4))   # = sigma*A, A = 1


def test_gap_gap_offset_interference_pressure():
    """-gapOffset -delta prescribes an interference fit on straight geometry. Two blocks in series
    with the penalty spring: sigma = delta/(2L/E + 1/epsN) exactly (constant stress, flat matched
    interface); held-load augmentation removes the penalty compliance -> sigma = delta*E/(2L)."""
    eps = 1.0e6
    pen = DELTA / (2.0 / E_BLK + 1.0 / eps)
    exact = DELTA * E_BLK / 2.0
    _blocks(("-gapOffset", -DELTA, "-augment", "request"), epsN=eps)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    assert _base_pressure() == pytest.approx(pen, rel=1e-10)
    # a second PHYSICAL step under `request` does not augment: still the penalty value.
    assert ops.analyze(1) == 0
    assert _base_pressure() == pytest.approx(pen, rel=1e-10)
    # inside the analyze_augmented bracket `request` augments -> the exact interference pressure.
    st, _, _ = _analyze_augmented(maxAug=40, augTol=1.0e-12)
    assert st == 0
    assert _base_pressure() == pytest.approx(exact, rel=1e-8)


def test_augment_modes_on_the_interference_blocks():
    """The three -augment modes on one closed form. `commit` (the default) augments on every
    physical commit (drifts toward delta*E/2L step by step); `never` stays on the penalty value
    even inside the analyze_augmented bracket."""
    eps = 1.0e6
    pen = DELTA / (2.0 / E_BLK + 1.0 / eps)
    _blocks(("-gapOffset", -DELTA), epsN=eps)                      # default = commit
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    p1 = _base_pressure()
    assert ops.analyze(1) == 0
    assert p1 == pytest.approx(pen, rel=1e-10)
    assert _base_pressure() > p1 * (1.0 + 1e-4)                    # augmented by the commit
    _blocks(("-gapOffset", -DELTA, "-augment", "never"), epsN=eps)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    st, _, _ = _analyze_augmented(maxAug=5, augTol=1.0e-12)
    assert st == 1                                                  # never converges: inert bracket
    assert _base_pressure() == pytest.approx(pen, rel=1e-10)


def test_gap_no_offset_is_inert():
    """A coincident interface with no shift carries no load (the baseline the offset acts on)."""
    _blocks(("-augment", "never"))
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    assert abs(_base_pressure()) < 1.0e-12


def test_augment_cohesive_bond_is_step_count_independent_under_never():
    """A cohesion-only mortar bond (mu = 0) with
    `-augment never` and an interference prestress, a shear ramp in 1 or 5 steps converges and gives
    the same answer (elastic stick is linear in the penalty)."""
    out = []
    for n in (1, 5):
        _blocks(("-gapOffset", -DELTA, "-augment", "never", "-cohesion", 100.0, "-epsT", 1.0e6),
                epsN=1.0e6, top_free_x=True, shear=5.0)
        ops.integrator("LoadControl", 1.0 / n)
        for _ in range(n):
            assert ops.analyze(1) == 0
        out.append(ops.nodeDisp(13, 1))
    assert out[0] != 0.0
    assert abs(out[1] - out[0]) <= 1.0e-10 * abs(out[0]), out


# --------------------------------------------------------------------------- tube in a ring
R_IN, R_IF, R_OUT = 0.3, 0.5, 1.2     # inner tube [R_IN, R_IF], ring [R_IF, R_OUT]
E_TUBE, NU_TUBE = 1.0e3, 0.3
EPS_C = 1.0e5                         # contact penalty ~ 20 E/h (200 E/h chatters on this
                                      # faceted interface)


def _polar(r, th, z):
    return (r * math.cos(th), r * math.sin(th), z)


def _tube_in_ring(nti=16, nto=24):
    """A tube (one radial layer, `nti` sectors) inside a ring (two radial layers, `nto` sectors),
    one element through z = [0, 0.5]; uz fixed everywhere (a plane-strain slice), the ring's outer
    face clamped. The interface r = R_IF is NON-matching (nti != nto): both are chord facets of the
    same circle, so as meshed there are spurious gaps AND penetrations. Returns the facet lists
    (master = tube outer face, slave = ring inner face) with each facet's centroid angle."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E_TUBE, NU_TUBE)
    H = 0.5

    def tag_t(i, j, k):          # tube: radial i in {0,1}, sector j, level k
        return 1 + i + 2 * (j % nti) + 2 * nti * k

    def tag_r(i, j, k):          # ring: radial i in {0,1,2}
        return 10001 + i + 3 * (j % nto) + 3 * nto * k

    for k in (0, 1):
        for j in range(nti):
            for i, r in enumerate((R_IN, R_IF)):
                ops.node(tag_t(i, j, k), *_polar(r, 2 * math.pi * j / nti, k * H))
                # uz fixed (plane-strain slice); uy fixed on the x axis (theta = 0, pi): the load is
                # along x, so this is the symmetry plane, and it removes the tube's free rigid
                # rotation about z (frictionless contact on a near-circle has no torsional stiffness).
                on_x_axis = (2 * j) % nti == 0
                ops.fix(tag_t(i, j, k), 0, 1 if on_x_axis else 0, 1)
        for j in range(nto):
            for i, r in enumerate((R_IF, 0.5 * (R_IF + R_OUT), R_OUT)):
                t = tag_r(i, j, k)
                ops.node(t, *_polar(r, 2 * math.pi * j / nto, k * H))
                ops.fix(t, 1, 1, 1) if i == 2 else ops.fix(t, 0, 0, 1)
    e = 1
    for j in range(nti):
        ops.element("stdBrick", e, tag_t(0, j, 0), tag_t(1, j, 0), tag_t(1, j + 1, 0), tag_t(0, j + 1, 0),
                    tag_t(0, j, 1), tag_t(1, j, 1), tag_t(1, j + 1, 1), tag_t(0, j + 1, 1), 1)
        e += 1
    for j in range(nto):
        for i in (0, 1):
            ops.element("stdBrick", e, tag_r(i, j, 0), tag_r(i + 1, j, 0), tag_r(i + 1, j + 1, 0),
                        tag_r(i, j + 1, 0), tag_r(i, j, 1), tag_r(i + 1, j, 1), tag_r(i + 1, j + 1, 1),
                        tag_r(i, j + 1, 1), 1)
            e += 1
    master = [((j + 0.5) * 2 * math.pi / nti,
               [tag_t(1, j, 0), tag_t(1, j + 1, 0), tag_t(1, j + 1, 1), tag_t(1, j, 1)])
              for j in range(nti)]
    slave = [((j + 0.5) * 2 * math.pi / nto,
              [tag_r(0, j, 0), tag_r(0, j + 1, 0), tag_r(0, j + 1, 1), tag_r(0, j, 1)])
             for j in range(nto)]
    loaded = [tag_t(0, j, k) for j in range(nti) for k in (0, 1)]
    tube_nodes = [tag_t(i, j, k) for i in (0, 1) for j in range(nti) for k in (0, 1)]
    return master, slave, loaded, tube_nodes


def _declare(master, slave, sectors, args, outward=False):
    """sectors=1: one contact over the whole closed interface; sectors=4: a quadrant split
    (facets binned by centroid angle), optionally with the quadrant bisector as -outward."""
    stag = 1
    for q in range(sectors):
        lo, hi = q * 2 * math.pi / sectors, (q + 1) * 2 * math.pi / sectors
        m = [n for th, f in master if lo <= th < hi for n in f]
        s = [n for th, f in slave if lo <= th < hi for n in f]
        ops.contactSurface(stag, "-master", 4, *m)
        ops.contactSurface(stag + 1, "-slave-segments", 4, *s)
        ow = ()
        if outward:
            b = 0.5 * (lo + hi)
            ow = ("-outward", math.cos(b), math.sin(b), 0.0)
        ops.contact(q + 1, stag, stag + 1, "-mortar", *args, *ow)
        stag += 2


def _lateral(loaded, F):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in loaded:
        ops.load(t, F / len(loaded), 0.0, 0.0)


def _mean_ux(nodes):
    return sum(ops.nodeDisp(t, 1) for t in nodes) / len(nodes)


def _tube_tie(sectors, extra):
    master, slave, loaded, tube = _tube_in_ring()
    _declare(master, slave, sectors, ("-tie", "-epsTie", 1.0e6) + tuple(extra))
    _lateral(loaded, 1.0)
    _static(tol=1e-12, system="UmfPack")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    return _mean_ux(tube)


def test_closed_cylinder_tie_matches_sector_split():
    """One tie around the whole cylinder equals the 4-sector split, with the default proximity
    condition and with -maxGap 0.1. A -maxGap wide enough to admit the far side of the ring ties
    every tube facet to the antipodal ring facets too (the clip accepts the anti-parallel facet),
    which stiffens the tube laterally: the bound the guard removes."""
    split = _tube_tie(4, ("-augment", "never"))
    default = _tube_tie(1, ("-augment", "never"))
    guarded = _tube_tie(1, ("-augment", "never", "-maxGap", 0.1))
    bare = _tube_tie(1, ("-augment", "never", "-maxGap", 10.0))
    assert default == pytest.approx(split, rel=1e-6)
    assert guarded == pytest.approx(split, rel=1e-6)
    assert bare < 0.99 * split, f"antipodal ties no longer over-stiff? bare {bare}, split {split}"


def test_gap_adjust_faceted_tube_starts_stress_free():
    """A faceted tube in a faceted ring (16 vs 24 chords of one circle): as meshed the interface has
    penetrations (tube nodes poke past the ring's chords) and gaps. Zero load: the unadjusted contact
    pushes the bodies apart; with -adjust nothing moves and nothing penetrates -- exactly."""
    res = {}
    for key, extra in (("raw", ()), ("adjust", ("-adjust",))):
        master, slave, loaded, tube = _tube_in_ring()
        _declare(master, slave, 4, ("-epsN", EPS_C, "-augment", "never") + extra, outward=True)
        _static(tol=1e-10, system="UmfPack", test="NormUnbalance")
        ops.integrator("LoadControl", 1.0)
        assert ops.analyze(1) == 0
        umax = max(abs(ops.nodeDisp(t, d)) for t in tube for d in (1, 2))
        res[key] = (umax, ops.ladrunoMortarPenetration())
    assert res["raw"][0] > 1.0e-6, res
    assert res["adjust"] == (0.0, 0.0), res


def test_gap_adjust_tolerance_branch():
    """`-adjust tol` adjusts only nodes whose reference nodal gap |gbar_ref| <= tol. On this symmetric
    faceted slice every load-bearing (penetrating) node has the same as-meshed |gbar_ref| in
    (6e-5, 7e-5] (measured), so the branch is exercised on both sides of that value: a tol below it
    adjusts none of them (the raw response, bit for bit), a tol above it adjusts all (exact zero)."""
    def umax(extra):
        master, slave, loaded, tube = _tube_in_ring()
        _declare(master, slave, 4, ("-epsN", EPS_C, "-augment", "never") + extra, outward=True)
        _static(tol=1e-10, system="UmfPack", test="NormUnbalance")
        ops.integrator("LoadControl", 1.0)
        assert ops.analyze(1) == 0
        return max(abs(ops.nodeDisp(t, d)) for t in tube for d in (1, 2))
    raw = umax(())
    assert raw > 1.0e-6
    assert umax(("-adjust", 1.0e-12)) == raw
    assert umax(("-adjust", 5.0e-5)) == raw          # just below the as-meshed penetration
    assert umax(("-adjust", 8.0e-5)) == 0.0          # just above it
    assert umax(("-adjust", 1.0)) == 0.0


def test_gap_adjust_then_load_and_interference():
    """With -adjust the faceted interface is usable: a lateral load converges (the tube bears on the
    front and opens at the back), and -adjust + -gapOffset gives a uniform radial shrink fit whose
    net lateral displacement is zero by symmetry."""
    master, slave, loaded, tube = _tube_in_ring()
    _declare(master, slave, 4, ("-epsN", EPS_C, "-augment", "never", "-adjust"), outward=True)
    _lateral(loaded, 1.0)
    _static(tol=1e-10, system="UmfPack", test="NormUnbalance")
    ops.integrator("LoadControl", 0.5)
    for _ in range(2):
        assert ops.analyze(1) == 0
    assert _mean_ux(tube) > 0.0
    # shrink fit: the tube moves radially INWARD everywhere, no net lateral motion.
    master, slave, loaded, tube = _tube_in_ring()
    _declare(master, slave, 4, ("-epsN", EPS_C, "-augment", "never", "-adjust",
                                "-gapOffset", -1.0e-4), outward=True)
    _static(tol=1e-10, system="UmfPack", test="NormUnbalance")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    assert abs(_mean_ux(tube)) < 1.0e-12
    for t in tube:
        x, y = ops.nodeCoord(t, 1), ops.nodeCoord(t, 2)
        ur = (x * ops.nodeDisp(t, 1) + y * ops.nodeDisp(t, 2)) / math.hypot(x, y)
        assert ur < 0.0


# ------------------------------------------------------------------------------ command surface
def _two_facets():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(1 + i, float(x), float(y), 0.0)
        ops.node(5 + i, float(x), float(y), 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)


def _refused(capfd, why, *args):
    """Refused for the named reason (a build without the options also raises, with "unexpected
    token")."""
    capfd.readouterr()
    with pytest.raises(ops.OpenSeesError):
        ops.contact(*args)
    err = capfd.readouterr()
    text = err.out + err.err
    assert why in text, text


@pytest.mark.parametrize("bad,why", [
    (("-augment", "sometimes"), "need commit|request|never"),
    (("-maxGap", 0.0), "need a distance > 0"),
    (("-maxGap", -1.0), "need a distance > 0"),
    (("-adjust", -1.0), "optional tolerance must be a value > 0"),
])
def test_refuses_bad_values(capfd, bad, why):
    _two_facets()
    _refused(capfd, why, 1, 1, 2, "-mortar", "-epsN", 1.0e5, *bad)


@pytest.mark.parametrize("flag", [("-augment", "never"), ("-maxGap", 0.1),
                                  ("-gapOffset", -1e-3), ("-adjust",)])
def test_flags_are_mortar_only(capfd, flag):
    _two_facets()
    ops.contactSurface(3, "-slave", 5, 6, 7, 8)
    _refused(capfd, "are -mortar options", 1, 1, 3, 1.0e5, 0.0, 0.0, *flag)


@pytest.mark.parametrize("flag", [("-gapOffset", -1e-3), ("-adjust",), ("-adjust", 1e-3)])
def test_gap_shift_refused_on_tie(capfd, flag):
    _two_facets()
    _refused(capfd, "do not apply to -tie", 1, 1, 2, "-mortar", "-tie", "-epsTie", 1.0e5, *flag)


def test_accepted_forms():
    """Every documented spelling parses (tie + augment/maxGap; contact + all four)."""
    _two_facets()
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", 1.0e5, "-augment", "request", "-maxGap", 0.1)
    _two_facets()
    ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e5, "-adjust", 1.0e-3, "-gapOffset", -1e-4,
                "-maxGap", 0.1, "-augment", "commit", "-outward", 0.0, 0.0, 1.0)
    _two_facets()
    ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e5, "-adjust", "-augment", "never")

# ---------------------------------------------------------------------- -maxGap on a tie
def _box_tie(extra=()):
    """Master: two parallel fixed facets, z = 0 and z = 1. Slave: one facet coincident with
    the top one, free in z, pulled up by 1. Returns the slave z-displacement."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for base, z in ((1, 0.0), (5, 1.0), (9, 1.0)):
        for t, x, y in zip(range(base, base + 4), (0.0, 1.0, 1.0, 0.0), (0.0, 0.0, 1.0, 1.0)):
            ops.node(t, x, y, z)
            ops.fix(t, 1, 1, 0 if base == 9 else 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave-segments", 4, 9, 10, 11, 12)
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", 1.0e4, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (9, 10, 11, 12):
        ops.load(t, 0.0, 0.0, 0.25)
    _static()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    return ops.nodeDisp(9, 3)


def test_maxgap_admitting_the_far_facet_doubles_the_bond():
    assert _box_tie() == pytest.approx(1.0e-4, rel=1e-10)
    assert _box_tie(("-maxGap", 0.5)) == pytest.approx(1.0e-4, rel=1e-10)
    assert _box_tie(("-maxGap", 2.0)) == pytest.approx(0.5e-4, rel=1e-10)


# ----------------------------------------------------- creased roof, shared multi-pair nodes
ALPHA = math.radians(20.0)


def _roof158(mu=0.0, coh=0.0, solid_e=0.0, flags=(), steps=4, pz=40.0, eps=1.0e6, delta=1.0e-3,
          kd=1.0e3, augment_never=False, mw=1.5, yw=0.5):
    """Roof z = -tan(a)|x|, ridge along y. Master: 2 fixed quads over |x| <= mw, -yw <= y <= 1+yw
    (the default oversizes the slave). Slave: 4 x 3 non-matching quads, ridge nodes shared by both
    flanks, lowered by
    delta; every slave node on kd springs (x, y, z). solid_e > 0 adds a 0.3-thick brick layer under
    the slave surface (stdBrick, E = solid_e). Load -pz on each surface node, force control.
    Returns (iterations per step, ok, surface node tags, all spring node tags)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ta = math.tan(ALPHA)
    xm, ym = [-mw, 0.0, mw], [-yw, 1.0 + yw]
    mt, t = {}, 1
    for j, y in enumerate(ym):
        for i, x in enumerate(xm):
            ops.node(t, x, y, -ta * abs(x))
            ops.fix(t, 1, 1, 1)
            mt[(i, j)] = t
            t += 1
    mq = []
    for i in range(2):
        mq += [mt[(i, 0)], mt[(i + 1, 0)], mt[(i + 1, 1)], mt[(i, 1)]]
    xs, ys = [-1.0, -0.45, 0.0, 0.45, 1.0], [0.0, 0.3, 0.75, 1.0]
    st, surf, t = {}, [], 101
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            ops.node(t, x, y, -ta * abs(x) - delta)
            st[(i, j)] = t
            surf.append(t)
            t += 1
    sq = []
    for j in range(len(ys) - 1):
        for i in range(len(xs) - 1):
            sq += [st[(i, j)], st[(i + 1, j)], st[(i + 1, j + 1)], st[(i, j + 1)]]
    ops.contactSurface(1, "-master", 4, *mq)
    ops.contactSurface(2, "-slave-segments", 4, *sq)
    opts = ["-mortar", "-epsN", eps, "-epsT", eps, "-outward", 0.0, 0.0, 1.0]
    if mu > 0.0:
        opts += ["-mu", mu]
    if coh > 0.0:
        opts += ["-cohesion", coh]
    if augment_never:
        opts += ["-augment", "never"]
    ops.contact(1, 1, 2, *opts, *flags)
    springs = list(surf)
    if solid_e > 0.0:
        ops.nDMaterial("ElasticIsotropic", 7, solid_e, 0.25)
        bot = {}
        for s in surf:
            c = ops.nodeCoord(s)
            bot[s] = 3000 + s
            ops.node(bot[s], c[0], c[1], c[2] - 0.3)
            springs.append(bot[s])
        for q in range(0, len(sq), 4):
            a, b, c, d = sq[q:q + 4]
            ops.element("stdBrick", 7000 + q, bot[a], bot[b], bot[c], bot[d], a, b, c, d, 7)
    ops.uniaxialMaterial("Elastic", 1, kd)
    for s in springs:
        g = 50000 + s
        ops.node(g, *ops.nodeCoord(s))
        ops.fix(g, 1, 1, 1)
        ops.element("zeroLength", 60000 + s, g, s, "-mat", 1, 1, 1, "-dir", 1, 2, 3)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for s in surf:
        ops.load(s, 0.0, 0.0, -pz)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1.0e-7, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / steps)
    ops.analysis("Static")
    its = []
    for _ in range(steps):
        if ops.analyze(1) != 0:
            return its, False, surf, springs
        its.append(ops.testIter())
    return its, True, surf, springs


# ---------------------------------------------------------------------------------------- (a)
SOLID = dict(mu=0.2, solid_e=1.0e11, pz=60.0, delta=1.0e-6, augment_never=True, steps=4)


def _descent_analytic(nsurf, nspring, pz, kd, kn, mu, delta, lam):
    """Rigid vertical descent w of the roof in full slip: the applied load lam*nsurf*pz is carried
    by the springs (kd*nspring*w) and the two flanks (area A = 1/cos(a) each), where the penalty
    pressure is p = kn*(delta + w)*cos(a) and the slip friction mu*p acts up-slope:
        lam*nsurf*pz = kd*nspring*w + 2*A*p*(cos(a) + mu*sin(a)).
    Returns (w, Fz) with Fz the vertical contact force on the slave."""
    ca, sa = math.cos(ALPHA), math.sin(ALPHA)
    c = 2.0 / ca * kn * ca * (ca + mu * sa)
    w = (lam * nsurf * pz - c * delta) / (kd * nspring + c)
    return w, c * (delta + w)


def test_roof_mu_multipair_roof_consistanttan_analytic():
    """(a) mu > 0 at shared multi-pair ridge nodes with the -consistanttan recipe: analytic contact
    force and <= 6 Newton iterations per step (measured 3, 4, 4, 4)."""
    its, ok, surf, springs = _roof158(flags=("-consistanttan",), **SOLID)
    assert ok, f"-consistanttan: Newton failed after {len(its)} step(s) {its}"
    assert max(its) <= 6, f"-consistanttan: iterations per step {its} (want <= 6)"
    kd = 1.0e3
    w, fz = _descent_analytic(len(surf), len(springs), SOLID["pz"], kd, 1.0e6, SOLID["mu"],
                              SOLID["delta"], 1.0)
    # contact force on the slave = applied load - spring forces (sum over all spring nodes)
    got = len(surf) * SOLID["pz"] + sum(kd * ops.nodeDisp(n, 3) for n in springs)
    fx = sum(kd * ops.nodeDisp(n, 1) for n in springs)
    assert got == pytest.approx(fz, rel=2e-4), f"Fz {got:.6f} != analytic {fz:.6f}"
    assert abs(fx) < 1e-8 * fz, f"spurious Fx = {fx:.3e}"


# ---------------------------------------------------------------------------------------- (c)
def _state(nodes):
    return [ops.nodeDisp(n) for n in nodes]


def test_roof_consistanttan_same_converged_state():
    """(c) -consistanttan changes the iterations only: when the symmetric default converges too, both
    tangents reach the same state (to the solver tolerance, not bitwise)."""
    _, ok1, surf, springs = _roof158(**SOLID)
    a = _state(springs)
    _, ok2, _, _ = _roof158(flags=("-consistanttan",), **SOLID)
    b = _state(springs)
    assert ok2, "-consistanttan failed"
    if not ok1:
        pytest.skip("the symmetric default did not converge on this platform: nothing to compare")
    scale = max(abs(v) for u in a for v in u)
    for ua, ub in zip(a, b):
        for va, vb in zip(ua, ub):
            assert abs(va - vb) <= 1e-6 * scale, f"{va} != {vb}"
