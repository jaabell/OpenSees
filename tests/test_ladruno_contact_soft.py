"""Soft-constraint contact penalty for explicit analysis (-soft, -edgeSoft).

Under an integrator whose element tangent assembles no stiffness (stock CentralDifference,
NewmarkExplicit, ExplicitDifference, ExplicitBathe), the penalty of a contact pair is replaced at
every evaluation by k = s 4 m_eff / dt^2, m_eff the effective mass of its gap mode, so that the
contact mode has omega dt = 2 sqrt(s) (0.63 for s = 0.1) and does not reduce the critical step.
At the structural step dt = 1e-3 a stiff penalty kn = 1e8 on a unit mass has omega dt = 10 and the
central-difference integration diverges.

Checks:

- Rigid plane, a unit mass pressed by a sustained load: the stiff penalty diverges and -soft stays
  bounded (|z| < 0.05), under NewmarkExplicit, stock CentralDifference and ExplicitDifference;
  the bare -soft (s = 0.1) is accepted.
- The mass of an element (a truss with density, no nodal mass) sizes the soft penalty.
- A free impact on the plane rebounds with restitution within 5 percent of 1 with -soft, while
  the stiff penalty blows up or injects energy (restitution > 2).
- Node-to-segment (3D and 2D) and the stick penalty of friction are stabilized the same way.
- Mortar facets (-mortar -soft, the segment-based soft penalty): a slave facet whose nodes all lie
  outside the master but which overlaps it is held, where node-to-segment contact lets it fall
  through; a matched facet under load is stable where the stiff penalty diverges.
- Edge-to-edge (-edgeSoft): a falling bar rebounds off a crossing edge at the structural step,
  the stiff edge penalty diverges.
- Under a static integrator -soft and -edgeSoft are bitwise inert.
- Refusals: -soft without a base penalty, -soft with -tie, -edgeSoft without -edgeedge, a
  non-positive scale, -soft on a 2D mortar contact.

Every test fails on the base with an unknown command.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


M = 1.0
DT = 1.0e-3
KN_STIFF = 1.0e8
RUNAWAY = 1.0


def _explicit(integrator="NewmarkExplicit"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    if integrator == "NewmarkExplicit":
        ops.integrator("NewmarkExplicit", 0.5)
    else:
        ops.integrator(integrator)
    ops.algorithm("Linear")
    ops.analysis("Transient")


def _run(node, dof, nsteps, dt=DT):
    """Integrate; returns (diverged, max |u|)."""
    maxabs = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return True, maxabs
        u = ops.nodeDisp(node, dof)
        if not math.isfinite(u):
            return True, maxabs
        maxabs = max(maxabs, abs(u))
        if maxabs > RUNAWAY:
            return True, maxabs
    return False, maxabs


def _static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


# ---------------------------------------------------------------- rigid plane
def _plane_loaded(kn, soft, integrator="NewmarkExplicit", v=-0.5, P=1.0e3, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.mass(1, M, M, M)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(1, "-slave", 1)
    extra = () if soft is None else (("-soft",) if soft == "bare" else ("-soft", soft))
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, kn, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    if v != 0.0:
        ops.setNodeVel(1, 3, v, "-commit")
    _explicit(integrator)
    return _run(1, 3, nsteps)


@pytest.mark.parametrize("integrator", ["NewmarkExplicit", "CentralDifference",
                                        "ExplicitDifference"])
def test_soft_stabilizes_the_plane(integrator):
    # stock CentralDifference starts with Ut-1 = Ut (no initial velocity), so the
    # loaded mass starts at rest there
    v = 0.0 if integrator != "NewmarkExplicit" else -0.5
    div_stiff, _ = _plane_loaded(KN_STIFF, None, integrator, v=v)
    div_soft, maxabs = _plane_loaded(KN_STIFF, 0.1, integrator, v=v)
    assert div_stiff
    assert not div_soft and maxabs < 0.05, maxabs


def test_bare_soft_uses_the_default_scale():
    div, maxabs = _plane_loaded(KN_STIFF, "bare", nsteps=2000)
    assert not div and maxabs < 0.05


def test_soft_uses_element_mass():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)               # no nodal mass
    ops.node(2, 0.0, 0.0, 1.0)
    ops.fix(2, 1, 1, 1)
    ops.fix(1, 1, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, 1.0e-3)
    ops.element("truss", 1, 1, 2, 1.0, 1, "-rho", 2.0)    # lumped mass 1 at node 1
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN_STIFF, "-soft", 0.1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)
    ops.setNodeVel(1, 3, -0.5, "-commit")
    _explicit()
    div, maxabs = _run(1, 3, 3000)
    assert not div, maxabs


def _plane_impact(soft, v_mag=2.0, gap0=0.02, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, gap0)
    ops.mass(1, M, M, M)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(1, "-slave", 1)
    extra = ("-soft", soft) if soft > 0.0 else ()
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, KN_STIFF, *extra)
    ops.setNodeVel(1, 3, -v_mag, "-commit")
    _explicit()
    v_out = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            return True, 0.0
        z = gap0 + ops.nodeDisp(1, 3)
        if not math.isfinite(z) or abs(z) > 1.0e6:
            return True, 0.0
        if z > 0.5 * gap0 and ops.nodeVel(1, 3) > 0.0:
            v_out = ops.nodeVel(1, 3)
    return False, abs(v_out / v_mag)


def test_soft_impact_restitution():
    blew, e_soft = _plane_impact(0.1)
    assert not blew and abs(1.0 - e_soft) < 0.05, e_soft
    blew, e_stiff = _plane_impact(0.0)
    assert blew or e_stiff > 2.0, e_stiff


def _plane_static(soft):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 0)
    ops.node(90, 0.0, 0.0, 0.0)
    ops.fix(90, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 90, 0.1)
    ops.element("zeroLength", 90, 90, 1, "-mat", 90, "-dir", 3)
    ops.contactSurface(1, "-slave", 1)
    extra = ("-soft", soft) if soft > 0.0 else ()
    ops.contactPlane(1, 1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0e6, *extra)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0e3)
    _static()
    return ops.nodeDisp(1, 3)


def test_soft_is_inert_under_statics():
    z_off = _plane_static(0.0)
    assert _plane_static(0.1) == z_off
    assert abs(z_off + 1.0e3 / 1.0e6) < 1.0e-9


# ------------------------------------------------------------ node-to-segment
def _segment_loaded(soft, kt=0.0, mu=0.0, Pt=0.0, free_x=False, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.node(5, 0.0, 0.0, 0.0)
    ops.mass(5, M, M, M)
    ops.fix(5, 0 if free_x else 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)
    extra = ("-outward", 0.0, 0.0, 1.0) + (("-soft", soft) if soft > 0.0 else ())
    ops.contact(1, 1, 2, KN_STIFF, kt, mu, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, Pt, 0.0, -1.0e3)
    if not free_x:
        ops.setNodeVel(5, 3, -0.5, "-commit")
    _explicit()
    return _run(5, 1 if free_x else 3, nsteps)


def test_soft_stabilizes_node_to_segment():
    assert _segment_loaded(0.0)[0]
    div, maxabs = _segment_loaded(0.1)
    assert not div and maxabs < 0.05, maxabs


def test_soft_stabilizes_the_stick_penalty():
    # cap = mu N = 5e3 > 50: continuous stick; the soft normal alone leaves a
    # stiff stick mode kt = 1e8
    assert _segment_loaded(0.0, kt=KN_STIFF, mu=5.0, Pt=50.0, free_x=True)[0]
    div, maxabs = _segment_loaded(0.1, kt=KN_STIFF, mu=5.0, Pt=50.0, free_x=True)
    assert not div, maxabs


def _segment_2d(soft, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, x in ((1, -1.0), (2, 1.0)):
        ops.node(t, x, 0.0)
        ops.fix(t, 1, 1)
    ops.node(5, 0.0, 0.0)
    ops.mass(5, M, M)
    ops.fix(5, 1, 0)
    ops.contactSurface(1, "-master", 2, 1, 2)
    ops.contactSurface(2, "-slave", 5)
    extra = ("-outward", 0.0, 1.0) + (("-soft", soft) if soft > 0.0 else ())
    ops.contact(1, 1, 2, KN_STIFF, 0.0, 0.0, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, -1.0e3)
    ops.setNodeVel(5, 2, -0.5, "-commit")
    _explicit()
    return _run(5, 2, nsteps)


def test_soft_stabilizes_2d_node_to_segment():
    assert _segment_2d(0.0)[0]
    div, maxabs = _segment_2d(0.1)
    assert not div and maxabs < 0.05, maxabs


# --------------------------------------------------------------------- mortar
def _master_quad():
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)


def _wide_facet_drop(nts, soft=0.1, z0=-0.005, v=-1.0, kn=1.0e6, nsteps=1500):
    """A slave facet [-0.5, 1.5]^2, every node outside the master [0, 1]^2, driven down."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_quad()
    for t, (x, y) in zip((5, 6, 7, 8), [(-0.5, -0.5), (1.5, -0.5), (1.5, 1.5), (-0.5, 1.5)]):
        ops.node(t, float(x), float(y), z0)
        ops.mass(t, M, M, M)
        ops.fix(t, 1, 1, 0)
    if nts:
        ops.contactSurface(2, "-slave", 5, 6, 7, 8)
        ops.contact(1, 1, 2, kn, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
        ops.contact(1, 1, 2, "-mortar", "-epsN", kn, "-soft", soft, "-outward", 0.0, 0.0, 1.0)
    for t in (5, 6, 7, 8):
        ops.setNodeVel(t, 3, v, "-commit")
    _explicit()
    min_z, last_z = z0, z0
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            break
        last_z = z0 + ops.nodeDisp(5, 3)
        min_z = min(min_z, last_z)
    return min_z, last_z


def test_soft_mortar_catches_an_overlap_without_nodes():
    nts_min, _ = _wide_facet_drop(nts=True)
    s2_min, s2_last = _wide_facet_drop(nts=False)
    assert nts_min < -0.5
    assert s2_min > -0.05
    assert s2_last > s2_min + 1.0e-4


def _matched_loaded(soft, nsteps=3000):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_quad()
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.mass(t, M, M, M)
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    extra = ("-soft", soft) if soft > 0.0 else ()
    ops.contact(1, 1, 2, "-mortar", "-epsN", KN_STIFF, *extra, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -25.0)
        ops.setNodeVel(t, 3, -0.5, "-commit")
    _explicit()
    return _run(5, 3, nsteps)


def test_soft_stabilizes_mortar():
    # 1000 steps: the four slave nodes are free masses joined only through the
    # facet, and over longer runs their repeated bounces accumulate the
    # engagement error of the non-smooth contact (rocking of a few cm by 3000
    # steps), which is not the instability tested here
    assert _matched_loaded(0.0, nsteps=1000)[0]
    div, maxabs = _matched_loaded(0.1, nsteps=1000)
    assert not div and maxabs < 0.05, maxabs


def _matched_static(soft):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master_quad()
    ops.uniaxialMaterial("Elastic", 90, 0.1)
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 0)
        ops.node(90 + t, float(x), float(y), 0.0)
        ops.fix(90 + t, 1, 1, 1)
        ops.element("zeroLength", 900 + t, 90 + t, t, "-mat", 90, "-dir", 3)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    extra = ("-soft", soft) if soft > 0.0 else ()
    ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e6, *extra, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -25.0)
    _static()
    return tuple(ops.nodeDisp(t, 3) for t in (5, 6, 7, 8))


def test_soft_mortar_is_inert_under_statics():
    assert _matched_static(0.1) == _matched_static(0.0)


# ----------------------------------------------------------------- edge-edge
EPS_EDGE = 1.0e8
ZC = 0.02
V0 = -2.0


def _bars(zc):
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    for t, x, y, z in [(5, -1.0, 0.0, zc), (6, 1.0, 0.0, zc),
                       (7, 1.0, 0.0, zc + 0.3), (8, -1.0, 0.0, zc + 0.3)]:
        ops.node(t, x, y, z)
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)


def _edge_impact(soft, nsteps=45):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _bars(ZC)
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)
        ops.mass(t, 1.0, 1.0, 1.0)
        ops.setNodeVel(t, 3, V0, "-commit")
    args = [1, 1, 2, "-mortar", "-epsN", EPS_EDGE, "-edgeedge", "-edgeBand", 0.15,
            "-outward", 0.0, 0.0, 1.0]
    if soft:
        args += ["-edgeSoft", 0.1]
    ops.contact(*args)
    _explicit()
    z = [ZC + ops.nodeDisp(5, 3)]
    for _ in range(nsteps):
        ops.analyze(1, DT)
        zi = ZC + ops.nodeDisp(5, 3)
        z.append(zi)
        if not math.isfinite(zi) or abs(zi) > 1.0e6:
            break
    return z


def test_edge_soft_rebounds_at_the_structural_step():
    z = _edge_impact(True)
    assert all(math.isfinite(zi) for zi in z)
    deepest = min(z)
    d_analytic = abs(V0) * math.sqrt(2.0 / (0.1 * 4.0 * 2.0 / (DT * DT)))
    assert deepest < 0.0
    assert -deepest < 5.0 * d_analytic, (deepest, d_analytic)
    assert z[-1] > deepest + 0.5 * (-deepest)
    assert ops.nodeVel(5, 3) > 0.0


def test_edge_stiff_penalty_diverges():
    z = _edge_impact(False)
    assert any((not math.isfinite(zi)) or abs(zi) > 1.0 for zi in z)


def _edge_static(soft):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _bars(-1.0e-4)
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, 500.0)
    for gnd, sn, x in [(305, 5, -1.0), (306, 6, 1.0)]:
        ops.node(gnd, x, 0.0, -1.0e-4)
        ops.fix(gnd, 1, 1, 1)
        ops.element("zeroLength", 900 + sn, gnd, sn, "-mat", 1, "-dir", 3)
    args = [1, 1, 2, "-mortar", "-epsN", 1.0e6, "-edgeedge", "-edgeBand", 0.15,
            "-outward", 0.0, 0.0, 1.0]
    if soft:
        args += ["-edgeSoft", 0.1]
    ops.contact(*args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 0.0, -200.0)
    _static()
    return tuple(ops.nodeDisp(t, 3) for t in (5, 6))


def test_edge_soft_is_inert_under_statics():
    assert _edge_static(True) == _edge_static(False)


# ------------------------------------------------------------------ refusals
@pytest.mark.parametrize("args", [
    (2, 1, 3, "-soft", 0.1),                                            # no base penalty
    (2, 1, 2, "-mortar", "-soft", 0.1),                                 # no base penalty
    (2, 1, 2, "-mortar", "-epsTie", 1.0e6, "-tie", "-soft", 0.1),       # tie
    (2, 1, 2, "-mortar", "-epsN", 1.0e6, "-edgeSoft", 0.1),             # no -edgeedge
    (2, 1, 3, 1.0e6, 0.0, 0.0, "-soft", -0.1),                          # non-positive scale
])
def test_soft_refusals(args):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.node(4 + t, float(x), float(y), 0.0)
    ops.node(9, 0.5, 0.5, 0.1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contactSurface(3, "-slave", 9)
    with pytest.raises(ops.OpenSeesError):
        ops.contact(*args)


def test_soft_on_2d_mortar_is_refused_at_analysis():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, x in ((1, -1.0), (2, 1.0)):
        ops.node(t, x, 0.0)
        ops.fix(t, 1, 1)
    for t, x in ((5, -1.0), (6, 1.0)):
        ops.node(t, x, -1.0e-8)
        ops.mass(t, M, M)
        ops.fix(t, 1, 0)
    ops.contactSurface(10, "-master", 2, 1, 2)
    ops.contactSurface(20, "-slave-segments", 2, 5, 6)
    ops.contact(1, 10, 20, "-mortar", "-epsN", 1.0e6, "-soft", 0.1, "-outward", 0.0, 1.0)
    _explicit()
    try:
        rc = ops.analyze(1, DT)
    except Exception:
        rc = -1
    assert rc < 0
