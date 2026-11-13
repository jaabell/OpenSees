"""Frictionless node-to-segment (NTS) penalty contact, 3D.

Checks (penalty contact has an O(load/kn) baseline, so the tolerances are
relative to the closed form):

- Static penetration of a slave pressed onto a fixed quad or tri facet equals
  P/kn to 1e-6; on a facet tilted by 30 degrees the free-z penetration equals
  P/(kn*nz^2) to 1e-5; reversing the facet winding changes nothing (1e-9).
- The NTS penetration matches a zeroLengthContactASDimplex node pair with the
  same kn and normal to 1e-6.
- Explicit impact on a fixed facet: restitution 0.9 < e <= 1.05 and peak
  penetration v*sqrt(m/kn) = 0.002 to 5 percent; a slave outside the facet
  passes through.
- Deformable stdBrick master: interface penetration P/kn (2 percent) in series
  with the brick shortening P*L/(E*A) (5 percent); block on block, P/(4*kn)
  per slave (5 percent).
- kn auto: penetration scales as 1/E (1 percent) and P/penetration equals
  0.10 * mean(n^T K_aa n) over the facet nodes of the owning element (2
  percent); a master face with no owning solid fails the analysis, an explicit
  kn on the same face runs.
- Bucket-sort broad phase: slave displacements on a 5x5 facet grid are
  bitwise identical to a single-bucket (brute force) search.
- -geomtan: identical result and iteration count on a flat facet, strictly
  fewer Newton iterations to the same root (1e-8) on a warped facet, identical
  result on ProfileSPD and FullGeneral.
- Rigid sphere (faceted master) pressed into an elastic stdBrick block whose
  top nodes are the slaves: contact radius equal to the Hertz value
  sqrt(R*delta) to 15 percent, centre-peaked pressure, load and radius
  growing with depth, radius ratio sqrt(2) to 0.30 and load ratio
  (delta2/delta1)^1.5 to 15 percent for a doubled depth.
- ladrunoContactForce returns the slave normal force (1e-6) while in contact
  and exactly 0.0 after the pair releases.
- Contact far from the origin (5e4 and 1e6) holds the slave and gives the
  same penetration as at the origin (1e-6).
- Handler: with no or a zero-force contact the result is bitwise identical to
  constraints Plain (explicit and Newmark); commits are counted; wipe clears
  the engine; nodes removed during the analysis are pruned and an emptied
  interaction is retired while the run continues.
- Nodes with 6 DOFs on either side of the interface: the summed contact force
  equals the applied load (1e-8), the column response equals an equalDOF
  bonded twin (1e-6) and the rotational DOFs follow their own analytic
  torsion.
- Refusals: duplicate contact tag, partial segment connectivity, missing node,
  and a non-trivial equationConstraint is reported as not enforced.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


KN = 1.0e6


# ---------------------------------------------------------------- helpers
def _quad_master(z=0.0, L=1.0, base=101, reverse=False):
    """Four fixed corner nodes of a flat square facet in the plane z, CCW."""
    h = L / 2.0
    tags = []
    for i, (x, y) in enumerate([(-h, -h), (h, -h), (h, h), (-h, h)]):
        t = base + i
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
        tags.append(t)
    return tags[::-1] if reverse else tags


def _static(tol=1.0e-12, iters=30, system="FullGeneral"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _explicit(system="Diagonal"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.integrator("NewmarkExplicit", 0.5)
    ops.algorithm("Linear")
    ops.analysis("Transient")


# ======================================================= handler basics
def _truss_chain(handler, integrator, with_contact, nsteps=40):
    """A 3D truss chain free along x, tip given an initial velocity."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    n = 5
    for i in range(n):
        ops.node(i + 1, float(i), 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    for i in range(2, n + 1):
        ops.fix(i, 0, 1, 1)
        ops.mass(i, 1.0, 1.0, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    for i in range(1, n):
        ops.element("Truss", i, i, i + 1, 1.0, 1)
    ops.setNodeVel(n, 1, 0.5, "-commit")
    if with_contact:
        # topology only: nothing ever comes into contact
        ops.contactSurface(1, "-master", 3, 1, 2, 3)
        ops.contactSurface(2, "-slave", 4, 5)
        ops.contact(1, 1, 2, 1.0e6, 0.0, 0.0)
    ops.constraints(handler)
    ops.numberer("Plain")
    if integrator == "Newmark":
        ops.system("FullGeneral")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1e-12, 25, 0)
        ops.algorithm("Newton")
    else:
        ops.system("Diagonal")
        ops.integrator("NewmarkExplicit", 0.5)
        ops.algorithm("Linear")
    ops.analysis("Transient")
    out = []
    for _ in range(nsteps):
        assert ops.analyze(1, 1.0e-3) == 0
        out.append(tuple(ops.nodeDisp(i + 1, 1) for i in range(n)))
    return out


@pytest.mark.parametrize("integrator", ["NewmarkExplicit", "Newmark"])
def test_handler_without_contact_matches_plain(integrator):
    ref = _truss_chain("Plain", integrator, with_contact=False)
    con = _truss_chain("LadrunoContact", integrator, with_contact=False)
    assert con == ref
    assert ops.ladrunoContactInfo()[:3] == [0, 0, 0]


@pytest.mark.parametrize("integrator", ["NewmarkExplicit", "Newmark"])
def test_zero_force_contact_matches_plain_and_commits(integrator):
    """A defined contact that never engages leaves the result bitwise
    unchanged, and the engine commit is called every step."""
    ref = _truss_chain("Plain", integrator, with_contact=False)
    con = _truss_chain("LadrunoContact", integrator, with_contact=True)
    assert con == ref
    info = ops.ladrunoContactInfo()      # [numContacts, numCommits, numReverts]
    assert info[0] == 1, info
    assert info[1] >= 40, info


def test_wipe_clears_contact_engine():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.node(3, 2.0, 0.0, 0.0)
    ops.contactSurface(1, "-master", 3, 1, 2, 3)
    ops.contactSurface(2, "-slave", 3)
    ops.contact(1, 1, 2, 1.0e6, 0.0, 0.0)
    assert ops.ladrunoContactInfo()[0] == 1
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    assert ops.ladrunoContactInfo()[:3] == [0, 0, 0]


def test_handler_rebuilds_after_element_removal():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i in range(4):
        ops.node(i + 1, float(i), 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    for i in range(2, 5):
        ops.fix(i, 0, 1, 1)
        ops.mass(i, 1.0, 1.0, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    for i in range(1, 4):
        ops.element("Truss", i, i, i + 1, 1.0, 1)
    _explicit()
    assert ops.analyze(10, 1.0e-3) == 0
    ops.remove("element", 3)
    assert ops.analyze(10, 1.0e-3) == 0


# ===================================================== bucket broad phase
def _grid_ntag(ix, iy, N):
    return 1 + ix * (N + 1) + iy


def _grid_model(N, cell=None, P=1.0e2):
    """Fixed N x N grid of unit facets at z=0, one slave just below the centre
    of each facet, free in z and loaded by -P. Returns the slave tags."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for ix in range(N + 1):
        for iy in range(N + 1):
            t = _grid_ntag(ix, iy, N)
            ops.node(t, float(ix), float(iy), 0.0)
            ops.fix(t, 1, 1, 1)
    master = []
    for cx in range(N):
        for cy in range(N):
            master += [_grid_ntag(cx, cy, N), _grid_ntag(cx + 1, cy, N),
                       _grid_ntag(cx + 1, cy + 1, N), _grid_ntag(cx, cy + 1, N)]
    slaves = []
    stag = 10000
    for cx in range(N):
        for cy in range(N):
            stag += 1
            ops.node(stag, cx + 0.5, cy + 0.5, -1.0e-8)
            ops.fix(stag, 1, 1, 0)
            slaves.append(stag)
    ops.contactSurface(10, "-master", 4, *master)
    ops.contactSurface(20, "-slave", *slaves)
    args = [1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0]
    if cell is not None:
        args += ["-cell", cell]
    ops.contact(*args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for s in slaves:
        ops.load(s, 0.0, 0.0, -P)
    _static(tol=1e-10, iters=50)
    assert ops.analyze(1) == 0
    return slaves


def test_bucket_sort_equals_brute_force():
    P = 1.0e2
    slaves = _grid_model(5, cell=None, P=P)
    d_bucket = [ops.nodeDisp(s, 3) for s in slaves]
    slaves2 = _grid_model(5, cell=1.0e12, P=P)     # one bucket: every facet
    d_brute = [ops.nodeDisp(s, 3) for s in slaves2]
    assert slaves == slaves2
    assert d_bucket == d_brute
    for d in d_bucket:
        assert abs(d - (-P / KN)) / (P / KN) < 0.02, d


@pytest.mark.parametrize("N", [1, 4])
def test_bucket_sort_every_slave_finds_its_facet(N):
    P = 1.0e2
    for s in _grid_model(N, P=P):
        dz = ops.nodeDisp(s, 3)
        assert abs(dz - (-P / KN)) / (P / KN) < 0.02, (s, dz)


# ========================================================= penalty NTS
def _slave_below(master_tags, outward=(0.0, 0.0, 1.0), P=1.0e3, tol=1e-12):
    """Slave 1e-8 below the facet at the origin, free z, loaded -P, static.
    Returns the penetration depth."""
    ops.node(1, 0.0, 0.0, -1.0e-8)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(10, "-master", len(master_tags), *master_tags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", *outward)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    _static(tol=tol, iters=40)
    assert ops.analyze(1) == 0
    return -(-1.0e-8 + ops.nodeDisp(1, 3))


def test_static_penetration_quad():
    P = 1.0e3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    pen = _slave_below(_quad_master(), P=P)
    assert pen > 0 and abs(pen - P / KN) / (P / KN) < 1e-6, pen


def test_static_penetration_tri():
    P = 5.0e2
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((101, 102, 103), ((-0.5, -0.5), (0.5, -0.5), (0.0, 0.6))):
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
    pen = _slave_below([101, 102, 103], P=P)
    assert pen > 0 and abs(pen - P / KN) / (P / KN) < 1e-6, pen


def test_static_penetration_oblique_facet():
    """Facet tilted 30 degrees about x, slave free in z only: the contact
    stiffness projected on z is kn*nz^2."""
    P = 1.0e3
    th = math.radians(30.0)
    c, s = math.cos(th), math.sin(th)
    n = (0.0, -s, c)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = []
    for i, (x, y0) in enumerate([(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]):
        ops.node(101 + i, x, y0 * c, y0 * s)
        ops.fix(101 + i, 1, 1, 1)
        mt.append(101 + i)
    pen = _slave_below(mt, outward=n, P=P, tol=1e-10)
    expect = P / (KN * n[2] * n[2])
    assert pen > 0 and abs(pen - expect) / expect < 1e-5, (pen, expect)


def test_winding_reversal_gives_identical_response():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    pen_fwd = _slave_below(_quad_master(reverse=False))
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    pen_rev = _slave_below(_quad_master(reverse=True))
    assert abs(pen_fwd - pen_rev) / pen_fwd < 1e-9, (pen_fwd, pen_rev)


def test_matches_zero_length_contact_asdimplex():
    P = 1.0e3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    pen_nts = _slave_below(_quad_master(), P=P)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(200, 0.0, 0.0, 0.0)
    ops.fix(200, 1, 1, 1)
    ops.node(1, 0.0, 0.0, -1.0e-8)
    ops.fix(1, 1, 1, 0)
    ops.element("zeroLengthContactASDimplex", 1, 200, 1, KN, 0.0, 0.0,
                "-orient", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-12, 30, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    pen_zlc = -(-1.0e-8 + ops.nodeDisp(1, 3))
    assert pen_nts > 0 and pen_zlc > 0
    assert abs(pen_nts - pen_zlc) / pen_zlc < 1e-6, (pen_nts, pen_zlc)


def _drop_on_quad(x0, y0, L, outward=True, gap0=0.01, v=2.0, nsteps=4000):
    """Unit mass gap0 above a fixed facet, thrown down at v, explicit.
    Returns (rc_ok, restitution, max penetration, final z)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _quad_master(L=L)
    ops.node(1, x0, y0, gap0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    if outward:
        ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contact(1, 10, 20, KN, 0.0, 0.0)
    ops.setNodeVel(1, 3, -v, "-commit")
    _explicit()
    dt = math.sqrt(1.0 / KN)
    max_pen, v_out = 0.0, 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return False, 0.0, max_pen, None
        z = gap0 + ops.nodeDisp(1, 3)
        vz = ops.nodeVel(1, 3)
        if z < 0:
            max_pen = max(max_pen, -z)
        if z > 0.5 * gap0 and vz > 0:
            v_out = vz
    return True, v_out / v, max_pen, gap0 + ops.nodeDisp(1, 3)


def test_explicit_impact_restitution():
    ok, e, max_pen, _ = _drop_on_quad(0.0, 0.0, 1.0)
    assert ok and 0.9 < e <= 1.05, e
    assert abs(max_pen - 0.002) / 0.002 < 0.05, max_pen


def test_explicit_impact_automatic_orientation():
    """Without -outward the orientation is taken from the slave side."""
    ok, e, max_pen, _ = _drop_on_quad(0.0, 0.0, 1.0, outward=False)
    assert ok and e > 0.9 and max_pen < 0.02, (e, max_pen)


def test_slave_outside_facet_passes_through():
    ok, _, _, zf = _drop_on_quad(5.0, 5.0, 0.2, nsteps=2000)
    assert ok and zf < -0.05, zf


# ===================================================== deformable master
def _cube(base, z0, L, E, nu, ele_tag, mat_tag, fix_bottom=False):
    """8-node stdBrick on [0,L]^2 x [z0,z0+L]. Returns (bottom4, top4)."""
    pts = [(0, 0, 0), (L, 0, 0), (L, L, 0), (0, L, 0),
           (0, 0, L), (L, 0, L), (L, L, L), (0, L, L)]
    tags = []
    for i, (x, y, z) in enumerate(pts):
        ops.node(base + i, float(x), float(y), z0 + z)
        tags.append(base + i)
    if fix_bottom:
        for t in tags[:4]:
            ops.fix(t, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", mat_tag, E, nu)
    ops.element("stdBrick", ele_tag, *tags, mat_tag)
    return tags[:4], tags[4:]


def test_slave_on_deformable_brick():
    E, nu, P, L = 2.0e4, 0.0, 1.0e2, 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _, top = _cube(1, 0.0, L, E, nu, 1, 1, fix_bottom=True)
    ops.node(99, L / 2, L / 2, L - 1.0e-8)
    ops.fix(99, 1, 1, 0)
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, 0.0, -P)
    _static(tol=1e-10, iters=50)
    assert ops.analyze(1) == 0
    pen = (L + ops.nodeDisp(top[0], 3)) - ((L - 1.0e-8) + ops.nodeDisp(99, 3))
    comp = -ops.nodeDisp(top[0], 3)
    assert abs(pen - P / KN) / (P / KN) < 0.02, pen
    expect = P * L / (E * L * L)
    assert abs(comp - expect) / expect < 0.05, comp


def test_block_on_block():
    E, nu, P, L = 5.0e4, 0.0, 2.0e2, 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _, top0 = _cube(1, 0.0, L, E, nu, 1, 1, fix_bottom=True)
    bot1, top1 = _cube(11, L - 1.0e-8, L, E, nu, 2, 2)
    for t in bot1 + top1:
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(10, "-master", 4, *top0)
    ops.contactSurface(20, "-slave", *bot1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top1:
        ops.load(t, 0.0, 0.0, -P / 4)
    _static(tol=1e-8, iters=60)
    assert ops.analyze(1) == 0
    pen = (L + ops.nodeDisp(top0[0], 3)) - ((L - 1.0e-8) + ops.nodeDisp(11, 3))
    expect = (P / 4) / KN
    assert pen > 0 and abs(pen - expect) / expect < 0.05, pen


# ================================================================ kn auto
def _cube_nodal_kzz(L, E, nu):
    """n^T K_aa n with n = +z for one node of an L-cube trilinear hexahedron
    (2x2x2 Gauss): sum over Gauss points of
    (bz^2 (lambda + 2 mu) + (bx^2 + by^2) mu) det J. Identical for every
    node by symmetry."""
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    g = 1.0 / math.sqrt(3.0)
    detJ = (L / 2.0) ** 3
    sx, sy, sz = 1.0, 1.0, 1.0           # node at (+1, +1, +1)
    k = 0.0
    for r in (-g, g):
        for s in (-g, g):
            for t in (-g, g):
                bx = 0.125 * sx * (1 + sy * s) * (1 + sz * t) * 2.0 / L
                by = 0.125 * sy * (1 + sx * r) * (1 + sz * t) * 2.0 / L
                bz = 0.125 * sz * (1 + sx * r) * (1 + sy * s) * 2.0 / L
                k += (bz * bz * (lam + 2 * mu) + (bx * bx + by * by) * mu) * detJ
    return k


def _auto_kn_press(L, E, nu, P):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _, top = _cube(1, 0.0, L, E, nu, 1, 1, fix_bottom=True)
    ops.node(99, L / 2, L / 2, L - 1.0e-8)
    ops.fix(99, 1, 1, 0)
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, "auto", "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, 0.0, -P)
    _static(tol=1e-10, iters=50)
    assert ops.analyze(1) == 0
    return (L + ops.nodeDisp(top[0], 3)) - ((L - 1.0e-8) + ops.nodeDisp(99, 3))


def test_kn_auto_converges_with_small_penetration():
    L, E = 1.0, 2.0e4
    P = 1.0e-3 * E * L * L
    pen = _auto_kn_press(L, E, 0.0, P)
    assert pen > 0.0
    assert pen / L < 0.10, pen


def test_kn_auto_scales_with_modulus():
    pen1 = _auto_kn_press(1.0, 1.0e4, 0.0, 1.0e2)
    pen3 = _auto_kn_press(1.0, 3.0e4, 0.0, 1.0e2)
    assert abs(pen1 / pen3 - 3.0) / 3.0 < 0.01, pen1 / pen3


def test_kn_auto_matches_element_stiffness_formula():
    L, E, nu, P = 1.0, 2.0e4, 0.0, 1.0e2
    kn_ref = 0.10 * _cube_nodal_kzz(L, E, nu)
    kn = P / _auto_kn_press(L, E, nu, P)
    assert abs(kn - kn_ref) / kn_ref < 0.02, (kn, kn_ref)


def _solidless_master(kn_args):
    """Fixed master face that belongs to no element; slave on a z spring."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    top = [1, 2, 3, 4]
    for t, (x, y) in zip(top, [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 1.0)
        ops.fix(t, 1, 1, 1)
    ops.node(99, 0.5, 0.5, 1.0 - 1.0e-8)
    ops.fix(99, 1, 1, 0)
    ops.node(98, 0.5, 0.5, 1.0 - 1.0e-8)
    ops.fix(98, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0e3)
    ops.element("zeroLength", 1, 98, 99, "-mat", 1, "-dir", 3)
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, *kn_args, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, 0.0, -1.0e2)


def test_kn_auto_without_owning_solid_fails():
    _solidless_master(["auto"])
    _static(tol=1e-10, iters=50)
    assert ops.analyze(1) != 0


def test_kn_auto_without_owning_solid_fails_explicit():
    _solidless_master(["auto"])
    ops.mass(99, 1.0, 1.0, 1.0)
    _explicit()
    assert ops.analyze(1, 1.0e-5) != 0


def test_explicit_kn_on_solidless_master_runs():
    _solidless_master([KN, 0.0, 0.0])
    _static(tol=1e-10, iters=50)
    assert ops.analyze(1) == 0
    free = -1.0e2 / 1.0e3                  # displacement with no contact
    assert abs(ops.nodeDisp(99, 3)) < 0.5 * abs(free)


# ======================================== -geomtan and the Hertz problem
KN_G = 2.0e6


def _warped_master(w, L=4.0, base=101):
    """Fixed quad with corners 0,2 raised by w and 1,3 lowered by w."""
    h = L / 2.0
    tags = []
    for i, (x, y, z) in enumerate([(-h, -h, w), (h, -h, -w), (h, h, w), (-h, h, -w)]):
        ops.node(base + i, x, y, z)
        ops.fix(base + i, 1, 1, 1)
        tags.append(base + i)
    return tags


def _slide_model(geom, w, solver="FullGeneral", Q=4.0e3, P=2.0e3, ka=2.0e4):
    """Slave on a (warped) facet, tangentially anchored by springs, loaded
    (Q, 0.7Q, -P). Returns (rc, iterations, ux, uy, uz)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _warped_master(w)
    ops.node(1, 0.0, 0.0, -1.0e-3)
    ops.node(2, 0.0, 0.0, -1.0e-3)
    ops.fix(2, 1, 1, 1)
    ops.node(3, 0.0, 0.0, -1.0e-3)
    ops.fix(3, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, ka)
    ops.element("zeroLength", 91, 2, 1, "-mat", 1, "-dir", 1)
    ops.element("zeroLength", 92, 3, 1, "-mat", 1, "-dir", 2)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    args = [1, 10, 20, KN_G, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0]
    if geom:
        args.append("-geomtan")
    ops.contact(*args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.7 * Q, -P)
    _static(tol=1e-11, iters=60, system=solver)
    rc = ops.analyze(1)
    return rc, ops.testIter(), ops.nodeDisp(1, 1), ops.nodeDisp(1, 2), ops.nodeDisp(1, 3)


def test_geomtan_flat_facet_unchanged():
    rc0, it0, x0, y0, z0 = _slide_model(False, 0.0)
    rc1, it1, x1, y1, z1 = _slide_model(True, 0.0)
    assert rc0 == 0 and rc1 == 0
    assert it0 == it1
    for a, b in ((x0, x1), (y0, y1), (z0, z1)):
        assert abs(b - a) <= 1e-14 + 1e-12 * abs(a)


def test_geomtan_curved_facet_fewer_iterations_same_root():
    rc_m, it_m, xm, ym, zm = _slide_model(False, 0.6)
    rc_g, it_g, xg, yg, zg = _slide_model(True, 0.6)
    assert rc_m == 0 and rc_g == 0
    assert it_g < it_m, (it_g, it_m)
    assert abs(xg - xm) < 1e-8 and abs(yg - ym) < 1e-8 and abs(zg - zm) < 1e-8


def test_geomtan_symmetric_solver_gives_same_result():
    rc_f, it_f, xf, yf, zf = _slide_model(True, 0.6, solver="FullGeneral")
    rc_s, it_s, xs, ys, zs = _slide_model(True, 0.6, solver="ProfileSPD")
    assert rc_f == 0 and rc_s == 0
    assert it_s == it_f
    assert abs(xs - xf) < 1e-9 and abs(ys - yf) < 1e-9 and abs(zs - zf) < 1e-9


def test_contact_force_query_equals_load():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _quad_master(L=20.0)
    ops.node(1, 0.0, 0.0, -1.0e-3)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, 1.0e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -250.0)
    _static()
    assert ops.analyze(1) == 0
    assert abs(ops.ladrunoContactForce(1) - 250.0) / 250.0 < 1e-6
    assert ops.ladrunoContactForce(99999) == 0.0


def test_contact_force_is_zero_after_release():
    """Slave on a soft z spring pressed into a fixed facet, then pulled back
    out: the reported normal force goes from the penalty value to 0.0."""
    KS, P, Q, g0 = 1.0e3, 1.0e2, 50.0, 1.0e-8
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _quad_master(L=2.0)
    ops.node(1, 0.0, 0.0, -g0)
    ops.fix(1, 1, 1, 0)
    ops.node(2, 0.0, 0.0, -g0)
    ops.fix(2, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 3)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    _static()
    assert ops.analyze(1) == 0
    # kn (g0 - u) + KS (-u) = P
    u = -(P - KN * g0) / (KN + KS)
    f_ref = KN * (g0 - u)
    f1 = ops.ladrunoContactForce(1)
    assert abs(f1 - f_ref) / f_ref < 1e-6, (f1, f_ref)

    ops.remove("loadPattern", 1)
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(1, 0.0, 0.0, Q)
    assert ops.analyze(1) == 0
    z = -g0 + ops.nodeDisp(1, 3)
    assert abs(z - (-g0 + Q / KS)) < 1e-9, z      # separated, spring only
    assert ops.ladrunoContactForce(1) == 0.0


def test_contact_force_is_zero_after_rebound():
    """Explicit drop on a fixed facet: once the slave has bounced clear the
    reported normal force is 0.0."""
    gap0 = 0.01
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _quad_master(L=1.0)
    ops.node(1, 0.0, 0.0, gap0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.setNodeVel(1, 3, -2.0, "-commit")
    _explicit()
    dt = math.sqrt(1.0 / KN)
    f_max = 0.0
    for _ in range(4000):
        assert ops.analyze(1, dt) == 0
        f_max = max(f_max, ops.ladrunoContactForce(1))
    assert f_max > 0.0
    assert gap0 + ops.nodeDisp(1, 3) > 0.5 * gap0
    assert ops.ladrunoContactForce(1) == 0.0


def _hertz_indent(nx=24, nz=8, W=2.0, H=1.0, R=20.0, delta=0.004, E=1000.0,
                  nu=0.3, kn=1.0e4, nm=41, geom=False, span=None):
    """Rigid sphere pressed into an elastic stdBrick block. The indenter is the
    master: fixed quad facets on the sphere surface, lowered by delta. The
    slaves are the top nodes of the block, so ladrunoContactForce gives the
    nodal contact force of the deformable body. Returns
    (rc, sorted [(r, pressure, force)] of the slaves in contact)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx, dz = W / nx, H / nz
    nid = {}
    t = 1
    for k in range(nz + 1):
        for j in range(nx + 1):
            for i in range(nx + 1):
                ops.node(t, -W / 2 + i * dx, -W / 2 + j * dx, -H + k * dz)
                nid[(i, j, k)] = t
                if k == 0:
                    ops.fix(t, 1, 1, 1)
                elif i in (0, nx) or j in (0, nx):
                    ops.fix(t, 1, 1, 0)
                t += 1
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    eid = 1
    for k in range(nz):
        for j in range(nx):
            for i in range(nx):
                ops.element("stdBrick", eid,
                            nid[(i, j, k)], nid[(i + 1, j, k)],
                            nid[(i + 1, j + 1, k)], nid[(i, j + 1, k)],
                            nid[(i, j, k + 1)], nid[(i + 1, j, k + 1)],
                            nid[(i + 1, j + 1, k + 1)], nid[(i, j + 1, k + 1)], 1)
                eid += 1
    if span is None:
        span = 2.0 * math.sqrt(R * delta)
    sdx = 2 * span / (nm - 1)
    mid = {}
    m = 400000
    for jj in range(nm):
        for ii in range(nm):
            x, y = -span + ii * sdx, -span + jj * sdx
            ops.node(m, x, y, (R - math.sqrt(R * R - x * x - y * y)) - delta)
            ops.fix(m, 1, 1, 1)
            mid[(ii, jj)] = m
            m += 1
    mt = []
    for jj in range(nm - 1):
        for ii in range(nm - 1):
            mt += [mid[(ii, jj)], mid[(ii + 1, jj)],
                   mid[(ii + 1, jj + 1)], mid[(ii, jj + 1)]]
    ops.contactSurface(10, "-master", 4, *mt)
    slaves, coord = [], {}
    for j in range(nx + 1):
        for i in range(nx + 1):
            x, y = -W / 2 + i * dx, -W / 2 + j * dx
            if abs(x) < span and abs(y) < span:
                s = nid[(i, j, nz)]
                slaves.append(s)
                coord[s] = math.hypot(x, y)
    ops.contactSurface(20, "-slave", *slaves)
    args = [1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, 0.0, -1.0]
    if geom:
        args.append("-geomtan")
    ops.contact(*args)
    ops.constraints("LadrunoContact")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-10, 60, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    rc = ops.analyze(1)
    data = []
    for s in slaves:
        f = ops.ladrunoContactForce(s)
        if f > 0.0:
            data.append((coord[s], f / (dx * dx), f))
    data.sort()
    return rc, data


def test_hertz_contact_radius():
    R, delta = 20.0, 0.004
    a_hertz = math.sqrt(R * delta)
    rc, data = _hertz_indent(R=R, delta=delta)
    assert rc == 0
    assert len(data) > 30, len(data)
    a_fe = max(r for r, _, _ in data)
    assert abs(a_fe - a_hertz) / a_hertz < 0.15, (a_fe, a_hertz)


def test_hertz_geomtan_same_patch():
    rc0, d0 = _hertz_indent(geom=False)
    rc1, d1 = _hertz_indent(geom=True)
    assert rc0 == 0 and rc1 == 0
    a0 = max(r for r, _, _ in d0)
    a1 = max(r for r, _, _ in d1)
    assert abs(a1 - a0) < 1e-9, (a0, a1)


def test_hertz_pressure_peaks_at_centre():
    rc, data = _hertz_indent()
    assert rc == 0
    a = max(r for r, _, _ in data)
    inner = [p for r, p, _ in data if r < a / 3.0]
    outer = [p for r, p, _ in data if r > 2.0 * a / 3.0]
    assert inner and outer
    assert sum(inner) / len(inner) > 1.5 * sum(outer) / len(outer)
    r_peak = max(data, key=lambda d: d[1])[0]
    assert r_peak < 0.5 * a, (r_peak, a)


def test_hertz_load_and_radius_grow_with_depth():
    R, d1v, d2v = 20.0, 0.003, 0.006
    span = 2.0 * math.sqrt(R * d2v)
    rc1, d1 = _hertz_indent(R=R, delta=d1v, span=span)
    rc2, d2 = _hertz_indent(R=R, delta=d2v, span=span)
    assert rc1 == 0 and rc2 == 0
    load1, load2 = sum(f for _, _, f in d1), sum(f for _, _, f in d2)
    a1, a2 = max(r for r, _, _ in d1), max(r for r, _, _ in d2)
    assert load2 > load1 and a2 > a1
    assert abs(a2 / a1 - math.sqrt(d2v / d1v)) < 0.30, a2 / a1
    # Hertz: P grows as delta^1.5
    assert abs(load2 / load1 - (d2v / d1v) ** 1.5) / (d2v / d1v) ** 1.5 < 0.15, load2 / load1


# ================================================ far-from-origin facets
TH_Z, TH_X = 0.37, 0.21


def _rot(x, y, z):
    x1 = x * math.cos(TH_Z) - y * math.sin(TH_Z)
    y1 = x * math.sin(TH_Z) + y * math.cos(TH_Z)
    return (x1, y1 * math.cos(TH_X) - z * math.sin(TH_X),
            y1 * math.sin(TH_X) + z * math.cos(TH_X))


def _press_rotated_strip(h, L, nseg=4, nsteps=3000):
    """Strip of nseg rotated facets of size h translated to (L, L, L); one
    slave over the strip pressed along -n by a constant load, explicit.
    Returns the minimum signed distance of the slave to the strip plane."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nx, ny, nz = _rot(0.0, 0.0, 1.0)
    mt = []
    for i in range(nseg + 1):
        xa, ya, za = _rot(i * h, -0.5 * h, 0.0)
        xb, yb, zb = _rot(i * h, 0.5 * h, 0.0)
        ops.node(1000 + i, xa + L, ya + L, za + L)
        ops.fix(1000 + i, 1, 1, 1)
        ops.node(2000 + i, xb + L, yb + L, zb + L)
        ops.fix(2000 + i, 1, 1, 1)
    for s in range(nseg):
        mt += [1000 + s, 1001 + s, 2001 + s, 2000 + s]
    gap0 = 1.0e-3 * h
    cx, cy, cz = _rot(1.5 * h, 0.0, 0.0)
    ops.node(1, cx + L + gap0 * nx, cy + L + gap0 * ny, cz + L + gap0 * nz)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", nx, ny, nz)
    F = 100.0
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, -F * nx, -F * ny, -F * nz)
    _explicit()
    dt = math.sqrt(1.0 / KN)
    min_d = gap0
    for k in range(nsteps):
        assert ops.analyze(1, dt) == 0, k
        d = gap0 + ops.nodeDisp(1, 1) * nx + ops.nodeDisp(1, 2) * ny + ops.nodeDisp(1, 3) * nz
        min_d = min(min_d, d)
    assert min_d < 0.0, min_d             # contact engaged
    return min_d


def test_contact_alive_far_from_origin():
    """500-unit facets at (5e4, 5e4, 5e4)."""
    assert _press_rotated_strip(500.0, 5.0e4) > -1.0


def test_penetration_independent_of_offset():
    d0 = _press_rotated_strip(500.0, 0.0)
    d1 = _press_rotated_strip(500.0, 5.0e4)
    assert abs(d1 - d0) < 1.0e-6, (d0, d1)


def test_unit_facets_at_one_million():
    assert _press_rotated_strip(1.0, 1.0e6) > -1.0e-2


# ============================================================ refusals
def _two_node_sets():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.node(5, 0.0, 0.0, 1.0e-4)
    ops.mass(5, 1.0, 1.0, 1.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)


def test_duplicate_contact_tag_refused():
    _two_node_sets()
    ops.contact(7, 1, 2, 1e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    with pytest.raises(ops.OpenSeesError):
        ops.contact(7, 1, 2, 1e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.contact(8, 1, 2, 1e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)


def test_partial_segment_connectivity_refused():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t in range(1, 8):
        ops.node(t, float(t), 0.0, 0.0)
    with pytest.raises(ops.OpenSeesError):
        ops.contactSurface(1, "-master", 4, 1, 2, 3, 4, 5, 6, 7)


def test_missing_node_refused_at_declaration():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
        ops.node(t, float(x), float(y), 0.0)
    with pytest.raises(ops.OpenSeesError):
        ops.contactSurface(1, "-master", 4, 1, 2, 3, 999)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)


def _contact_with_equation_constraint():
    _two_node_sets()
    ops.node(6, 5.0, 0.0, 0.0)
    ops.mass(6, 1.0, 1.0, 1.0)
    ops.node(7, 6.0, 0.0, 0.0)
    ops.mass(7, 1.0, 1.0, 1.0)
    ops.equationConstraint(6, 1, 1.0, 7, 1, -0.5)     # u6x = 0.5 u7x
    ops.contact(1, 1, 2, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, 0.0, -100.0)
    _explicit()


def test_nontrivial_equation_constraint_reported(capfd):
    _contact_with_equation_constraint()
    assert ops.analyze(1, 1.0e-4) == 0
    err = capfd.readouterr().err
    assert "EQ_Constraint at node 6" in err and "NOT ENFORCED" in err, err[-300:]


def test_nontrivial_equation_constraint_ignored_contact_works():
    _contact_with_equation_constraint()
    ops.load(7, 50.0, 0.0, 0.0)
    for _ in range(200):
        assert ops.analyze(1, 1.0e-3) == 0
    u6, u7 = ops.nodeDisp(6, 1), ops.nodeDisp(7, 1)
    assert abs(u7) > 1.0e-6
    assert abs(u6 - 0.5 * u7) > 1.0e-9
    assert abs(1.0e-4 + ops.nodeDisp(5, 3) + 100.0 / KN) < 5.0e-4


# ============================================== nodes removed at run time
def _brick_master():
    """Fixed-bottom stdBrick; its top face (5, 6, 7, 8) is the master."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.0)
    pts = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
           (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    for t, (x, y, z) in zip(range(1, 9), pts):
        ops.node(t, float(x), float(y), float(z))
    ops.element("stdBrick", 101, *range(1, 9), 1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    for n in (5, 6, 7, 8):
        ops.fix(n, 1, 1, 0)
        ops.mass(n, 1.0, 1.0, 1.0)


def _slave_set(grounded):
    """Four slaves just below the master face. Grounded slaves have a z
    spring to an anchor; the others are fixed and carry no element or load,
    so they can be removed."""
    xy = {11: (0, 0), 12: (1, 0), 13: (1, 1), 14: (0, 1)}
    ops.uniaxialMaterial("Elastic", 90, 1.0e4)
    for t in sorted(xy):
        ops.node(t, float(xy[t][0]), float(xy[t][1]), 0.995)
        ops.fix(t, 1, 1, 0 if t in grounded else 1)
        ops.mass(t, 1.0, 1.0, 1.0)
    for t in grounded:
        ops.node(900 + t, float(xy[t][0]), float(xy[t][1]), 0.995)
        ops.fix(900 + t, 1, 1, 1)
        ops.element("zeroLength", 900 + t, 900 + t, t, "-mat", 90, "-dir", 3)
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave", *sorted(xy))
    ops.contact(1, 1, 2, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)


def test_removed_slave_node_is_pruned():
    _brick_master()
    grounded = (11, 12, 13)
    _slave_set(grounded)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in grounded:
        ops.load(n, 0.0, 0.0, -10.0)
    _explicit()
    assert ops.analyze(1, 1.0e-5) == 0
    ops.remove("node", 14)
    assert ops.analyze(1, 1.0e-5) == 0
    for n in grounded:                    # the remaining pairs still act
        assert ops.ladrunoContactForce(n) > 0.0, n
    assert ops.ladrunoContactInfo()[0] == 1


def test_contact_with_emptied_surface_is_retired(capfd):
    _brick_master()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.0, 0.0, -10.0)
    _slave_set(grounded=())
    _explicit()
    assert ops.analyze(1, 1.0e-5) == 0
    for n in (11, 12, 13, 14):
        ops.remove("node", n)
    assert ops.analyze(1, 1.0e-5) == 0
    out = capfd.readouterr()
    assert "retir" in (out.out + out.err).lower(), (out.out + out.err)[-2000:]


# ============================================ nodes with rotational DOFs
COL_E, COL_A, COL_G, COL_J, COL_I = 1.0e4, 0.25, 4.0e3, 0.02, 0.01
COL_L = 4.0
COL_P = 100.0
COL_M = 2.0                               # torque on each column top
COL_KN = 1.0e3 * COL_E * COL_A / COL_L
COL_GAP0 = 1.0e-9
CORNERS = [(0, 0), (1, 0), (1, 1), (0, 1)]


def _columns():
    """Four elastic 3D beam columns (ndf 6), bases fixed, tops held
    laterally, each top loaded by a torque COL_M. Returns the top tags."""
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.geomTransf("Linear", 1, 1.0, 0.0, 0.0)
    tops = []
    for i, (x, y) in enumerate(CORNERS):
        b, t = 1 + i, 11 + i
        ops.node(b, float(x), float(y), 0.0)
        ops.fix(b, 1, 1, 1, 1, 1, 1)
        ops.node(t, float(x), float(y), COL_L)
        ops.fix(t, 1, 1, 0, 0, 0, 0)
        ops.element("elasticBeamColumn", 1 + i, b, t, COL_A, COL_E, COL_G,
                    COL_J, COL_I, COL_I, 1)
        tops.append(t)
    return tops


def _platen(z, fix_lateral):
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    tags = []
    for i, (x, y) in enumerate(CORNERS):
        ops.node(9001 + i, float(x), float(y), z)
        if fix_lateral:
            ops.fix(9001 + i, 1, 1, 0)
        tags.append(9001 + i)
    return tags


def _column_loads(tops, platen):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in platen:
        ops.load(t, 0.0, 0.0, -COL_P / 4.0)
    for t in tops:
        ops.load(t, 0.0, 0.0, 0.0, 0.0, 0.0, COL_M)


def _column_contact(reverse):
    ops.wipe()
    tops = _columns()
    platen = _platen(COL_L - COL_GAP0, fix_lateral=True)
    if reverse:
        ops.contactSurface(10, "-master", 4, *platen)
        ops.contactSurface(20, "-slave", *tops)
        ops.contact(1, 10, 20, COL_KN, 0.0, 0.0, "-outward", 0.0, 0.0, -1.0)
    else:
        ops.contactSurface(10, "-master", 4, *tops)
        ops.contactSurface(20, "-slave", *platen)
        ops.contact(1, 10, 20, COL_KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    _column_loads(tops, platen)
    _static(iters=50)
    rc = ops.analyze(1)
    slaves = tops if reverse else platen
    return rc, tops, platen, [ops.ladrunoContactForce(t) for t in slaves]


def _column_bonded():
    ops.wipe()
    tops = _columns()
    platen = _platen(COL_L, fix_lateral=False)
    for t, p in zip(tops, platen):
        ops.equalDOF(t, p, 1, 2, 3)
    _column_loads(tops, platen)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return {t: [ops.nodeDisp(t, d) for d in range(1, 7)] for t in tops}


@pytest.mark.parametrize("reverse", [False, True])
def test_six_dof_nodes_in_contact(reverse):
    """Platen (ndf 3) on four beam-column tops (ndf 6); with reverse=True the
    column tops are the slaves. Contact acts on the translations only."""
    rc, tops, platen, fn = _column_contact(reverse)
    assert rc == 0
    assert abs(sum(fn) - COL_P) <= 1.0e-8 * COL_P, fn
    assert min(fn) > 0.0, fn
    d = [ops.nodeDisp(t, 3) - ops.nodeDisp(p, 3) for t, p in zip(tops, platen)]
    dref = COL_P / (4.0 * COL_KN) - COL_GAP0
    assert max(abs(x - dref) for x in d) <= 1.0e-6 * dref, (d, dref)
    shortening = (COL_P / 4.0) * COL_L / (COL_E * COL_A)
    twist = COL_M * COL_L / (COL_G * COL_J)
    got = {t: [ops.nodeDisp(t, k) for k in range(1, 7)] for t in tops}
    for t in tops:
        assert abs(-got[t][2] - shortening) <= 1.0e-6 * shortening, got[t]
        assert abs(got[t][5] - twist) <= 1.0e-9 * twist, got[t]
        assert abs(got[t][3]) < 1.0e-12 and abs(got[t][4]) < 1.0e-12, got[t]
    bonded = _column_bonded()
    for t in tops:
        assert abs(got[t][2] - bonded[t][2]) <= 1.0e-6 * shortening
        assert abs(got[t][5] - bonded[t][5]) <= 1.0e-9 * twist
