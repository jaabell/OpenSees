"""2D node-to-segment contact against its 3D twin (a one-layer stdBrick slice).

A plane-strain quad block (E = 1e4, nu = 0.25) rests on a wider quad slab through a 2D
node-to-segment contact, loaded by P on its top with an asymmetric split (1/3, 2/3), x fixed on
every node. The 3D twin is the same pair as one layer of stdBrick elements with u_z fixed on every
node (so the solution does not depend on z and the brick reproduces the plane-strain quad): the
slab overhangs the block in x and z, the block is one layer of thickness 1, its four bottom nodes
are the slaves of a quad master facet, and kn is halved because two 3D slave nodes stand for one
2D slave node.

Checks: the summed contact force equals P on both (1e-9), and the 2D and 3D contact forces and
gaps at the two slave stations agree (1e-5).

The file fails on the base with an unknown command.
"""
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E = 1.0e4
NU = 0.25
TZ = 1.0
W = 1.0
H = 1.0
KN = 1.0e6
SEED = 1.0e-8
P = 1.0e2
LOADW = (1.0 / 3.0, 2.0 / 3.0)
XS = (0.0, W)
MX0, MX1 = -1.5 * W, 2.5 * W
MZ0, MZ1 = -1.5 * TZ, 2.5 * TZ
THICK_SLAB = MZ1 - MZ0


def _static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _run_2d():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    for t, x in zip((1, 2), (MX0, MX1)):
        ops.node(t, x, 0.0)
        ops.fix(t, 1, 1)
    for t, x in zip((3, 4), (MX0, MX1)):
        ops.node(t, x, H)
        ops.fix(t, 1, 0)
    ops.element("quad", 1, 1, 2, 4, 3, THICK_SLAB, "PlaneStrain", 1)
    for t, x in zip((5, 6), XS):
        ops.node(t, x, H - SEED)
        ops.fix(t, 1, 0)
    for t, x in zip((7, 8), XS):
        ops.node(t, x, 2.0 * H)
        ops.fix(t, 1, 0)
    ops.element("quad", 2, 5, 6, 8, 7, 1.0, "PlaneStrain", 1)
    ops.contactSurface(10, "-master", 2, 3, 4)
    ops.contactSurface(20, "-slave", 5, 6)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, wgt in zip((7, 8), LOADW):
        ops.load(t, 0.0, -P * wgt)
    _static()
    out = []
    for t, x in zip((5, 6), XS):
        s = (x - MX0) / (MX1 - MX0)
        um = (1.0 - s) * ops.nodeDisp(3, 2) + s * ops.nodeDisp(4, 2)
        out.append((ops.ladrunoContactForce(t), ops.nodeDisp(t, 2) - um - SEED))
    return out


def _run_3d():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    slab = tuple(range(11, 19))
    block = tuple(range(21, 29))
    slab_pts = [(MX0, 0.0, MZ0), (MX1, 0.0, MZ0), (MX1, H, MZ0), (MX0, H, MZ0),
                (MX0, 0.0, MZ1), (MX1, 0.0, MZ1), (MX1, H, MZ1), (MX0, H, MZ1)]
    for t, (x, y, z) in zip(slab, slab_pts):
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1 if y == 0.0 else 0, 1)
    ops.element("stdBrick", 1, *slab, 1)
    y0, y1 = H - SEED, 2.0 * H
    blk_pts = [(XS[0], y0, 0.0), (XS[1], y0, 0.0), (XS[1], y1, 0.0), (XS[0], y1, 0.0),
               (XS[0], y0, TZ), (XS[1], y0, TZ), (XS[1], y1, TZ), (XS[0], y1, TZ)]
    for t, (x, y, z) in zip(block, blk_pts):
        ops.node(t, x, y, z)
        ops.fix(t, 1, 0, 1)
    ops.element("stdBrick", 2, *block, 1)
    mast = (slab[3], slab[2], slab[6], slab[7])          # the y = H face
    slave = ((block[0], block[4]), (block[1], block[5]))
    top = ((block[3], block[7]), (block[2], block[6]))
    ops.contactSurface(10, "-master", 4, *mast)
    ops.contactSurface(20, "-slave", *[t for pair in slave for t in pair])
    ops.contact(1, 10, 20, KN / 2.0, 0.0, 0.0, "-outward", 0.0, 1.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for pair, wgt in zip(top, LOADW):
        for t in pair:
            ops.load(t, 0.0, -P * wgt / 2.0, 0.0)
    _static()
    out = []
    for pair, x in zip(slave, XS):
        f = sum(ops.ladrunoContactForce(t) for t in pair)
        s = (x - MX0) / (MX1 - MX0)
        gaps = []
        for t in pair:
            z = ops.nodeCoord(t, 3)
            e = (z - MZ0) / (MZ1 - MZ0)
            wts = ((1.0 - s) * (1.0 - e), s * (1.0 - e), s * e, (1.0 - s) * e)
            um = sum(wi * ops.nodeDisp(ti, 2) for wi, ti in zip(wts, mast))
            gaps.append(ops.nodeDisp(t, 2) - um - SEED)
        out.append((f, 0.5 * (gaps[0] + gaps[1])))
    return out


def test_slice_twin():
    r2 = _run_2d()
    r3 = _run_3d()
    assert abs(sum(f for f, _ in r2) - P) < 1.0e-9 * P
    assert abs(sum(f for f, _ in r3) - P) < 1.0e-9 * P
    for (f2, g2), (f3, g3) in zip(r2, r3):
        assert abs(f2 - f3) < 1.0e-5 * P, (r2, r3)
        assert abs(g2 - g3) < 1.0e-5 * abs(g2), (r2, r3)
