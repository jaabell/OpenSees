"""ASDEmbeddedNodeElement whole-host families (tri6, quad4/8/9, hexa8/20/27).

The embedded node must follow the host's own field, u_C = sum_a N_a(xi*) u_a,
over all host nodes, with xi* from the inverse isoparametric map. A field is
prescribed on the host nodes and the embedded node is free:

* a complete quadratic field on an affine host is reproduced exactly by the
  quadratic families;
* a linear field is reproduced exactly by every isoparametric host, also on a
  distorted (curved) host, which checks the inverse map;
* tri/quad hosts also work as faces in 3D.

Node orders are STKO / MpcCore's (see the element source).
"""
import pytest

ops = pytest.importorskip("opensees")

QX = [-1, 1, 1, -1, 0, 1, 0, -1, 0]
QY = [-1, -1, 1, 1, -1, 0, 1, 0, 0]
HX = [-1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0]
HY = [-1, -1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, -1, 1, 1, -1, 0, 1, 0, 0, -1, 0, 1, 0, 0, 0]
HZ = [-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, 0, 0, 0, 0, 1, 0]
T6 = [(0, 0), (1, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0, 0.5)]


def natural_nodes(fam):
    if fam == "tri6":
        return T6
    if fam.startswith("quad"):
        n = int(fam[4:])
        return [(QX[i], QY[i]) for i in range(n)]
    n = int(fam[4:])
    return [(HX[i], HY[i], HZ[i]) for i in range(n)]


def geom(xi, ndm, distort):
    a, b = xi[0], xi[1]
    c = xi[2] if len(xi) > 2 else 0.0
    x = 1.3 * a + 0.4 * b + 0.2 * c + 2.0
    y = -0.3 * a + 0.9 * b + 0.1 * c + 1.0
    z = 0.2 * a - 0.25 * b + 1.1 * c + 0.5
    if distort:
        x += 0.15 * a * b + 0.1 * b * b
        y += 0.12 * a * a - 0.1 * a * c
        z += 0.1 * b * c
    return (x, y) if ndm == 2 else (x, y, z)


def field(p, deg):
    X, Y = p[0], p[1]
    Z = p[2] if len(p) > 2 else 0.0
    u = 1e-3 * (1 + 2 * X - Y + 0.5 * Z)
    v = 1e-3 * (-2 + 0.3 * X + Y - Z)
    w = 1e-3 * (0.5 - X + 0.2 * Y + 3 * Z)
    if deg == 2:
        u += 1e-3 * (0.7 * X * X - 0.4 * X * Y + 0.3 * Y * Z)
        v += 1e-3 * (0.2 * Y * Y + 0.5 * X * Z - 0.6 * Z * Z)
        w += 1e-3 * (-0.3 * X * Y + 0.8 * Z * Z + 0.1 * X * X)
    return (u, v, w)


def embedded_error(fam, ndm, xi_p, distort, deg, shape):
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndm)
    nat = natural_nodes(fam)
    p = geom(xi_p, ndm, distort)
    ops.node(1, *p)
    rn = []
    for k, xi in enumerate(nat):
        ops.node(k + 2, *geom(xi, ndm, distort))
        rn.append(k + 2)
    args = ["ASDEmbeddedNodeElement", 1, 1, *rn, "-K", 1e14]
    if shape:
        args += ["-shape", shape]
    ops.element(*args)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in rn:
        u = field(ops.nodeCoord(n), deg)
        for d in range(ndm):
            ops.sp(n, d + 1, u[d])
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 5)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ue = field(p, deg)
    err = max(abs(ops.nodeDisp(1, d + 1) - ue[d]) for d in range(ndm))
    return err / max(abs(ue[d]) for d in range(ndm))


QUADRATIC = ("tri6", "quad8", "quad9", "hexa20", "hexa27")

CASES = [
    # family, ndm, -shape, natural coordinates of the embedded node
    ("tri6", 2, "", (0.2, 0.3)),
    ("quad4", 2, "quad", (0.3, -0.6)),
    ("quad8", 2, "quad", (0.3, -0.6)),
    ("quad9", 2, "", (0.55, 0.8)),
    ("hexa8", 3, "hexa", (0.3, -0.6, 0.45)),
    ("hexa20", 3, "", (0.95, -0.2, -0.5)),
    ("hexa27", 3, "hexa", (0.1, 0.7, -0.9)),
    ("tri6", 3, "tri", (0.25, 0.15)),
    ("quad8", 3, "quad", (0.3, -0.6)),
    ("quad9", 3, "quad", (-0.4, 0.2)),
]


@pytest.mark.parametrize("fam,ndm,shape,xi", CASES)
def test_linear_field_on_distorted_host(fam, ndm, shape, xi):
    assert embedded_error(fam, ndm, xi, True, 1, shape) < 1e-8


@pytest.mark.parametrize("fam,ndm,shape,xi", [c for c in CASES if c[0] in QUADRATIC])
def test_quadratic_field_on_affine_host(fam, ndm, shape, xi):
    assert embedded_error(fam, ndm, xi, False, 2, shape) < 1e-8


@pytest.mark.parametrize("ndm,nret,shape", [
    (2, 8, "hexa"),     # a hexahedron needs 3D
    (3, 9, "hexa"),     # wrong count for the shape
    (3, 5, ""),         # no host has 5 nodes
    (2, 4, "blob"),     # unknown shape
])
def test_refused_hosts(ndm, nret, shape):
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndm)
    for i in range(1, nret + 2):
        ops.node(i, *([0.1 * i, 0.13 * (i % 7), 0.17 * (i % 5)][:ndm]))
    args = ["ASDEmbeddedNodeElement", 1, *range(1, nret + 2)]
    if shape:
        args += ["-shape", shape]
    with pytest.raises(Exception):
        ops.element(*args)
