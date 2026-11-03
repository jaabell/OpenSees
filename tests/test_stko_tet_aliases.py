"""Interim aliases so models exported by STKO run unchanged.

  element ASDTet   $tag $n1..$n4  $matTag   -> FourNodeTetrahedron, same node order
  element ASDTet10 $tag $n1..$n10 $matTag   -> TenNodeTetrahedron, nodes 9 and 10 swapped
  element ASDEmbeddedNodeElement $tag $C $R1..$R10 ...  -> 10-node tet host, same order
                                                           as ASDTet10

STKO writes the tet10 mid-edge nodes as 5=(1,2) 6=(2,3) 7=(3,1) 8=(1,4) 9=(2,4) 10=(3,4);
TenNodeTetrahedron (gmsh) has 9=(3,4) 10=(2,4).
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


CORNERS = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0),
           3: (0.0, 1.0, 0.0), 4: (0.0, 0.0, 1.0)}
STKO_EDGES = [(1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4)]
STKO_TO_GMSH = [1, 2, 3, 4, 5, 6, 7, 8, 10, 9]


def field(x, y, z):
    # quadratic, and different on edges (2,4) and (3,4), so a wrong swap shows
    return (x * x + 2.0 * y - z, 0.5 * x * y, y * z - 0.3 * x)


def build_nodes(n=10):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    xyz = dict(CORNERS)
    if n == 10:
        for k, (a, b) in enumerate(STKO_EDGES):
            xyz[5 + k] = tuple((CORNERS[a][d] + CORNERS[b][d]) / 2.0 for d in range(3))
    for tag, x in xyz.items():
        ops.node(tag, *x)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    return xyz


def solve_prescribed(xyz, tags):
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for tag in tags:
        u = field(*xyz[tag])
        for d in range(3):
            ops.sp(tag, d + 1, u[d])
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormUnbalance", 1.0e-8, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.reactions()


def tet10_reactions(cmd, conn):
    xyz = build_nodes()
    ops.element(cmd, 1, *conn, 1)
    solve_prescribed(xyz, range(1, 11))
    return [ops.nodeReaction(t) for t in range(1, 11)]


def test_asdtet10_matches_tennodetetrahedron_with_swap():
    stko = list(range(1, 11))
    gmsh = [stko[i - 1] for i in STKO_TO_GMSH]
    r_alias = tet10_reactions("ASDTet10", stko)
    r_ref = tet10_reactions("TenNodeTetrahedron", gmsh)
    r_noswap = tet10_reactions("TenNodeTetrahedron", stko)
    for a, b in zip(r_alias, r_ref):
        assert a == pytest.approx(b, rel=1e-10, abs=1e-10)
    # the swap matters for this field
    assert max(abs(a - b) for ra, rb in zip(r_noswap, r_ref) for a, b in zip(ra, rb)) > 1e-3


def test_asdtet_matches_fournodetetrahedron():
    res = []
    for cmd in ("ASDTet", "FourNodeTetrahedron"):
        build_nodes(4)
        ops.element(cmd, 1, 1, 2, 3, 4, 1)
        solve_prescribed(CORNERS, range(1, 5))
        res.append([ops.nodeReaction(t) for t in range(1, 5)])
    for a, b in zip(*res):
        assert a == pytest.approx(b, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("cmd,n", [("ASDTet", 4), ("ASDTet10", 10)])
def test_corotational_rejected(cmd, n):
    build_nodes()
    with pytest.raises(Exception):
        ops.element(cmd, 1, *range(1, n + 1), 1, "-corotational")


def embedded_disp(retained, *opts):
    xyz = build_nodes(len(retained))
    c = (0.2, 0.1, 0.15)
    ops.node(11, *c)
    ops.element("ASDEmbeddedNodeElement", 1, 11, *retained, "-K", 1.0e14, *opts)
    solve_prescribed(xyz, retained)
    return [ops.nodeDisp(11, d + 1) for d in range(3)], field(*c)


@pytest.mark.parametrize("opts", [(), ("-shape", "tet")])
def test_embedded_tet10_stko_order(opts):
    u, exact = embedded_disp(list(range(1, 11)), *opts)
    assert u == pytest.approx(exact, rel=1e-8, abs=1e-10)


def test_embedded_tet4_shape_tet():
    u, _ = embedded_disp([1, 2, 3, 4], "-shape", "tet")
    # linear interpolation of the corner values
    L = (0.55, 0.2, 0.1, 0.15)
    lin = [sum(L[i] * field(*CORNERS[i + 1])[d] for i in range(4)) for d in range(3)]
    assert u == pytest.approx(lin, rel=1e-8, abs=1e-10)


@pytest.mark.parametrize("retained,opts", [
    (list(range(1, 11)), ("-shape", "quad")),
    (list(range(1, 11)), ("-rot",)),
    (list(range(1, 8)), ()),
])
def test_embedded_bad_input_rejected(retained, opts):
    build_nodes()
    ops.node(11, 0.2, 0.1, 0.15)
    with pytest.raises(Exception):
        ops.element("ASDEmbeddedNodeElement", 1, 11, *retained, "-K", 1.0e14, *opts)
