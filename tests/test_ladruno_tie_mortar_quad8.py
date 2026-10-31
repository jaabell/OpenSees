"""Tests for LadrunoTie integral mortar ties with quadratic (serendipity)
facets: a 20NodeBrick face (quad8) tied onto a stdBrick face (quad4), and
quad8 / tri6 master facets.

The kernel keeps the geometry on the corner polygon and evaluates the full
serendipity basis at the mapped Gauss points; this is exact for flat facets
with the midside nodes at the edge midpoints, which is checked.

Mesh: one stdBrick hexahedron x in [0, 1] (master face quad4 at x = 1) and one
20NodeBrick x in [1, 2] whose x = 1 face is the quad8 slave facet. The corner
nodes coincide with the master ones; the midside nodes have no master
counterpart.

- Patch: u_x = gamma z imposed on every node except the tied ones reaches the
  tied quad8 face exactly (corners and midsides), standard and dual basis.
- Split column: interface at U/2 and master stress E U / 2, for two DOF
  numberings.
- Master basis: a linear field prescribed on tri6 or quad8 master facets is
  transferred exactly to a quad4 slave held by soft springs.
- Refusals: tri6 slave facets, a curved quad8 (midside off its edge midpoint),
  quadratic facets in the collocation mode, a quad8 slave off the master
  plane, a two-facet slave tilted about a ridge on the master plane, and a
  doubly-listed master facet. A quad4 slave on tri6 master facets is accepted.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# a refused LadrunoTie command raises the interpreter error (an unknown
# command would raise AttributeError instead and fail these tests)
OPS_ERROR = getattr(ops, "OpenSeesError", RuntimeError)


_MASTER = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
           5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
MASTER_FACET = [2, 3, 7, 6]

# 20-node hexahedron x in [1, 2]: corners 1-8, bottom edge midsides 9-12, top
# edge midsides 13-16, vertical edge midsides 17-20 (tags 101..120)
_S = {101: (1, 0, 0), 102: (2, 0, 0), 103: (2, 1, 0), 104: (1, 1, 0),
      105: (1, 0, 1), 106: (2, 0, 1), 107: (2, 1, 1), 108: (1, 1, 1),
      109: (1.5, 0, 0), 110: (2, 0.5, 0), 111: (1.5, 1, 0), 112: (1, 0.5, 0),
      113: (1.5, 0, 1), 114: (2, 0.5, 1), 115: (1.5, 1, 1), 116: (1, 0.5, 1),
      117: (1, 0, 0.5), 118: (2, 0, 0.5), 119: (2, 1, 0.5), 120: (1, 1, 0.5)}

# the slave x = 1 face as a quad8 facet: corners first, then the midsides in edge order
QUAD8_FACET = [101, 104, 108, 105, 112, 120, 116, 117]
SLAVE_IFACE = set(QUAD8_FACET)
_ALL = {**_MASTER, **_S}


def _mesh(E=1000.0, nu=0.0, dual=True):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, 0.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("20NodeBrick", 2, *range(101, 121), 1)
    args = ["-mortar", "-slaveFacets", 8, 1, *QUAD8_FACET,
            "-masterFacets", 4, 1, *MASTER_FACET]
    if dual:
        args.append("-dual")
    ops.LadrunoTie(*args)


def _static_solve(handler="Lagrange", numberer="Plain"):
    ops.constraints(handler)
    ops.numberer(numberer)
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _gp_stresses(ele):
    ops.eleResponse(ele, "forces")
    s = list(ops.eleResponse(ele, "stresses"))
    return [s[i * 6:(i + 1) * 6] for i in range(len(s) // 6)]


@pytest.mark.parametrize("dual", [True, False])
def test_quad8_patch_linear_exact(dual):
    E, gamma = 1000.0, 1e-3
    _mesh(E=E, nu=0.0, dual=dual)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_IFACE:
            continue
        ops.fix(t, 0, 1, 1)
        ops.sp(t, 1, gamma * float(z))
    _static_solve("Lagrange")
    for t in QUAD8_FACET:
        z = _ALL[t][2]
        got = ops.nodeDisp(t, 1)
        assert abs(got - gamma * z) <= 1e-9, f"tied node {t} (z = {z}): u_x = {got} != {gamma * z}"
    gps = _gp_stresses(1)
    ref = gps[0]
    target = E * gamma / 2.0
    atol = 1e-7 * E * gamma
    for v in gps:
        for si, ri in zip(v, ref):
            assert abs(si - ri) <= atol, f"non-constant master stress: {v} vs {ref}"
    assert abs(max(abs(c) for c in ref) - target) <= 1e-5 * target + atol


@pytest.mark.parametrize("numberer", ["Plain", "RCM"])
def test_quad8_split_column(numberer):
    """Series column through the quad8 tie. The tied nodes get y and z from the
    tie only: fixing them as well would add linearly dependent Lagrange rows."""
    E, U = 1000.0, 0.02
    _mesh(E=E, nu=0.0, dual=True)
    x2_face = [t for t, (x, y, z) in _S.items() if abs(x - 2.0) < 1e-12]
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t not in SLAVE_IFACE:
            ops.fix(t, 0, 1, 1)
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)
        elif t in x2_face:
            ops.sp(t, 1, U)
    _static_solve("Lagrange", numberer)
    assert abs(ops.nodeDisp(2, 1) - 0.5 * U) <= 1e-6 * U
    for t in QUAD8_FACET:
        assert abs(ops.nodeDisp(t, 1) - 0.5 * U) <= 1e-6 * U, \
            f"tied node {t}: {ops.nodeDisp(t, 1)} != {0.5 * U}"
    got = _gp_stresses(1)[0]
    assert abs(got[0] - E * U / 2.0) <= 1e-5 * E * U, f"master sigma_xx {got[0]}"


def _spring_rig_solve(master_nodes, master_facets, npsM, field):
    """Masters prescribed to a linear u_x field; each slave node of one quad4
    facet held by a soft spring to the ground, so u_s = P u_m."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in master_nodes.items():
        ops.node(t, float(x), float(y), float(z))
        ops.fix(t, 0, 1, 1)
    sv = {501: (1, 0.2, 0.2), 502: (1, 0.8, 0.2), 503: (1, 0.8, 0.8), 504: (1, 0.2, 0.8)}
    ops.uniaxialMaterial("Elastic", 91, 10.0)
    for t, (x, y, z) in sv.items():
        ops.node(t, x, y, z)
        ops.node(t + 100, x, y, z)
        ops.fix(t + 100, 1, 1, 1)
        ops.element("zeroLength", t, t + 100, t, "-mat", 91, 91, 91, "-dir", 1, 2, 3)
    ops.LadrunoTie("-mortar", "-dual", "-slaveFacets", 4, 1, 501, 502, 503, 504,
                   "-masterFacets", npsM, len(master_facets), *sum(master_facets, []),
                   "-outward", 1.0, 0.0, 0.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in master_nodes.items():
        ops.sp(t, 1, field(y, z))
    _static_solve("Lagrange")
    return sv


_TRI6_NODES = {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1), 4: (1, 0, 1),
               5: (1, 0.5, 0), 6: (1, 1, 0.5), 7: (1, 0.5, 0.5),
               8: (1, 0.5, 1), 9: (1, 0, 0.5)}
_TRI6_A = [1, 2, 3, 5, 6, 7]
_TRI6_B = [1, 3, 4, 7, 8, 9]


def test_tri6_master_linear_field_exact():
    field = lambda y, z: 0.3 + 1.7 * y - 0.9 * z
    sv = _spring_rig_solve(_TRI6_NODES, [_TRI6_A, _TRI6_B], 6, field)
    for t, (x, y, z) in sv.items():
        assert abs(ops.nodeDisp(t, 1) - field(y, z)) <= 1e-9, f"tri6-master transfer wrong at {t}"


def test_quad8_master_linear_field_exact():
    mn = {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1), 4: (1, 0, 1),
          5: (1, 0.5, 0), 6: (1, 1, 0.5), 7: (1, 0.5, 1), 8: (1, 0, 0.5)}
    field = lambda y, z: -0.2 + 0.8 * y + 1.3 * z
    sv = _spring_rig_solve(mn, [[1, 2, 3, 4, 5, 6, 7, 8]], 8, field)
    for t, (x, y, z) in sv.items():
        assert abs(ops.nodeDisp(t, 1) - field(y, z)) <= 1e-9, f"quad8-master transfer wrong at {t}"


def _bare_master():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _MASTER.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.0, 0.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)


def _nodes(d):
    for t, (x, y, z) in d.items():
        ops.node(t, float(x), float(y), float(z))


def test_refuse_tri6_slave():
    _bare_master()
    tri6 = {201: (1, 0, 0), 202: (1, 1, 0), 203: (1, 1, 1),
            204: (1, 0.5, 0), 205: (1, 1, 0.5), 206: (1, 0.5, 0.5)}
    _nodes(tri6)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-mortar", "-slaveFacets", 6, 1, *sorted(tri6),
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_curved_quad8():
    _bare_master()
    q8 = {201: (1, 0, 0), 202: (1, 1, 0), 203: (1, 1, 1), 204: (1, 0, 1),
          205: (1, 0.5, 0.1), 206: (1, 1, 0.5), 207: (1, 0.5, 1), 208: (1, 0, 0.5)}
    _nodes(q8)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-mortar", "-slaveFacets", 8, 1, *sorted(q8),
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_quad8_off_surface():
    _bare_master()
    q8 = {201: (1.4, 0, 0), 202: (1.4, 1, 0), 203: (1.4, 1, 1), 204: (1.4, 0, 1),
          205: (1.4, 0.5, 0), 206: (1.4, 1, 0.5), 207: (1.4, 0.5, 1), 208: (1.4, 0, 0.5)}
    _nodes(q8)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-mortar", "-slaveFacets", 8, 1, *sorted(q8),
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_collocation_quadratic_master():
    _bare_master()
    q8 = {201: (1, 0, 0), 202: (1, 1, 0), 203: (1, 1, 1), 204: (1, 0, 1),
          205: (1, 0.5, 0), 206: (1, 1, 0.5), 207: (1, 0.5, 1), 208: (1, 0, 0.5)}
    _nodes(q8)
    ops.node(300, 1.0, 0.4, 0.4)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-slaveNodes", 1, 300, "-masterFacets", 8, 1, *sorted(q8))


def test_refuse_tilted_ridge_slave():
    """Two slave facets tilted in opposite senses with their shared ridge on the
    master plane: the signed gap cancels at the ridge nodes, the L1 gap does
    not, and the tie is refused."""
    _bare_master()
    d = 0.2
    acc = {201: (1 - d, 0, 0), 202: (1, 0.5, 0), 203: (1, 0.5, 1), 204: (1 - d, 0, 1),
           205: (1 + d, 1, 0), 206: (1 + d, 1, 1)}
    _nodes(acc)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-mortar",
                       "-slaveFacets", 4, 2, 201, 202, 203, 204, 202, 205, 206, 203,
                       "-masterFacets", 4, 1, *MASTER_FACET, "-outward", 1.0, 0.0, 0.0)


def test_refuse_duplicate_master_facets():
    _bare_master()
    sn = {201: (1, 0.2, 0.2), 202: (1, 0.8, 0.2), 203: (1, 0.8, 0.8), 204: (1, 0.2, 0.8)}
    _nodes(sn)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 201, 202, 203, 204,
                       "-masterFacets", 4, 2, *MASTER_FACET, *MASTER_FACET)


def test_accept_tri6_master():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _nodes(_TRI6_NODES)
    sv = {501: (1, 0.2, 0.2), 502: (1, 0.8, 0.2), 503: (1, 0.8, 0.8), 504: (1, 0.2, 0.8)}
    _nodes(sv)
    ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 501, 502, 503, 504,
                   "-masterFacets", 6, 2, *_TRI6_A, *_TRI6_B, "-outward", 1.0, 0.0, 0.0)
