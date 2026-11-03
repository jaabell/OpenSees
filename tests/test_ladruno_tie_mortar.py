"""Tests for the LadrunoTie integral mortar tie (-mortar, -dual).

The mortar tie assembles D_IJ = int N_I N_J and M_IK = int N_I phi_K over the
clipped slave/master overlap and writes u_s = P u_m with P = D^-1 M (dense),
or with the dual basis (-dual) that makes D diagonal (P local).

Mesh: a non-matching solid interface at x = 1. The master block is one
stdBrick hexahedron, x in [0, 1], whose x = 1 face is one quad facet. The
slave block, x in [1, 2], is meshed 2 x 2 in (y, z); its x = 1 face is a 3 x 3
node grid (4 quad facets). The slave nodes at y or z = 0.5 have no master
counterpart.

- Patch test: u_x = gamma z imposed on every node except the tied ones gives
  the constant shear stress E gamma / 2 at every Gauss point (standard and
  dual basis; Lagrange and Penalty).
- Split column: the tied column reproduces the conforming column.
- Tie residual: the slave facets refine the master facet, so the mortar
  projection of the bilinear master trace is that trace itself, and under a
  general load every tied node satisfies u_s = sum N_k(y_s, z_s) u_k exactly.
- Closed surface: a tie around a tube (16 master and 24 slave facets) under a
  rigid rotation about the axis. The far side of the tube projects onto the
  auxiliary plane of every near facet; those pairs are excluded by their
  separation, so the slaves follow the rotation instead of being averaged
  with the opposite side of the tube.
- Refusals: uncovered slave, slave surface off the master plane, a node on
  both surfaces, -dual without -mortar. A non-affine slave facet spanning
  several master facets is accepted.
"""
import math

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
MASTER_FACET = [2, 3, 7, 6]       # (y, z) = (0,0) (1,0) (1,1) (0,1) on x = 1

_YS = [0.0, 0.5, 1.0]
_ZS = [0.0, 0.5, 1.0]
_XS = [1.0, 2.0]


def _stag(ix, iy, iz):
    return 100 + ix * 9 + iy * 3 + iz


def _slave_nodes():
    nodes = {}
    for ix, x in enumerate(_XS):
        for iy, y in enumerate(_YS):
            for iz, z in enumerate(_ZS):
                nodes[_stag(ix, iy, iz)] = (x, y, z)
    return nodes


def _slave_hexes():
    hexes = []
    for jy in range(2):
        for jz in range(2):
            n = lambda ix, dy, dz: _stag(ix, jy + dy, jz + dz)
            hexes.append([n(0, 0, 0), n(1, 0, 0), n(1, 1, 0), n(0, 1, 0),
                          n(0, 0, 1), n(1, 0, 1), n(1, 1, 1), n(0, 1, 1)])
    return hexes


def _slave_facets():
    fac = []
    for jy in range(2):
        for jz in range(2):
            fac += [_stag(0, jy, jz), _stag(0, jy + 1, jz),
                    _stag(0, jy + 1, jz + 1), _stag(0, jy, jz + 1)]
    return fac


SLAVE_NODES = _slave_nodes()
SLAVE_IFACE = [_stag(0, iy, iz) for iy in range(3) for iz in range(3)]
_ALL = {**_MASTER, **SLAVE_NODES}


def _mesh(E=1000.0, nu=0.0, dual=False):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, 0.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for e, conn in enumerate(_slave_hexes(), start=2):
        ops.element("stdBrick", e, *conn, 1)
    sf = _slave_facets()
    args = ["-mortar", "-slaveFacets", 4, len(sf) // 4, *sf,
            "-masterFacets", 4, 1, *MASTER_FACET]
    if dual:
        args.append("-dual")
    ops.LadrunoTie(*args)


def _gp_stresses(ele):
    ops.eleResponse(ele, "forces")
    s = list(ops.eleResponse(ele, "stresses"))
    return [s[i * 6:(i + 1) * 6] for i in range(len(s) // 6)]


def _static_solve(handler="Lagrange", system="FullGeneral"):
    if handler == "Penalty":
        ops.constraints("Penalty", 1.0e12, 1.0e12)
    else:
        ops.constraints(handler)
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", 1e-12, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


@pytest.mark.parametrize("dual", [False, True])
@pytest.mark.parametrize("handler,rtol", [("Lagrange", 1e-6), ("Penalty", 1e-5)])
def test_mortar_constant_stress_patch(dual, handler, rtol):
    E, gamma = 1000.0, 1e-3
    _mesh(E=E, nu=0.0, dual=dual)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_IFACE:
            continue
        ops.fix(t, 0, 1, 1)
        ops.sp(t, 1, gamma * float(z))
    _static_solve(handler)

    allgp = []
    for e in range(1, 6):
        allgp += _gp_stresses(e)
    ref = allgp[0]
    target = E * gamma / 2.0
    atol = rtol * target
    for v in allgp:
        for si, ri in zip(v, ref):
            assert abs(si - ri) <= atol, f"non-constant stress across the tie: {v} vs {ref}"
    comps = [abs(c) for c in ref]
    assert abs(max(comps) - target) <= atol, f"stress {ref} != {target}"
    assert sum(1 for c in comps if c > 1e-3 * target) == 1, f"not pure shear: {ref}"


def _monolithic_axial(E, U):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nodes, tag = {}, 1
    for xi in (0.0, 1.0, 2.0):
        for (y, z) in [(0, 0), (1, 0), (1, 1), (0, 1)]:
            nodes[tag] = (xi, y, z)
            ops.node(tag, float(xi), float(y), float(z))
            tag += 1
    ops.nDMaterial("ElasticIsotropic", 1, E, 0.0, 0.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 2, 5, 6, 7, 8, 9, 10, 11, 12, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in nodes.items():
        ops.fix(t, 0, 1, 1)
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, 1, U)
    _static_solve("Transformation")
    return _gp_stresses(1)[0]


@pytest.mark.parametrize("dual", [False, True])
def test_mortar_split_column_matches_conforming(dual):
    E, U = 1000.0, 0.02
    ref_stress = _monolithic_axial(E, U)
    _mesh(E=E, nu=0.0, dual=dual)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_IFACE:
            continue
        ops.fix(t, 0, 1, 1)
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, 1, U)
    _static_solve("Lagrange")
    mid = ops.nodeDisp(2, 1)
    assert abs(mid - 0.5 * U) <= 1e-6 * U, f"mid-plane {mid} != {0.5 * U}"
    got = _gp_stresses(1)[0]
    for gi, ri in zip(got, ref_stress):
        assert abs(gi - ri) <= 1e-5 * E * U, f"stress {got} != {ref_stress}"


@pytest.mark.parametrize("dual", [False, True])
@pytest.mark.parametrize("handler,rtol", [("Lagrange", 1e-9), ("Penalty", 1e-6)])
def test_mortar_tie_residual_general_load(dual, handler, rtol):
    """Clamped master base, an unsymmetric load on the slave end (nu = 0.3):
    u_s = sum N_k(y_s, z_s) u_k on every tied node and DOF."""
    _mesh(E=1000.0, nu=0.3, dual=dual)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _MASTER.items():
        if abs(x) < 1e-12:
            ops.fix(t, 1, 1, 1)
    ops.load(_stag(1, 0, 2), 4.0, 1.0, -2.0)
    ops.load(_stag(1, 2, 0), -1.0, 3.0, 2.0)
    ops.load(_stag(1, 1, 1), 2.0, -2.0, 1.0)
    _static_solve(handler)

    corner = {(0, 0): 2, (1, 0): 3, (1, 1): 7, (0, 1): 6}
    umax = max(abs(ops.nodeDisp(t, d)) for t in SLAVE_NODES for d in (1, 2, 3))
    assert umax > 1e-4, "the load did not deform the model"
    worst = 0.0
    for t in SLAVE_IFACE:
        _, y, z = SLAVE_NODES[t]
        w = {corner[(0, 0)]: (1 - y) * (1 - z), corner[(1, 0)]: y * (1 - z),
             corner[(1, 1)]: y * z, corner[(0, 1)]: (1 - y) * z}
        for d in (1, 2, 3):
            um = sum(wk * ops.nodeDisp(m, d) for m, wk in w.items())
            worst = max(worst, abs(ops.nodeDisp(t, d) - um))
    assert worst <= rtol * umax, f"mortar tie residual {worst} (umax {umax})"


# ---------------------------------------------------------------- closed surface
R_TUBE, L_TUBE, W_ROT = 0.5, 1.0, 0.01


def _ring(n, tag0):
    nodes, fac = {}, []
    for k in range(n):
        a = 2.0 * math.pi * k / n
        for j, z in enumerate((0.0, L_TUBE)):
            nodes[tag0 + 2 * k + j] = (R_TUBE * math.cos(a), R_TUBE * math.sin(a), z)
    for k in range(n):
        k1 = (k + 1) % n
        fac += [tag0 + 2 * k, tag0 + 2 * k1, tag0 + 2 * k1 + 1, tag0 + 2 * k + 1]
    return nodes, fac


@pytest.mark.parametrize("tol", [0.1, 3.5])
def test_mortar_closed_tube_follows_rigid_rotation(tol):
    """One tie around a closed tube, masters prescribed to a rigid rotation
    about the tube axis, each slave held by a soft spring to the ground. The
    slaves must follow the rotation (u_t = W R); binding the antipodal facets
    as well gives an error of 97% of W R (the opposite side moves the other
    way), and the tie is refused at tol <= 1 by the per-facet gap check.

    tol = 3.5 is a value at which, without the separation test, the far-side
    pairs pass every per-facet check."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mn, mf = _ring(16, 1)
    sn, sf = _ring(24, 1001)
    ops.uniaxialMaterial("Elastic", 9, 1.0e-3)
    for t, (x, y, z) in mn.items():
        ops.node(t, x, y, z)
    for t, (x, y, z) in sn.items():
        ops.node(t, x, y, z)
        ops.node(t + 5000, x, y, z)
        ops.fix(t + 5000, 1, 1, 1)
        ops.element("zeroLength", t, t + 5000, t, "-mat", 9, 9, 9, "-dir", 1, 2, 3)
    ops.LadrunoTie("-mortar", "-slaveFacets", 4, len(sf) // 4, *sf,
                   "-masterFacets", 4, len(mf) // 4, *mf,
                   "-tol", tol, "-outward", 1.0, 0.3, 0.2)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in mn.items():
        ops.sp(t, 1, -W_ROT * y)
        ops.sp(t, 2, W_ROT * x)
        ops.sp(t, 3, 0.0)
    _static_solve("Lagrange", system="UmfPack")
    err = 0.0
    for t, (x, y, z) in sn.items():
        err = max(err, abs(ops.nodeDisp(t, 1) + W_ROT * y),
                  abs(ops.nodeDisp(t, 2) - W_ROT * x), abs(ops.nodeDisp(t, 3)))
    # the faceting offset between the 16- and 24-gon (<= 0.0096 radially)
    # leaves an error of 0.9% of W R; bound at 2%
    assert err <= 0.02 * W_ROT * R_TUBE, f"slaves do not follow the rotation: err/(W R) = {err / (W_ROT * R_TUBE)}"


# ---------------------------------------------------------------- refusals
def _bare_master():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _MASTER.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.0, 0.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)


def test_mortar_refuse_uncovered_slave():
    _bare_master()
    sn = {201: (1, 0, 0), 202: (1, 2, 0), 203: (1, 2, 1), 204: (1, 0, 1)}
    for t, (x, y, z) in sn.items():
        ops.node(t, float(x), float(y), float(z))
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 201, 202, 203, 204,
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_mortar_refuse_off_surface():
    _bare_master()
    sn = {201: (1.4, 0, 0), 202: (1.4, 1, 0), 203: (1.4, 1, 1), 204: (1.4, 0, 1)}
    for t, (x, y, z) in sn.items():
        ops.node(t, float(x), float(y), float(z))
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 201, 202, 203, 204,
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_mortar_refuse_non_disjoint():
    _bare_master()
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, 2, 3, 7, 6,
                       "-masterFacets", 4, 1, *MASTER_FACET)


def test_dual_requires_mortar():
    _bare_master()
    ops.node(99, 1.0, 0.5, 0.5)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-dual", "-slaveNodes", 1, 99, "-masterFacets", 4, 1, *MASTER_FACET)


def test_mortar_nonaffine_slave_accepted():
    """A trapezoidal slave facet spanning four master facets: the multi-piece
    overlap and the self-clip area differ in quadrature only, inside the
    coverage slack, so the tie is accepted. Slave nodes carry no mass."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mtag, t = {}, 1
    for iy, y in enumerate([0.0, 0.5, 1.0]):
        for iz, z in enumerate([0.0, 0.5, 1.0]):
            mtag[(iy, iz)] = t
            ops.node(t, 1.0, y, z)
            t += 1
    mfac = []
    for iy in range(2):
        for iz in range(2):
            mfac += [mtag[(iy, iz)], mtag[(iy + 1, iz)], mtag[(iy + 1, iz + 1)], mtag[(iy, iz + 1)]]
    sverts = [(1.0, 0.10, 0.15), (1.0, 0.90, 0.20), (1.0, 0.78, 0.85), (1.0, 0.18, 0.92)]
    stags = []
    for i, (x, y, z) in enumerate(sverts):
        ops.node(500 + i, x, y, z)
        stags.append(500 + i)
    ops.LadrunoTie("-mortar", "-slaveFacets", 4, 1, *stags,
                   "-masterFacets", 4, len(mfac) // 4, *mfac)
