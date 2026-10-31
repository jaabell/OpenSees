"""Tests for the LadrunoTie shell-to-solid plane-section tie (-shellSolid),
static Lagrange analyses with SSPbrick and ShellMITC4.

An ndf-6 shell EDGE (the master polyline) is tied to an ndf-3 solid FACE (the
slaves) by the linearized rigid plane-section map

    u_s = sum_{j=A,B} N_j(xi) [u_j + theta_j x d],   d = x_s - xbar(xi) (frozen)

three rows per solid face node; the drilling coefficient is ~zero and is
dropped, so the shell drilling DOF stays free.

Mesh: a solid block x in [0, 1] (4 x 4 x 2 SSPbrick, thickness t = 0.1 in z)
butting a shell strip x in [1, 2] (two ShellMITC4, split at y = 0.5, mid-surface
z = 0). The 15 solid face nodes at x = 1 (y in {0, .25, .5, .75, 1}, z in
{-0.05, 0, 0.05}) are the slaves; the solid y grid nests 2:1 in the shell one.

- Membrane: a uniform pull (nu = 0, equal EA) gives the interface at Delta/2
  and a constant sigma_xx in the solid.
- Rigid motion: all shell DOFs prescribed to a rigid rotation and
  translation; every solid node follows the rigid field.
- Section state: a prescribed constant-moment section state (w0, theta_y0)
  puts the face nodes on u = N (u_j + theta_j x d), linear in z.
- Bending: an imposed far-edge rotation on members of equal EI crosses the tie
  (interface rotation ~ theta/2) and the face satisfies the tie rows.
- Hinge: the same load through translation-only rows (equationConstraint) is
  a hinge: the shell rotates freely and the solid face stays unbent.
- Refusals: -mortar, -hermite, -dof, -outward, an ndf-3 master edge, a slave
  outside the polyline, the -thickness arm check, a node on both sides,
  -masterEdge without -shellSolid. Without -thickness the arm check is skipped.
- A slave without mass is accepted.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# a refused LadrunoTie command raises the interpreter error (an unknown
# command would raise AttributeError instead and fail these tests)
OPS_ERROR = getattr(ops, "OpenSeesError", RuntimeError)


E, NU, T = 1000.0, 0.0, 0.1
C = T / 2.0

XS = [0.0, 0.25, 0.5, 0.75, 1.0]
YS = [0.0, 0.25, 0.5, 0.75, 1.0]
ZS = [-C, 0.0, C]

_SHELL = {201: (1.0, 0.0), 202: (1.0, 0.5), 203: (1.0, 1.0),
          211: (2.0, 0.0), 212: (2.0, 0.5), 213: (2.0, 1.0)}
EDGE_MASTERS = [201, 202, 203]
FAR_EDGE = [211, 212, 213]
_SHELL_QUADS = [[201, 211, 212, 202], [202, 212, 213, 203]]
MASTER_EDGE = [2, 201, 202, 202, 203]      # -masterEdge <nseg = 2> a1 b1 a2 b2


def _solid_tag(i, j, k):
    return 1000 + i * 100 + j * 10 + k


def _face_nodes():
    return [_solid_tag(len(XS) - 1, j, k) for j in range(len(YS)) for k in range(3)]


def _face_coords():
    return {_solid_tag(len(XS) - 1, j, k): (YS[j], ZS[k])
            for j in range(len(YS)) for k in range(3)}


def _mesh(tie=True):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, 0.0)
    for i, x in enumerate(XS):
        for j, y in enumerate(YS):
            for k, z in enumerate(ZS):
                ops.node(_solid_tag(i, j, k), x, y, z)
    eid = 1
    for i in range(len(XS) - 1):
        for j in range(len(YS) - 1):
            for k in range(len(ZS) - 1):
                n = [_solid_tag(i, j, k), _solid_tag(i + 1, j, k),
                     _solid_tag(i + 1, j + 1, k), _solid_tag(i, j + 1, k),
                     _solid_tag(i, j, k + 1), _solid_tag(i + 1, j, k + 1),
                     _solid_tag(i + 1, j + 1, k + 1), _solid_tag(i, j + 1, k + 1)]
                ops.element("SSPbrick", eid, *n, 1)
                eid += 1
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y) in _SHELL.items():
        ops.node(t, x, y, 0.0)
    ops.section("ElasticMembranePlateSection", 1, E, NU, T, 0.0)
    for e, conn in enumerate(_SHELL_QUADS, start=100):
        ops.element("ShellMITC4", e, *conn, 1)
    if tie:
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 15, *_face_nodes(),
                       "-masterEdge", *MASTER_EDGE, "-thickness", T)


def _plane_section_pred(y, z):
    """Right-hand side of the tie rows for a face slave at (1, y, z):
    u = sum N_j [u_j + theta_j x d], d = (0, 0, z), theta x d = (theta_y z, -theta_x z, 0)."""
    if y <= 0.5:
        (mA, mB), xi = (201, 202), y / 0.5
    else:
        (mA, mB), xi = (202, 203), (y - 0.5) / 0.5
    pred = [0.0, 0.0, 0.0]
    for m, N in ((mA, 1.0 - xi), (mB, xi)):
        u = [ops.nodeDisp(m, d) for d in (1, 2, 3)]
        th = [ops.nodeDisp(m, d) for d in (4, 5, 6)]
        pred[0] += N * (u[0] + th[1] * z)
        pred[1] += N * (u[1] - th[0] * z)
        pred[2] += N * u[2]
    return pred


def _assert_face_satisfies_rows(tol):
    for s, (y, z) in _face_coords().items():
        pred = _plane_section_pred(y, z)
        for d in (1, 2, 3):
            err = abs(ops.nodeDisp(s, d) - pred[d - 1])
            assert err <= tol, f"face node {s} dof {d} violates the tie by {err}"


def _static_solve():
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 30)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def test_membrane_patch_crosses_tie():
    DELTA = 1e-3
    _mesh()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(len(YS)):
        for k in range(3):
            ops.fix(_solid_tag(0, j, k), 1, 1, 1)
    for t in FAR_EDGE:
        ops.sp(t, 1, DELTA)
    _static_solve()
    for s in _face_nodes():
        assert abs(ops.nodeDisp(s, 1) - 0.5 * DELTA) <= 1e-8 * DELTA, f"face node {s} ux != Delta/2"
    for m in EDGE_MASTERS:
        assert abs(ops.nodeDisp(m, 1) - 0.5 * DELTA) <= 1e-8 * DELTA
    sig = E * 0.5 * DELTA
    for e in range(1, 33):
        ops.eleResponse(e, "forces")
        s = list(ops.eleResponse(e, "stresses"))
        assert abs(s[0] - sig) <= 1e-6 * sig, f"element {e} sigma_xx {s[0]} != {sig}"
    _assert_face_satisfies_rows(1e-12 * DELTA + 1e-15)


def test_rigid_motion_carries_the_solid():
    U0 = (0.003, -0.002, 0.001)
    OM = (0.010, -0.020, 0.015)
    P0 = (0.3, -0.2, 0.05)

    def cross(a, b):
        return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

    _mesh()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for m, (x, y) in _SHELL.items():
        u = cross(OM, (x - P0[0], y - P0[1], 0.0 - P0[2]))
        vals = [U0[0]+u[0], U0[1]+u[1], U0[2]+u[2], OM[0], OM[1], OM[2]]
        for d in range(6):
            ops.sp(m, d + 1, vals[d])
    _static_solve()
    err = 0.0
    for i, x in enumerate(XS):
        for j, y in enumerate(YS):
            for k, z in enumerate(ZS):
                u = cross(OM, (x - P0[0], y - P0[1], z - P0[2]))
                for d in range(3):
                    err = max(err, abs(ops.nodeDisp(_solid_tag(i, j, k), d + 1) - (U0[d] + u[d])))
    assert err < 1e-10, f"solid did not follow the rigid field (err = {err})"


def test_plane_section_state_on_face():
    W0, TH0 = 2.0e-3, -4.0e-3
    _mesh()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for m in _SHELL:
        vals = [0.0, 0.0, W0, 0.0, TH0, 0.0]
        for d in range(6):
            ops.sp(m, d + 1, vals[d])
    _static_solve()
    for s, (y, z) in _face_coords().items():
        exp = (TH0 * z, 0.0, W0)
        for d in (1, 2, 3):
            err = abs(ops.nodeDisp(s, d) - exp[d - 1])
            assert err <= 1e-10, f"face node {s} dof {d} off the plane section by {err}"


THETA = 0.01


def _bend_bcs_and_solve():
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(len(YS)):
        for k in range(3):
            ops.fix(_solid_tag(0, j, k), 1, 1, 1)
    for t in FAR_EDGE:
        ops.sp(t, 5, THETA)
    _static_solve()


def test_bending_crosses_tie():
    """Equal-EI members under an imposed end rotation with a free tip: a
    constant moment, theta(x) ~ THETA x / 2, so the shell edge rotates
    ~THETA/2 and the face carries a through-thickness couple."""
    _mesh()
    _bend_bcs_and_solve()
    th_edge = sum(ops.nodeDisp(m, 5) for m in EDGE_MASTERS) / 3.0
    assert 0.3 * THETA < th_edge < 0.7 * THETA, f"interface rotation {th_edge} not ~THETA/2"
    _assert_face_satisfies_rows(1e-10)
    g = {}
    for s, (y, z) in _face_coords().items():
        g.setdefault(y, {})[z] = ops.nodeDisp(s, 1)
    for y, col in g.items():
        grad = (col[C] - col[-C]) / (2.0 * C)
        assert abs(grad) > 0.2 * THETA, f"no through-thickness couple at y = {y}"


def test_translations_only_rows_are_a_hinge():
    _mesh(tie=False)
    for s, (y, z) in _face_coords().items():
        if y <= 0.5:
            (mA, mB), xi = (201, 202), y / 0.5
        else:
            (mA, mB), xi = (202, 203), (y - 0.5) / 0.5
        for d in (1, 2, 3):
            args = [s, d, 1.0]
            for m, N in ((mA, 1.0 - xi), (mB, xi)):
                if abs(N) > 1e-12:
                    args += [m, d, -N]
            ops.equationConstraint(*args)
    _bend_bcs_and_solve()
    th_edge = sum(ops.nodeDisp(m, 5) for m in EDGE_MASTERS) / 3.0
    assert abs(th_edge - THETA) < 0.05 * THETA, f"hinge should rotate freely (theta_edge = {th_edge})"
    fc = _face_coords()
    for y in YS:
        top = ops.nodeDisp([s for s, (yy, zz) in fc.items() if yy == y and zz == C][0], 1)
        bot = ops.nodeDisp([s for s, (yy, zz) in fc.items() if yy == y and zz == -C][0], 1)
        assert abs((top - bot) / (2.0 * C)) < 0.05 * THETA, f"hinge carried a couple at y = {y}"


@pytest.mark.parametrize("extra", [("-mortar",), ("-hermite",), ("-dof", 3, 1, 2, 3),
                                   ("-outward", 0.0, 0.0, 1.0)])
def test_refuse_incompatible_options(extra):
    _mesh(tie=False)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 15, *_face_nodes(),
                       "-masterEdge", *MASTER_EDGE, *extra)


def test_refuse_low_ndf_master_edge():
    _mesh(tie=False)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 15, *_face_nodes(),
                       "-masterEdge", 1, _solid_tag(0, 0, 0), _solid_tag(0, 2, 0))


def test_refuse_outside_polyline():
    _mesh(tie=False)
    ops.node(999, 1.0, 1.5, 0.0, "-ndf", 3)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 999, "-masterEdge", *MASTER_EDGE)


def test_arm_check_with_and_without_thickness():
    """A slave 2t off the mid-surface is refused with -thickness and accepted
    without it (the check is then skipped)."""
    _mesh(tie=False)
    ops.node(999, 1.0, 0.5, 2.0 * T, "-ndf", 3)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 999,
                       "-masterEdge", *MASTER_EDGE, "-thickness", T)
    ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 999, "-masterEdge", *MASTER_EDGE)


def test_refuse_non_disjoint():
    _mesh(tie=False)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 201,
                       "-masterEdge", *MASTER_EDGE, "-thickness", T)


def test_refuse_master_edge_without_shell_solid():
    _mesh(tie=False)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-slaveNodes", 15, *_face_nodes(), "-masterEdge", *MASTER_EDGE)


def test_massless_slave_accepted():
    _mesh(tie=False)
    ops.node(999, 1.0, 0.5, 0.02, "-ndf", 3)       # on the face, no element, no mass
    ops.LadrunoTie("-shellSolid", "-slaveNodes", 1, 999,
                   "-masterEdge", *MASTER_EDGE, "-thickness", T)
