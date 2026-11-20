"""Tests for LadrunoTie on ndf-6 shell nodes (rotations, -hermite), static
Lagrange analyses with ShellMITC4.

For a 3D ndf-6 node the tied DOFs default to 1..6, so the rotations are tied
with the same transfer as the translations, theta_s = sum P_sk theta_k.

- Bending: a constant moment crosses a non-conforming butt-joint tie (master
  shell x in [0, 1], slave shells x in [1, 2] split at y = 0.5): the tied strip
  reproduces the conforming strip (interface rotation and transverse
  displacement, tip displacement). With -hermite the result is the same.
- Hinge: tying translations only (-dof 3 1 2 3) leaves a hinge; the moment
  does not cross.
- Mortar: a rigid rotation about y (uz = -phi x, theta_y = phi) crosses a
  non-matching mortar shell tie (2 x 2 master, 3 x 3 slave) on all 6 DOFs,
  standard and dual basis.
- Hermite: bending that varies along the interface (w = kappa y^2 / 2,
  theta_x = kappa y) crosses exactly with -hermite at a mid-edge slave node;
  the linear transfer gives the linear interpolant kappa/4 instead of the
  exact kappa/8 there.
- Refusals: -hermite with -mortar, -hermite with a partial -dof, -hermite with
  a slave inside a master facet.
- ShellMITC4 lumps no rotational mass; the 6-DOF tie is accepted without
  nodal rotary mass.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# a refused LadrunoTie command raises the interpreter error (an unknown
# command would raise AttributeError instead and fail these tests)
OPS_ERROR = getattr(ops, "OpenSeesError", RuntimeError)


E, NU, H = 1000.0, 0.0, 0.1
THETA = 0.01           # imposed tip rotation (rad)
UZ, THY = 3, 5         # 1-based DOF indices


def _shell_section():
    ops.section("ElasticMembranePlateSection", 1, E, NU, H, 0.0)


_MASTER = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0)}
MASTER_FACET = [1, 2, 3, 4]
_SLAVE = {11: (1, 0, 0), 12: (2, 0, 0), 13: (2, 0.5, 0), 14: (1, 0.5, 0),
          15: (2, 1, 0), 16: (1, 1, 0)}
_SLAVE_QUADS = [[11, 12, 13, 14], [14, 13, 15, 16]]
SLAVE_IFACE = [11, 14, 16]                     # x = 1 slave nodes (14 non-conforming)
_ALL = {**_MASTER, **_SLAVE}


def _mesh_butt():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    _shell_section()
    ops.element("ShellMITC4", 1, *MASTER_FACET, 1)
    for e, conn in enumerate(_SLAVE_QUADS, start=2):
        ops.element("ShellMITC4", e, *conn, 1)


def _static_solve(handler="Lagrange"):
    ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 30)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _apply_bending_bcs(nodes):
    """Clamp the x = 0 root, impose theta_y = THETA at the x = 2 tip."""
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in nodes.items():
        if abs(x) < 1e-12:
            ops.fix(t, 1, 1, 1, 1, 1, 1)
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, THY, THETA)


def _monolithic_bending():
    """Conforming strip x in [0, 2]. Returns (theta_y, uz) at x = 1 and uz at the tip."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nd = {1: (0, 0, 0), 2: (0, 1, 0), 3: (1, 0, 0), 4: (1, 1, 0), 5: (2, 0, 0), 6: (2, 1, 0)}
    for t, (x, y, z) in nd.items():
        ops.node(t, float(x), float(y), float(z))
    _shell_section()
    ops.element("ShellMITC4", 1, 1, 3, 4, 2, 1)
    ops.element("ShellMITC4", 2, 3, 5, 6, 4, 1)
    _apply_bending_bcs(nd)
    _static_solve("Transformation")
    return ops.nodeDisp(3, THY), ops.nodeDisp(3, UZ), ops.nodeDisp(5, UZ)


@pytest.mark.parametrize("extra", [(), ("-hermite",)])
def test_collocation_bending_crosses_tie(extra):
    th_ref, uz_ref, tip_ref = _monolithic_bending()
    assert abs(th_ref) > 1e-6, "reference did not bend"

    _mesh_butt()
    ops.LadrunoTie("-slaveNodes", 3, *SLAVE_IFACE, "-masterFacets", 4, 1, *MASTER_FACET, *extra)
    _apply_bending_bcs(_ALL)
    _static_solve("Lagrange")

    tol = 1e-6 * max(abs(th_ref), abs(uz_ref), abs(tip_ref))
    assert abs(ops.nodeDisp(2, THY) - th_ref) <= tol, "interface rotation did not cross"
    assert abs(ops.nodeDisp(2, UZ) - uz_ref) <= tol, "interface transverse displacement did not cross"
    assert abs(ops.nodeDisp(14, THY) - th_ref) <= tol, "non-conforming node rotation wrong"
    assert abs(ops.nodeDisp(12, UZ) - tip_ref) <= tol, "slave tip displacement wrong"


def test_collocation_translations_only_is_a_hinge():
    th_ref, _, _ = _monolithic_bending()
    _mesh_butt()
    ops.LadrunoTie("-slaveNodes", 3, *SLAVE_IFACE, "-masterFacets", 4, 1, *MASTER_FACET,
                   "-dof", 3, 1, 2, 3)
    _apply_bending_bcs(_ALL)
    _static_solve("Lagrange")
    assert abs(ops.nodeDisp(2, THY)) < 1e-3 * abs(th_ref), \
        "moment crossed a translations-only tie"


def _shell_grid(n, tag0):
    idx = lambda i, j: tag0 + j * (n + 1) + i
    nodes = {idx(i, j): (i / n, j / n, 0.0) for j in range(n + 1) for i in range(n + 1)}
    quads = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)]
             for j in range(n) for i in range(n)]
    return nodes, quads


def _overlap_shells():
    m_nodes, m_quads = _shell_grid(2, 1)               # nodes 1..9, centre 5
    s_nodes, s_quads = _shell_grid(3, 101)             # nodes 101..116 (non-matching)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in {**m_nodes, **s_nodes}.items():
        ops.node(t, x, y, z)
    _shell_section()
    eid = 1
    for conn in m_quads + s_quads:
        ops.element("ShellMITC4", eid, *conn, 1)
        eid += 1
    sf = [n for q in s_quads for n in q]
    mf = [n for q in m_quads for n in q]
    return m_nodes, s_nodes, sf, mf, len(s_quads), len(m_quads)


@pytest.mark.parametrize("extra", [(), ("-dual",)])
def test_mortar_rigid_rotation_crosses_tie(extra):
    phi = 0.02
    m_nodes, s_nodes, sf, mf, nsq, nmq = _overlap_shells()
    ops.LadrunoTie("-mortar", "-slaveFacets", 4, nsq, *sf, "-masterFacets", 4, nmq, *mf,
                   "-outward", 0.0, 0.0, 1.0, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in m_nodes.items():
        if t == 5:
            continue                                   # free interior master node
        ops.fix(t, 1, 1, 0, 1, 0, 1)
        ops.sp(t, UZ, -phi * x)
        ops.sp(t, THY, phi)
    _static_solve("Lagrange")
    tol = 1e-6
    for t, (x, y, z) in s_nodes.items():
        assert abs(ops.nodeDisp(t, UZ) - (-phi * x)) <= tol, f"slave {t} uz wrong"
        assert abs(ops.nodeDisp(t, THY) - phi) <= tol, f"slave {t} theta_y wrong"
        for d in (1, 2, 4, 6):
            assert abs(ops.nodeDisp(t, d)) <= tol, f"slave {t} spurious DOF {d}"


KAPPA = 0.04           # curvature along the interface (y): w = kappa y^2 / 2, theta_x = kappa y


def _mesh_hermite(*extra):
    """Master quad [0, 1]^2 and three bare slave nodes on its x = 1 edge (y = 0,
    0.5, 1). The master carries the exact bending field; the tie carries it to
    the slaves (all slave DOFs are tied)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in _MASTER.items():
        ops.node(t, float(x), float(y), float(z))
    slaves = {21: (1.0, 0.0), 22: (1.0, 0.5), 23: (1.0, 1.0)}
    for t, (x, y) in slaves.items():
        ops.node(t, x, y, 0.0)
    _shell_section()
    ops.element("ShellMITC4", 1, *MASTER_FACET, 1)
    ops.LadrunoTie("-slaveNodes", 3, 21, 22, 23, "-masterFacets", 4, 1, *MASTER_FACET, *extra)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _MASTER.items():
        if abs(y) < 1e-12:
            ops.fix(t, 1, 1, 1, 1, 1, 1)
        else:
            ops.fix(t, 1, 1, 0, 0, 1, 1)
            ops.sp(t, 3, 0.5 * KAPPA * y * y)
            ops.sp(t, 4, KAPPA * y)
    _static_solve("Lagrange")
    return slaves


def test_hermite_along_interface_bending_exact():
    slaves = _mesh_hermite("-hermite")
    tol = 1e-9 * KAPPA
    for t, (x, y) in slaves.items():
        assert abs(ops.nodeDisp(t, 3) - 0.5 * KAPPA * y * y) <= tol, f"slave {t} w wrong"
        assert abs(ops.nodeDisp(t, 4) - KAPPA * y) <= tol, f"slave {t} theta_x wrong"
        for d in (1, 2, 5, 6):
            assert abs(ops.nodeDisp(t, d)) <= tol, f"slave {t} spurious DOF {d}"


def test_linear_tie_gives_interpolant_along_interface():
    _mesh_hermite()
    w_linear = 0.5 * (0.0 + 0.5 * KAPPA)               # kappa / 4
    w_exact = 0.5 * KAPPA * 0.25                        # kappa / 8
    w_tie = ops.nodeDisp(22, 3)
    assert abs(w_tie - w_linear) <= 1e-9 * KAPPA
    assert abs(w_tie - w_exact) > 0.9 * abs(w_linear - w_exact)


def test_hermite_refused_with_mortar():
    m_nodes, s_nodes, sf, mf, nsq, nmq = _overlap_shells()
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-mortar", "-hermite", "-slaveFacets", 4, nsq, *sf,
                       "-masterFacets", 4, nmq, *mf, "-outward", 0.0, 0.0, 1.0)


def test_hermite_refused_with_partial_dofs():
    _mesh_butt()
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-slaveNodes", 3, *SLAVE_IFACE, "-masterFacets", 4, 1, *MASTER_FACET,
                       "-hermite", "-dof", 3, 1, 2, 3)


def test_hermite_refused_for_interior_projection():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y, z) in _MASTER.items():
        ops.node(t, float(x), float(y), float(z))
    ops.node(99, 0.5, 0.5, 0.0)                        # centre of the master facet
    _shell_section()
    ops.element("ShellMITC4", 1, *MASTER_FACET, 1)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-slaveNodes", 1, 99, "-masterFacets", 4, 1, *MASTER_FACET, "-hermite")


def test_rotational_tie_accepted_without_rotary_mass():
    """Default 6-DOF collocation and mortar ties on ShellMITC4 nodes with no
    rotational mass are accepted (Lagrange and Penalty need no mass)."""
    _mesh_butt()
    ops.LadrunoTie("-slaveNodes", 3, *SLAVE_IFACE, "-masterFacets", 4, 1, *MASTER_FACET)
    m_nodes, s_nodes, sf, mf, nsq, nmq = _overlap_shells()
    ops.LadrunoTie("-mortar", "-slaveFacets", 4, nsq, *sf, "-masterFacets", 4, nmq, *mf,
                   "-outward", 0.0, 0.0, 1.0)
