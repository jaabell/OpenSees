"""Tests for the KinematicCoupling (RBE2-type) and DistributingCoupling
(RBE3-type) elements.

KinematicCoupling: the slaves follow the rigid-body motion of the reference
node, u_i = u_R + theta_R x d_i, d_i = x_i - x_R.
DistributingCoupling: the reference node follows the weighted least-squares
rigid-body fit of the independent nodes, and a load at the reference node is
distributed to them as a statically equivalent force set.

Kinematic checks use a free dependent side: its only stiffness is the
penalty tie, so equilibrium drives the gap to zero for any penalty value.
Force checks fix one side, load the other and read the reactions.

Each group is skipped when the build does not provide the element.
"""
import math
import os
import tempfile

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _has_element(name):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.node(3, 0.0, 1.0, 0.0)
    ops.node(4, 0.0, 0.0, 1.0)
    try:
        ops.element(name, 1, 1, 3, 2, 3, 4, "-k", 1.0e6)
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


HAS_KC = _has_element("KinematicCoupling")
HAS_DC = _has_element("DistributingCoupling")
needs_kc = pytest.mark.skipif(not HAS_KC, reason="this build has no KinematicCoupling element")
needs_dc = pytest.mark.skipif(not HAS_DC, reason="this build has no DistributingCoupling element")

CORNERS = [(1.0, 1.0, 0.0), (-1.0, 1.0, 0.0), (-1.0, -1.0, 0.0), (1.0, -1.0, 0.0)]
HANDLERS = {
    "Transformation": ("Transformation",),
    "Penalty": ("Penalty", 1.0e14, 1.0e14),
}


# --------------------------------------------------------------- helpers
def _face(kind, ref_xyz=(0.0, 0.0, 0.0), node_ndf=3, extra=(), kt=1.0e7):
    """Reference node 1 (6 DOF) coupled to the four corners of a square face
    (nodes 2..5, node_ndf DOFs each)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, *ref_xyz, "-ndf", 6)
    for i, c in enumerate(CORNERS):
        ops.node(2 + i, *c, "-ndf", node_ndf)
    ops.element(kind, 1, 1, 4, 2, 3, 4, 5, *extra, "-k", kt)


def _solve_static(handler="Transformation", nsteps=1):
    ops.constraints(*HANDLERS[handler])
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0


def _pattern():
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)


def _prescribe(node, values):
    for d, v in enumerate(values):
        ops.sp(node, d + 1, v)


def _reactions(nodes):
    ops.reactions()
    return {nd: ops.nodeReaction(nd) for nd in nodes}


# ===================================================== KinematicCoupling
@needs_kc
@pytest.mark.parametrize("handler", sorted(HANDLERS))
def test_kc_rigid_translation(handler):
    _face("KinematicCoupling")
    c = (0.003, -0.002, 0.001)
    _pattern()
    _prescribe(1, (c[0], c[1], c[2], 0.0, 0.0, 0.0))
    _solve_static(handler)
    for i in range(4):
        d = ops.nodeDisp(2 + i)
        for k in range(3):
            assert d[k] == pytest.approx(c[k], rel=1e-7, abs=1e-12)
    assert max(abs(v) for v in ops.eleResponse(1, "gap")) < 1e-9


@needs_kc
@pytest.mark.parametrize("handler", sorted(HANDLERS))
def test_kc_rotation_transport(handler):
    """Offset reference, prescribed rotation phi about z: each slave moves by
    theta x d_i = phi (-d_iy, d_ix, 0)."""
    xoff, phi = 2.0, 1.0e-3
    _face("KinematicCoupling", ref_xyz=(xoff, 0.0, 0.0))
    _pattern()
    _prescribe(1, (0.0, 0.0, 0.0, 0.0, 0.0, phi))
    _solve_static(handler)
    for i, c in enumerate(CORNERS):
        dx, dy = c[0] - xoff, c[1]
        u = ops.nodeDisp(2 + i)
        assert u[0] == pytest.approx(-phi * dy, rel=1e-6, abs=1e-12)
        assert u[1] == pytest.approx(phi * dx, rel=1e-6, abs=1e-12)
        assert abs(u[2]) < 1e-9


@needs_kc
def test_kc_slave_rotations_tied_by_default():
    _face("KinematicCoupling", node_ndf=6)
    phi = 1.0e-3
    _pattern()
    _prescribe(1, (0.0, 0.0, 0.0, 0.0, 0.0, phi))
    _solve_static()
    for i in range(4):
        assert ops.nodeDisp(2 + i)[5] == pytest.approx(phi, rel=1e-6)
    assert ops.eleResponse(1, "tiedDOFs")[0] == pytest.approx(24.0)


@needs_kc
def test_kc_dof_translations_only():
    """-dof 1 2 3 on 6-DOF slaves: translations follow the transport term,
    slave rotations are not tied (fixed here to keep the system regular)."""
    _face("KinematicCoupling", node_ndf=6, extra=("-dof", 1, 2, 3))
    for i in range(4):
        ops.fix(2 + i, 0, 0, 0, 1, 1, 1)
    phi = 1.0e-3
    _pattern()
    _prescribe(1, (0.0, 0.0, 0.0, 0.0, 0.0, phi))
    _solve_static()
    for i, c in enumerate(CORNERS):
        u = ops.nodeDisp(2 + i)
        assert u[0] == pytest.approx(-phi * c[1], rel=1e-6, abs=1e-12)
        assert u[1] == pytest.approx(phi * c[0], rel=1e-6, abs=1e-12)
        assert u[5] == 0.0
    assert ops.eleResponse(1, "tiedDOFs")[0] == pytest.approx(12.0)


@needs_kc
def test_kc_force_transfer():
    _face("KinematicCoupling")
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    Fx = 1000.0
    _pattern()
    ops.load(1, Fx, 0.0, 0.0, 0.0, 0.0, 0.0)
    _solve_static()
    R = _reactions((2, 3, 4, 5))
    assert sum(R[nd][0] for nd in R) == pytest.approx(-Fx, rel=1e-5)


@needs_kc
def test_kc_moment_into_solid_face():
    """A moment at the 6-DOF reference node enters the 3-DOF face as a force
    couple: zero net force, moment balance about the reference node."""
    _face("KinematicCoupling")
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    Mz = 1000.0
    _pattern()
    ops.load(1, 0.0, 0.0, 0.0, 0.0, 0.0, Mz)
    _solve_static()
    R = _reactions((2, 3, 4, 5))
    assert abs(sum(R[nd][0] for nd in R)) < 1e-3
    assert abs(sum(R[nd][1] for nd in R)) < 1e-3
    mz = sum(CORNERS[i][0] * R[2 + i][1] - CORNERS[i][1] * R[2 + i][0] for i in range(4))
    assert mz == pytest.approx(-Mz, rel=1e-5)


@needs_kc
def test_kc_single_rigid_link_3d():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, 1.5, 0.5, 0.0)
    ops.element("KinematicCoupling", 1, 1, 1, 2, "-k", 1.0e7)
    ux, phi = 0.002, 1.0e-3
    _pattern()
    _prescribe(1, (ux, 0.0, 0.0, 0.0, 0.0, phi))
    _solve_static()
    u = ops.nodeDisp(2)
    assert u[0] == pytest.approx(ux - phi * 0.5, rel=1e-6, abs=1e-12)
    assert u[1] == pytest.approx(phi * 1.5, rel=1e-6, abs=1e-12)


@needs_kc
def test_kc_2d_rigid_link():
    """2D: a 3-DOF reference node drives a 3-DOF slave, u = u_R + theta x d."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 2.0, 1.0)
    ops.element("KinematicCoupling", 1, 1, 1, 2, "-k", 1.0e7)
    ux, uy, phi = 0.001, -0.002, 1.0e-3
    _pattern()
    _prescribe(1, (ux, uy, phi))
    _solve_static()
    u = ops.nodeDisp(2)
    assert u[0] == pytest.approx(ux - phi * 1.0, rel=1e-6, abs=1e-12)
    assert u[1] == pytest.approx(uy + phi * 2.0, rel=1e-6, abs=1e-12)
    assert u[2] == pytest.approx(phi, rel=1e-6)


@needs_kc
def test_kc_derived_rotational_penalty():
    """Default K_r = K_t l^2; corners of the unit face have |d_i|^2 = 2."""
    _face("KinematicCoupling", node_ndf=6)
    assert ops.eleResponse(1, "kt")[0] == pytest.approx(1.0e7, rel=1e-12)
    assert ops.eleResponse(1, "kr")[0] == pytest.approx(2.0e7, rel=1e-12)


@needs_kc
def test_kc_self_tie_is_inert():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.element("KinematicCoupling", 1, 1, 2, 1, 2, "-k", 1.0e7)
    assert ops.eleResponse(1, "tiedDOFs")[0] == pytest.approx(0.0)


@needs_kc
def test_kc_extra_dof_slave_needs_dof_list():
    """A slave with ndf 4 (e.g. a u-p node) is refused without -dof and
    accepted with an explicit translation list."""
    def build(with_dof):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
        for i, c in enumerate(CORNERS):
            ops.node(2 + i, *c, "-ndf", 4)
        args = ["KinematicCoupling", 1, 1, 4, 2, 3, 4, 5]
        if with_dof:
            args += ["-dof", 1, 2, 3]
        try:
            ops.element(*args, "-k", 1.0e7)
        except Exception:
            pass
        return 1 in ops.getEleTags()

    assert not build(with_dof=False)
    assert build(with_dof=True)
    assert ops.eleResponse(1, "tiedDOFs")[0] == pytest.approx(12.0)


@needs_kc
def test_kc_augmented_lagrangian_closes_gap():
    """Slave on a ground spring k_s, reference prescribed. With penalty only
    the gap is -k_s/(K_t+k_s) u_R; the per-step multiplier update of
    -enforce al reduces it by k_s/(K_t+k_s) per committed step."""
    def run(enforce, nsteps):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 1.0, 0.0)
        ops.node(3, 1.0, 0.0)
        ops.fix(3, 1, 1, 1)
        ops.uniaxialMaterial("Elastic", 1, 1.0e3)
        ops.element("zeroLength", 2, 3, 2, "-mat", 1, "-dir", 1)
        ops.element("KinematicCoupling", 1, 1, 1, 2, "-k", 1.0e4, "-enforce", enforce)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        _prescribe(1, (0.01, 0.0, 0.0))
        _solve_static(nsteps=nsteps)
        return ops.nodeDisp(2, 1) - 0.01

    g_pen = run("penalty", 8)
    assert g_pen == pytest.approx(-1.0e3 / 1.1e4 * 0.01, rel=1e-8)
    g_al = run("al", 8)
    assert abs(g_al) < 1.0e-6 * abs(g_pen)


# ===================================================== DistributingCoupling
def _prescribe_face(field):
    """Prescribe every translation of every independent node (zeros too), so
    that the reference node carries the only free DOFs."""
    _pattern()
    for nd, u in field.items():
        _prescribe(nd, u)


@needs_dc
def test_dc_weighted_average_translation():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, -1.0, 0.0, 0.0)
    ops.node(3, 1.0, 0.0, 0.0)
    ops.fix(1, 0, 0, 0, 1, 1, 1)
    ops.element("DistributingCoupling", 1, 1, 2, 2, 3, "-w", 0.25, 0.75, "-k", 1.0e7)
    _prescribe_face({2: (0.004, 0.0, 0.0), 3: (0.0, 0.0, 0.0)})
    _solve_static()
    d = ops.nodeDisp(1)
    assert d[0] == pytest.approx(0.25 * 0.004, rel=1e-8)
    assert abs(d[1]) < 1e-12 and abs(d[2]) < 1e-12


@needs_dc
@pytest.mark.parametrize("handler", sorted(HANDLERS))
def test_dc_rotation_fit(handler):
    """A rotation field u_i = phi e_z x r_i on the face: the free reference
    node at the centroid recovers theta_z = phi and no translation."""
    _face("DistributingCoupling")
    phi = 1.0e-3
    _prescribe_face({2 + i: (-phi * c[1], phi * c[0], 0.0) for i, c in enumerate(CORNERS)})
    _solve_static(handler)
    d = ops.nodeDisp(1)
    assert d[5] == pytest.approx(phi, rel=1e-7)
    for k in range(5):
        assert abs(d[k]) < 1e-9


@needs_dc
def test_dc_offset_reference_transport():
    """Reference offset by xoff: under the same rotation field it translates
    by theta x (x_R - x_c) = (0, phi xoff, 0)."""
    xoff, phi = 2.0, 1.0e-3
    _face("DistributingCoupling", ref_xyz=(xoff, 0.0, 0.0))
    _prescribe_face({2 + i: (-phi * c[1], phi * c[0], 0.0) for i, c in enumerate(CORNERS)})
    _solve_static()
    d = ops.nodeDisp(1)
    assert d[5] == pytest.approx(phi, rel=1e-7)
    assert d[1] == pytest.approx(phi * xoff, rel=1e-6)
    assert abs(d[0]) < 1e-9 and abs(d[2]) < 1e-9


@needs_dc
@pytest.mark.parametrize("handler", sorted(HANDLERS))
def test_dc_rigid_translation(handler):
    _face("DistributingCoupling")
    c = (0.003, -0.002, 0.001)
    _prescribe_face({2 + i: c for i in range(4)})
    _solve_static(handler)
    d = ops.nodeDisp(1)
    for k in range(3):
        assert d[k] == pytest.approx(c[k], rel=1e-7, abs=1e-12)
    for k in range(3, 6):
        assert abs(d[k]) < 1e-9
    assert max(abs(v) for v in ops.eleResponse(1, "gap")) < 1e-9


@needs_dc
def test_dc_equal_force_split():
    _face("DistributingCoupling")
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    Fx = 1000.0
    _pattern()
    ops.load(1, Fx, 0.0, 0.0, 0.0, 0.0, 0.0)
    _solve_static()
    R = _reactions((2, 3, 4, 5))
    assert sum(R[nd][0] for nd in R) == pytest.approx(-Fx, rel=1e-5)
    for nd in R:
        assert R[nd][0] == pytest.approx(-Fx / 4.0, rel=1e-4)


@needs_dc
def test_dc_weighted_force_split():
    """Reference at the weighted centroid, weights 1..4: a force F at the
    reference is split as w_i/W F."""
    w = [1.0, 2.0, 3.0, 4.0]
    W = sum(w)
    xc = [sum(wi * c[k] for wi, c in zip(w, CORNERS)) / W for k in range(3)]
    _face("DistributingCoupling", ref_xyz=tuple(xc), extra=("-w", *w))
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    Fx = 1000.0
    _pattern()
    ops.load(1, Fx, 0.0, 0.0, 0.0, 0.0, 0.0)
    _solve_static()
    R = _reactions((2, 3, 4, 5))
    for i in range(4):
        assert R[2 + i][0] == pytest.approx(-Fx * w[i] / W, rel=1e-5)
        assert abs(R[2 + i][1]) < 1e-6 * Fx


@needs_dc
def test_dc_moment_into_solid_face():
    _face("DistributingCoupling")
    for nd in (2, 3, 4, 5):
        ops.fix(nd, 1, 1, 1)
    Mz = 1000.0
    _pattern()
    ops.load(1, 0.0, 0.0, 0.0, 0.0, 0.0, Mz)
    _solve_static()
    R = _reactions((2, 3, 4, 5))
    assert abs(sum(R[nd][0] for nd in R)) < 1e-3
    assert abs(sum(R[nd][1] for nd in R)) < 1e-3
    mz = sum(CORNERS[i][0] * R[2 + i][1] - CORNERS[i][1] * R[2 + i][0] for i in range(4))
    assert mz == pytest.approx(-Mz, rel=1e-5)


@needs_dc
def test_dc_rotation_axes_flat_and_collinear():
    _face("DistributingCoupling")
    assert ops.eleResponse(1, "rotationAxes")[0] == pytest.approx(3.0)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    ops.node(2, -1.0, 0.0, 0.0)
    ops.node(3, 0.0, 0.0, 0.0)
    ops.node(4, 1.0, 0.0, 0.0)
    ops.element("DistributingCoupling", 1, 1, 3, 2, 3, 4, "-k", 1.0e7)
    assert ops.eleResponse(1, "rotationAxes")[0] == pytest.approx(2.0)


@needs_dc
def test_dc_derived_rotational_penalty():
    """Default K_r = K_t sum w_i |r_i|^2 / W = 2 K_t on the unit face."""
    _face("DistributingCoupling")
    assert ops.eleResponse(1, "kt")[0] == pytest.approx(1.0e7, rel=1e-12)
    assert ops.eleResponse(1, "kr")[0] == pytest.approx(2.0e7, rel=1e-12)


@needs_dc
def test_dc_unequal_weights_weighted_mean():
    _face("DistributingCoupling", extra=("-w", 1.0, 2.0, 3.0, 4.0))
    ops.fix(1, 0, 0, 0, 1, 1, 1)
    ux = [0.001, 0.002, 0.003, 0.004]
    _prescribe_face({2 + i: (ux[i], 0.0, 0.0) for i in range(4)})
    _solve_static()
    expected = sum(w * u for w, u in zip([1, 2, 3, 4], ux)) / 10.0
    assert ops.nodeDisp(1)[0] == pytest.approx(expected, rel=1e-7)


@needs_dc
def test_dc_2d_rotation_fit():
    """2D: a 3-DOF reference node on three 2-DOF independents recovers the
    in-plane rotation of a rigid rotation field."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0, "-ndf", 3)
    tri = [(1.0, 0.0), (-0.5, 0.8660254037844386), (-0.5, -0.8660254037844386)]
    for i, c in enumerate(tri):
        ops.node(2 + i, *c)
    ops.element("DistributingCoupling", 1, 1, 3, 2, 3, 4, "-k", 1.0e7)
    phi = 1.0e-3
    _prescribe_face({2 + i: (-phi * c[1], phi * c[0]) for i, c in enumerate(tri)})
    _solve_static()
    d = ops.nodeDisp(1)
    assert d[2] == pytest.approx(phi, rel=1e-6)
    assert abs(d[0]) < 1e-9 and abs(d[1]) < 1e-9
    assert ops.eleResponse(1, "rotationAxes")[0] == pytest.approx(1.0)


@needs_dc
def test_dc_augmented_lagrangian_closes_gap():
    """Reference node on a ground spring, independents prescribed. Penalty
    leaves u_R = K_t/(K_t+k_s) u; -enforce al converges to u."""
    def run(enforce, nsteps):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(10, 0.0, 0.0)
        ops.fix(10, 1, 1, 1)
        tri = [(1.0, 0.0), (-0.5, 0.8660254037844386), (-0.5, -0.8660254037844386)]
        for i, c in enumerate(tri):
            ops.node(2 + i, *c)
        ops.uniaxialMaterial("Elastic", 1, 1.0e3)
        ops.element("zeroLength", 2, 10, 1, "-mat", 1, "-dir", 1)
        ops.element("DistributingCoupling", 1, 1, 3, 2, 3, 4, "-k", 1.0e4,
                    "-enforce", enforce)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for i in range(3):
            _prescribe(2 + i, (0.01, 0.0, 0.0))
        _solve_static(nsteps=nsteps)
        return ops.nodeDisp(1, 1) - 0.01

    g_pen = run("penalty", 8)
    assert g_pen == pytest.approx(-1.0e3 / 1.1e4 * 0.01, rel=1e-8)
    g_al = run("al", 8)
    assert abs(g_al) < 1.0e-6 * abs(g_pen)


# ===================================================== shared behaviour
def _skip_if_missing(name):
    if name == "KinematicCoupling" and not HAS_KC:
        pytest.skip("this build has no KinematicCoupling element")
    if name == "DistributingCoupling" and not HAS_DC:
        pytest.skip("this build has no DistributingCoupling element")


COUPLINGS = ["KinematicCoupling", "DistributingCoupling"]
RAYLEIGH_ON = (0.0, 1.0e-3, 1.0e-3, 1.0e-3)   # no alphaM: the nodes carry the mass
RAYLEIGH_OFF = (0.0, 0.0, 0.0, 0.0)


def _transient(name, rayleigh, nsteps=20, dt=1.0e-3):
    """Massed 6-DOF reference node coupled to a fixed face; Newmark."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.mass(1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    for i, c in enumerate(CORNERS):
        ops.node(2 + i, *c)
        ops.fix(2 + i, 1, 1, 1, 1, 1, 1)
    ops.element(name, 1, 1, 4, 2, 3, 4, 5, "-k", 1.0e4)
    _pattern()
    ops.load(1, 100.0, -50.0, 25.0, 10.0, -20.0, 30.0)
    ops.rayleigh(*rayleigh)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 30)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    hist = []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        hist.append(tuple(ops.nodeVel(1)))
    return hist


@pytest.mark.parametrize("name", COUPLINGS)
def test_rayleigh_factors_are_not_applied(name):
    """The coupling has no mass and no damping: a run with stiffness
    proportional Rayleigh factors equals the undamped run exactly, and the
    element damping force is zero while moving."""
    _skip_if_missing(name)
    off = _transient(name, RAYLEIGH_OFF)
    on = _transient(name, RAYLEIGH_ON)
    assert max(abs(v) for step in off for v in step) > 0.0
    assert on == off
    damp = ops.eleResponse(1, "dampingForce")
    assert len(damp) > 0 and all(x == 0.0 for x in damp)
    assert all(math.isfinite(v) for v in on[-1])


@pytest.mark.parametrize("name", COUPLINGS)
def test_database_roundtrip(name):
    """Save and restore through the database: the broker recreates the
    element and recvSelf restores the augmented Lagrangian multipliers."""
    _skip_if_missing(name)

    def skeleton():
        _face(name, extra=("-enforce", "al"), kt=1.0e4)
        for nd in (2, 3, 4, 5):
            ops.fix(nd, 1, 1, 1)

    def probe():
        return list(ops.eleResponse(1, "lambda")) + [ops.eleResponse(1, "kr")[0]]

    skeleton()
    _pattern()
    ops.load(1, 100.0, 50.0, 0.0, 0.0, 0.0, 20.0)
    _solve_static(nsteps=2)
    before = probe()
    assert any(abs(v) > 0.0 for v in before[:-1])

    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "coupling_db")
        try:
            ops.database("File", path)
        except Exception as exc:
            pytest.skip(f"database not available in this build: {exc}")
        ops.save(1)
        ops.wipe()
        skeleton()
        ops.database("File", path)
        ops.restore(1)
        after = probe()
        ops.wipe()
    assert after == pytest.approx(before, rel=1e-12, abs=1e-12)
