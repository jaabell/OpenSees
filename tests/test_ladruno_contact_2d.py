"""2D (-ndm 2) contact: rigid line, node-to-segment, corners, friction, passenger DOFs.

Checks:

- Rigid line (contactPlane nx ny px py kn): a node seeded below y = 0 and loaded by P settles at
  -P/kn (1e-12); the 9-value form with nz = pz = 0 gives the same answer bitwise; nz != 0 is refused.
- Node-to-segment on a 2-node master line: penetration -P/kn (1e-12); a slave node with ndf 3 (the
  third DOF fixed, carried as a passenger) gives the bitwise same answer as ndf 2.
- Orientation: a flush interface (slave on the master line) without -outward is refused at
  analysis time; a slave 0.5 above the master closes onto it with the analytic spring-contact
  equilibrium (1e-9), and the orientation vote, -outward 0 1 and -outward winding give bitwise
  equal answers.
- A master listed out of chain order (a node shared by two segments that are not consecutive head
  to tail) is refused; -outward winding with a disjoint master is refused; -outward winding with
  -mortar or together with -outward ox oy, three -outward components on 2D surfaces and
  -thickness without -mortar are refused.
- Concave corner: a slave pushed diagonally into the re-entrant corner of a floor and a wall,
  where both segments refuse it, is held by the vertex pair at the radial penetration
  (P + ks d0)/(kn + ks) (1e-9) along the diagonal (1e-9).
- Friction (scalar return map along the segment tangent): a sub-cone tangential load sticks at
  Q/(k + kt) and a load above the cone slips to (Q - mu N)/k (1e-9), with the same answer on
  ProfileSPD and FullGeneral.

Every test fails on the base with an unknown command.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


KN = 1.0e6
SEED = 1.0e-8    # seeded penetration: the pair is active from the first iteration
P = 10.0


def _static(system="FullGeneral", tol=1.0e-12, maxit=30):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", tol, maxit, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _analyze_rc():
    try:
        return ops.analyze(1)
    except Exception:
        return -1


# ----------------------------------------------------------------- rigid line
def _rigid_line(nine=False, nz=0.0):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, -SEED)
    ops.fix(1, 1, 0)
    ops.contactSurface(20, "-slave", 1)
    if nine:
        ops.contactPlane(1, 20, 0.0, 1.0, nz, 0.0, 0.0, 0.0, KN)
    else:
        ops.contactPlane(1, 20, 0.0, 1.0, 0.0, 0.0, KN)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P)
    _static()
    assert ops.analyze(1) == 0
    return ops.nodeDisp(1, 2) - SEED


def test_rigid_line_penetration():
    y = _rigid_line()
    assert abs(y + P / KN) < 1.0e-12


def test_rigid_line_nine_value_form():
    assert _rigid_line(nine=True) == _rigid_line()


def test_rigid_line_nonzero_nz_is_refused():
    with pytest.raises(ops.OpenSeesError):
        _rigid_line(nine=True, nz=1.0)


# ------------------------------------------------------- node-to-segment
def _flat_master():
    for t, x in ((1, -1.0), (2, 1.0)):
        ops.node(t, x, 0.0)
        ops.fix(t, 1, 1)
    ops.contactSurface(10, "-master", 2, 1, 2)


def _seeded_slave(ndf=2):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _flat_master()
    ops.model("basic", "-ndm", 2, "-ndf", ndf)
    ops.node(5, 0.0, -SEED)
    if ndf == 2:
        ops.fix(5, 1, 0)
    else:
        ops.fix(5, 1, 0, 1)
    ops.contactSurface(20, "-slave", 5)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, *([0.0, -P] + [0.0] * (ndf - 2)))
    _static()
    assert ops.analyze(1) == 0
    return ops.nodeDisp(5, 2) - SEED


def test_segment_penetration():
    y = _seeded_slave()
    assert abs(y + P / KN) < 1.0e-12
    assert abs(ops.ladrunoContactForce(5) - P) < 1.0e-9 * P


def test_passenger_dof_on_ndf3_slave():
    assert _seeded_slave(ndf=3) == _seeded_slave(ndf=2)


def test_flush_interface_without_outward_is_refused():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _flat_master()
    ops.node(5, 0.0, -SEED)
    ops.fix(5, 1, 0)
    ops.contactSurface(20, "-slave", 5)
    ops.contact(1, 10, 20, KN, 0.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, -P)
    _static()
    assert _analyze_rc() < 0


K_SP = 1.0e3
P_CLOSE = 1000.0


def _closing(orient):
    """Slave 0.5 above the master on a spring k, pushed down by P: it closes the gap and
    settles at y = (k 0.5 - P)/(k + kn)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _flat_master()
    ops.node(5, 0.0, 0.5)
    ops.fix(5, 1, 0)
    ops.node(6, 0.0, 0.5)
    ops.fix(6, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, K_SP)
    ops.element("zeroLength", 1, 6, 5, "-mat", 1, "-dir", 2)
    ops.contactSurface(20, "-slave", 5)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, *orient)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, -P_CLOSE)
    _static(tol=1.0e-12, maxit=50)
    assert ops.analyze(1) == 0
    return 0.5 + ops.nodeDisp(5, 2)


def test_orientation_vote_closes_the_gap():
    y = _closing(())
    assert abs(y - (K_SP * 0.5 - P_CLOSE) / (K_SP + KN)) < 1.0e-9


def test_vote_outward_and_winding_agree():
    y_vote = _closing(())
    assert _closing(("-outward", 0.0, 1.0)) == y_vote
    assert _closing(("-outward", "winding")) == y_vote


def _two_segments(chain, winding=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, x in ((1, -1.0), (2, 0.0), (3, 1.0), (4, 2.0)):
        ops.node(t, x, 0.0)
        ops.fix(t, 1, 1)
    ops.node(5, 0.5, -SEED)
    ops.fix(5, 1, 0)
    ops.contactSurface(10, "-master", 2, *chain)
    ops.contactSurface(20, "-slave", 5)
    orient = ("-outward", "winding") if winding else ("-outward", 0.0, 1.0)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, *orient)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, -P)
    _static()
    return _analyze_rc()


def test_chained_master_runs():
    assert _two_segments([1, 2, 2, 3]) == 0
    assert abs(ops.nodeDisp(5, 2) - SEED + P / KN) < 1.0e-12


def test_master_out_of_chain_order_is_refused():
    assert _two_segments([1, 2, 3, 2]) < 0


def test_winding_with_disjoint_master_is_refused():
    assert _two_segments([1, 2, 3, 4], winding=True) < 0


@pytest.mark.parametrize("args", [
    (1, 10, 20, KN, 0.0, 0.0, "-outward", "winding", "-outward", 0.0, 1.0),
    (1, 10, 21, "-mortar", "-epsN", KN, "-outward", "winding"),
    (1, 10, 21, "-mortar", "-epsN", KN, "-outward", 0.0, 1.0, 0.0),   # 3 components in 2D
    (1, 10, 20, KN, 0.0, 0.0, "-thickness", 2.0),                    # -thickness without -mortar
])
def test_2d_option_refusals(args):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _flat_master()
    ops.node(5, -1.0, -SEED)
    ops.node(6, 1.0, -SEED)
    ops.contactSurface(20, "-slave", 5, 6)
    ops.contactSurface(21, "-slave-segments", 2, 5, 6)
    with pytest.raises(ops.OpenSeesError):
        ops.contact(*args)


# ----------------------------------------------------------- concave corner
def test_concave_corner_vertex_pair():
    """Floor (-1,0)->(0,0) and wall (0,0)->(0,1), the slave on the free side x < 0, y > 0
    (-outward -1 1). A slave at (d, -d), behind the corner, is refused by both segments; the
    vertex pair pushes it back along the diagonal."""
    d = 1.0e-6
    ks = 1.0e-3 * KN
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, x, y in ((1, -1.0, 0.0), (2, 0.0, 0.0), (3, 0.0, 1.0)):
        ops.node(t, x, y)
        ops.fix(t, 1, 1)
    ops.node(5, d, -d)
    ops.node(6, d, -d)
    ops.fix(6, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, ks)
    ops.element("zeroLength", 1, 6, 5, "-mat", 1, 1, "-dir", 1, 2)
    ops.contactSurface(10, "-master", 2, 1, 2, 2, 3)
    ops.contactSurface(20, "-slave", 5)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", -1.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, P / math.sqrt(2.0), -P / math.sqrt(2.0))
    _static()
    assert ops.analyze(1) == 0
    x = d + ops.nodeDisp(5, 1)
    y = -d + ops.nodeDisp(5, 2)
    r = math.hypot(x, y)
    r_exp = (P + ks * d * math.sqrt(2.0)) / (KN + ks)
    assert abs(r - r_exp) < 1.0e-9 * r_exp
    assert abs(x + y) < 1.0e-9 * r             # along the diagonal (round-off of d + u)
    assert abs(ops.ladrunoContactForce(5) - KN * r) < 1.0e-9 * KN * r


# ------------------------------------------------------------------ friction
KT = 1.0e4
MU = 0.3
N_LOAD = 100.0
K_X = 100.0


def _friction(Q, system="FullGeneral"):
    """Slave pressed by N on a rigid segment, a spring K_X in x and a tangential load Q."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _flat_master()
    ops.node(5, 0.0, -SEED)
    ops.node(6, 0.0, -SEED)
    ops.fix(6, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, K_X)
    ops.element("zeroLength", 1, 6, 5, "-mat", 1, "-dir", 1)
    ops.contactSurface(20, "-slave", 5)
    ops.contact(1, 10, 20, KN, KT, MU, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, Q, -N_LOAD)
    _static(system=system, tol=1.0e-12, maxit=50)
    assert ops.analyze(1) == 0
    return ops.nodeDisp(5, 1)


def test_friction_stick():
    Q = 10.0                      # < mu N = 30
    assert abs(_friction(Q) - Q / (K_X + KT)) < 1.0e-9 * Q / (K_X + KT)


def test_friction_slip():
    Q = 50.0                      # > mu N = 30
    u = _friction(Q)
    assert abs(u - (Q - MU * N_LOAD) / K_X) < 1.0e-9


def test_friction_symmetric_tangent_is_solver_safe():
    assert abs(_friction(50.0, "ProfileSPD") - _friction(50.0, "FullGeneral")) < 1.0e-12
