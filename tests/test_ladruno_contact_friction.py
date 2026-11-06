"""Coulomb and Tresca friction on node-to-segment (NTS) penalty contact, 3D.

Checks:

- Block on a 30 degree incline under gravity slides at g(sin - mu cos)
  (3 percent; mu = 0 gives g sin); a driven block slides at (Q - mu N)/m and
  the friction force equals mu N for two drives (2 percent).
- Below the cone the block is held within the undamped overshoot 2Q/kt; with
  friction the exit speed is lower and Q x = m v^2/2 + mu N x (5 percent).
- mu does not perturb the static normal penetration P/kn (1e-6, and 1e-9
  between mu = 0 and mu = 0.5).
- Implicit: a static stick converges to Q/kt (1e-4) on FullGeneral and
  ProfileSPD; a static slip on an anchor spring converges to
  (Q - mu N)/ka (1e-3) with friction mu N, also with -consistanttan (no more
  iterations than the symmetric tangent, same root to 1e-6) and on
  ProfileSPD; Newmark sliding acceleration (Q - mu N)/m (5 percent).
- The explicit sliding acceleration is (Q - mu N)/m (2 percent).
- A slave that drifts tangentially while separated feels no friction at
  engagement (acceleration below 5 percent of kt * drift / m).
- Commits are counted per step; a failed step fires the revert; a rejected
  first-engagement step leaves no trace (1e-12); ops.reset() drops the
  friction history (1e-12); wipe and re-run reproduces the result (1e-6).
- Friction drags a deformable stdBrick master in the slide direction.
- Tresca cap: with -mu 0 -cohesion c the slip force equals c; with
  -mu m -cohesion c -tauMax t it equals min(m N + c, t); the anchor reaction
  balances the applied load (1e-3). Below c the interface sticks.
- Refusals: kt < 0, and -tauMax without -mu or -cohesion.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


KN = 1.0e6
KT = 1.0e5


# ---------------------------------------------------------------- helpers
def _flat_master(L=20.0, base=101):
    """Four fixed corner nodes of a flat square facet at z = 0 (normal +z)."""
    h = L / 2.0
    tags = []
    for i, (x, y) in enumerate([(-h, -h), (h, -h), (h, h), (-h, h)]):
        ops.node(base + i, x, y, 0.0)
        ops.fix(base + i, 1, 1, 1)
        tags.append(base + i)
    return tags


def _tilted_master(theta, L=1.0, base=101):
    """Fixed facet rotated about x by theta; normal (0, -sin, cos)."""
    c, s = math.cos(theta), math.sin(theta)
    h = L / 2.0
    tags = []
    for i, (x, y0) in enumerate([(-h, -h), (h, -h), (h, h), (-h, h)]):
        ops.node(base + i, x, y0 * c, y0 * s)
        ops.fix(base + i, 1, 1, 1)
        tags.append(base + i)
    return tags


def _explicit(system="Diagonal"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.integrator("NewmarkExplicit", 0.5)
    ops.algorithm("Linear")
    ops.analysis("Transient")


def _static(tol, iters, system="FullGeneral"):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _contact_on_slave(contact_args, slave=1):
    ops.contactSurface(10, "-master", 4, *_flat_master())
    ops.contactSurface(20, "-slave", slave)
    ops.contact(1, 10, 20, *contact_args)


# ======================================================= explicit Coulomb
def test_incline_sliding_acceleration():
    g, m, theta = 9.81, 2.0, math.radians(30.0)
    c, s = math.cos(theta), math.sin(theta)
    n = (0.0, -s, c)
    down = (0.0, -c, -s)
    W = m * g
    dt, k1, k2 = 1.0e-3, 80, 180

    def slope_accel(mu):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        mt = _tilted_master(theta)
        dn = W * c / KN                    # normal pre-penetration
        ops.node(1, 0.0, dn * s, -dn * c)
        ops.mass(1, m, m, m)
        ops.contactSurface(10, "-master", 4, *mt)
        ops.contactSurface(20, "-slave", 1)
        ops.contact(1, 10, 20, KN, KT, mu, "-outward", *n)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, 0.0, 0.0, -W)
        _explicit()
        v = []
        for _ in range(k2 + 1):
            assert ops.analyze(1, dt) == 0
            v.append(sum(ops.nodeVel(1, d + 1) * down[d] for d in range(3)))
        return (v[k2] - v[k1]) / ((k2 - k1) * dt)

    a_fric, a_free = slope_accel(0.3), slope_accel(0.0)
    a_th_fric, a_th_free = g * (s - 0.3 * c), g * s
    assert abs(a_fric - a_th_fric) / a_th_fric < 0.03, (a_fric, a_th_fric)
    assert abs(a_free - a_th_free) / a_th_free < 0.03, (a_free, a_th_free)


def _driven_block(mu, Q, P=1000.0, m=1.0, nsteps=150, k1=40, k2=140):
    """Slave pre-penetrated to z = -P/kn, held by -P, driven by +Q in x,
    explicit. Returns (a_x, final vx, final x, max |x|)."""
    dt = 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, -P / KN)
    ops.mass(1, m, m, m)
    _contact_on_slave([KN, KT, mu, "-outward", 0.0, 0.0, 1.0])
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    _explicit()
    vx, xmax = [], 0.0
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0
        vx.append(ops.nodeVel(1, 1))
        xmax = max(xmax, abs(ops.nodeDisp(1, 1)))
    return (vx[k2] - vx[k1]) / ((k2 - k1) * dt), vx[-1], ops.nodeDisp(1, 1), xmax


def test_slip_force_capped_at_mu_n():
    P, m, mu = 1000.0, 1.0, 0.3
    for Q, a_th in [(600.0, 300.0), (900.0, 600.0)]:
        a_x, _, _, _ = _driven_block(mu, Q, P=P, m=m)
        assert abs(a_x - a_th) / a_th < 0.02, (Q, a_x)
        assert abs((Q - m * a_x) - mu * P) / (mu * P) < 0.02, (Q, a_x)


def test_stick_below_cone():
    Q = 100.0                              # mu N = 300 > 2Q
    _, _, _, xmax = _driven_block(0.3, Q, nsteps=300, k1=100, k2=200)
    assert xmax < 1.2 * (2.0 * Q / KT), xmax


def test_friction_dissipates_energy():
    P, m, mu, Q = 1000.0, 1.0, 0.3, 600.0
    _, v_fric, x_fric, _ = _driven_block(mu, Q, P=P, m=m)
    _, v_free, _, _ = _driven_block(0.0, Q, P=P, m=m)
    assert v_fric < v_free
    w_diss = Q * x_fric - 0.5 * m * v_fric ** 2
    assert abs(w_diss - mu * P * x_fric) / (mu * P * x_fric) < 0.05, w_diss


def test_explicit_sliding_acceleration():
    a_x, _, _, _ = _driven_block(0.3, 600.0, nsteps=150, k1=40, k2=140)
    assert abs(a_x - 300.0) / 300.0 < 0.02, a_x


def test_wipe_and_rerun_reproduces():
    a1, _, _, _ = _driven_block(0.3, 600.0)
    a2, _, _, _ = _driven_block(0.3, 600.0)
    assert abs(a1 - 300.0) / 300.0 < 0.02 and abs(a2 - 300.0) / 300.0 < 0.02
    assert abs(a1 - a2) < 1.0e-6, (a1, a2)


def test_commit_counted_every_step():
    nsteps = 50
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, -1000.0 / KN)
    ops.mass(1, 1.0, 1.0, 1.0)
    _contact_on_slave([KN, KT, 0.3, "-outward", 0.0, 0.0, 1.0])
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 600.0, 0.0, -1000.0)
    _explicit()
    for _ in range(nsteps):
        assert ops.analyze(1, 1.0e-3) == 0
    info = ops.ladrunoContactInfo()        # [numContacts, numCommits, numReverts]
    assert info[0] == 1 and info[1] == nsteps, info


def test_no_friction_at_engagement_after_separated_drift():
    """The slave flies in +x while separated, then lands. At the step whose
    residual first sees the penetration the tangential acceleration is ~0;
    measuring slip from the original position instead of the engagement point
    would give about kt * drift / m."""
    m, gap0, V, Vz, dt = 1.0, 0.01, 0.5, 2.0, 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, gap0)
    ops.mass(1, m, m, m)
    _contact_on_slave([KN, KT, 0.3, "-outward", 0.0, 0.0, 1.0])
    ops.setNodeVel(1, 1, V, "-commit")
    ops.setNodeVel(1, 3, -Vz, "-commit")
    _explicit()
    ax, x_eng = None, None
    for _ in range(60):
        assert ops.analyze(1, dt) == 0
        if gap0 + ops.nodeDisp(1, 3) < 0.0:
            ax = ops.nodeAccel(1, 1)
            x_eng = ops.nodeDisp(1, 1)
            break
    assert ax is not None
    assert x_eng > 1.0e-4, x_eng
    assert abs(ax) < 0.05 * KT * x_eng / m, (ax, KT * x_eng / m)


def test_failed_step_fires_revert():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, -1.0e-3)
    ops.mass(1, 1.0, 1.0, 1.0)
    _contact_on_slave([KN, KT, 0.3, "-outward", 0.0, 0.0, 1.0])
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 600.0, 0.0, -1000.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.test("NormDispIncr", 1e-30, 1, 0)  # cannot converge
    ops.algorithm("Newton")
    ops.analysis("Transient")
    assert ops.analyze(1, 1.0e-3) != 0
    assert ops.ladrunoContactInfo()[2] >= 1


def test_static_mu_does_not_change_normal_penetration():
    def pen(mu):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.node(1, 0.0, 0.0, -1.0e-8)
        ops.fix(1, 1, 1, 0)
        ops.contactSurface(10, "-master", 4, *_flat_master(L=2.0))
        ops.contactSurface(20, "-slave", 1)
        ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, 0.0, 0.0, -1.0e3)
        _static(1e-12, 40)
        assert ops.analyze(1) == 0
        return -(-1.0e-8 + ops.nodeDisp(1, 3))

    pen_mu0, pen_mu5 = pen(0.0), pen(0.5)
    expect = 1.0e3 / KN
    for p in (pen_mu0, pen_mu5):
        assert p > 0 and abs(p - expect) / expect < 1e-6, p
    assert abs(pen_mu0 - pen_mu5) / expect < 1e-9


def test_friction_drags_deformable_master():
    E, nu, rho, L = 2.0e4, 0.0, 1.0, 1.0
    mu, P, Q, m = 0.3, 1.0e3, 6.0e2, 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    pts = [(0, 0, 0), (L, 0, 0), (L, L, 0), (0, L, 0),
           (0, 0, L), (L, 0, L), (L, L, L), (0, L, L)]
    tags = []
    for i, (x, y, z) in enumerate(pts):
        ops.node(101 + i, float(x), float(y), float(z))
        tags.append(101 + i)
    for t in tags[:4]:
        ops.fix(t, 1, 1, 1)
    top = tags[4:]
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    ops.element("stdBrick", 1, *tags, 1)
    ops.node(1, L / 2, L / 2, L - P / KN)
    ops.mass(1, m, m, m)
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    ops.rayleigh(0.0, 0.0, 1.0e-4, 0.0)
    _explicit(system="ProfileSPD")
    for _ in range(1500):
        assert ops.analyze(1, 2.0e-4) == 0
    top_x = sum(ops.nodeDisp(t, 1) for t in top) / 4.0
    assert top_x > 1.0e-5, top_x
    assert ops.nodeDisp(1, 1) > top_x


# ======================================================= implicit Coulomb
def _stick_model(mu, Q, P=1000.0, solver="FullGeneral"):
    """Slave at the normal equilibrium z = -P/kn, free x and z, loaded
    (Q, 0, -P) below the cone. Static. Returns the analyze code."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, -P / KN)
    ops.fix(1, 0, 1, 0)
    _contact_on_slave([KN, KT, mu, "-outward", 0.0, 0.0, 1.0])
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    _static(1e-10, 40, system=solver)
    return ops.analyze(1)


@pytest.mark.parametrize("solver", ["FullGeneral", "ProfileSPD"])
def test_static_stick(solver):
    Q = 200.0
    assert _stick_model(0.5, Q, solver=solver) == 0
    assert abs(ops.nodeDisp(1, 1) - Q / KT) / (Q / KT) < 1e-4
    assert abs(ops.nodeDisp(1, 3)) < 1e-6


def _anchored_slave(contact_args, Q, P, ka=1.0e4, z0=-1.0e-8,
                    solver="FullGeneral", tol=1e-9):
    """Slave free in x and z with an x anchor spring ka, loaded (Q, 0, -P).
    Static. Returns (rc, x, friction force Q - ka x, iterations, anchor
    reaction)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, z0)
    ops.fix(1, 0, 1, 0)
    ops.node(2, 0.0, 0.0, z0)
    ops.fix(2, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, ka)
    ops.element("zeroLength", 99, 2, 1, "-mat", 1, "-dir", 1)
    _contact_on_slave(contact_args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    _static(tol, 60, system=solver)
    rc = ops.analyze(1)
    it = ops.testIter()
    x = ops.nodeDisp(1, 1)
    ops.reactions()
    return rc, x, Q - ka * x, it, ops.nodeReaction(2, 1)


@pytest.mark.parametrize("solver", ["FullGeneral", "ProfileSPD"])
def test_static_slip_with_anchor(solver):
    mu, P, Q, ka = 0.5, 1000.0, 800.0, 1.0e4
    rc, x, fric, _, _ = _anchored_slave(
        [KN, KT, mu, "-outward", 0.0, 0.0, 1.0], Q, P, ka, solver=solver)
    assert rc == 0
    expect = (Q - mu * P) / ka
    assert abs(x - expect) / expect < 1e-3, x
    assert abs(fric - mu * P) / (mu * P) < 1e-3, fric


def test_consistent_tangent_slip():
    mu, P, Q, ka = 0.5, 1000.0, 800.0, 1.0e4
    base = [KN, KT, mu, "-outward", 0.0, 0.0, 1.0]
    rc_c, x_c, f_c, it_c, _ = _anchored_slave(base + ["-consistanttan"], Q, P, ka)
    rc_s, x_s, _, it_s, _ = _anchored_slave(base, Q, P, ka)
    assert rc_c == 0 and rc_s == 0
    expect = (Q - mu * P) / ka
    assert abs(x_c - expect) / expect < 1e-3
    assert abs(x_c - x_s) / expect < 1e-6
    assert abs(f_c - mu * P) / (mu * P) < 1e-3
    assert it_c <= it_s, (it_c, it_s)


def test_newmark_sliding_acceleration():
    mu, P, Q, m, dt = 0.3, 1000.0, 600.0, 1.0, 1.0e-3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, -P / KN)
    ops.fix(1, 0, 1, 0)
    ops.mass(1, m, m, m)
    _contact_on_slave([KN, KT, mu, "-outward", 0.0, 0.0, 1.0])
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.test("NormDispIncr", 1e-8, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Transient")
    vx = []
    for _ in range(120):
        assert ops.analyze(1, dt) == 0
        vx.append(ops.nodeVel(1, 1))
    a = (vx[110] - vx[40]) / (70 * dt)
    assert abs(a - 300.0) / 300.0 < 0.05, a


# ============================================== revert and reset history
def _implicit_drag(fail_first):
    """Slave starts separated so friction first engages inside a step. With
    fail_first, one step of 1e-3 is attempted under an unreachable
    tolerance before the 5e-4 steps. Returns the x trajectory."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 1.0e-4)
    ops.mass(1, 1.0, 1.0, 1.0)
    _contact_on_slave([KN, KT, 0.3, "-outward", 0.0, 0.0, 1.0])
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 600.0, 0.0, -1000.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Newton")
    ops.analysis("Transient")
    if fail_first:
        ops.test("NormDispIncr", 1e-30, 1, 0)
        assert ops.analyze(1, 1.0e-3) != 0
        assert ops.ladrunoContactInfo()[2] >= 1
    ops.test("NormDispIncr", 1e-10, 50, 0)
    traj = []
    for _ in range(20):
        assert ops.analyze(1, 0.5e-3) == 0
        traj.append(ops.nodeDisp(1, 1))
    assert 1.0e-4 + ops.nodeDisp(1, 3) < 0.0       # in contact
    return traj


def test_rejected_step_leaves_no_engagement_trace():
    ref = _implicit_drag(False)
    retry = _implicit_drag(True)
    assert max(abs(a - b) for a, b in zip(retry, ref)) < 1.0e-12


def test_reset_drops_friction_history():
    dt = math.sqrt(1.0 / KN)

    def run():
        traj = []
        for _ in range(500):
            assert ops.analyze(1, dt) == 0
            traj.append(ops.nodeDisp(1, 1))
        return traj

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, -1.0e-3)
    ops.mass(1, 1.0, 1.0, 1.0)
    _contact_on_slave([KN, KT, 0.3, "-outward", 0.0, 0.0, 1.0])
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 600.0, 0.0, -1000.0)
    _explicit()
    t1 = run()
    assert t1[-1] > 1.0e-3, t1[-1]                 # it slid
    ops.reset()
    ops.wipeAnalysis()
    _explicit()
    t2 = run()
    assert max(abs(a - b) for a, b in zip(t1, t2)) < 1.0e-12


# ============================================================ Tresca cap
P_N = 1000.0
Q_T = 800.0
K_A = 1.0e4


def _tresca(extra, mu=0.0, Q=Q_T):
    rc, x, fric, _, r_anchor = _anchored_slave(
        [KN, KT, mu, "-outward", 0.0, 0.0, 1.0] + extra, Q, P_N, K_A,
        z0=-P_N / KN)
    assert rc == 0
    assert abs(Q + r_anchor - fric) <= 1.0e-9 * Q    # anchor balance
    return x, fric


def test_cohesion_only_slip_force():
    c = 150.0
    x, fric = _tresca(["-mu", 0.0, "-cohesion", c])
    assert abs(fric - c) / c < 1e-3, fric
    assert abs(x - (Q_T - c) / K_A) / ((Q_T - c) / K_A) < 1e-3, x


def test_cohesion_holds_below_cap():
    c, Q = 150.0, 100.0
    x, fric = _tresca(["-cohesion", c], Q=Q)
    assert abs(x - Q / (KT + K_A)) / (Q / (KT + K_A)) < 1e-3, x
    assert fric < c


@pytest.mark.parametrize("mu, c, tau, cap", [
    (0.5, 100.0, 400.0, 400.0),           # tauMax governs
    (0.5, 100.0, 1.0e4, 600.0),           # mu N + c governs
])
def test_tresca_cap_on_coulomb_cone(mu, c, tau, cap):
    x, fric = _tresca(["-cohesion", c, "-tauMax", tau], mu=mu)
    assert abs(fric - cap) / cap < 1e-3, (fric, cap)
    assert abs(x - (Q_T - cap) / K_A) / ((Q_T - cap) / K_A) < 1e-3, x


def test_mu_option_equals_positional_mu():
    x_opt, f_opt = _tresca(["-mu", 0.5, "-cohesion", 100.0, "-tauMax", 1.0e4])
    x_pos, f_pos = _tresca(["-cohesion", 100.0, "-tauMax", 1.0e4], mu=0.5)
    assert abs(f_opt - 600.0) / 600.0 < 1e-3, f_opt
    assert abs(x_opt - x_pos) <= 1e-12 + 1e-9 * abs(x_pos)


# ============================================================ refusals
def _refused(declare, run):
    """True if `declare` raises, or if it is accepted and `run` fails."""
    try:
        declare()
    except Exception:
        return True
    return run() != 0


def _refusal_rig():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in zip((1, 2, 3, 4), [(-1, -1), (1, -1), (1, 1), (-1, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.node(5, 0.0, 0.0, -1.0e-4)
    ops.mass(5, 1.0, 1.0, 1.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)

    def run():
        _explicit()
        return ops.analyze(1, 1.0e-4)
    return run


def test_negative_kt_refused():
    _refusal_rig()
    with pytest.raises(ops.OpenSeesError):
        ops.contact(7, 1, 2, KN, -1.0, 0.4, "-outward", 0.0, 0.0, 1.0)


def test_tau_max_alone_refused():
    run = _refusal_rig()
    assert _refused(
        lambda: ops.contact(7, 1, 2, KN, KT, 0.0, "-outward", 0.0, 0.0, 1.0,
                            "-tauMax", 50.0),
        run)
