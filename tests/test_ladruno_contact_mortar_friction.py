"""Coulomb/Tresca friction on mortar segment-to-segment contact.

Checks:

- Static stick (Tresca cap 50, Q = 20) converges and moves Q/(epsT A) (10 percent); static
  slip on tangential springs gives (Q - cap)/(4 kx) (5 percent); the stick solve is the same
  on ProfileSPD and FullGeneral (1e-9); ModifiedNewton (initial stiffness) converges in
  stick; Coulomb stick converges.
- Tangential augmentation: held-load commits drive the penalty stick creep, which is
  O(1/epsT), below 1e-6 for epsT = 1e4 and 1e6 with the same result (1e-7).
- -consistanttan on a Coulomb slip reaches the root of the symmetric tangent (1e-7) in no
  more iterations; -consistentTan is the same option (bitwise).
- Refusals: -tie with friction, -tauMax without -mu or -cohesion, -epsT without -mortar.

Every test fails on the base with an unknown command.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


EPS = 1.0e4
DELTA = 1.0e-4    # initial normal penetration, contact active from the first iteration
COH = 50.0        # Tresca cap COH * A, independent of the pressure


def _master():
    for t, x, y in [(1, 0, 0), (2, 1, 0), (3, 1, 1), (4, 0, 1)]:
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)


def _slave(z, fixes=(0, 0, 0)):
    for t, x, y in [(5, 0, 0), (6, 1, 0), (7, 1, 1), (8, 0, 1)]:
        ops.node(t, float(x), float(y), z)
        ops.fix(t, *fixes)


def _surfaces():
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)


def _x_springs(kx):
    ops.uniaxialMaterial("Elastic", 1, kx)
    for i, t in enumerate((5, 6, 7, 8)):
        gt = 200 + t
        ops.node(gt, *[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)][i])
        ops.fix(gt, 1, 1, 1)
        ops.element("zeroLength", 300 + t, gt, t, "-mat", 1, "-dir", 1)


def _loads(qx, pz):
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, qx / 4.0, 0.0, -pz / 4.0)


def _solve(system="FullGeneral", algo="Newton", tol=1.0e-10, maxit=60):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", tol, maxit, 0)
    ops.algorithm(algo)
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    return ops.analyze(1)


def _frictional_static(qx, pz=200.0, mu=0.0, cohesion=COH, system="FullGeneral"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master()
    _slave(-DELTA)
    _surfaces()
    opts = ["-mortar", "-epsN", EPS, "-epsT", EPS, "-outward", 0.0, 0.0, 1.0]
    if mu > 0.0:
        opts += ["-mu", mu]
    if cohesion > 0.0:
        opts += ["-cohesion", cohesion]
    ops.contact(1, 1, 2, *opts)
    _loads(qx, pz)
    rc = _solve(system)
    return rc, ops.nodeDisp(5, 1)


def test_static_stick():
    qx = 20.0                                    # below the cap 50
    rc, ux = _frictional_static(qx)
    assert rc == 0
    assert abs(ux) == pytest.approx(qx / EPS, rel=0.10)


def test_static_slip_on_springs():
    qx, pz, kx = 120.0, 200.0, 1.0e3            # above the cap 50
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master()
    _slave(-DELTA)
    _x_springs(kx)
    _surfaces()
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-epsT", EPS, "-cohesion", COH,
                "-outward", 0.0, 0.0, 1.0)
    _loads(qx, pz)
    assert _solve() == 0
    ux = ops.nodeDisp(5, 1)
    assert ux == pytest.approx((qx - COH) / (4.0 * kx), rel=0.05)


def test_symmetric_solver():
    rc_sym, ux_sym = _frictional_static(20.0, system="ProfileSPD")
    rc_full, ux_full = _frictional_static(20.0, system="FullGeneral")
    assert rc_sym == 0 and rc_full == 0
    assert abs(ux_sym - ux_full) < 1.0e-9


def test_modified_newton_stick():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master()
    _slave(-DELTA)
    _surfaces()
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-epsT", EPS, "-cohesion", COH,
                "-outward", 0.0, 0.0, 1.0)
    _loads(20.0, 200.0)
    assert _solve(algo="ModifiedNewton") == 0


def test_coulomb_stick():
    rc, ux = _frictional_static(5.0, mu=0.5, cohesion=0.0)
    assert rc == 0 and abs(ux) < 1.0


def _stick_creep(epsN, augment):
    q = 20.0
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master()
    _slave(-DELTA, fixes=(0, 1, 1))              # x free, normal held at -DELTA
    _surfaces()
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-epsT", epsN, "-cohesion", COH,
                "-outward", 0.0, 0.0, 1.0)
    _loads(q, 0.0)
    assert _solve(tol=1.0e-11, maxit=50) == 0
    u_penalty = ops.nodeDisp(5, 1)
    if not augment:
        return u_penalty, None
    ops.integrator("LoadControl", 0.0)          # each held commit is one tangential update
    for _ in range(15):
        assert ops.analyze(1) == 0
    return u_penalty, ops.nodeDisp(5, 1)


def test_tangential_augmentation_removes_stick_creep():
    pen_lo, aug_lo = _stick_creep(1.0e4, True)
    pen_hi, aug_hi = _stick_creep(1.0e6, True)
    assert pen_lo > 20.0 * pen_hi > 0.0
    assert abs(aug_lo) < 1.0e-6 and abs(aug_hi) < 1.0e-6
    assert abs(aug_lo - aug_hi) < 1.0e-7
    assert abs(aug_lo) < 1.0e-2 * abs(pen_lo)


def _coulomb_slip(flag):
    q, pz, kx, mu, epsN = 120.0, 200.0, 1.0e3, 0.5, 1.0e4
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _master()
    _slave(-DELTA, fixes=(0, 1, 1))
    _x_springs(kx)
    _surfaces()
    opts = ["-mortar", "-epsN", epsN, "-mu", mu, "-epsT", epsN, "-outward", 0.0, 0.0, 1.0]
    if flag:
        opts.append(flag)
    ops.contact(1, 1, 2, *opts)
    _loads(q, pz)
    rc = _solve()
    return rc, ops.nodeDisp(5, 1), ops.testIter()


def test_consistent_tangent_same_root_fewer_iterations():
    rc_sym, ux_sym, it_sym = _coulomb_slip(None)
    rc_con, ux_con, it_con = _coulomb_slip("-consistanttan")
    assert rc_sym == 0 and rc_con == 0
    assert abs(ux_sym - ux_con) < 1.0e-7
    assert it_con <= it_sym


def test_consistentTan_alias():
    a = _coulomb_slip("-consistanttan")
    b = _coulomb_slip("-consistentTan")
    assert a == b


def _two_facets():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(1 + i, float(x), float(y), 0.0)
        ops.node(5 + i, float(x), float(y), 0.0)
    _surfaces()
    ops.contactSurface(3, "-slave", 5, 6, 7, 8)


@pytest.mark.parametrize("args", [
    (1, 1, 2, "-mortar", "-tie", "-epsTie", 1.0e5, "-mu", 0.3),
    (1, 1, 2, "-mortar", "-tie", "-epsTie", 1.0e5, "-cohesion", 1.0),
    (1, 1, 2, "-mortar", "-tie", "-epsTie", 1.0e5, "-tauMax", 2.0),
    (1, 1, 2, "-mortar", "-epsN", 1.0e5, "-tauMax", 2.0),
    (1, 1, 3, 1.0e5, 1.0e4, 0.3, "-epsT", 1.0e5),
])
def test_refusals(args):
    _two_facets()
    with pytest.raises(ops.OpenSeesError):
        ops.contact(*args)
