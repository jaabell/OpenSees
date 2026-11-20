"""CentralDifferenceLadruno explicit leap-frog integrator.

Checks, on SDOF systems and truss bars:
  * second-order convergence (log-log slope 1.7 to 2.3);
  * undamped stability boundary omega dt = 2 and the damped limit
    (2/omega)(sqrt(1 + xi^2) - xi) for xi = 0.02, 0.1, 0.5;
  * agreement with Newmark (average acceleration) and CentralDifference within
    3 % at dt = T/200;
  * 1-D wave speed sqrt(E/rho) within 5 %, exact rigid-body motion;
  * criticalTimeStep = le/c on a one-element bar (2 % tolerance);
  * first-step accuracy with nonzero u0 and v0 (error below 1e-3 and not
    larger than CentralDifference);
  * stiffness-proportional damping reduces dt_cr about quadratically while
    mass-proportional damping leaves it unchanged;
  * the mass tangent cache gives the same trajectory bit for bit as
    -noMassCache, through an element removal and a failed-and-retried step;
  * -commitSolveState is bit-identical for rate-independent materials and a
    bounded change for a Maxwell damper;
  * NaN and Inf accelerations stop the step; the -divergence check does not
    trip at kinetic-energy troughs and does trip on growth;
  * revertToLastStep: a failed step followed by a retry at the same dt or at
    dt/2 matches a clean run (1e-9).
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

CDL = "CentralDifferenceLadruno"
TWO_PI = 2.0 * math.pi


# --------------------------------------------------------------------------
# pure-python numeric helpers (numpy-free)
# --------------------------------------------------------------------------
def _rms(xs):
    return math.sqrt(sum(v * v for v in xs) / len(xs))


def _linfit_slope(xs, ys):
    """Least-squares slope of ys vs xs (both 1-D, equal length)."""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / sxx


def _all_finite(xs):
    return all(math.isfinite(v) for v in xs)


def _max_abs(xs):
    return max(abs(v) for v in xs)


# --------------------------------------------------------------------------
# shared models
# --------------------------------------------------------------------------
def _sdof(k, m, integrator_args, dt, nsteps, d0=0.0, v0=0.0):
    """Free-vibration SDOF; returns (t[], u[]). integrator_args = tuple."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, m, 0.0)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    if d0 != 0.0:
        ops.setNodeDisp(2, 1, d0, "-commit")
    if v0 != 0.0:
        ops.setNodeVel(2, 1, v0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    t = [ops.getTime()]
    u = [ops.nodeDisp(2, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        t.append(ops.getTime())
        u.append(ops.nodeDisp(2, 1))
    return t, u


def _stable_at(integrator_args, omega_dt, m=1.0, w=TWO_PI):
    """Run an SDOF at Omega = w*dt; True if it does not blow up."""
    k = m * w * w
    dt = omega_dt / w
    _, u = _sdof(k, m, integrator_args, dt, 200, d0=1.0)
    if not _all_finite(u):
        return False
    early = _max_abs(u[:40])
    late = _max_abs(u[-40:])
    return late <= 5.0 * max(early, 1e-12)


def _build_bar(N, L, E, A, rho):
    """Fixed-free chain of N truss elements along x; returns element length."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    le = L / N
    for i in range(N + 1):
        ops.node(i + 1, i * le, 0.0)
    ops.fix(1, 1, 1)
    for i in range(1, N + 1):
        ops.fix(i + 1, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(N):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", rho)
    return le


# ==========================================================================
# order of accuracy ~ 2.0
# ==========================================================================
def test_cdl_order():
    m = 1.0
    w = TWO_PI
    k = m * w * w
    v0 = 1.0
    Tend = 1.0
    log_dts, log_errs = [], []
    for N in (40, 80, 160, 320, 640):
        dt = Tend / N
        t, u = _sdof(k, m, (CDL,), dt, N, v0=v0)
        err = [ui - (v0 / w) * math.sin(w * ti) for ti, ui in zip(t, u)]
        log_dts.append(math.log(dt))
        log_errs.append(math.log(_rms(err)))
    slope = _linfit_slope(log_dts, log_errs)
    assert 1.7 <= slope <= 2.3, "convergence slope = %.3f (expect ~2.0)" % slope


# ==========================================================================
# undamped stability boundary: stable for w*dt < 2, unstable just above
# ==========================================================================
def test_cdl_stability():
    assert _stable_at((CDL,), 1.9), "should be stable at Omega=1.9 (< 2.0)"
    assert not _stable_at((CDL,), 2.1), "should be unstable at Omega=2.1 (> 2.0)"


# ==========================================================================
# damped stability limit (2/w)(sqrt(1+xi^2)-xi), xi = 0.02/0.1/0.5
# ==========================================================================
def _sdof_damped(k, m, alphaM, dt, nsteps, d0=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, m, 0.0)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.setNodeDisp(2, 1, d0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.rayleigh(alphaM, 0.0, 0.0, 0.0)
    ops.integrator(CDL)
    ops.analysis("Transient")
    u = [ops.nodeDisp(2, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None
        u.append(ops.nodeDisp(2, 1))
    return u


@pytest.mark.parametrize("xi", [0.02, 0.1, 0.5])
def test_cdl_damped(xi):
    m = 1.0
    w = TWO_PI
    k = m * w * w
    alphaM = 2.0 * xi * w
    dt_cr = (2.0 / w) * (math.sqrt(1.0 + xi * xi) - xi)  # damped central-difference limit
    u_lo = _sdof_damped(k, m, alphaM, 0.9 * dt_cr, 200)  # under -> stable
    u_hi = _sdof_damped(k, m, alphaM, 1.1 * dt_cr, 200)  # over  -> blows up
    stable_lo = (
        u_lo is not None
        and _all_finite(u_lo)
        and _max_abs(u_lo[-20:]) <= _max_abs(u_lo[:20])
    )
    unstable_hi = (
        u_hi is None
        or not _all_finite(u_hi)
        or _max_abs(u_hi[-20:]) > 5.0 * _max_abs(u_hi[:20])
    )
    assert stable_lo, "xi=%.2f: should be stable below dt_cr=%.4f" % (xi, dt_cr)
    assert unstable_hi, "xi=%.2f: should be unstable above dt_cr=%.4f" % (xi, dt_cr)


# ==========================================================================
# cross-check vs Newmark (avg accel) and CentralDifference at small dt
# ==========================================================================
def test_cdl_crosscheck():
    m = 1.0
    w = TWO_PI
    k = m * w * w
    d0 = 1.0
    dt = (TWO_PI / w) / 200
    N = 400
    _, ucdl = _sdof(k, m, (CDL,), dt, N, d0=d0)
    _, unm = _sdof(k, m, ("Newmark", 0.5, 0.25), dt, N, d0=d0)
    _, ucd = _sdof(k, m, ("CentralDifference",), dt, N, d0=d0)
    e_nm = max(abs(a - b) for a, b in zip(ucdl, unm)) / d0
    e_cd = max(abs(a - b) for a, b in zip(ucdl, ucd)) / d0
    assert e_nm < 0.03, "max|CDL-Newmark| = %.3f%% (expect < 3%%)" % (100 * e_nm)
    assert e_cd < 0.03, "max|CDL-CD| = %.3f%% (expect < 3%%)" % (100 * e_cd)


# ==========================================================================
# 1-D wave speed exact = sqrt(E/rho)
# ==========================================================================
def test_cdl_wave():
    N = 100
    L = 20.0
    E = 100.0
    A = 1.0
    rho = 1.0
    P = 1.0
    c = math.sqrt(E / rho)
    le = _build_bar(N, L, E, A, rho)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(N + 1, P, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    dt = 0.5 * le / c
    v_particle = P / (A * rho * c)
    vth = 0.5 * v_particle
    dists = [5.0, 8.0, 11.0, 14.0]
    nodes = [int(round((L - d) / le)) + 1 for d in dists]
    arrival = {nd: None for nd in nodes}
    tcur = 0.0
    for _ in range(int(2.2 * L / c / dt)):
        ops.analyze(1, dt)
        tcur += dt
        for nd in nodes:
            if arrival[nd] is None and abs(ops.nodeVel(nd, 1)) > vth:
                arrival[nd] = tcur
    pts = [(d, arrival[nd]) for d, nd in zip(dists, nodes) if arrival[nd]]
    assert len(pts) >= 2, "wave front not detected at >= 2 stations"
    slope = _linfit_slope([p[0] for p in pts], [p[1] for p in pts])  # = 1/speed
    c_meas = 1.0 / slope
    err = abs(c_meas - c) / c
    assert err < 0.05, "measured wave speed %.3f vs sqrt(E/rho)=%.3f (err %.1f%%)" % (
        c_meas,
        c,
        100 * err,
    )


# ==========================================================================
# rigid-body momentum conserved exactly (zero stiffness, nonzero v0)
# ==========================================================================
def test_cdl_rigid_body():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 0, 1)
    m = 2.0
    v0 = 3.0
    ops.mass(1, m, 0.0)
    ops.setNodeVel(1, 1, v0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    dt = 0.01
    N = 100
    for _ in range(N):
        ops.analyze(1, dt)
    u = ops.nodeDisp(1, 1)
    v = ops.nodeVel(1, 1)
    u_exact = v0 * (N * dt)  # free particle: u = v0 t, v = v0
    assert abs(v - v0) < 1e-9, "velocity drift: v=%.9f vs v0=%.1f" % (v, v0)
    assert abs(u - u_exact) < 1e-6, "position: u=%.6f vs v0*t=%.6f" % (u, u_exact)


# ==========================================================================
# criticalTimeStep() = le/c exact on a 1-element bar (2/omega_max)
# ==========================================================================
def test_cdl_dtcr_bar():
    N = 1
    L = 1.0
    E = 100.0
    A = 1.0
    rho = 1.0
    le = _build_bar(N, L, E, A, rho)
    c = math.sqrt(E / rho)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")
    ops.analyze(1, 0.1 * le / c)  # one priming step triggers the dt_cr compute
    dtcr = ops.criticalTimeStep()
    target = le / c  # single free DOF: omega_max = 2c/le -> dt_cr = le/c
    err = abs(dtcr - target) / target
    assert err < 0.02, "criticalTimeStep()=%.5f vs le/c=%.5f (err %.1f%%)" % (
        dtcr,
        target,
        100 * err,
    )


# ==========================================================================
# first-step correctness : nonzero v0 AND nonzero
#        a0 (= -w^2 u0). The starter reproduces u(t)=u0 cos wt + (v0/w) sin wt
#        from step 1, where the legacy CentralDifference carries a 1st-step error.
# ==========================================================================
def test_cdl_first_step():
    m = 1.0
    w = TWO_PI
    k = m * w * w
    u0 = 1.0
    v0 = 0.7  # both nonzero -> a0 = -w^2 u0 != 0
    dt = (TWO_PI / w) / 200
    N = 5  # look at the FIRST few steps only
    t, u = _sdof(k, m, (CDL,), dt, N, d0=u0, v0=v0)
    ue = [u0 * math.cos(w * ti) + (v0 / w) * math.sin(w * ti) for ti in t]
    err_cdl = max(abs(a - b) for a, b in zip(u, ue))
    # legacy class on the same problem: expected larger first-step error
    _, ucd = _sdof(k, m, ("CentralDifference",), dt, N, d0=u0, v0=v0)
    err_cd = max(abs(a - b) for a, b in zip(ucd, ue))
    assert err_cdl < 1e-3, "max|CDL-exact| = %.2e (expect < 1e-3)" % err_cdl
    assert err_cdl <= err_cd, (
        "CDL first-step error %.2e should be <= legacy CD's %.2e" % (err_cdl, err_cd)
    )


# ==========================================================================
# stiffness-proportional damping: on a stiff bar dt_cr collapses ~quadratically
#        as Rayleigh stiffness damping betaK rises, while the same nominal
#        damping applied as mass-proportional alphaM stays benign.
# ==========================================================================
def _bar_dtcr_cdl(extra_rayleigh):
    N = 20
    L = 1.0
    E = 1.0e6
    A = 1.0
    rho = 1.0
    le = _build_bar(N, L, E, A, rho)
    c = math.sqrt(E / rho)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    if extra_rayleigh is not None:
        ops.rayleigh(*extra_rayleigh)
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")
    ops.analyze(1, 0.1 * le / c)  # one priming step triggers the dt_cr compute
    return ops.criticalTimeStep()


def test_cdl_betak_trap():
    dt0 = _bar_dtcr_cdl(None)                      # undamped reference
    dt_bk1 = _bar_dtcr_cdl((0.0, 1e-5, 0.0, 0.0))  # small betaK
    dt_bk2 = _bar_dtcr_cdl((0.0, 4e-5, 0.0, 0.0))  # 4x betaK -> ~1/4 the step
    dt_am = _bar_dtcr_cdl((50.0, 0.0, 0.0, 0.0))   # alphaM -> benign
    # betaK collapses dt_cr ~quadratically (4x betaK -> ~4x smaller dt_cr);
    # alphaM leaves it essentially unchanged.
    assert dt_bk1 < 0.95 * dt0, "betaK should shrink dt_cr (%.3e vs %.3e)" % (dt_bk1, dt0)
    assert dt_bk2 < 0.6 * dt_bk1, (
        "4x betaK should ~quarter dt_cr (%.3e vs %.3e)" % (dt_bk2, dt_bk1)
    )
    assert dt_am > 0.9 * dt0, "alphaM should leave dt_cr ~unchanged (%.3e vs %.3e)" % (
        dt_am,
        dt0,
    )




# ==========================================================================
# mass tangent cache: bit-identical to -noMassCache
# ==========================================================================
_MC_N = 20
_MC_E, _MC_A, _MC_RHO = 1.0e4, 1.0, 1.0
_MC_TIP = _MC_N + 1
_MC_TWIN = 1000             # parallel twin of element N//2, removable
_MC_DT = 0.004              # below dt_cr of the bar


def _mc_build(int_args=()):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.uniaxialMaterial("Elastic", 1, _MC_E)
    for n in range(1, _MC_N + 2):
        ops.node(n, float(n - 1))
    ops.fix(1, 1)
    for e in range(1, _MC_N + 1):
        ops.element("Truss", e, e, e + 1, _MC_A, 1, "-rho", _MC_RHO)
    mid = _MC_N // 2
    ops.element("Truss", _MC_TWIN, mid, mid + 1, _MC_A, 1, "-rho", _MC_RHO)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_MC_TIP, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl", *int_args)
    ops.analysis("Transient")


def _mc_steps(n, traj):
    for _ in range(n):
        assert ops.analyze(1, _MC_DT) == 0
        traj.append(ops.nodeDisp(_MC_TIP, 1))


def _mc_fail_once():
    # the -divergence check stops an oversized step within a few attempts;
    # the failed analyze() calls Domain::revertToLastCommit and
    # revertToLastStep
    for _ in range(25):
        if ops.analyze(1, 60.0 * _MC_DT) != 0:
            return True
    return False


def _assert_identical(a, b, label):
    assert len(a) == len(b)
    bad = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    assert not bad, "%s: trajectories differ at steps %s" % (label, bad[:5])


def test_mass_cache_plain_bit_identity():
    runs = []
    for args in ((), ("-noMassCache",)):
        _mc_build(args)
        traj = []
        _mc_steps(300, traj)
        runs.append(traj)
    _assert_identical(runs[0], runs[1], "plain")
    assert max(abs(u) for u in runs[0]) > 0.0


def test_mass_cache_element_removal():
    runs = []
    for args in ((), ("-noMassCache",)):
        _mc_build(args)
        traj = []
        _mc_steps(100, traj)
        ops.remove("element", _MC_TWIN)       # domain change: mass reassembled
        _mc_steps(200, traj)
        runs.append(traj)
    _assert_identical(runs[0], runs[1], "element removal")
    assert runs[0][-1] != runs[0][50]


def test_mass_cache_revert_and_retry():
    runs = []
    for args in ((), ("-noMassCache",)):
        _mc_build(("-divergence", 5.0) + args)
        traj = []
        _mc_steps(100, traj)
        assert _mc_fail_once(), "the oversized step must fail"
        _mc_steps(200, traj)
        runs.append(traj)
    _assert_identical(runs[0], runs[1], "revert and retry")


# ==========================================================================
# -commitSolveState
# ==========================================================================
def _cs_build(int_args, viscous=None):
    """Fixed-free truss bar. viscous: None (elastic), 'zl-damper' or
    'truss-damper' (Maxwell ViscousDamper, a rate material with internal
    state), 'zl-viscous' (memoryless Viscous dashpot)."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.uniaxialMaterial("Elastic", 1, _MC_E)
    ops.uniaxialMaterial("Viscous", 2, 5.0, 1.0)
    ops.uniaxialMaterial("ViscousDamper", 3, 1.0e4, 5.0, 1.0)
    for n in range(1, _MC_N + 2):
        ops.node(n, float(n - 1))
    ops.fix(1, 1)
    n_ele = _MC_N - 1 if viscous == "truss-damper" else _MC_N
    for e in range(1, n_ele + 1):
        ops.element("Truss", e, e, e + 1, _MC_A, 1, "-rho", _MC_RHO)
    if viscous == "truss-damper":
        ops.element("Truss", _MC_N, _MC_N, _MC_N + 1, _MC_A, 3, "-rho", _MC_RHO)
    elif viscous in ("zl-damper", "zl-viscous"):
        ops.node(900, float(_MC_N))
        ops.fix(900, 1)
        mat = 3 if viscous == "zl-damper" else 2
        ops.element("zeroLength", 901, 900, _MC_TIP, "-mat", mat, "-dir", 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_MC_TIP, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl", *int_args)
    ops.analysis("Transient")


def _cs_run(int_args, viscous=None, nsteps=300):
    _cs_build(int_args, viscous)
    traj = []
    _mc_steps(nsteps, traj)
    return traj


def test_commit_solve_state_rate_independent_bit_identity():
    t_on = _cs_run(("-commitSolveState",))
    t_off = _cs_run(())
    _assert_identical(t_on, t_off, "elastic bar")
    assert max(abs(u) for u in t_on) > 0.0


@pytest.mark.parametrize("path,expect_delta", [
    ("zl-damper", True),
    ("truss-damper", True),
    ("zl-viscous", False),
])
def test_commit_solve_state_rate_dependent(path, expect_delta):
    t_on = _cs_run(("-commitSolveState",), viscous=path)
    t_off = _cs_run((), viscous=path)
    differs = any(x != y for x, y in zip(t_on, t_off))
    if expect_delta:
        # a rate material with internal state commits a v_{n+1/2}-rate state:
        # a change, bounded at 5 % of the peak displacement
        assert differs, "%s: expected a change for a rate material with state" % path
        peak = max(abs(u) for u in t_off)
        delta = max(abs(x - y) for x, y in zip(t_on, t_off))
        assert delta < 0.05 * peak, "%s: change %.3e above 5%% of %.3e" % (path, delta, peak)
    else:
        # a memoryless dashpot never reads its committed state back
        assert not differs, "%s: memoryless dashpot must be identical" % path


def test_commit_solve_state_with_removal_and_revert():
    def run(int_args):
        _cs_build(("-divergence", 5.0) + tuple(int_args))
        ops.element("Truss", 1000, _MC_N // 2, _MC_N // 2 + 1, _MC_A, 1, "-rho", _MC_RHO)
        ops.domainChange()
        traj = []
        _mc_steps(80, traj)
        ops.remove("element", 1000)
        _mc_steps(80, traj)
        assert _mc_fail_once()
        _mc_steps(80, traj)
        return traj

    _assert_identical(run(("-commitSolveState",)), run(()), "removal and revert")


# ==========================================================================
# NaN / Inf accelerations and the -divergence check
# ==========================================================================
def _sdof_seeded(integrator_args, u0=None, v0=None):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, 1.0, 0.0)
    ops.uniaxialMaterial("Elastic", 1, 100.0)     # omega = 10, dt_cr = 0.2
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    if u0 is not None:
        ops.setNodeDisp(2, 1, u0, "-commit")
    if v0 is not None:
        ops.setNodeVel(2, 1, v0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")


def test_nan_acceleration_stops_the_step(capfd):
    # a NaN committed displacement makes F_int and the solved acceleration NaN
    _sdof_seeded((CDL,), u0=float("nan"))
    capfd.readouterr()
    rc = ops.analyze(1, 1.0e-3)
    err = capfd.readouterr().err
    assert rc != 0, "a NaN acceleration must fail the step"
    assert "non-finite" in err
    assert ops.getTime() == pytest.approx(0.0, abs=1e-15)


def test_inf_acceleration_stops_the_step(capfd):
    # dt = 0.21 > dt_cr = 0.2: the amplitude overflows to Inf
    _sdof_seeded((CDL,), u0=0.1)
    capfd.readouterr()
    rc = 0
    for _ in range(3000):
        rc = ops.analyze(1, 0.21)
        if rc != 0:
            break
    err = capfd.readouterr().err
    assert rc != 0, "a divergent run must stop"
    assert "non-finite" in err


def test_divergence_check_no_trip_at_kinetic_energy_troughs(capfd):
    _sdof_seeded((CDL, "-divergence", 10.0), v0=1.0)
    capfd.readouterr()
    rc = ops.analyze(400, 0.01)
    err = capfd.readouterr().err
    assert rc == 0 and "ABORT" not in err, err


def test_divergence_check_trips_on_growth(capfd):
    _sdof_seeded((CDL, "-divergence", 10.0), v0=1.0)
    capfd.readouterr()
    rc = ops.analyze(200, 0.25)                   # above dt_cr = 0.2
    err = capfd.readouterr().err
    assert rc != 0 and "kinetic-energy proxy" in err, err


# ==========================================================================
# revertToLastStep
# ==========================================================================
_RV_K, _RV_M = 100.0, 1.0              # omega = 10, dt_cr = 0.2
_RV_OMEGA = math.sqrt(_RV_K / _RV_M)
_RV_DT_OK = 0.005
_RV_DT_BAD = 0.5                       # one step trips the -divergence check
_RV_BEFORE, _RV_AFTER = 200, 400
# Circular motion: equal springs in x and y with u_y0 = A, v_x0 = A omega keep
# the total kinetic energy constant, so only the faulty step trips the check.
_RV_A0 = 0.1


def _rv_build(integrator_args, seed=None):
    if seed is None:
        seed = (0.0, _RV_A0, _RV_A0 * _RV_OMEGA, 0.0, 0.0, -_RV_A0 * _RV_OMEGA ** 2)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.mass(2, _RV_M, _RV_M)
    ops.uniaxialMaterial("Elastic", 1, _RV_K)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, 1, "-dir", 1, 2)
    ux, uy, vx, vy, ax, ay = seed
    for dof, u, v, a in ((1, ux, vx, ax), (2, uy, vy, ay)):
        if u != 0.0:
            ops.setNodeDisp(2, dof, u, "-commit")
        if v != 0.0:
            ops.setNodeVel(2, dof, v, "-commit")
        if a != 0.0:
            ops.setNodeAccel(2, dof, a, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")


def _rv_state():
    return (ops.nodeDisp(2, 1), ops.nodeDisp(2, 2),
            ops.nodeVel(2, 1), ops.nodeVel(2, 2),
            ops.nodeAccel(2, 1), ops.nodeAccel(2, 2))


def _rv_run(nsteps, dt):
    out = []
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        out.append((ops.getTime(), _rv_state()))
    return out


def _rv_match(cont, ref, tol, what):
    assert len(cont) == len(ref), "%s: %d steps vs %d" % (what, len(cont), len(ref))
    for i, ((_, sf), (_, sr)) in enumerate(zip(cont, ref)):
        for got, exp in zip(sf, sr):
            assert abs(got - exp) < tol, (
                "%s: step %d got %.12g expected %.12g" % (what, i, got, exp))


def test_revert_then_retry_at_same_dt():
    args = (CDL, "-divergence", 10.0)
    _rv_build(args)
    before = _rv_run(_RV_BEFORE, _RV_DT_OK)
    assert len(before) == _RV_BEFORE
    t_committed = before[-1][0]
    assert ops.analyze(1, _RV_DT_BAD) != 0
    assert ops.getTime() == pytest.approx(t_committed, abs=1e-12)
    cont = _rv_run(_RV_AFTER, _RV_DT_OK)

    _rv_build(args)
    clean = _rv_run(_RV_BEFORE + _RV_AFTER, _RV_DT_OK)
    assert len(clean) == _RV_BEFORE + _RV_AFTER
    _rv_match(cont, clean[_RV_BEFORE:], 1e-9, "retry at dt")


def test_revert_then_retry_at_half_dt():
    # the retry restarts from the committed state, so the reference is a fresh
    # analysis seeded with that state
    args = (CDL, "-divergence", 10.0)
    dt_retry = _RV_DT_OK / 2.0
    _rv_build(args)
    assert len(_rv_run(_RV_BEFORE, _RV_DT_OK)) == _RV_BEFORE
    committed = _rv_state()
    assert ops.analyze(1, _RV_DT_BAD) != 0
    cont = _rv_run(_RV_AFTER, dt_retry)

    _rv_build(args, seed=committed)
    ref = _rv_run(_RV_AFTER, dt_retry)
    assert len(ref) == _RV_AFTER
    _rv_match(cont, ref, 1e-9, "retry at dt/2")
