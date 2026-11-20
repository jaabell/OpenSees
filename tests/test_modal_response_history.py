"""modalResponseHistory: piecewise-linear modal-superposition transient.

After `eigen` + `modalProperties` each retained mode is a decoupled SDOF
advanced by a constant per-mode recurrence that is exact for a load varying
linearly within each step (Nigam and Jennings 1969); the response is
recovered as u = sum_a (phi_a Vscale) q_a and one domain step is committed
per time station, so ordinary recorders capture the history.

  test_sdof_branches[under/critical/over]
        single-DOF oscillator under base motion against an independent numpy
        exact piecewise-linear integration of u'' + 2 xi w u' + w^2 u = -u_g''
        to 1e-9 (Gamma and the eigenvector scale cancel for one DOF, so this
        checks participation, recovery and recurrence together).
  test_mdof_vs_opensees_newmark
        2-DOF chain against a direct OpenSees Newmark transient with
        mass-proportional damping, 2e-3 (Newmark truncation).
  test_mdof_stiffness_damping_vs_full_matrix
        Rayleigh a0 + a1 damping against an independent numpy full-matrix
        Newmark of (M, a0 M + a1 K, K).
  test_mode_subset_plumbing / test_rayleigh_matches_modal_ratio
        internal consistency.
  -load channel (nodal forces P(t) = s(t) P, absolute response):
        SDOF force history exact to 1e-9, 2-DOF against a direct Newmark run
        driving the same pattern, -load with P = -M R equal to -baseAccel to
        1e-10, and invariance under modalProperties -unorm.
  test_guard_*: refusals.

Element-level stiffness-proportional Rayleigh damping on Truss elements does
not reproduce the assembled C = a1 K (a direct Newmark run with betaK differs
from the exact classical-modal solution by several percent, independent of
dt), so the comparisons with OpenSees Newmark use mass-proportional damping
only and the stiffness-damped case uses a numpy full-matrix Newmark.

Small models use `eigen -fullGenLapack` (the default ARPACK solver requires
NEV < N). Every test fails on a build without modalResponseHistory (unknown
command).
"""
import math

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


# ---------------------------------------------------------------------------
# Independent numpy exact piecewise-linear SDOF oracle
# ---------------------------------------------------------------------------
def _blocks(w, d, dt):
    CRIT = 1e-8
    hd = 0.5 * d
    if w == 0.0:
        if d == 0.0:
            A = np.array([[1.0, dt], [0.0, 1.0]])
        else:
            e = math.exp(-d * dt)
            A = np.array([[1.0, (1.0 - e) / d], [0.0, e]])
    else:
        e = math.exp(-hd * dt)
        Delta = w * w - hd * hd
        if abs(Delta) <= CRIT * (w * w):
            c, sd = 1.0, dt
        elif Delta > 0.0:
            wd = math.sqrt(Delta); c = math.cos(wd * dt); sd = math.sin(wd * dt) / wd
        else:
            wh = math.sqrt(-Delta); c = math.cosh(wh * dt); sd = math.sinh(wh * dt) / wh
        A = e * np.array([[c + hd * sd, sd], [-w * w * sd, c - hd * sd]])
    if w == 0.0 and d == 0.0:
        B = np.array([[dt * dt / 3.0, dt * dt / 6.0], [dt / 2.0, dt / 2.0]])
        return A, B
    B = np.zeros((2, 2))
    for j, (f0, f1) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        r = (f1 - f0) / dt
        if w == 0.0:
            invd, invd2 = 1.0 / d, 1.0 / (d * d)
            qp0 = 0.0
            qpd0 = f0 * invd - r * invd2
            qpDt = (f0 * dt + r * dt * dt * 0.5) * invd - r * dt * invd2
            qpdDt = (f0 + r * dt) * invd - r * invd2
        else:
            w2 = w * w; iw2 = 1.0 / w2; iw4 = iw2 * iw2
            qp0 = f0 * iw2 - d * r * iw4
            qpd0 = r * iw2
            qpDt = (f0 + r * dt) * iw2 - d * r * iw4
            qpdDt = r * iw2
        qh0, qhd0 = -qp0, -qpd0
        B[0, j] = qpDt + A[0, 0] * qh0 + A[0, 1] * qhd0
        B[1, j] = qpdDt + A[1, 0] * qh0 + A[1, 1] * qhd0
    return A, B


def sdof_pwl(omega, xi, dt, ag):
    """Exact PWL response of u'' + 2 xi w u' + w^2 u = -ag(t), from rest."""
    A, B = _blocks(omega, 2.0 * xi * omega, dt)
    n = len(ag)
    u = np.zeros(n); v = np.zeros(n)
    for k in range(n - 1):
        f0, f1 = -ag[k], -ag[k + 1]
        u[k + 1] = A[0, 0] * u[k] + A[0, 1] * v[k] + B[0, 0] * f0 + B[0, 1] * f1
        v[k + 1] = A[1, 0] * u[k] + A[1, 1] * v[k] + B[1, 0] * f0 + B[1, 1] * f1
    return u


def newmark_full(M, C, K, R, ag, dt, dof):
    """Independent full-matrix Newmark (avg-accel 0.5,0.25), relative form,
    returns the history of response `dof` incl the initial rest state."""
    ndof = M.shape[0]
    g, b = 0.5, 0.25
    u = np.zeros(ndof); v = np.zeros(ndof)
    a = np.linalg.solve(M, -M @ R * ag[0] - C @ v - K @ u)
    Keff = K + (g / (b * dt)) * C + (1.0 / (b * dt * dt)) * M
    out = [u[dof]]
    for n in range(len(ag) - 1):
        p1 = -M @ R * ag[n + 1]
        rhs = (p1
               + M @ ((1 / (b * dt * dt)) * u + (1 / (b * dt)) * v + (1 / (2 * b) - 1) * a)
               + C @ ((g / (b * dt)) * u + (g / b - 1) * v + dt * (g / (2 * b) - 1) * a))
        un = np.linalg.solve(Keff, rhs)
        an = (1 / (b * dt * dt)) * (un - u) - (1 / (b * dt)) * v - (1 / (2 * b) - 1) * a
        vn = v + dt * ((1 - g) * a + g * an)
        u, v, a = un, vn, an
        out.append(u[dof])
    return np.array(out)


# ---------------------------------------------------------------------------
# model builders (1-D axial chain: shear-building analogue, ndf=1)
# ---------------------------------------------------------------------------
def _ground_motion(n, dt, seed=3):
    rng = np.random.default_rng(seed)
    t = np.arange(n + 1) * dt
    ag = (0.8 * np.sin(2 * math.pi * 1.3 * t)
          + 0.4 * np.sin(2 * math.pi * 3.7 * t + 0.5)
          + 0.15 * np.cumsum(rng.standard_normal(n + 1)) * math.sqrt(dt))
    ag[0] = 0.0
    return ag


def _timeseries(tag, dt, ag):
    # pad 2 trailing samples so the last analysis station (t=N*dt) samples the
    # record's interior, not its exact end (getFactor at the end returns 0 due
    # to n*dt floating-point round-off past the final abscissa).
    vals = ag.tolist() + [float(ag[-1])] * 2
    ops.timeSeries('Path', tag, '-dt', dt, '-values', *vals)


def _build_sdof(m, k):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.fix(1, 1)
    ops.node(2, 1.0); ops.mass(2, m)
    ops.uniaxialMaterial('Elastic', 1, k)   # E=k
    ops.element('Truss', 1, 1, 2, 1.0, 1)    # A=1, L=1 -> stiffness = E*A/L = k
    return math.sqrt(k / m)


def _build_chain(masses, ks):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.fix(1, 1)
    for i, mm in enumerate(masses, start=2):
        ops.node(i, float(i - 1)); ops.mass(i, mm)
    for e, kk in enumerate(ks, start=1):
        ops.uniaxialMaterial('Elastic', e, kk)
        ops.element('Truss', e, e, e + 1, 1.0, e)


def _chain_MK(masses, ks):
    n = len(masses)
    M = np.diag(masses)
    K = np.zeros((n, n))
    for i, kk in enumerate(ks):
        # dof i-1 (below) and i (above); base dof i=0 is fixed (node 2 -> dof 0)
        lo = i - 1
        if lo >= 0:
            K[lo, lo] += kk
        K[i, i] += kk
        if lo >= 0:
            K[lo, i] -= kk; K[i, lo] -= kk
    return M, K


def _record_mrh(node, ndof, dt, nsteps, damp_args, base_dir=1, ts_tag=1,
                modes=None, path=None):
    ops.recorder('Node', '-file', str(path), '-time', '-precision', 16,
                 '-node', node, '-dof', ndof, 'disp')
    args = ['-dt', dt, '-nsteps', nsteps, *damp_args,
            '-baseAccel', ts_tag, '-dir', base_dir]
    if modes is not None:
        args += ['-modes', *modes]
    ops.modalResponseHistory(*args)
    ops.remove('recorders')
    data = np.loadtxt(str(path))
    return data[:, 0], data[:, 1]


# ---------------------------------------------------------------------------
# SDOF exactness across damping branches
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("xi,tag", [(0.05, "under"), (1.0, "critical"), (1.6, "over")])
def test_sdof_branches(tmp_path, xi, tag):
    m, k = 2.0, 500.0            # w = sqrt(250) ~ 15.8 rad/s
    dt, N = 0.004, 500
    ag = _ground_motion(N, dt)

    w = _build_sdof(m, k)
    _timeseries(1, dt, ag)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _, u = _record_mrh(2, 1, dt, N, ['-damp', xi], path=tmp_path / f"sdof_{tag}.out")

    u_ref = sdof_pwl(w, xi, dt, ag)
    err = np.max(np.abs(u - u_ref))
    peak = np.max(np.abs(u_ref))
    assert peak > 1e-6, "degenerate oracle"
    assert err < 1e-9, f"{tag}: max|u-u_ref|={err:.2e} (peak {peak:.3e})"


# ---------------------------------------------------------------------------
# end-to-end vs OpenSees's own Newmark (mass-proportional damping)
# ---------------------------------------------------------------------------
def test_mdof_vs_opensees_newmark(tmp_path):
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    dt, N = 0.0015, 700
    a0 = 0.30                      # mass-proportional only (see the module docstring)
    ag = _ground_motion(N, dt, seed=7)

    # direct OpenSees Newmark
    _build_chain(masses, ks)
    _timeseries(1, dt, ag)
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
    ops.rayleigh(a0, 0.0, 0.0, 0.0)
    nk_path = tmp_path / "newmark.out"
    ops.recorder('Node', '-file', str(nk_path), '-time', '-precision', 16,
                 '-node', 3, '-dof', 1, 'disp')
    ops.constraints('Plain'); ops.numberer('Plain'); ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-12, 20); ops.algorithm('Linear')
    ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')
    for _ in range(N):
        ops.analyze(1, dt)
    ops.remove('recorders')
    nk = np.loadtxt(str(nk_path))

    # modal transient (all modes)
    _build_chain(masses, ks)
    _timeseries(1, dt, ag)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    mt, mu = _record_mrh(3, 1, dt, N, ['-rayleigh', a0, 0.0],
                         path=tmp_path / "modal.out")

    # align by time (modal has the extra t=0 station)
    ui = np.interp(nk[:, 0], mt, mu)
    peak = np.max(np.abs(nk[:, 1]))
    rel = np.max(np.abs(ui - nk[:, 1])) / peak
    assert peak > 1e-6
    assert rel < 2e-3, f"modal vs OpenSees-Newmark rel err {rel:.2e}"


# ---------------------------------------------------------------------------
# damped multi-mode case: modal vs an independent full-matrix Newmark
# ---------------------------------------------------------------------------
def test_mdof_stiffness_damping_vs_full_matrix(tmp_path):
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    dt, N = 0.0015, 700
    a0, a1 = 0.30, 0.002
    ag = _ground_motion(N, dt, seed=7)

    _build_chain(masses, ks)
    _timeseries(1, dt, ag)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    mt, mu = _record_mrh(3, 1, dt, N, ['-rayleigh', a0, a1], path=tmp_path / "modal.out")

    M, K = _chain_MK(masses, ks)
    C = a0 * M + a1 * K
    ref = newmark_full(M, C, K, np.ones(len(masses)), ag, dt, dof=1)  # node3 = dof 1

    # modal station j (t=j*dt) vs full-matrix step j
    n = min(len(mu), len(ref))
    peak = np.max(np.abs(ref[:n]))
    rel = np.max(np.abs(mu[:n] - ref[:n])) / peak
    assert peak > 1e-6
    assert rel < 1e-3, f"modal vs full-matrix Newmark rel err {rel:.2e}"


# ---------------------------------------------------------------------------
# plumbing / consistency
# ---------------------------------------------------------------------------
def test_mode_subset_plumbing(tmp_path):
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    dt, N = 0.004, 400
    ag = _ground_motion(N, dt, seed=11)

    def run(modes):
        _build_chain(masses, ks)
        _timeseries(1, dt, ag)
        ops.eigen('-fullGenLapack', 2); ops.modalProperties()
        _, u = _record_mrh(3, 1, dt, N, ['-damp', 0.03], modes=modes,
                           path=tmp_path / f"sub_{'_'.join(map(str, modes)) if modes else 'all'}.out")
        return u

    all_modes = run(None)
    both = run([1, 2])
    one = run([1])
    assert np.max(np.abs(both - all_modes)) < 1e-10, "explicit {1,2} != default(all)"
    assert np.max(np.abs(one - all_modes)) > 1e-6, "single mode indistinguishable from all"


def test_rayleigh_matches_modal_ratio(tmp_path):
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    dt, N = 0.004, 400
    a0, a1 = 0.25, 0.0015
    ag = _ground_motion(N, dt, seed=5)

    _build_chain(masses, ks)
    _timeseries(1, dt, ag)
    lam = ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    xis = [a0 / (2.0 * math.sqrt(l)) + a1 * math.sqrt(l) / 2.0 for l in lam]
    _, u_ray = _record_mrh(3, 1, dt, N, ['-rayleigh', a0, a1], path=tmp_path / "ray.out")

    _build_chain(masses, ks)
    _timeseries(1, dt, ag)
    ops.eigen('-fullGenLapack', 2); ops.modalProperties()
    _, u_mod = _record_mrh(3, 1, dt, N, ['-modalDamp', *xis], path=tmp_path / "mod.out")

    assert np.max(np.abs(u_ray - u_mod)) < 1e-9


# ---------------------------------------------------------------------------
# refuse-loudly guards
# ---------------------------------------------------------------------------
def test_guard_no_modalproperties():
    _build_sdof(2.0, 500.0)
    _timeseries(1, 0.01, np.array([0.0, 1.0, 0.0]))
    ops.eigen('-fullGenLapack', 1)
    with pytest.raises(Exception):
        ops.modalResponseHistory('-dt', 0.01, '-nsteps', 2, '-damp', 0.05,
                                 '-baseAccel', 1, '-dir', 1)


def test_guard_bad_dir():
    _build_sdof(2.0, 500.0)
    _timeseries(1, 0.01, np.array([0.0, 1.0, 0.0]))
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.modalResponseHistory('-dt', 0.01, '-nsteps', 2, '-damp', 0.05,
                                 '-baseAccel', 1, '-dir', 5)   # ndf==1


def test_guard_no_baseaccel():
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.modalResponseHistory('-dt', 0.01, '-nsteps', 2, '-damp', 0.05,
                                 '-dir', 1)   # missing -baseAccel


def test_guard_duplicate_modes():
    # a repeated mode index would silently double-count that mode -> must refuse
    _build_chain([2.0, 1.5], [800.0, 500.0])
    _timeseries(1, 0.01, np.array([0.0, 1.0, 0.0]))
    ops.eigen('-fullGenLapack', 2); ops.modalProperties()
    with pytest.raises(Exception):
        ops.modalResponseHistory('-dt', 0.01, '-nsteps', 2, '-damp', 0.05,
                                 '-baseAccel', 1, '-dir', 1, '-modes', 1, 1)


# ===========================================================================
# -load channel: nodal-force transient P(t) = s(t)*P
#
#   f_a(t) = (psi_a^T P / m~_a) s(t) through the same exact PWL recurrence.
#   Response is ABSOLUTE.
# ===========================================================================
def _add_force_pattern(tag, ts_tag, loads):
    ops.pattern('Plain', tag, ts_tag)
    for n, v in loads.items():
        ops.load(n, v)


def _record_mrh_load(node, dt, nsteps, damp_args, pattern_tag, series_tag,
                     path):
    ops.recorder('Node', '-file', str(path), '-time', '-precision', 16,
                 '-node', node, '-dof', 1, 'disp')
    ops.modalResponseHistory('-dt', dt, '-nsteps', nsteps, *damp_args,
                             '-load', pattern_tag, '-series', series_tag)
    ops.remove('recorders')
    data = np.loadtxt(str(path))
    return data[:, 0], data[:, 1]


def test_load_sdof_force_exact(tmp_path):
    # SDOF under F(t) = P*s(t): m u'' + c u' + k u = P s  <=>  the sdof_pwl
    # oracle with ag = -(P/m) s. Pins fcoef = psi^T P/m~ end-to-end at 1e-9.
    m, k, xi, P = 2.0, 500.0, 0.05, 3.0
    dt, N = 0.004, 500
    s = _ground_motion(N, dt, seed=11)

    w = _build_sdof(m, k)
    _timeseries(1, dt, s)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _add_force_pattern(2, 1, {2: P})
    _, u = _record_mrh_load(2, dt, N, ['-damp', xi], 2, 1,
                            tmp_path / "lsdof.out")

    u_ref = sdof_pwl(w, xi, dt, -(P / m) * s)
    err = np.max(np.abs(u - u_ref))
    peak = np.max(np.abs(u_ref))
    assert peak > 1e-6
    assert err < 1e-9, f"-load SDOF max|u-u_ref|={err:.2e} (peak {peak:.3e})"


def test_load_mdof_vs_opensees_newmark(tmp_path):
    # multi-node force history vs a direct
    # OpenSees Newmark run driving the SAME pattern (independent of the modal
    # machinery end-to-end).
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    dt, N = 0.0015, 700
    a0 = 0.30                      # mass-proportional only (see the module docstring)
    s = _ground_motion(N, dt, seed=13)
    loads = {2: 4.0, 3: -1.5}

    # direct Newmark with the pattern active
    _build_chain(masses, ks)
    _timeseries(1, dt, s)
    _add_force_pattern(2, 1, loads)
    ops.rayleigh(a0, 0.0, 0.0, 0.0)
    nk_path = tmp_path / "lnewmark.out"
    ops.recorder('Node', '-file', str(nk_path), '-time', '-precision', 16,
                 '-node', 3, '-dof', 1, 'disp')
    ops.constraints('Plain'); ops.numberer('Plain'); ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-12, 20); ops.algorithm('Linear')
    ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')
    for _ in range(N):
        ops.analyze(1, dt)
    ops.remove('recorders')
    nk = np.loadtxt(str(nk_path))

    # modal -load transient (all modes; same pattern shape + series)
    _build_chain(masses, ks)
    _timeseries(1, dt, s)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    _add_force_pattern(2, 1, loads)
    mt, mu = _record_mrh_load(3, dt, N, ['-rayleigh', a0, 0.0], 2, 1,
                              tmp_path / "lmodal.out")

    ui = np.interp(nk[:, 0], mt, mu)
    peak = np.max(np.abs(nk[:, 1]))
    rel = np.max(np.abs(ui - nk[:, 1])) / peak
    assert peak > 1e-6
    assert rel < 2e-3, f"-load modal vs OpenSees-Newmark rel err {rel:.2e}"


def test_load_equals_baseaccel_channel(tmp_path):
    # -baseAccel -dir 1 == -load with P = -M R and the SAME series: the modal
    # loads are identical (Gamma = psi^T M R/m~), so the histories must match
    # to machine precision.
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    dt, N = 0.004, 300
    ag = _ground_motion(N, dt, seed=17)

    _build_chain(masses, ks)
    _timeseries(1, dt, ag)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    _, ub = _record_mrh(3, 1, dt, N, ['-damp', 0.04],
                        path=tmp_path / "chan_base.out")

    _build_chain(masses, ks)
    _timeseries(1, dt, ag)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    _add_force_pattern(2, 1, {2: -2.0, 3: -1.5})   # P = -M R
    _, ul = _record_mrh_load(3, dt, N, ['-damp', 0.04], 2, 1,
                             tmp_path / "chan_load.out")
    assert np.allclose(ub, ul, rtol=1e-10, atol=1e-14)


def test_load_unorm_invariance_transient(tmp_path):
    # default vs -unorm modalProperties: psi^T P/m~ is normalization-invariant,
    # so the committed history must be identical (the sole pin on the second
    # Vscale factor; default normalization has Vscale=1 and cannot see it).
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    dt, N = 0.004, 300
    s = _ground_motion(N, dt, seed=19)
    loads = {2: 1.0, 3: 3.0}

    _build_chain(masses, ks)
    _timeseries(1, dt, s)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    _add_force_pattern(2, 1, loads)
    _, ud = _record_mrh_load(3, dt, N, ['-damp', 0.03], 2, 1,
                             tmp_path / "un_def.out")

    _build_chain(masses, ks)
    _timeseries(1, dt, s)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties('-unorm')
    _add_force_pattern(2, 1, loads)
    _, uu = _record_mrh_load(3, dt, N, ['-damp', 0.03], 2, 1,
                             tmp_path / "un_un.out")
    assert np.allclose(ud, uu, rtol=1e-10, atol=1e-14)


def test_guard_load_without_series():
    _build_sdof(2.0, 500.0)
    _timeseries(1, 0.01, np.array([0.0, 1.0, 0.0]))
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    _add_force_pattern(2, 1, {2: 1.0})
    with pytest.raises(Exception):
        ops.modalResponseHistory('-dt', 0.01, '-nsteps', 2, '-damp', 0.05,
                                 '-load', 2)


def test_guard_load_and_baseaccel_exclusive_transient():
    _build_sdof(2.0, 500.0)
    _timeseries(1, 0.01, np.array([0.0, 1.0, 0.0]))
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    _add_force_pattern(2, 1, {2: 1.0})
    with pytest.raises(Exception):
        ops.modalResponseHistory('-dt', 0.01, '-nsteps', 2, '-damp', 0.05,
                                 '-baseAccel', 1, '-dir', 1,
                                 '-load', 2, '-series', 1)


def test_guard_series_without_load():
    _build_sdof(2.0, 500.0)
    _timeseries(1, 0.01, np.array([0.0, 1.0, 0.0]))
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.modalResponseHistory('-dt', 0.01, '-nsteps', 2, '-damp', 0.05,
                                 '-baseAccel', 1, '-dir', 1, '-series', 1)
