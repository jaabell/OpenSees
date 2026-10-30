"""randomResponse: stationary random response, PSD -> RMS.

For a base-acceleration excitation given as a one-sided PSD G(f) in Hz
((accel)^2/Hz; sigma_ug^2 = int_0^inf G df, the convention of engineering
specifications), the response auto-PSD is

    G_xx(f) = |H_x(f)|^2 G(f),

with H_x the frequencyResponse transfer function, and the band moments and
statistics on the sweep grid (trapezoid rule):

    m_k = int f^k G_xx df,  RMS = sqrt(m0),  nu0 = sqrt(m2/m0) [Hz],
    E[peak] = (sqrt(2 ln(nu0 T)) + 0.5772/sqrt(2 ln(nu0 T))) * RMS.

Against the two-sided rad/s PSD of random-vibration texts, G(f) = 4 pi S(Om);
the white-noise SDOF variance in this convention is sigma_x^2 = G0/(8 xi w^3).

  test_white_noise_sdof_analytic
        flat G0 through an SDOF: RMS against sqrt(G0/(8 xi w^3)).
  test_mdof_psd_vs_direct[disp/vel/accel]
        the {f, Gin, Gxx} table against an independent numpy direct complex
        solve |(K - Om^2 M + i Om C)^-1 (-M R)|^2 G(f), row for row, plus the
        RMS integral (does not use the modal formula); G_xx >= 0.
  test_monte_carlo_time_domain
        RMS against a synthetic time-history realization of the target PSD
        marched through the exact SDOF recurrence in numpy (fixed seed, 8 %).
        A one-sided/two-sided factor error reads ~41 %, a Hz/rad error ~150 %.
  test_stats_narrowband
        nu0 = sqrt(m2/m0) at fn for a lightly damped mode; Davenport peak.
  -load channel, damping-channel equivalence, mode subset, -out/return
        consistency and the refusals (rigid mode with f = 0 in band, zero-damped
        mode resonant in band, single-point grid, negative PSD).

Small models use `eigen -fullGenLapack` (ARPACK needs NEV < N). Every
non-guard test fails on a build without randomResponse (unknown command).
"""
import math

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

TWO_PI = 2.0 * math.pi


# ---------------------------------------------------------------------------
# numpy references (independent of the modal formula where it matters)
# ---------------------------------------------------------------------------
def frf_direct(M, K, C, R, freqs, dof, resp="disp"):
    """Direct complex solve (K - Om^2 M + i Om C)^{-1}(-M R) per frequency."""
    out = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        Om = TWO_PI * f
        u = np.linalg.solve(K - Om * Om * M + 1j * Om * C, -(M @ R))[dof]
        if resp == "vel":
            u *= 1j * Om
        elif resp == "accel":
            u *= -(Om * Om)
        out[i] = u
    return out


def trapz_moments(f, g):
    m0 = np.trapezoid(g, f)
    m2 = np.trapezoid(f * f * g, f)
    return m0, m2


def _build_sdof(m, k):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.fix(1, 1)
    ops.node(2, 1.0); ops.mass(2, m)
    ops.uniaxialMaterial('Elastic', 1, k)
    ops.element('Truss', 1, 1, 2, 1.0, 1)


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
    M = np.diag(np.asarray(masses, float))
    K = np.zeros((n, n))
    for i, kk in enumerate(ks):
        lo = i - 1
        if lo >= 0:
            K[lo, lo] += kk
            K[lo, i] -= kk; K[i, lo] -= kk
        K[i, i] += kk
    return M, K


def _flat_psd(tag, G0):
    """One-sided PSD G(f) = G0 for all f (Constant series sampled at f)."""
    ops.timeSeries('Constant', tag, '-factor', G0)


def _band_psd(tag, f1, f2, G0):
    """Band-limited flat PSD: G0 on [f1, f2], 0 outside (Path sampled at f)."""
    ops.timeSeries('Path', tag, '-time', 0.0, 0.999 * f1, f1, f2,
                   '-values', 0.0, 0.0, G0, G0)


# ---------------------------------------------------------------------------
# white-noise SDOF vs the analytic variance
# ---------------------------------------------------------------------------
def test_white_noise_sdof_analytic(tmp_path):
    m, k, xi = 2.0, 800.0, 0.05
    w = math.sqrt(k / m)
    fn = w / TWO_PI
    G0 = 0.13

    _build_sdof(m, k)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, G0)
    rms = ops.randomResponse('-freq', 1e-3, 30.0 * fn, 2500, '-biased',
                             '-baseAccel', '-dir', 1, '-inputPSD', 1,
                             '-damp', xi, '-node', 2, '-dof', 1)
    ref = math.sqrt(G0 / (8.0 * xi * w ** 3))
    # exact analytic is the infinite band; [0, 30 fn] truncates ~nothing
    assert abs(rms - ref) / ref < 2e-3, f"RMS {rms:.6e} vs analytic {ref:.6e}"


# ---------------------------------------------------------------------------
# PSD table + RMS vs an independent direct complex solve
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("resp", ["disp", "vel", "accel"])
def test_mdof_psd_vs_direct(tmp_path, resp):
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    a0, a1 = 0.30, 2.0e-3
    f1, f2, G0 = 0.2, 6.0, 0.05

    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    _band_psd(1, f1, f2, G0)
    path = tmp_path / f"psd_{resp}.out"
    ret = ops.randomResponse('-freq', f1, f2, 400, '-biased',
                             '-baseAccel', '-dir', 1, '-inputPSD', 1,
                             '-rayleigh', a0, a1, '-node', 4, '-dof', 1,
                             '-resp', resp, '-stats', '-out', str(path))
    rows = np.loadtxt(str(path))
    f, gin, gxx = rows[:, 0], rows[:, 1], rows[:, 2]

    # a PSD is non-negative by construction
    assert np.all(gxx >= 0.0) and np.all(gin >= 0.0)

    # the sampled input PSD is the Path series (linear interp, 0 outside)
    gin_ref = np.interp(f, [0.0, 0.999 * f1, f1, f2], [0.0, 0.0, G0, G0],
                        left=0.0, right=0.0)
    assert np.allclose(gin, gin_ref, rtol=1e-12, atol=1e-15)

    # row-for-row: Gxx == |H_direct|^2 * Gin  (independent of the modal sum)
    M, K = _chain_MK(masses, ks)
    C = a0 * M + a1 * K
    H = frf_direct(M, K, C, np.ones(3), f, dof=2, resp=resp)
    gxx_ref = np.abs(H) ** 2 * gin
    scale = gxx_ref.max()
    assert scale > 0
    assert np.max(np.abs(gxx - gxx_ref)) / scale < 1e-8

    # and the integral: returned stats {rms, nu0, m0, m2} match the same trapz
    rms, nu0, m0, m2 = ret
    m0_ref, m2_ref = trapz_moments(f, gxx_ref)
    assert abs(m0 - m0_ref) / m0_ref < 1e-8
    assert abs(m2 - m2_ref) / m2_ref < 1e-8
    assert abs(rms - math.sqrt(m0_ref)) / math.sqrt(m0_ref) < 1e-8
    assert abs(nu0 - math.sqrt(m2_ref / m0_ref)) / nu0 < 1e-8


# ---------------------------------------------------------------------------
# RMS vs a Monte-Carlo synthetic time history (fixed seed)
# ---------------------------------------------------------------------------
def test_monte_carlo_time_domain():
    rng = np.random.default_rng(12345)
    m, k, xi = 2.0, 800.0, 0.05
    w = math.sqrt(k / m); fn = w / TWO_PI; d = 2.0 * xi * w
    G0, f1, f2 = 0.02, 0.2, 12.0

    # ---- numpy leg: synthesize ug''(t) with the target ONE-SIDED-Hz PSD ----
    df = 1.0 / 128.0                     # full period T = 128 s
    fk = np.arange(f1, f2, df) + 0.5 * df
    amp = math.sqrt(2.0 * G0 * df)
    phi = rng.uniform(0.0, TWO_PI, len(fk))
    dt = 0.01
    t = np.arange(0.0, 128.0, dt)
    ug = (amp * np.cos(TWO_PI * np.outer(t, fk) + phi)).sum(axis=1)

    # exact PWL SDOF recurrence (underdamped closed form, pure numpy)
    wd = w * math.sqrt(1.0 - xi * xi)
    e = math.exp(-0.5 * d * dt)
    c, s = math.cos(wd * dt), math.sin(wd * dt) / wd
    A = np.array([[e * (c + 0.5 * d * s), e * s],
                  [e * (-w * w * s), e * (c - 0.5 * d * s)]])
    B = np.empty((2, 2))
    for j in range(2):                    # ramp-response load block
        f0, fq = (1.0, 0.0) if j == 0 else (0.0, 1.0)
        r = (fq - f0) / dt
        qp0 = f0 / w ** 2 - d * r / w ** 4
        qpd0 = r / w ** 2
        qpDt = (f0 + r * dt) / w ** 2 - d * r / w ** 4
        h = A @ np.array([-qp0, -qpd0])
        B[0, j] = qpDt + h[0]
        B[1, j] = qpd0 + h[1]
    x = np.zeros(2)
    q = np.empty(len(t))
    fl = -ug                              # Gamma = 1 for the SDOF
    for i in range(len(t)):
        q[i] = x[0]
        if i + 1 < len(t):
            x = A @ x + B @ np.array([fl[i], fl[i + 1]])
    n0 = int(10.0 / (xi * fn) / dt)       # discard ~10 damped periods
    sig_mc = q[n0:].std()

    # ---- OpenSees leg: randomResponse on the same band-limited PSD ----
    _build_sdof(m, k)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _band_psd(1, f1, f2, G0)
    rms = ops.randomResponse('-freq', f1, f2, 1500, '-biased',
                             '-baseAccel', '-dir', 1, '-inputPSD', 1,
                             '-damp', xi, '-node', 2, '-dof', 1)

    err = abs(sig_mc - rms) / rms
    # one realization of a narrow-band process: a few % sampling error; the
    # factor bugs this pins read 41% (one/two-sided) and 150% (Hz/rad).
    assert err < 0.08, f"MC {sig_mc:.5e} vs spectral {rms:.5e} ({err * 100:.1f}%)"


# ---------------------------------------------------------------------------
# stats: narrow-band nu0 ~ fn; Davenport peak; -duration plumbing
# ---------------------------------------------------------------------------
def test_stats_narrowband():
    m, k, xi = 2.0, 800.0, 0.02
    w = math.sqrt(k / m); fn = w / TWO_PI
    _build_sdof(m, k)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, 0.04)
    T = 600.0
    ret = ops.randomResponse('-freq', 1e-3, 25.0 * fn, 2000, '-biased',
                             '-baseAccel', '-dir', 1, '-inputPSD', 1,
                             '-damp', xi, '-node', 2, '-dof', 1,
                             '-stats', '-duration', T)
    assert len(ret) == 5, "stats with -duration = {rms, nu0, m0, m2, peak}"
    rms, nu0, m0, m2, peak = ret

    assert abs(rms - math.sqrt(m0)) / rms < 1e-12
    assert abs(nu0 - math.sqrt(m2 / m0)) / nu0 < 1e-12
    # narrow-band process: zero-upcrossing rate sits at the natural frequency
    assert abs(nu0 - fn) / fn < 0.02
    # Davenport peak factor recomputed from the returned nu0
    lnt = math.sqrt(2.0 * math.log(nu0 * T))
    assert abs(peak - (lnt + 0.5772 / lnt) * rms) / peak < 1e-12

    # without -duration: 4 stats entries
    ret4 = ops.randomResponse('-freq', 1e-3, 25.0 * fn, 2000, '-biased',
                              '-baseAccel', '-dir', 1, '-inputPSD', 1,
                              '-damp', xi, '-node', 2, '-dof', 1, '-stats')
    assert len(ret4) == 4
    assert abs(ret4[0] - rms) / rms < 1e-12


# ---------------------------------------------------------------------------
# velocity RMS: sigma_v == w * sigma_x for white noise through an SDOF
# ---------------------------------------------------------------------------
def test_vel_rms_white_noise():
    m, k, xi = 2.0, 800.0, 0.05
    w = math.sqrt(k / m); fn = w / TWO_PI
    _build_sdof(m, k)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, 0.13)
    args = ('-freq', 1e-3, 40.0 * fn, 3000, '-biased', '-baseAccel', '-dir', 1,
            '-inputPSD', 1, '-damp', xi, '-node', 2, '-dof', 1)
    rms_d = ops.randomResponse(*args)
    rms_v = ops.randomResponse(*args, '-resp', 'vel')
    assert abs(rms_v - w * rms_d) / (w * rms_d) < 5e-3


# ---------------------------------------------------------------------------
# damping-channel equivalence: rayleigh == the equivalent per-mode list
# ---------------------------------------------------------------------------
def test_damping_channel_equiv():
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    a0, a1 = 0.25, 1.5e-3
    _build_chain(masses, ks)
    lam = ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    xis = [a0 / (2.0 * math.sqrt(l)) + a1 * math.sqrt(l) / 2.0 for l in lam]
    _band_psd(1, 0.2, 5.0, 0.03)
    args = ('-freq', 0.2, 5.0, 300, '-biased', '-baseAccel', '-dir', 1,
            '-inputPSD', 1, '-node', 3, '-dof', 1)
    r_ray = ops.randomResponse(*args, '-rayleigh', a0, a1)
    r_mod = ops.randomResponse(*args, '-modalDamp', *xis)
    assert abs(r_ray - r_mod) / r_ray < 1e-10


# ---------------------------------------------------------------------------
# mode subset: single retained mode vs the numpy single-mode reference
# ---------------------------------------------------------------------------
def test_mode_subset(tmp_path):
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    xi = 0.04
    f1, f2, G0 = 0.2, 6.0, 0.05
    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    _band_psd(1, f1, f2, G0)
    path = tmp_path / "m1.out"
    rms = ops.randomResponse('-freq', f1, f2, 400, '-biased',
                             '-baseAccel', '-dir', 1, '-inputPSD', 1,
                             '-damp', xi, '-node', 4, '-dof', 1,
                             '-modes', 1, '-out', str(path))
    rows = np.loadtxt(str(path))
    f, gin = rows[:, 0], rows[:, 1]

    # numpy single-mode modal reference (mass-normalized phi)
    M, K = _chain_MK(masses, ks)
    Lam, Q = np.linalg.eigh(M)
    Mis = Q @ np.diag(1.0 / np.sqrt(Lam)) @ Q.T
    w2, U = np.linalg.eigh(Mis @ K @ Mis)
    Phi = Mis @ U
    w1 = math.sqrt(w2[0])
    Gam = Phi[:, 0] @ (M @ np.ones(3))
    Om = TWO_PI * f
    H1 = Phi[2, 0] * (-Gam) / (w1 * w1 - Om * Om + 1j * Om * 2.0 * xi * w1)
    m0_ref = np.trapezoid(np.abs(H1) ** 2 * gin, f)
    assert abs(rms - math.sqrt(m0_ref)) / rms < 1e-8


# ---------------------------------------------------------------------------
# -out file == returned values (trapz consistency), scalar return form
# ---------------------------------------------------------------------------
def test_out_file_matches_scalar_return(tmp_path):
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _band_psd(1, 0.5, 8.0, 0.02)
    path = tmp_path / "tab.out"
    rms = ops.randomResponse('-freq', 0.5, 8.0, 250, '-lin',
                             '-baseAccel', '-dir', 1, '-inputPSD', 1,
                             '-damp', 0.03, '-node', 2, '-dof', 1,
                             '-out', str(path))
    assert isinstance(rms, float)
    rows = np.loadtxt(str(path))
    m0, _ = trapz_moments(rows[:, 0], rows[:, 2])
    assert abs(rms - math.sqrt(m0)) / rms < 1e-12


# ---------------------------------------------------------------------------
# zero-damped mode OUTSIDE the band is legal (finite in-band integrand)
# ---------------------------------------------------------------------------
def test_zero_damping_out_of_band_is_legal(tmp_path):
    m, k = 2.0, 800.0
    fn = math.sqrt(k / m) / TWO_PI
    _build_sdof(m, k)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, 0.02)
    # band strictly above the (undamped) resonance
    path = tmp_path / "oob.out"
    rms = ops.randomResponse('-freq', 3.0 * fn, 6.0 * fn, 200, '-lin',
                             '-baseAccel', '-dir', 1, '-inputPSD', 1,
                             '-damp', 0.0, '-node', 2, '-dof', 1,
                             '-out', str(path))
    rows = np.loadtxt(str(path))
    w = math.sqrt(k / m)
    Om = TWO_PI * rows[:, 0]
    gxx_ref = np.abs(-1.0 / (w * w - Om * Om)) ** 2 * rows[:, 1]
    assert np.max(np.abs(rows[:, 2] - gxx_ref)) / gxx_ref.max() < 1e-8
    m0, _ = trapz_moments(rows[:, 0], gxx_ref)
    assert abs(rms - math.sqrt(m0)) / rms < 1e-8


# ---------------------------------------------------------------------------
# refuse-loudly guards
# ---------------------------------------------------------------------------
def test_guard_zero_damping_in_band():
    # the resonance inside the band with d==0: variance integral diverges
    m, k = 2.0, 800.0
    fn = math.sqrt(k / m) / TWO_PI
    _build_sdof(m, k)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, 0.02)
    with pytest.raises(Exception):
        ops.randomResponse('-freq', 0.1, 3.0 * fn, 100, '-lin',
                           '-baseAccel', '-dir', 1, '-inputPSD', 1,
                           '-damp', 0.0, '-node', 2, '-dof', 1)


def test_guard_missing_input_psd():
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    with pytest.raises(Exception):
        ops.randomResponse('-freq', 0.1, 5.0, 100, '-baseAccel', '-dir', 1,
                           '-damp', 0.05, '-node', 2, '-dof', 1)


def test_guard_negative_psd():
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    ops.timeSeries('Path', 1, '-time', 0.0, 2.0, 5.0,
                   '-values', 0.1, -0.05, 0.1)     # dips negative in-band
    with pytest.raises(Exception):
        ops.randomResponse('-freq', 0.1, 5.0, 100, '-baseAccel', '-dir', 1,
                           '-inputPSD', 1, '-damp', 0.05, '-node', 2, '-dof', 1)


def test_guard_single_point_grid():
    # RMS is a band integral: nf==1 has zero measure -> refuse
    # (frequencyResponse allows it)
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, 0.02)
    with pytest.raises(Exception):
        ops.randomResponse('-freq', 2.0, 2.0, 1, '-baseAccel', '-dir', 1,
                           '-inputPSD', 1, '-damp', 0.05, '-node', 2, '-dof', 1)


def test_guard_negative_damp():
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, 0.02)
    with pytest.raises(Exception):
        ops.randomResponse('-freq', 0.1, 5.0, 100, '-baseAccel', '-dir', 1,
                           '-inputPSD', 1, '-damp', -0.05, '-node', 2, '-dof', 1)


def test_guard_bad_duration():
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, 0.02)
    with pytest.raises(Exception):
        ops.randomResponse('-freq', 0.1, 5.0, 100, '-baseAccel', '-dir', 1,
                           '-inputPSD', 1, '-damp', 0.05, '-node', 2, '-dof', 1,
                           '-duration', -10.0)


def test_guard_rigid_mode_with_fmin_zero():
    # a RIGID mode's FRF diverges at its own f=0 for ANY damping
    # (denom -Om^2 + iOm*d -> 0), so Rayleigh a0>0 must NOT rescue it; a d<=0
    # test alone would let it through to an inf/NaN RMS.
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.mass(1, 2.0)      # unrestrained: 1 rigid + 1 elastic
    ops.node(2, 1.0); ops.mass(2, 1.5)
    ops.uniaxialMaterial('Elastic', 1, 600.0)
    ops.element('Truss', 1, 1, 2, 1.0, 1)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    _flat_psd(1, 0.02)
    with pytest.raises(Exception):
        ops.randomResponse('-freq', 0.0, 5.0, 100, '-baseAccel', '-dir', 1,
                           '-inputPSD', 1, '-rayleigh', 0.3, 1e-3,
                           '-node', 2, '-dof', 1)
    # with fmin > 0 the rigid mode's in-band integrand is finite -> legal
    rms = ops.randomResponse('-freq', 0.5, 5.0, 200, '-baseAccel', '-dir', 1,
                             '-inputPSD', 1, '-rayleigh', 0.3, 1e-3,
                             '-node', 2, '-dof', 1)
    assert math.isfinite(rms) and rms > 0.0


def test_duration_return_shape_is_stable():
    # -duration ALWAYS returns 5 entries; a peak
    # that cannot be estimated (nu0*T <= 1) is NaN, never silently dropped.
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, 0.02)
    ret = ops.randomResponse('-freq', 0.5, 15.0, 400, '-biased',
                             '-baseAccel', '-dir', 1, '-inputPSD', 1,
                             '-damp', 0.03, '-node', 2, '-dof', 1,
                             '-duration', 1e-9)   # nu0*T << 1
    assert len(ret) == 5
    assert math.isnan(ret[4])
    assert math.isfinite(ret[0]) and ret[0] > 0.0


def test_guard_no_modalproperties():
    # NOTE: DomainModalProperties survives wipe(), so this model must NOT
    # reproduce the spectrum of any earlier test in this file: an identical
    # eigenvalue would make the stale
    # snapshot legitimately indistinguishable from a fresh one. k=512 is unique.
    _build_sdof(2.0, 512.0)
    ops.eigen('-fullGenLapack', 1)
    _flat_psd(1, 0.02)
    with pytest.raises(Exception):
        ops.randomResponse('-freq', 0.1, 5.0, 100, '-baseAccel', '-dir', 1,
                           '-inputPSD', 1, '-damp', 0.05, '-node', 2, '-dof', 1)


# ===========================================================================
# -load channel: stationary nodal-force excitation P(t) = s(t)*P, G_s(f) input
# ===========================================================================
def test_force_psd_white_noise_sdof():
    # white force PSD G_F through an SDOF: sigma_u = sqrt(G_F/(8 xi w^3 m^2))
    # (load_frf_oracle.py check 4).
    m, k, xi = 2.0, 800.0, 0.05
    w = math.sqrt(k / m); fn = w / TWO_PI
    GF = 0.4
    _build_sdof(m, k)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    ops.timeSeries('Constant', 1, '-factor', GF)
    ops.timeSeries('Constant', 111)
    ops.pattern('Plain', 11, 111)
    ops.load(2, 1.0)                     # unit force shape; PSD carries scale
    rms = ops.randomResponse('-freq', 1e-3, 30.0 * fn, 2500, '-biased',
                             '-load', 11, '-inputPSD', 1,
                             '-damp', xi, '-node', 2, '-dof', 1)
    ref = math.sqrt(GF / (8.0 * xi * w ** 3 * m * m))
    assert abs(rms - ref) / ref < 2e-3, f"RMS {rms:.6e} vs analytic {ref:.6e}"


def test_load_random_vs_direct(tmp_path):
    # {f, Gin, Gxx} rows vs the independent direct solve |(K-Om^2 M+iOm C)^-1 P|^2 G
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    a0, a1 = 0.30, 2.0e-3
    f1, f2, G0 = 0.2, 6.0, 0.05
    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    _band_psd(1, f1, f2, G0)
    ops.timeSeries('Constant', 112)
    ops.pattern('Plain', 12, 112)
    ops.load(3, 5.0); ops.load(4, -2.0)
    path = tmp_path / "lpsd.out"
    rms = ops.randomResponse('-freq', f1, f2, 400, '-biased',
                             '-load', 12, '-inputPSD', 1,
                             '-rayleigh', a0, a1, '-node', 4, '-dof', 1,
                             '-out', str(path))
    rows = np.loadtxt(str(path))
    f, gin, gxx = rows[:, 0], rows[:, 1], rows[:, 2]

    M, K = _chain_MK(masses, ks)
    C = a0 * M + a1 * K
    P = np.array([0.0, 5.0, -2.0])
    H = np.array([np.linalg.solve(K - (TWO_PI * ff) ** 2 * M
                                  + 1j * (TWO_PI * ff) * C, P)[2] for ff in f])
    gxx_ref = np.abs(H) ** 2 * gin
    assert np.max(np.abs(gxx - gxx_ref)) / gxx_ref.max() < 1e-8
    m0_ref = np.trapezoid(gxx_ref, f)
    assert abs(rms - math.sqrt(m0_ref)) / rms < 1e-8


def test_guard_random_load_and_baseaccel_exclusive():
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _flat_psd(1, 0.02)
    ops.timeSeries('Constant', 113)
    ops.pattern('Plain', 13, 113)
    ops.load(2, 1.0)
    with pytest.raises(Exception):
        ops.randomResponse('-freq', 0.5, 5.0, 100, '-baseAccel', '-dir', 1,
                           '-load', 13, '-inputPSD', 1, '-damp', 0.05,
                           '-node', 2, '-dof', 1)
