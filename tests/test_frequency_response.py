"""frequencyResponse / steadyStateDynamics: modal frequency response.

After `eigen` + `modalProperties`, the steady harmonic response of a linear
model to a uniform base acceleration u_g''(t) = amp e^{+iOm t} (relative
formulation) is the modal sum

    u_hat(Om) = amp * sum_a (phi_a Vscale_a)(-Gamma_a) H_a(Om),
    H_a(Om)   = 1/(w_a^2 - Om^2 + i Om d_a),   d_a = c_a/m_a.

`frequencyResponse` returns the complex transfer function {f, Re, Im};
`steadyStateDynamics` the amplitude {f, |.|}. Frequencies in Hz (Om = 2 pi f).

  test_mdof_frf_vs_direct[disp/vel/accel]
        Rayleigh-damped 3-DOF chain against an independent numpy direct
        complex solve (K - Om^2 M + i Om C)^{-1}(-M R), which does not use the
        modal formula, so sign, normalization and the H_a form are all
        checked (magnitude and phase, 1e-8).
  test_sdof_frf_closed_form
        1-DOF |H| and phase against the closed form, including the resonant
        peak 1/(2 xi w^2) and the +90 deg response phase at resonance.
  test_frf_peaks_at_modal_freqs
        -biased sweep: |FRF| local maxima within 1 % of the eigen frequencies.
  test_ssd_static_limit
        steadyStateDynamics at Om -> 0 equals the static response |K^-1 M R|.
  -load channel (nodal forces amp P e^{+iOm t}): against the direct solve
        (K - Om^2 M + i Om C)^{-1} P for translational and rotational loads,
        -load with P = -M R equal to -baseAccel, invariance under
        modalProperties -unorm, the static limit K^-1 P.
  test_ssd_equals_frf_magnitude, damping-channel equivalence, grid shapes
        and the refusals.

Small models use `eigen -fullGenLapack` (ARPACK needs NEV < N). Every
non-guard test fails on a build without these commands (unknown command).
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
# numpy references
# ---------------------------------------------------------------------------
def modal_eig(M, K):
    """Mass-normalized undamped modes: phi^T M phi = I, K phi = M phi diag(w^2)."""
    Lam, Q = np.linalg.eigh(M)
    Mis = Q @ np.diag(1.0 / np.sqrt(Lam)) @ Q.T
    w2, U = np.linalg.eigh(Mis @ K @ Mis)
    Phi = Mis @ U
    return np.sqrt(np.clip(w2, 0.0, None)), Phi


def frf_direct(M, K, C, R, freqs, dof, resp="disp", amp=1.0):
    """Physical relative-response FRF via the DIRECT complex solve."""
    out = []
    for f in freqs:
        Om = TWO_PI * f
        D = K - Om * Om * M + 1j * Om * C
        u = np.linalg.solve(D, -(M @ R))[dof]
        if resp == "vel":
            u *= 1j * Om
        elif resp == "accel":
            u *= -Om * Om
        out.append(amp * u)
    return np.array(out)


def frf_modal(M, K, dfun, R, freqs, dof, resp="disp", amp=1.0):
    """Modal-sum FRF with an arbitrary per-mode d_a=dfun(w_a) (covers uniform xi)."""
    w, Phi = modal_eig(M, K)
    out = []
    for f in freqs:
        Om = TWO_PI * f
        uc = 0.0 + 0.0j
        for a in range(Phi.shape[1]):
            phi = Phi[:, a]
            Gam = phi @ (M @ R)
            H = 1.0 / (w[a] * w[a] - Om * Om + 1j * Om * dfun(w[a]))
            uc += phi[dof] * (-Gam) * H
        if resp == "vel":
            uc *= 1j * Om
        elif resp == "accel":
            uc *= -Om * Om
        out.append(amp * uc)
    return np.array(out)


# ---------------------------------------------------------------------------
# model builders (1-D axial chain, ndf=1)
# ---------------------------------------------------------------------------
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
        K[i, i] += kk
        if lo >= 0:
            K[lo, i] -= kk; K[i, lo] -= kk
    return M, K


def _run(cmd, path, **_):
    """Call frequencyResponse/steadyStateDynamics, return the {f, ...} table as an (nf, ncol) array.

    Reads the -out file (source of truth) and cross-checks the returned list.
    """
    res = cmd()
    arr = np.loadtxt(str(path))
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if res:  # returned table should carry the same values as the -out file
        r = np.array(res, dtype=float).ravel()
        assert r.size == arr.size, (r.size, arr.size)
        assert np.allclose(r, arr.ravel(), rtol=1e-9, atol=0), "return != -out file"
    return arr


# ---------------------------------------------------------------------------
# FRF vs the direct complex solve
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("resp", ["disp", "vel", "accel"])
def test_mdof_frf_vs_direct(tmp_path, resp):
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    a0, a1 = 0.30, 2.0e-3
    fmin, fmax, nf = 0.05, 5.0, 220

    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    path = tmp_path / f"frf_{resp}.out"
    arr = _run(lambda: ops.frequencyResponse(
        '-freq', fmin, fmax, nf, '-lin', '-baseAccel', '-dir', 1,
        '-rayleigh', a0, a1, '-node', 4, '-dof', 1, '-resp', resp,
        '-out', str(path)), path)

    freqs = np.linspace(fmin, fmax, nf)
    assert np.allclose(arr[:, 0], freqs, rtol=0, atol=1e-9)

    M, K = _chain_MK(masses, ks)
    C = a0 * M + a1 * K
    ref = frf_direct(M, K, C, np.ones(3), freqs, dof=2, resp=resp)  # node4 = dof 2

    got = arr[:, 1] + 1j * arr[:, 2]
    scale = np.max(np.abs(ref))
    err = np.max(np.abs(got - ref)) / scale
    assert scale > 1e-8
    assert err < 1e-8, f"{resp}: FRF vs direct rel err {err:.2e}"


# ---------------------------------------------------------------------------
# SDOF closed form: resonant peak height + phase
# ---------------------------------------------------------------------------
def test_sdof_frf_closed_form(tmp_path):
    m, k, xi = 2.0, 800.0, 0.05
    w0 = math.sqrt(k / m)
    f0 = w0 / TWO_PI

    _build_sdof(m, k)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    # sample exactly at resonance and at DC
    path = tmp_path / "sdof.out"
    arr = _run(lambda: ops.frequencyResponse(
        '-freq', f0, f0, 1, '-baseAccel', '-dir', 1, '-damp', xi,
        '-node', 2, '-dof', 1, '-out', str(path)), path)
    u_res = arr[0, 1] + 1j * arr[0, 2]

    peak_analytic = 1.0 / (2.0 * xi * w0 * w0)        # |H(w0)| = 1/(2 xi w0^2)
    assert abs(abs(u_res) - peak_analytic) / peak_analytic < 1e-6
    # response phase at resonance is +90 deg (u = +i/(w0 d))
    assert abs(math.degrees(math.atan2(u_res.imag, u_res.real)) - 90.0) < 1e-4

    # static limit u(0) = -m/k
    p0 = tmp_path / "sdof0.out"
    a0 = _run(lambda: ops.frequencyResponse(
        '-freq', 1e-6, 1e-6, 1, '-baseAccel', '-dir', 1, '-damp', xi,
        '-node', 2, '-dof', 1, '-out', str(p0)), p0)
    assert abs((a0[0, 1] + 1j * a0[0, 2]).real - (-m / k)) < 1e-8


# ---------------------------------------------------------------------------
# FRF magnitude peaks at the modal frequencies (resonance-biased grid)
# ---------------------------------------------------------------------------
def test_frf_peaks_at_modal_freqs(tmp_path):
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    _build_chain(masses, ks)
    lam = ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    modal_hz = np.sort([math.sqrt(l) / TWO_PI for l in lam])

    path = tmp_path / "biased.out"
    arr = _run(lambda: ops.steadyStateDynamics(
        '-freq', 0.05, 6.0, 300, '-biased', '-baseAccel', '-dir', 1,
        '-damp', 0.02, '-node', 4, '-dof', 1, '-out', str(path)), path)
    f, mag = arr[:, 0], arr[:, 1]

    peaks = [f[i] for i in range(1, len(f) - 1)
             if mag[i] > mag[i - 1] and mag[i] > mag[i + 1]]
    assert len(peaks) >= 3, f"expected >=3 resonances, found {len(peaks)}"
    for mf in modal_hz:
        assert min(abs(np.array(peaks) - mf)) / mf < 0.01, f"no peak near {mf:.3f} Hz"


# ---------------------------------------------------------------------------
# SSD static limit
# ---------------------------------------------------------------------------
def test_ssd_static_limit(tmp_path):
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    path = tmp_path / "ssd0.out"
    arr = _run(lambda: ops.steadyStateDynamics(
        '-freq', 1e-6, 1e-6, 1, '-baseAccel', '-dir', 1, '-damp', 0.03,
        '-node', 4, '-dof', 1, '-out', str(path)), path)

    M, K = _chain_MK(masses, ks)
    ustat = np.linalg.solve(K, -(M @ np.ones(3)))
    assert abs(arr[0, 1] - abs(ustat[2])) / abs(ustat[2]) < 1e-6


# ---------------------------------------------------------------------------
# SSD amplitude == |frequencyResponse|
# ---------------------------------------------------------------------------
def test_ssd_equals_frf_magnitude(tmp_path):
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    common = dict(node=3)
    pf = tmp_path / "frf.out"; ps = tmp_path / "ssd.out"
    frf = _run(lambda: ops.frequencyResponse(
        '-freq', 0.1, 4.0, 60, '-baseAccel', '-dir', 1, '-damp', 0.04,
        '-node', 3, '-dof', 1, '-out', str(pf)), pf)
    ssd = _run(lambda: ops.steadyStateDynamics(
        '-freq', 0.1, 4.0, 60, '-baseAccel', '-dir', 1, '-damp', 0.04,
        '-node', 3, '-dof', 1, '-out', str(ps)), ps)
    mag = np.hypot(frf[:, 1], frf[:, 2])
    assert np.allclose(ssd[:, 1], mag, rtol=1e-12, atol=0)


# ---------------------------------------------------------------------------
# damping-channel equivalence: rayleigh == equivalent per-mode modalDamp
# ---------------------------------------------------------------------------
def test_damping_channel_equiv(tmp_path):
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    a0, a1 = 0.25, 1.5e-3
    _build_chain(masses, ks)
    lam = ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    xis = [a0 / (2.0 * math.sqrt(l)) + a1 * math.sqrt(l) / 2.0 for l in lam]

    pr = tmp_path / "ray.out"; pm = tmp_path / "mod.out"
    r = _run(lambda: ops.frequencyResponse(
        '-freq', 0.1, 4.0, 80, '-baseAccel', '-dir', 1, '-rayleigh', a0, a1,
        '-node', 3, '-dof', 1, '-out', str(pr)), pr)
    m = _run(lambda: ops.frequencyResponse(
        '-freq', 0.1, 4.0, 80, '-baseAccel', '-dir', 1, '-modalDamp', *xis,
        '-node', 3, '-dof', 1, '-out', str(pm)), pm)
    assert np.allclose(r, m, rtol=1e-10, atol=1e-14)


# ---------------------------------------------------------------------------
# uniform-xi channel vs the modal reference (covers -damp against numpy)
# ---------------------------------------------------------------------------
def test_uniform_damping_vs_modal_reference(tmp_path):
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    xi = 0.05
    fmin, fmax, nf = 0.05, 5.0, 150
    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    path = tmp_path / "uni.out"
    arr = _run(lambda: ops.frequencyResponse(
        '-freq', fmin, fmax, nf, '-baseAccel', '-dir', 1, '-damp', xi,
        '-node', 4, '-dof', 1, '-out', str(path)), path)

    M, K = _chain_MK(masses, ks)
    freqs = np.linspace(fmin, fmax, nf)
    ref = frf_modal(M, K, lambda w: 2.0 * xi * w, np.ones(3), freqs, dof=2)
    got = arr[:, 1] + 1j * arr[:, 2]
    assert np.max(np.abs(got - ref)) / np.max(np.abs(ref)) < 1e-8


# ---------------------------------------------------------------------------
# grid checks: log spacing + biased superset
# ---------------------------------------------------------------------------
def test_log_grid_and_biased_superset(tmp_path):
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()

    pl = tmp_path / "log.out"
    lg = _run(lambda: ops.frequencyResponse(
        '-freq', 0.1, 10.0, 25, '-log', '-baseAccel', '-dir', 1, '-damp', 0.03,
        '-node', 3, '-dof', 1, '-out', str(pl)), pl)
    f = lg[:, 0]
    assert abs(f[0] - 0.1) < 1e-12 and abs(f[-1] - 10.0) < 1e-9
    ratios = f[1:] / f[:-1]
    assert np.allclose(ratios, ratios[0], rtol=1e-9)  # geometric

    pb = tmp_path / "bias.out"
    bi = _run(lambda: ops.frequencyResponse(
        '-freq', 0.1, 10.0, 25, '-biased', '-baseAccel', '-dir', 1, '-damp', 0.03,
        '-node', 3, '-dof', 1, '-out', str(pb)), pb)
    assert bi.shape[0] > 25                      # clustered extras added
    assert np.all(np.diff(bi[:, 0]) > 0)         # strictly sorted, de-duped


# ---------------------------------------------------------------------------
# refuse-loudly guards
# ---------------------------------------------------------------------------
def test_guard_no_modalproperties():
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1)
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-baseAccel', '-dir', 1,
                              '-damp', 0.05, '-node', 2, '-dof', 1)


def test_guard_bad_dof():
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-baseAccel', '-dir', 1,
                              '-damp', 0.05, '-node', 2, '-dof', 3)  # ndf==1


def test_guard_bad_node():
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-baseAccel', '-dir', 1,
                              '-damp', 0.05, '-node', 99, '-dof', 1)


def test_guard_missing_freq():
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-baseAccel', '-dir', 1, '-damp', 0.05,
                              '-node', 2, '-dof', 1)


def test_guard_missing_damping():
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-baseAccel', '-dir', 1,
                              '-node', 2, '-dof', 1)


def test_guard_duplicate_modes():
    _build_chain([2.0, 1.5], [800.0, 500.0])
    ops.eigen('-fullGenLapack', 2); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-baseAccel', '-dir', 1,
                              '-damp', 0.05, '-node', 3, '-dof', 1, '-modes', 1, 1)


def test_guard_bad_resp():
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-baseAccel', '-dir', 1,
                              '-damp', 0.05, '-node', 2, '-dof', 1, '-resp', 'bogus')


# NOTE: the element-wise stale-`modalProperties` guard in run() (the same as
# modalResponseHistory's check()) is intentionally NOT unit-tested here: `modalProperties()` installs a
# DirectIntegrationAnalysis that shadows `eigen`, so a second same-NEV `eigen`
# fails ("no EigenSOE") before a stale snapshot can be constructed via openseespy.
# The guard never false-positives on the in-sync path (all compute tests above run
# eigen -> modalProperties -> frequencyResponse and pass), and is the same code as
# the modalResponseHistory staleness check.


def test_guard_negative_damp():
    # negative -damp xi is refused: d_a=2 xi w<0 CONJUGATES the FRF (same |H|,
    # flipped phase), so a '-0.05' typo would pass every magnitude check. (xi=0 is
    # legal: the undamped resonance is a genuine singularity.)
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-baseAccel', '-dir', 1,
                              '-damp', -0.05, '-node', 2, '-dof', 1)


def test_guard_negative_modaldamp_entry():
    # a single negative entry in the -modalDamp list is refused (same trap).
    _build_chain([2.0, 1.5], [800.0, 500.0])
    ops.eigen('-fullGenLapack', 2); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-baseAccel', '-dir', 1,
                              '-modalDamp', 0.03, -0.01, '-node', 3, '-dof', 1)


def test_guard_degenerate_sweep():
    # fmax==fmin with nf>1 would emit nf identical rows -> refuse (use nf==1).
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 3.0, 3.0, 5, '-baseAccel', '-dir', 1,
                              '-damp', 0.05, '-node', 2, '-dof', 1)


def test_fixed_response_node_is_zero(tmp_path):
    # A fully fixed node still gets a (zero) eigenvector row, so getNumEigenvectors
    # is >=1 and the non-exiting probe passes; the FRF is then identically zero
    # (the node cannot move). The point: Node::getEigenvectors() exit(0) must NOT be
    # reached; that path is reserved for a genuine orphan node (no eigenvector at
    # all), which the probe would refuse. Here the run completes and returns 0.
    _build_sdof(2.0, 500.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    path = tmp_path / "fixed.out"
    arr = _run(lambda: ops.frequencyResponse(
        '-freq', 0.1, 5.0, 5, '-baseAccel', '-dir', 1, '-damp', 0.05,
        '-node', 1, '-dof', 1, '-out', str(path)), path)
    assert np.allclose(arr[:, 1:], 0.0)


# ===========================================================================
# -load channel (nodal-force FRF)
#
#   u(Om) = amp * sum_a psi_a H_a(Om) (psi_a^T P)/m~_a,  psi_a = evec*Vscale,
#   m~_a = generalizedMasses()(a), normalization-invariant (any Vscale).
# ===========================================================================
def frf_direct_load(M, K, C, P, freqs, dof, resp="disp", amp=1.0):
    """Independent direct complex solve (K - Om^2 M + i Om C)^{-1} P."""
    out = []
    for f in freqs:
        Om = TWO_PI * f
        u = np.linalg.solve(K - Om * Om * M + 1j * Om * C, P)[dof]
        if resp == "vel":
            u *= 1j * Om
        elif resp == "accel":
            u *= -Om * Om
        out.append(amp * u)
    return np.array(out)


def _add_load_pattern(tag, loads):
    """Plain pattern holding nodal loads {nodeTag: value} (ndf=1 chains)."""
    ops.timeSeries('Constant', 100 + tag)
    ops.pattern('Plain', tag, 100 + tag)
    for n, v in loads.items():
        ops.load(n, v)


@pytest.mark.parametrize("resp", ["disp", "vel", "accel"])
def test_load_frf_vs_direct(tmp_path, resp):
    # multi-node force shape vs the independent
    # direct complex solve (does NOT reuse the modal formula).
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    a0, a1 = 0.30, 2.0e-3
    fmin, fmax, nf = 0.05, 5.0, 220

    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    _add_load_pattern(1, {3: 5.0, 4: -2.0})   # nodes 3,4 = dofs 1,2
    path = tmp_path / f"lfrf_{resp}.out"
    arr = _run(lambda: ops.frequencyResponse(
        '-freq', fmin, fmax, nf, '-lin', '-load', 1,
        '-rayleigh', a0, a1, '-node', 4, '-dof', 1, '-resp', resp,
        '-out', str(path)), path)

    freqs = np.linspace(fmin, fmax, nf)
    M, K = _chain_MK(masses, ks)
    C = a0 * M + a1 * K
    P = np.array([0.0, 5.0, -2.0])
    ref = frf_direct_load(M, K, C, P, freqs, dof=2, resp=resp)

    got = arr[:, 1] + 1j * arr[:, 2]
    scale = np.max(np.abs(ref))
    assert scale > 1e-12
    err = np.max(np.abs(got - ref)) / scale
    assert err < 1e-8, f"{resp}: -load FRF vs direct rel err {err:.2e}"


def test_load_frf_unorm_invariance(tmp_path):
    # modalProperties -unorm rescales the eigenvectors (Vscale != 1); the
    # psi(psi^T P)/m~ weight must be EXACTLY invariant. This is the pin that a
    # missing/extra Vscale factor (a Vscale^2 mistake) cannot survive.
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    _add_load_pattern(1, {2: 1.0, 4: 3.0})
    args = ('-freq', 0.05, 5.0, 120, '-lin', '-load', 1, '-damp', 0.04,
            '-node', 4, '-dof', 1)
    p1 = tmp_path / "def.out"
    a_def = _run(lambda: ops.frequencyResponse(*args, '-out', str(p1)), p1)

    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties('-unorm')
    _add_load_pattern(1, {2: 1.0, 4: 3.0})
    p2 = tmp_path / "unorm.out"
    a_un = _run(lambda: ops.frequencyResponse(*args, '-out', str(p2)), p2)

    assert np.allclose(a_def, a_un, rtol=1e-10, atol=1e-14)


def test_load_static_limit(tmp_path):
    # SSD at Om->0 equals the static solution |K^-1 P| (spectral identity,
    # exact with all modes retained).
    masses = [2.0, 3.0, 1.5]; ks = [1200.0, 900.0, 600.0]
    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    _add_load_pattern(1, {2: 4.0, 3: -1.0})
    path = tmp_path / "lstat.out"
    arr = _run(lambda: ops.steadyStateDynamics(
        '-freq', 1e-6, 1e-6, 1, '-load', 1, '-damp', 0.03,
        '-node', 4, '-dof', 1, '-out', str(path)), path)
    M, K = _chain_MK(masses, ks)
    ustat = np.linalg.solve(K, np.array([4.0, -1.0, 0.0]))
    assert abs(arr[0, 1] - abs(ustat[2])) / abs(ustat[2]) < 1e-6


def test_load_equals_baseaccel_for_P_minus_MR(tmp_path):
    # -baseAccel -dir 1 is the special case P = -M R: the two channels must
    # agree row-for-row on the same model.
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    _build_chain(masses, ks)
    ops.eigen('-fullGenLapack', 2)
    ops.modalProperties()
    _add_load_pattern(1, {2: -2.0, 3: -1.5})  # P = -M R (unit base accel)
    args = ('-freq', 0.1, 4.0, 80, '-damp', 0.04, '-node', 3, '-dof', 1)
    pb = tmp_path / "ba.out"; pl = tmp_path / "ld.out"
    a_b = _run(lambda: ops.frequencyResponse(
        *args, '-baseAccel', '-dir', 1, '-out', str(pb)), pb)
    a_l = _run(lambda: ops.frequencyResponse(
        *args, '-load', 1, '-out', str(pl)), pl)
    assert np.allclose(a_b, a_l, rtol=1e-10, atol=1e-14)


def test_load_amp_scaling(tmp_path):
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1)
    ops.modalProperties()
    _add_load_pattern(1, {2: 3.0})
    args = ('-freq', 0.5, 4.0, 40, '-load', 1, '-damp', 0.05,
            '-node', 2, '-dof', 1)
    p1 = tmp_path / "a1.out"; p2 = tmp_path / "a2.out"
    a1 = _run(lambda: ops.frequencyResponse(*args, '-out', str(p1)), p1)
    a2 = _run(lambda: ops.frequencyResponse(
        *args, '-amp', 2.5, '-out', str(p2)), p2)
    assert np.allclose(a2[:, 1:], 2.5 * a1[:, 1:], rtol=1e-12, atol=0)


def test_guard_load_and_baseaccel_exclusive():
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    _add_load_pattern(1, {2: 1.0})
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-baseAccel', '-dir', 1,
                              '-load', 1, '-damp', 0.05, '-node', 2, '-dof', 1)


def test_guard_load_missing_pattern():
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-load', 99,
                              '-damp', 0.05, '-node', 2, '-dof', 1)


def test_guard_load_pattern_without_nodal_loads():
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    ops.timeSeries('Constant', 105)
    ops.pattern('Plain', 5, 105)              # empty pattern
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-load', 5,
                              '-damp', 0.05, '-node', 2, '-dof', 1)


def test_load_moment_on_rotational_dof_vs_direct(tmp_path):
    # a MOMENT load on a rotational DOF (ndf=3 beam node) vs the independent
    # direct solve; the ndf=1 chains
    # above never exercise a rotational load entry.
    E, A, Iz, L = 30000.0, 20.0, 800.0, 100.0
    m, J = 3.0, 40.0                      # tip mass + rotary inertia
    xi = 0.04
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1, 1)
    ops.node(2, L, 0.0);   ops.mass(2, m, m, J)
    ops.geomTransf('Linear', 1)
    ops.element('elasticBeamColumn', 1, 1, 2, A, E, Iz, 1)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()
    ops.timeSeries('Constant', 120)
    ops.pattern('Plain', 20, 120)
    ops.load(2, 0.0, 2.0, 5.0)            # shear + tip MOMENT
    path = tmp_path / "moment.out"
    fmin, fmax, nf = 0.05, 4.0, 150
    arr = _run(lambda: ops.frequencyResponse(
        '-freq', fmin, fmax, nf, '-lin', '-load', 20, '-damp', xi,
        '-node', 2, '-dof', 3, '-out', str(path)), path)

    # cantilever free-node operator (Ux, Uy, Rz at the tip)
    K = np.array([[E * A / L, 0.0, 0.0],
                  [0.0, 12 * E * Iz / L ** 3, -6 * E * Iz / L ** 2],
                  [0.0, -6 * E * Iz / L ** 2, 4 * E * Iz / L]])
    M = np.diag([m, m, J])
    P = np.array([0.0, 2.0, 5.0])
    freqs = np.linspace(fmin, fmax, nf)
    # modal-load reference with uniform xi (same channel as -damp)
    Lam, Q = np.linalg.eigh(M)
    Mis = Q @ np.diag(1.0 / np.sqrt(Lam)) @ Q.T
    w2, U = np.linalg.eigh(Mis @ K @ Mis)
    Phi = Mis @ U
    w = np.sqrt(np.clip(w2, 0.0, None))
    ref = np.array([
        sum(Phi[2, a] * (Phi[:, a] @ P)
            / (w[a] ** 2 - (TWO_PI * f) ** 2 + 1j * (TWO_PI * f) * 2 * xi * w[a])
            for a in range(3))
        for f in freqs])
    got = arr[:, 1] + 1j * arr[:, 2]
    err = np.max(np.abs(got - ref)) / np.max(np.abs(ref))
    assert err < 1e-8, f"moment-load FRF vs direct rel err {err:.2e}"


def test_guard_load_pattern_with_sp():
    # a pattern carrying sp constraints does not define a pure force shape ->
    # refused loudly rather than silently dropping the sp part.
    _build_sdof(2.0, 800.0)
    ops.eigen('-fullGenLapack', 1); ops.modalProperties()
    ops.timeSeries('Constant', 106)
    ops.pattern('Plain', 6, 106)
    ops.load(2, 1.0)
    ops.sp(2, 1, 0.01)
    with pytest.raises(Exception):
        ops.frequencyResponse('-freq', 0.1, 5.0, 10, '-load', 6,
                              '-damp', 0.05, '-node', 2, '-dof', 1)
