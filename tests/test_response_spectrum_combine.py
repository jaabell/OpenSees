"""responseSpectrumAnalysis -combine: SRSS, CQC, ABS and TenPercent.

Without -combine, responseSpectrumAnalysis computes the per-mode modal
displacements and commits one domain step per mode, leaving the modal
combination to the user. -combine {SRSS|CQC|ABS|TenPercent} combines the
per-mode peak nodal displacements into a single committed field (CQC with the
correlation coefficient of Der Kiureghian 1981; TenPercent per U.S. NRC
Regulatory Guide 1.92 Rev. 1).

Reference: the per-mode fields of the unchanged `-mode k` path (one run per
mode, nodeDisp) combined in numpy by the textbook rules, compared with the
-combine field. Both use u = V MPF Sa / lambda, so the check isolates the
combination.

  test_combine_matches_permode[SRSS/CQC/ABS/TenPercent]  to 1e-9; fails on a
        build without -combine (the option is ignored there and the last
        mode's field is committed).
  test_cqc_modalDamp_list   per-mode damping ratios through -modalDamp.
  test_cqc_bounded_by_abs   0 <= CQC <= ABS.
  test_default_unchanged    without -combine, `-mode k` still yields the
        analytic per-mode field.
  test_cqc_needs_damping / test_combine_mode_exclusive / test_bad_rule  refusals.

Small models use `eigen -fullGenLapack` (ARPACK needs NEV < N).
"""
import math

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


# --- numpy reference combination (Der Kiureghian CQC) ----------------------
def _rho(wi, wj, xi_i, xi_j):
    r = wj / wi
    num = 8.0 * math.sqrt(xi_i * xi_j) * (xi_i + r * xi_j) * r ** 1.5
    den = (1 - r * r) ** 2 + 4 * xi_i * xi_j * r * (1 + r * r) \
        + 4 * (xi_i ** 2 + xi_j ** 2) * r * r
    return num / den


def combine(rule, R, w, xi, T):
    n = len(R); R = np.asarray(R, float)
    if rule == "ABS":
        return float(np.sum(np.abs(R)))
    if rule == "SRSS":
        return float(np.sqrt(np.sum(R * R)))
    if rule == "TenPercent":
        s = float(np.sum(R * R))
        for i in range(n):
            for j in range(i + 1, n):
                if abs(T[i] - T[j]) <= 0.10 * min(T[i], T[j]):
                    s += 2.0 * abs(R[i] * R[j])
        return math.sqrt(max(s, 0.0))
    # CQC
    s = 0.0
    for i in range(n):
        for j in range(n):
            rho = 1.0 if i == j else _rho(w[i], w[j], xi[i], xi[j])
            s += rho * R[i] * R[j]
    return math.sqrt(max(s, 0.0))


# --- model: 1-D axial chain (ndf=1) ----------------------------------------
def _build_chain(masses, ks):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0); ops.fix(1, 1)
    for i, mm in enumerate(masses, start=2):
        ops.node(i, float(i - 1)); ops.mass(i, mm)
    for e, kk in enumerate(ks, start=1):
        ops.uniaxialMaterial('Elastic', e, kk)
        ops.element('Truss', e, e, e + 1, 1.0, e)


TN = [0.02, 0.1, 0.3, 0.6, 1.0, 3.0]
SA = [0.40, 1.00, 1.20, 0.90, 0.60, 0.20]


def _setup(masses=(2.0, 1.5, 1.0), ks=(1200.0, 800.0, 500.0)):
    _build_chain(list(masses), list(ks))
    lam = ops.eigen('-fullGenLapack', len(masses))
    ops.modalProperties()
    return lam


def _permode_fields(nmodes, nodes, dir=1):
    """per-mode nodal disp R[mode][node] from the existing -mode k runs."""
    R = np.zeros((nmodes, len(nodes)))
    for k in range(1, nmodes + 1):
        ops.responseSpectrumAnalysis(dir, '-Tn', *TN, '-Sa', *SA, '-mode', k)
        for jn, nd in enumerate(nodes):
            R[k - 1, jn] = ops.nodeDisp(nd, 1)
    return R


@pytest.mark.parametrize("rule", ["SRSS", "CQC", "ABS", "TenPercent"])
def test_combine_matches_permode(rule):
    masses = [2.0, 1.5, 1.0]; ks = [1200.0, 800.0, 500.0]
    xi0 = 0.05
    nodes = [2, 3, 4]
    lam = _setup(masses, ks)
    nm = len(lam)
    w = [math.sqrt(l) for l in lam]
    T = [2 * math.pi / wi for wi in w]
    xi = [xi0] * nm

    # per-mode reference from the existing -mode runs
    R = _permode_fields(nm, nodes)
    ref = [combine(rule, R[:, jn], w, xi, T) for jn in range(len(nodes))]

    # the -combine field
    args = [1, '-Tn', *TN, '-Sa', *SA, '-combine', rule]
    if rule == "CQC":
        args += ['-damp', xi0]
    ops.responseSpectrumAnalysis(*args)
    got = [ops.nodeDisp(nd, 1) for nd in nodes]

    for jn, nd in enumerate(nodes):
        assert abs(got[jn] - ref[jn]) < 1e-9, \
            f"{rule} node {nd}: got {got[jn]:.12e} ref {ref[jn]:.12e}"
        # combined magnitude is non-negative
        assert got[jn] >= -1e-14


def test_cqc_bounded_by_abs():
    # CQC (with 0<=rho<=1) is always bounded above by ABS and is non-negative.
    # NOTE: CQC is NOT always >= SRSS: with mixed-sign modal contributions the
    # positive-rho cross terms are negative and pull CQC below SRSS (correct).
    masses = [2.0, 1.5, 1.0]; ks = [1200.0, 800.0, 500.0]
    _setup(masses, ks)
    def field(rule, extra=()):
        ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-combine', rule, *extra)
        return abs(ops.nodeDisp(4, 1))
    cqc = field("CQC", ('-damp', 0.05)); ab = field("ABS")
    assert 0.0 <= cqc <= ab + 1e-9, (cqc, ab)


def test_default_unchanged():
    # without -combine, the existing per-mode path must still produce the analytic
    # mode-1 field: u = V*Vscale*MPF*Sa/lambda
    masses = [2.0, 1.5]; ks = [800.0, 500.0]
    lam = _setup(masses, ks)
    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-mode', 1)
    u3 = ops.nodeDisp(3, 1)
    # analytic mode-1 at node 3
    w1 = math.sqrt(lam[0]); T1 = 2 * math.pi / w1
    Sa1 = np.interp(T1, TN, SA)
    phi = ops.nodeEigenvector(3, 1, 1)
    # participation factor / generalized mass from the model
    M = np.diag(masses); R = np.ones(2)
    phivec = np.array([ops.nodeEigenvector(2, 1, 1), ops.nodeEigenvector(3, 1, 1)])
    Gamma = (phivec @ M @ R) / (phivec @ M @ phivec)
    u_ref = phi * Gamma * Sa1 / lam[0]
    assert abs(u3 - u_ref) < 1e-9, (u3, u_ref)


def test_cqc_needs_damping():
    _setup([2.0, 1.5], [800.0, 500.0])
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-combine', 'CQC')


def test_cqc_modalDamp_list():
    # exercise a FULL per-mode -modalDamp list (distinct ratios) and check the
    # CQC field matches the numpy reference using those per-mode ratios. This is the
    # only test that drives the multi-value -modalDamp list parser.
    masses = [2.0, 1.5, 1.0]; ks = [1200.0, 800.0, 500.0]
    nodes = [2, 3, 4]
    lam = _setup(masses, ks)
    nm = len(lam)
    w = [math.sqrt(l) for l in lam]
    T = [2 * math.pi / wi for wi in w]
    xis = [0.03, 0.05, 0.07][:nm]

    R = _permode_fields(nm, nodes)
    ref = [combine("CQC", R[:, jn], w, xis, T) for jn in range(len(nodes))]

    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA,
                                 '-combine', 'CQC', '-modalDamp', *xis)
    got = [ops.nodeDisp(nd, 1) for nd in nodes]
    for jn, nd in enumerate(nodes):
        assert math.isfinite(got[jn]), f"node {nd}: non-finite {got[jn]}"
        assert abs(got[jn] - ref[jn]) < 1e-9, \
            f"CQC modalDamp node {nd}: got {got[jn]:.12e} ref {ref[jn]:.12e}"


def test_combine_mode_exclusive():
    _setup([2.0, 1.5], [800.0, 500.0])
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA,
                                     '-combine', 'SRSS', '-mode', 1)


def test_bad_rule():
    _setup([2.0, 1.5], [800.0, 500.0])
    with pytest.raises(Exception):
        ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-combine', 'BOGUS')
