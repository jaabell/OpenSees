"""complexEigen: complex (state-space) modal analysis of non-classically
damped models.

Kernel entry `complexEigen -qz p Mt Ct Kt`: the quadratic eigenproblem
(lambda^2 Mt + lambda Ct + Kt) z = 0 is solved through the symmetric-block
state-space pencil and LAPACK dggev. Checked against an independent numpy
oracle that uses a different linearization (the companion form
d/dt {u; v} = [[0, I], [-M^-1 K, -M^-1 C]] {u; v} with numpy.linalg.eig) to
1e-10 on hand-built and seeded random pencils, and against the analytic SDOF
roots (underdamped, overdamped, rigid plus decay).

Domain entry `complexEigen [numModes]`: after `eigen`, M and C are projected
element by element onto the real modes (getMass/getDamp plus the nodal mass
and alphaM terms) and Kt = sym(Mt diag(w^2)).

  closed form, pure global Rayleigh: zeta_a = aM/(2 w_a) + bK w_a/2 and
      w_d = w_a sqrt(1 - zeta^2) against the domain's own eigen to 1e-8; the
      default projection equals the -closedForm path to 1e-10.
  non-classical: a dashpot on one DOF of a 2-DOF chain (p = N, so the
      projection is exact) against the full-order state-space eigenvalues to
      1e-8; the same for global betaK on Truss elements with a zeroLength
      dashpot (zeroLength -doRayleigh defaults to 0), region -rayleigh on one
      element, betaKinit, and a dashpot on an equalDOF retained-node
      dependent.
  base-isolated two-mass stick: isolation-mode zeta/omegaD against the exact
      eigenvalues to 1e-8 and against the logarithmic decrement and damped
      period of a Newmark free decay to 5 %.
  Node recorder complexEigenRe<N>/complexEigenIm<N>: recorded shapes carry
      the phase convention (largest entry real-positive), a phased shape for
      the non-classical chain, and the real eigenvector for classical damping.
  guards: no prior eigen (also after wipe), more modes than eigen produced,
      betaKinit/betaKcomm under -closedForm, malformed -qz input.

Every test fails on a build without complexEigen (unknown command).
"""
import math

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


# ---------------------------------------------------------------------------
# independent numpy oracle (companion form, not the symmetric-block pencil)
# ---------------------------------------------------------------------------
def raw_roots(M, C, K):
    """All 2p eigenvalues of (lambda^2 M + lambda C + K) z = 0, unsorted."""
    M = np.asarray(M, dtype=float)
    C = np.asarray(C, dtype=float)
    K = np.asarray(K, dtype=float)
    p = M.shape[0]
    Minv = np.linalg.inv(M)
    comp = np.zeros((2 * p, 2 * p))
    comp[:p, p:] = np.eye(p)
    comp[p:, :p] = -Minv @ K
    comp[p:, p:] = -Minv @ C
    return np.linalg.eigvals(comp)


def sorted_roots(M, C, K):
    """Raw roots sorted by (|lambda|, Re, Im) for multiset comparison."""
    lam = raw_roots(M, C, K)
    order = np.lexsort((lam.imag, lam.real, np.abs(lam)))
    return lam[order]


STRIDE = 7  # omega0, omegaD, zeta, Re, Im, kind, resid


def expand_cpp_modes(flat):
    """Expand the complexEigen return (STRIDE doubles per reported mode)
    back into the raw-root multiset: an underdamped entry (kind 0)
    contributes lambda and conj(lambda); overdamped/rigid entries (kind 1/2)
    contribute their single real root."""
    flat = list(flat)
    assert len(flat) % STRIDE == 0, \
        f"expected {STRIDE} doubles per mode, got {len(flat)}"
    roots = []
    for i in range(0, len(flat), STRIDE):
        re, im, kind = flat[i + 3], flat[i + 4], int(round(flat[i + 5]))
        if kind == 0:
            roots.append(complex(re, im))
            roots.append(complex(re, -im))
        else:
            roots.append(complex(re, 0.0))
    lam = np.array(roots)
    order = np.lexsort((lam.imag, lam.real, np.abs(lam)))
    return lam[order]


def modes_table(flat):
    """The complexEigen return as a list of dicts, in the reporting order."""
    kinds = {0: "underdamped", 1: "overdamped", 2: "rigid"}
    out = []
    for i in range(0, len(flat), STRIDE):
        out.append({
            "omega0": flat[i],
            "omegaD": flat[i + 1],
            "zeta": flat[i + 2],
            "lambda": complex(flat[i + 3], flat[i + 4]),
            "kind": kinds[int(round(flat[i + 5]))],
            "resid": flat[i + 6],
        })
    return out

def qz(M, C, K, tol=None):
    """Call the C++ kernel; returns the flat -qz output list."""
    M = np.atleast_2d(np.asarray(M, dtype=float))
    C = np.atleast_2d(np.asarray(C, dtype=float))
    K = np.atleast_2d(np.asarray(K, dtype=float))
    p = M.shape[0]
    ops.wipe()
    args = ["-qz", p]
    args += [float(x) for x in M.flatten()]
    args += [float(x) for x in C.flatten()]
    args += [float(x) for x in K.flatten()]
    if tol is not None:
        args += ["-tol", float(tol)]
    out = ops.complexEigen(*args)
    if isinstance(out, float):
        out = (out,)
    return out


def assert_roots_match(flat, M, C, K, rtol=1e-10):
    got = expand_cpp_modes(flat)
    ref = sorted_roots(M, C, K)
    assert len(got) == len(ref)
    scale = max(np.max(np.abs(ref)), 1.0)
    np.testing.assert_allclose(got, ref, rtol=0.0, atol=rtol * scale)


def assert_residuals_small(flat, M, C, K, rtol=1e-9):
    """Every reported mode must satisfy its own quadratic residual; this
    ties the reported lambda to the reported eigenvector (catches
    conjugate-pairing mistakes the eigenvalue-only comparison cannot see)."""
    K = np.atleast_2d(np.asarray(K, dtype=float))
    scale = max(np.linalg.norm(K, ord=np.inf), 1.0)
    for m in modes_table(flat):
        assert m["resid"] <= rtol * scale, \
            f"mode lambda={m['lambda']}: resid {m['resid']} > {rtol * scale}"


# ------------------------------------------------------------------ analytic
def test_sdof_underdamped():
    w = 2.0 * math.pi * 1.7
    zeta = 0.05
    flat = qz([[1.0]], [[2.0 * zeta * w]], [[w * w]])
    modes = modes_table(flat)
    assert len(modes) == 1
    m = modes[0]
    assert m["kind"] == "underdamped"
    assert m["omega0"] == pytest.approx(w, rel=1e-10)
    assert m["zeta"] == pytest.approx(zeta, rel=1e-10)
    assert m["omegaD"] == pytest.approx(w * math.sqrt(1 - zeta**2), rel=1e-10)
    assert m["lambda"].real == pytest.approx(-zeta * w, rel=1e-10)


def test_sdof_overdamped():
    w = 3.0
    zeta = 1.5
    flat = qz([[1.0]], [[2.0 * zeta * w]], [[w * w]])
    modes = modes_table(flat)
    assert len(modes) == 2
    lam_exact = sorted([-zeta * w + w * math.sqrt(zeta**2 - 1),
                        -zeta * w - w * math.sqrt(zeta**2 - 1)])
    lam_got = sorted(m["lambda"].real for m in modes)
    for got, exact in zip(lam_got, lam_exact):
        assert got == pytest.approx(exact, rel=1e-10)
    for m in modes:
        assert m["kind"] == "overdamped"
        assert m["omegaD"] == 0.0
        assert m["zeta"] == pytest.approx(1.0, rel=1e-12)  # -Re/|lambda|


def test_rigid_plus_decay():
    c = 0.3
    flat = qz([[1.0]], [[c]], [[0.0]])
    modes = modes_table(flat)
    assert len(modes) == 2
    kinds = sorted(m["kind"] for m in modes)
    assert kinds == ["overdamped", "rigid"]
    rigid = next(m for m in modes if m["kind"] == "rigid")
    decay = next(m for m in modes if m["kind"] == "overdamped")
    assert rigid["omega0"] == 0.0
    assert abs(rigid["lambda"]) <= 1e-10 * c
    assert decay["lambda"].real == pytest.approx(-c, rel=1e-10)


# ----------------------------------------------------- classical / oracle
def test_classical_rayleigh():
    """Rayleigh-only reduced matrices are the classical case: every mode
    stays 'real' (a plain damped oscillator per mode) with the textbook
    zeta_k = aM/(2 w_k) + bK w_k / 2, checked directly on the kernel."""
    w = np.array([2.1, 5.3, 9.8, 17.2])
    aM, bK = 0.35, 0.004
    Mt = np.eye(4)
    Kt = np.diag(w**2)
    Ct = aM * Mt + bK * Kt
    flat = qz(Mt, Ct, Kt)
    modes = modes_table(flat)
    assert len(modes) == 4
    for m, wk in zip(modes, w):
        zk = aM / (2 * wk) + bK * wk / 2
        assert m["kind"] == "underdamped"
        assert m["omega0"] == pytest.approx(wk, rel=1e-10)
        assert m["zeta"] == pytest.approx(zk, rel=1e-10)
    assert_roots_match(flat, Mt, Ct, Kt)
    assert_residuals_small(flat, Mt, Ct, Kt)


def test_2dof_nonclassical():
    """Two-mass chain with a dashpot on ONE dof, the canonical
    non-classical problem (C not simultaneously diagonalizable)."""
    M = np.eye(2)
    K = np.array([[5.0, -1.0], [-1.0, 1.0]])
    C = np.array([[0.0, 0.0], [0.0, 0.4]])
    flat = qz(M, C, K)
    modes = modes_table(flat)
    assert len(modes) == 2
    assert all(m["kind"] == "underdamped" for m in modes)
    # non-classical: the two damping ratios are NOT a Rayleigh pair's
    # signature; just require both positive and distinct
    z1, z2 = modes[0]["zeta"], modes[1]["zeta"]
    assert z1 > 0 and z2 > 0 and abs(z1 - z2) > 1e-4
    assert_roots_match(flat, M, C, K)
    assert_residuals_small(flat, M, C, K)


def test_random_pencils_vs_numpy():
    rng = np.random.default_rng(4600)
    for p in (3, 8):
        for _ in range(3):
            w = np.sort(rng.uniform(1.0, 30.0, size=p))
            Mt = np.eye(p)
            Kt = np.diag(w**2)
            # symmetric PSD damping with off-diagonal (non-classical) coupling
            R = rng.uniform(-1.0, 1.0, size=(p, p))
            Ct = 0.05 * (R @ R.T) + 0.02 * np.diag(w)
            flat = qz(Mt, Ct, Kt)
            assert_roots_match(flat, Mt, Ct, Kt)
            # the residual assertion ties each reported lambda to ITS
            # eigenvector, the check that catches beta<0 conjugate
            # mis-pairing
            assert_residuals_small(flat, Mt, Ct, Kt)
            # Ct is PSD by construction, so all modes decay
            for m in modes_table(flat):
                if m["kind"] == "underdamped":
                    assert -1e-12 < m["zeta"] < 1.0


# ------------------------------------------------------------------- guards
def test_bad_input_rejected():
    with pytest.raises(Exception):
        ops.complexEigen("-qz", 2, 1.0, 2.0, 3.0)  # far too few entries
    with pytest.raises(Exception):
        ops.complexEigen("-notAMode", 1, 1.0, 1.0, 1.0)


# ------------------------------------------------- domain path, closed form
def _cantilever(n_ele=8, n_eigen=4):
    """2D elastic cantilever with lumped translational masses; returns the
    undamped circular frequencies from the domain's own `eigen` run.
    8 elements (24 equations): ARPACK needs nev + 8 below the number of
    equations."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    L, E, A, Iz, m = 3.0, 30.0e9, 0.09, 6.75e-4, 8.0e3
    for i in range(n_ele + 1):
        ops.node(i + 1, 0.0, i * L)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    for i in range(n_ele):
        ops.element("elasticBeamColumn", i + 1, i + 1, i + 2, A, E, Iz, 1)
    for i in range(2, n_ele + 2):
        ops.mass(i, m, m, 0.0)
    lam = ops.eigen(n_eigen)
    return np.sqrt(np.asarray(lam))


def test_classical_frame():
    """On a pure global-Rayleigh model the complex solve
    must reproduce the classical zeta_a = aM/(2 w_a) + bK w_a / 2 and
    w_d = w_a sqrt(1 - zeta^2) against the domain's own eigen frequencies."""
    w = _cantilever()
    aM, bK = 0.25, 0.0008
    ops.rayleigh(aM, bK, 0.0, 0.0)
    modes = modes_table(ops.complexEigen())
    assert len(modes) == len(w)
    for m, wk in zip(modes, w):
        zk = aM / (2 * wk) + bK * wk / 2
        assert m["kind"] == "underdamped"
        assert m["omega0"] == pytest.approx(wk, rel=1e-8)
        assert m["zeta"] == pytest.approx(zk, rel=1e-8)
        assert m["omegaD"] == pytest.approx(wk * math.sqrt(1 - zk**2),
                                            rel=1e-8)
        assert m["lambda"].real == pytest.approx(-zk * wk, rel=1e-8)
        assert m["resid"] <= 1e-9 * max(w) ** 2


def test_undamped_domain():
    """No rayleigh call: factors are zero, so lambda = +- i w exactly.

    Doubles as the wipe regression: the preceding test set
    rayleigh(0.25, ...); if wipe() did not reset the domain Rayleigh copy,
    zeta here would be ~0.06, not 0."""
    w = _cantilever()
    modes = modes_table(ops.complexEigen())
    assert len(modes) == len(w)
    for m, wk in zip(modes, w):
        assert m["kind"] == "underdamped"
        assert m["omega0"] == pytest.approx(wk, rel=1e-10)
        assert abs(m["zeta"]) <= 1e-12
        assert m["omegaD"] == pytest.approx(wk, rel=1e-10)


def test_mode_subset():
    w = _cantilever(n_eigen=4)
    ops.rayleigh(0.1, 0.001, 0.0, 0.0)
    got2 = modes_table(ops.complexEigen(2))
    assert len(got2) == 2
    got3 = modes_table(ops.complexEigen("-numModes", 3, "-tol", 1.0e-9))
    assert len(got3) == 3
    for m, wk in zip(got3, w[:3]):
        assert m["omega0"] == pytest.approx(wk, rel=1e-8)
    # requesting more modes than eigen produced must be refused, not padded
    with pytest.raises(Exception):
        ops.complexEigen(5)


def test_closed_form_scoped_rayleigh_warns(capfd):
    """region -rayleigh writes factors straight onto elements and bypasses
    the domain copy; the closed form must WARN, not silently
    report a zeta that ignores the scoped damping."""
    _cantilever()
    ops.rayleigh(0.1, 0.001, 0.0, 0.0)
    modes_table(ops.complexEigen("-closedForm"))
    assert "differ from the last global" not in capfd.readouterr().err
    ops.region(1, "-ele", 1, 2, "-rayleigh", 0.5, 0.0, 0.0, 0.0)
    modes = modes_table(ops.complexEigen("-closedForm"))
    assert len(modes) == 4  # still answers (global factors), but flagged:
    assert "differ from the last global" in capfd.readouterr().err
    # the DEFAULT assembled path handles scoped damping exactly -> no warning
    modes_table(ops.complexEigen())
    assert "differ from the last global" not in capfd.readouterr().err


# ------------------------------------- domain path, element-by-element projection
def _chain2dof(k1=400.0, k2=150.0, m2=2.0, m3=1.0, c=1.2, dashpot=True):
    """1D two-DOF spring-mass chain, optional grounded dashpot on DOF 2
    (node 3), the canonical non-classical textbook system. Returns the
    exact (M, K, C_dashpot_only) numpy operators for oracle use."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.node(3, 2.0)
    ops.fix(1, 1)
    ops.mass(2, m2)
    ops.mass(3, m3)
    # truss springs: k = EA/L with A=L=1 -> E=k
    ops.uniaxialMaterial("Elastic", 1, k1)
    ops.uniaxialMaterial("Elastic", 2, k2)
    # -doRayleigh 1 explicit: Truss (like the zeroLength family) defaults the
    # flag to 0, so by default trusses feel no global Rayleigh damping. The
    # assembled projection reproduces that (C is exactly getDamp()), so the
    # oracle must model elements that do carry the factors. With default
    # trusses the assembled zeta is zero while the closed form claims
    # bK*w/2: the closed form is the one that differs from the transient.
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-doRayleigh", 1)
    ops.element("Truss", 2, 2, 3, 1.0, 2, "-doRayleigh", 1)
    if dashpot:
        ops.node(4, 2.0)
        ops.fix(4, 1)
        ops.uniaxialMaterial("Viscous", 3, c, 1.0)
        ops.element("zeroLength", 3, 4, 3, "-mat", 3, "-dir", 1)
    M = np.diag([m2, m3])
    K = np.array([[k1 + k2, -k2], [-k2, k2]])
    C = np.zeros((2, 2))
    if dashpot:
        C[1, 1] = c
    return M, K, C


def _domain_lambdas(flat):
    """Reported physical lambdas (+Im representative) sorted by omega0."""
    return [m["lambda"] for m in modes_table(flat)]


def _oracle_lambdas(M, C, K):
    """Independent full-order state-space eig; one root per physical mode
    (the +Im representative of each conjugate pair), sorted by |lambda|."""
    n = M.shape[0]
    A = np.block([[np.zeros((n, n)), np.eye(n)],
                  [-np.linalg.solve(M, K), -np.linalg.solve(M, C)]])
    lams = np.linalg.eigvals(A)
    keep = [l for l in lams if l.imag > 1e-12]
    keep += [l.real + 0j for l in lams if abs(l.imag) <= 1e-12]
    return sorted(keep, key=lambda l: (abs(l), l.real))


def test_assembled_equals_closed_form():
    """On a pure GLOBAL-Rayleigh model the assembled projection must
    reproduce the closed form to near machine precision; pins the element
    AND nodal-mass/nodal-alphaM loops at once (cantilever mass is nodal)."""
    _cantilever()
    ops.rayleigh(0.3, 0.0006, 0.0, 0.0)
    assembled = modes_table(ops.complexEigen())
    closed = modes_table(ops.complexEigen("-closedForm"))
    assert len(assembled) == len(closed) == 4
    for b, a in zip(assembled, closed):
        assert b["lambda"].real == pytest.approx(a["lambda"].real, rel=1e-10)
        assert b["lambda"].imag == pytest.approx(a["lambda"].imag, rel=1e-10)
        assert b["zeta"] == pytest.approx(a["zeta"], rel=1e-10)
        assert b["resid"] <= 1e-8


def test_2dof_nonclassical_domain():
    """A dashpot on ONE dof of a 2-DOF chain:
    C is not simultaneously diagonalizable with M and K. With p = N = 2 the
    projection is EXACT, so the domain solve must match the full-order
    state-space oracle to 1e-8."""
    M, K, C = _chain2dof()
    ops.eigen("-fullGenLapack", 2)
    got = _domain_lambdas(ops.complexEigen())
    ref = _oracle_lambdas(M, C, K)
    assert len(got) == 2
    scale = max(abs(r) for r in ref)
    for g, r in zip(got, ref):
        assert abs(g - r) <= 1e-8 * scale
    # non-classical fingerprint: the two zetas are distinct and both nonzero
    z = [m["zeta"] for m in modes_table(ops.complexEigen())]
    assert z[0] > 1e-4 and z[1] > 1e-4 and abs(z[0] - z[1]) > 1e-4


def test_zerolength_doRayleigh_default():
    """zeroLength -doRayleigh defaults to 0: global betaK adds
    stiffness damping to the TRUSSES but NOT to the dashpot's zeroLength.
    The assembled C must equal bK*K_trusses + C_dashpot exactly."""
    bK = 0.002
    M, K, C = _chain2dof()
    ops.rayleigh(0.0, bK, 0.0, 0.0)
    ops.eigen("-fullGenLapack", 2)
    got = _domain_lambdas(ops.complexEigen())
    # oracle: trusses carry bK*K (zeroLength has NO stiffness and doRayleigh=0
    # so it contributes only the material dashpot)
    ref = _oracle_lambdas(M, C + bK * K, K)
    scale = max(abs(r) for r in ref)
    for g, r in zip(got, ref):
        assert abs(g - r) <= 1e-8 * scale


def test_betaK0_assembled():
    """betaKinit is refused by the closed form but handled by the assembled
    projection; on a
    LINEAR elastic model K0 == Kt so zeta = bK0*w/2 classically."""
    bK0 = 0.0012
    w = _cantilever()
    ops.rayleigh(0.0, 0.0, bK0, 0.0)
    modes = modes_table(ops.complexEigen())
    for m, wk in zip(modes, w):
        assert m["zeta"] == pytest.approx(bK0 * wk / 2, rel=1e-8)
    with pytest.raises(Exception):
        ops.complexEigen("-closedForm")


def test_scoped_region_exact():
    """region -rayleigh on ONE element: silently wrong under -closedForm
    (it warns), exact under the assembled default."""
    bKa = 0.004
    M, K, C = _chain2dof(dashpot=False)
    k1 = 400.0
    ops.region(1, "-ele", 1, "-rayleigh", 0.0, bKa, 0.0, 0.0)
    ops.eigen("-fullGenLapack", 2)
    got = _domain_lambdas(ops.complexEigen())
    # element 1 stiffness matrix restricted to free DOFs: k1 on DOF (2,2)...
    # truss 1 connects fixed node 1 to node 2 -> contributes k1 at (0,0)
    Kel1 = np.array([[k1, 0.0], [0.0, 0.0]])
    ref = _oracle_lambdas(M, bKa * Kel1, K)
    scale = max(abs(r) for r in ref)
    for g, r in zip(got, ref):
        assert abs(g - r) <= 1e-8 * scale


def test_equaldof_dependent_damper():
    """MP constraint: a dashpot attached to an equalDOF dependent of node 3
    must damp exactly like one on node 3 itself;
    pins that the (default Transformation) handler distributes eigenvectors
    to dependent nodes, so the assembler projects with their physical shape.
    ('constraints Plain' with MP constraints would write zero rows there and
    lose the damper, as documented in LadrunoDampingAssembler.h.)"""
    c = 1.2
    M, K, C = _chain2dof(dashpot=False)
    ops.node(5, 2.0)
    ops.equalDOF(3, 5, 1)
    ops.node(4, 2.0)
    ops.fix(4, 1)
    ops.uniaxialMaterial("Viscous", 3, c, 1.0)
    ops.element("zeroLength", 3, 4, 5, "-mat", 3, "-dir", 1)
    ops.eigen("-fullGenLapack", 2)
    C[1, 1] = c  # physically identical to the dashpot on the retained node
    got = _domain_lambdas(ops.complexEigen())
    ref = _oracle_lambdas(M, C, K)
    scale = max(abs(r) for r in ref)
    for g, r in zip(got, ref):
        assert abs(g - r) <= 1e-8 * scale


# ------------------------------------------- base isolation, mode-shape output
def test_base_isolated_stick(tmp_path):
    """A stiff superstructure mass on a
    soft, heavily-damped isolation layer. The isolation-mode zeta/omegaD
    from complexEigen must match (a) the exact full-order state-space
    eigenvalues at 1e-8 and (b) a logarithmic-decrement fit of the actual
    Newmark free-decay time history (the same getDamp() C the transient
    feels) to a few percent."""
    m_iso, m_sup = 1.0, 1.0
    k_iso, k_sup = 4.0, 400.0
    c_iso = 0.6
    M, K, C = _chain2dof(k1=k_iso, k2=k_sup, m2=m_iso, m3=m_sup,
                         c=0.0, dashpot=False)
    # dashpot at the ISOLATION level (node 2), not the superstructure
    ops.node(4, 1.0)
    ops.fix(4, 1)
    ops.uniaxialMaterial("Viscous", 3, c_iso, 1.0)
    ops.element("zeroLength", 3, 4, 2, "-mat", 3, "-dir", 1)
    C[0, 0] = c_iso
    ops.eigen("-fullGenLapack", 2)

    # (a) exact oracle
    modes = modes_table(ops.complexEigen())
    ref = _oracle_lambdas(M, C, K)
    scale = max(abs(r) for r in ref)
    for m, r in zip(modes, ref):
        assert abs(m["lambda"] - r) <= 1e-8 * scale
    iso = modes[0]  # isolation mode = lowest omega0
    assert iso["zeta"] > 0.05  # heavily damped, the mode the design hinges on

    # (b) physical decay: release from a static push at the isolation level
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    ops.analyze(1)
    ops.setTime(0.0)
    ops.remove("loadPattern", 1)
    ops.wipeAnalysis()
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    T_iso = 2 * math.pi / iso["omegaD"]
    dt = T_iso / 200.0
    hist = []
    for _ in range(int(6 * T_iso / dt)):
        ops.analyze(1, dt)
        hist.append(ops.nodeDisp(2, 1))
    u = np.asarray(hist)
    # peak picking -> log decrement + damped period
    pk = [i for i in range(1, len(u) - 1)
          if u[i] > u[i - 1] and u[i] >= u[i + 1] and u[i] > 0]
    assert len(pk) >= 3
    ratios = [u[pk[i]] / u[pk[i + 1]] for i in range(len(pk) - 1)]
    delta = float(np.mean(np.log(ratios)))
    zeta_th = delta / math.sqrt(4 * math.pi**2 + delta**2)
    omega_d_th = 2 * math.pi / (float(np.mean(np.diff(pk))) * dt)
    # the free decay mixes both modes slightly: 5 % band
    assert zeta_th == pytest.approx(iso["zeta"], rel=0.05)
    assert omega_d_th == pytest.approx(iso["omegaD"], rel=0.05)


def test_recorder_roundtrip(tmp_path):
    """recorder Node ... complexEigenRe1/Im1 writes exactly the pushed
    psi_1, and the non-classical chain has a genuinely PHASED shape
    (nonzero Im after phase normalization)."""
    _chain2dof()
    ops.eigen("-fullGenLapack", 2)
    f_re = str(tmp_path / "re1.out")
    f_im = str(tmp_path / "im1.out")
    ops.recorder("Node", "-file", f_re, "-precision", 15, "-node", 2, 3,
                 "-dof", 1, "complexEigenRe1")
    ops.recorder("Node", "-file", f_im, "-precision", 15, "-node", 2, 3,
                 "-dof", 1, "complexEigenIm1")
    ops.complexEigen()
    ops.record()
    ops.remove("recorders")
    re = np.loadtxt(f_re, ndmin=2)[-1]
    im = np.loadtxt(f_im, ndmin=2)[-1]
    assert re.shape == (2,) and im.shape == (2,)
    # phase normalization: the largest-|psi| entry is real-positive
    mags = np.hypot(re, im)
    jmax = int(np.argmax(mags))
    assert re[jmax] > 0 and abs(im[jmax]) <= 1e-10 * mags[jmax]
    # non-classical fingerprint: the OTHER entry carries genuine phase
    jo = 1 - jmax
    assert abs(im[jo]) > 1e-6 * mags[jmax]


def test_classical_zero_phase(tmp_path):
    """Classical degeneracy: with pure global Rayleigh the complex shapes
    collapse to the REAL eigenvectors: Im ~ 0 and Re equals
    nodeEigenvector up to the documented phase convention."""
    _cantilever(n_eigen=3)
    ops.rayleigh(0.2, 0.001, 0.0, 0.0)
    f_re = str(tmp_path / "re2.out")
    f_im = str(tmp_path / "im2.out")
    nodes = [2, 3, 4, 5]
    ops.recorder("Node", "-file", f_re, "-precision", 15, "-node", *nodes,
                 "-dof", 1, "complexEigenRe2")
    ops.recorder("Node", "-file", f_im, "-precision", 15, "-node", *nodes,
                 "-dof", 1, "complexEigenIm2")
    ops.complexEigen()
    ops.record()
    ops.remove("recorders")
    re = np.loadtxt(f_re, ndmin=2)[-1]
    im = np.loadtxt(f_im, ndmin=2)[-1]
    phi = np.array([ops.nodeEigenvector(n, 2)[0] for n in nodes])
    scale = float(np.max(np.abs(phi)))
    assert np.max(np.abs(im)) <= 1e-8 * scale          # zero phase spread
    # Re == phi up to one common factor (mode-scale convention)
    j = int(np.argmax(np.abs(phi)))
    np.testing.assert_allclose(re / re[j], phi / phi[j], atol=1e-8)


def test_guards():
    # no prior eigen result in the domain
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    with pytest.raises(Exception):
        ops.complexEigen()
    # betaKinit / betaKcomm are not covered by the CLOSED FORM -> refused
    # there (the default assembled path handles them: test_betaK0_assembled)
    _cantilever()
    ops.rayleigh(0.1, 0.0, 0.001, 0.0)
    with pytest.raises(Exception):
        ops.complexEigen("-closedForm")
    _cantilever()
    ops.rayleigh(0.1, 0.0, 0.0, 0.001)
    with pytest.raises(Exception):
        ops.complexEigen("-closedForm")
    # ...and setting them back to zero un-refuses (domain copy tracks the
    # LAST rayleigh call, matching what the elements actually carry)
    ops.rayleigh(0.1, 0.002, 0.0, 0.0)
    assert len(modes_table(ops.complexEigen("-closedForm"))) == 4
