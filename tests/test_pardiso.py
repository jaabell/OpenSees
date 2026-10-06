"""Tests for `system Pardiso` (MKL PARDISO, PARDISOGenLinSOE/Solver).

The whole module is skipped when the build has no `system Pardiso` (built
without MKL, or a non-Windows build without -DMKL_PARDISO_LINUX=ON).

Threaded PARDISO is not bit-reproducible run to run, so results are compared
with a relative tolerance, never bitwise.
"""
import os
import re

import pytest

# MKL reads this at its first initialization.
os.environ.setdefault("MKL_NUM_THREADS", "1")

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E, NU = 200000.0, 0.3
KMOD = E / (3.0 * (1.0 - 2.0 * NU))
GMOD = E / (2.0 * (1.0 + NU))
S0, HISO = 250.0, 2000.0
NX = 4          # 4x4x4 bricks, 300 free DOF
LEL = 100.0
NSTEPS = 5
TOL = 1.0e-8


def _pardiso_available():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    try:
        ops.system("Pardiso")
        ok = True
    except Exception:
        ok = False
    ops.wipe()
    return ok


pytestmark = pytest.mark.skipif(not _pardiso_available(),
                                reason="this build has no system Pardiso")


def _nid(i, j, k):
    return 1 + i + (NX + 1) * (j + (NX + 1) * k)


def _run(system_args, nonlinear=True, algorithm="Newton"):
    """Cantilever block of stdBrick elements, fixed at the base and loaded
    laterally and vertically at the top. Returns the top-centre (ux, uz).

    nonlinear=True uses J2Plasticity; the load takes the base past yield, so
    the tangent changes between Newton iterations."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    if nonlinear:
        ops.nDMaterial("J2Plasticity", 1, KMOD, GMOD, S0, S0, 0.0, HISO)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    for k in range(NX + 1):
        for j in range(NX + 1):
            for i in range(NX + 1):
                ops.node(_nid(i, j, k), i * LEL, j * LEL, k * LEL)
    for j in range(NX + 1):
        for i in range(NX + 1):
            ops.fix(_nid(i, j, 0), 1, 1, 1)
    tag = 1
    for k in range(NX):
        for j in range(NX):
            for i in range(NX):
                ops.element(
                    "stdBrick", tag,
                    _nid(i, j, k), _nid(i + 1, j, k),
                    _nid(i + 1, j + 1, k), _nid(i, j + 1, k),
                    _nid(i, j, k + 1), _nid(i + 1, j, k + 1),
                    _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1),
                    1)
                tag += 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(NX + 1):
        for i in range(NX + 1):
            ops.load(_nid(i, j, NX), 4.0e5, 0.0, -4.0e5)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", TOL, 50)
    ops.algorithm(algorithm)
    ops.integrator("LoadControl", 1.0 / NSTEPS)
    ops.analysis("Static")
    for s in range(NSTEPS):
        assert ops.analyze(1) == 0, (
            f"analyze failed at step {s + 1} with system {system_args}")
    tip = _nid(NX // 2, NX // 2, NX)
    out = (ops.nodeDisp(tip, 1), ops.nodeDisp(tip, 3))
    ops.wipe()
    return out


def _close(got, ref, rtol):
    return all(abs(g - r) <= rtol * max(abs(r), 1e-30) for g, r in zip(got, ref))


def _text(capfd):
    cap = capfd.readouterr()
    return cap.out + cap.err


@pytest.mark.parametrize("matrix_type", [0, 1, 2])
@pytest.mark.parametrize("nonlinear", [False, True])
def test_matches_umfpack(matrix_type, nonlinear):
    """Every -matrixType agrees with UmfPack on a linear and a plastic model.
    The J2 tangent is symmetric and stays positive definite, so 1 and 2 are
    both valid here."""
    ref = _run(["UmfPack"], nonlinear)
    got = _run(["Pardiso", "-matrixType", matrix_type], nonlinear)
    assert _close(got, ref, 1e-9), f"Pardiso {got} vs UmfPack {ref}"


@pytest.mark.parametrize("matrix_type", [0, 2])
def test_factorization_reuse_modified_newton(matrix_type):
    """ModifiedNewton keeps the tangent, so later solves in a step reuse the
    stored factorization (phase 33 only). Stale factors would show up as a
    different converged path from UmfPack under the same algorithm."""
    ref = _run(["UmfPack"], algorithm="ModifiedNewton")
    got = _run(["Pardiso", "-matrixType", matrix_type],
               algorithm="ModifiedNewton")
    assert _close(got, ref, 1e-9), f"Pardiso {got} vs UmfPack {ref}"


def test_krylov_matches_direct(capfd):
    """-krylov reuses the previous factors as a CGS preconditioner for a
    changed tangent; the result matches the direct solve to the CGS
    tolerance."""
    ref = _run(["Pardiso"])
    got = _run(["Pardiso", "-krylov", 6])
    text = _text(capfd)
    assert "CGS solve (phase 23)" not in text
    assert _close(got, ref, 1e-6), f"-krylov 6 {got} vs direct {ref}"


def test_krylov_refused_for_symmetric_indefinite(capfd):
    """PARDISO documents no CGS mode for mtype -2: warn and solve directly."""
    ref = _run(["Pardiso", "-matrixType", 2])
    got = _run(["Pardiso", "-matrixType", 2, "-krylov", 6])
    text = _text(capfd)
    assert "-krylov is not available" in text
    assert _close(got, ref, 1e-9)


def test_stats_printed_per_factorization(capfd):
    _run(["Pardiso", "-matrixType", 2, "-stats"])
    text = _text(capfd)
    assert text.count("PARDISO stats:") >= NSTEPS, text
    assert "matrixType=-2" in text
    assert "fact memory KB iparm(17)" in text


def test_bad_option_values_degrade(capfd):
    """A string -matrixType and an unknown flag warn and fall back to the
    unsymmetric default; the analysis still runs on PARDISO."""
    ref = _run(["Pardiso"])
    got = _run(["Pardiso", "-matrixType", "2", "-bogus"])
    text = _text(capfd)
    assert "failed to get -matrixType" in text
    assert "unknown option -bogus" in text
    assert _close(got, ref, 1e-9)


# ---- half-storage asymmetry check -------------------------------------------
# A non-associated Drucker-Prager block (rho_bar != rho) is symmetric while
# elastic and turns unsymmetric after yield, long after the first sampling
# window, so the warning must come from a later sampling window.

_DP_K, _DP_G = 27777.78, 9259.26


def _dp_ramp(system_args, elastic=False, nsteps=48):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    if elastic:
        e = 9.0 * _DP_K * _DP_G / (3.0 * _DP_K + _DP_G)
        nu = (3.0 * _DP_K - 2.0 * _DP_G) / (2.0 * (3.0 * _DP_K + _DP_G))
        ops.nDMaterial("ElasticIsotropic", 1, e, nu)
    else:
        # K, G, sigma_y, rho, rho_bar, Kinf, Ko, delta1, delta2, H, theta, density
        ops.nDMaterial("DruckerPrager", 1, _DP_K, _DP_G, 5.0, 0.398, 0.1,
                       0.0, 0.0, 0.0, 0.0, 1000.0, 1.0, 0.0)
    n = 2
    h = 1.0 / n

    def nid(i, j, k):
        return 1 + i + (n + 1) * (j + (n + 1) * k)

    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                ops.node(nid(i, j, k), i * h, j * h, k * h)
    for j in range(n + 1):
        for i in range(n + 1):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    tag = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                ops.element("stdBrick", tag,
                            nid(i, j, k), nid(i + 1, j, k),
                            nid(i + 1, j + 1, k), nid(i, j + 1, k),
                            nid(i, j, k + 1), nid(i + 1, j, k + 1),
                            nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1), 1)
                tag += 1
    ntop = (n + 1) * (n + 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(n + 1):
        for i in range(n + 1):
            ops.load(nid(i, j, n), 2.025 / ntop, 0.0, -3.0 / ntop)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1.0e-8, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    done = 0
    # Half storage of an unsymmetric tangent solves the reflected system, so
    # Newton may stop converging after yield; the warning is what is tested.
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        done += 1
    ops.wipe()
    return done


def test_late_asymmetry_is_reported_once(capfd):
    _dp_ramp(["Pardiso", "-matrixType", 2])
    text = _text(capfd).lower()
    assert text.count("is unsymmetric") == 1, text
    m = re.search(r"tangent assembly (\d+)", text)
    assert m is not None and int(m.group(1)) >= 64, text


def test_symmetric_tangent_does_not_warn(capfd):
    assert _dp_ramp(["Pardiso", "-matrixType", 1], elastic=True) == 48
    assert "is unsymmetric" not in _text(capfd).lower()


# ---- -deterministic ----------------------------------------------------------
# Last in the module: the MKL reproducibility mode is process-wide and stays
# on for the rest of the process.

def test_deterministic_matches_default(capfd):
    ref = _run(["Pardiso"])
    capfd.readouterr()
    got = _run(["Pardiso", "-deterministic"])
    text = _text(capfd)
    if "unknown option -deterministic" in text:
        pytest.skip("this build has no -deterministic option")
    # MKL may refuse the mode once it has been used in this process (on Linux
    # after any earlier solve); the notice is printed either way.
    assert "PARDISO deterministic mode" in text, text
    assert _close(got, ref, 1e-9)

