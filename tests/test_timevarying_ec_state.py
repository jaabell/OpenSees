"""Native state-operation tests of nDMaterial TimeVaryingEC (experimental).

Phase B of the EC plan: every check runs on the native material through
OpenSeesPy (NDTest), not on a Python prototype. Equality is bitwise unless a
tolerance is stated. Physical closed forms used:
  * isotropic elasticity: sigma = C(K, G) eps_m exactly, whatever the history;
  * shear ordering: Voigt xx yy zz xy yz zx, engineering shear strain in,
    tensor shear stress out.
"""
import math

import pytest

ops = pytest.importorskip("opensees")

# reference material (MPa)
K0, G0 = 175000.0, 80769.230769230766
FY = 355.0
J2_TOL = 1e-8   # J2Plasticity local tolerance, 1e-8 sigma_0
REF_ELASTIC, REF_J2, EC = 1, 2, 10


def C_iso(K, G):
    lam = K - 2.0 * G / 3.0
    C = [[0.0] * 6 for _ in range(6)]
    for i in range(3):
        for j in range(3):
            C[i][j] = lam
        C[i][i] = lam + 2.0 * G
        C[i + 3][i + 3] = G
    return C


def matvec(C, v):
    return [sum(C[i][j] * v[j] for j in range(6)) for i in range(6)]


def E_of(K, G):
    return 9.0 * K * G / (3.0 * K + G)


# law: driver 0..3, stiffness and strength change (A = strength map)
LAW_T = [0.0, 1.0, 2.0, 3.0]
LAW_K = [K0, 0.8 * K0, 0.6 * K0, 0.5 * K0]
LAW_G = [G0, 0.7 * G0, 0.5 * G0, 0.4 * G0]
LAW_E = [E_of(k, g) for k, g in zip(LAW_K, LAW_G)]
LAW_A = [1.0, 1.5, 2.5, 4.0]


def build(ref="J2", law=None, tag=EC):
    t, E, K, A = law or (LAW_T, LAW_E, LAW_K, LAW_A)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    if ref == "J2":
        # perfect J2: sigInf = sig0, delta = H = eta = 0
        ops.nDMaterial("J2Plasticity", REF_J2, K0, G0, FY, FY, 0.0, 0.0, 0.0)
        r = REF_J2
    else:
        ops.nDMaterial("ElasticIsotropic", REF_ELASTIC, E_of(K0, G0), (3 * K0 - 2 * G0) / (2 * (3 * K0 + G0)))
        r = REF_ELASTIC
    ops.nDMaterial("TimeVaryingEC", tag, r, len(t), *t, *E, *K, *A)
    return tag


def trial(tag, eps):
    return ops.NDTest("TrySetStrain", tag, *eps)


def commit(tag):
    return ops.NDTest("TryCommit", tag)


def stress(tag):
    return list(ops.NDTest("GetStress", tag))


def strain(tag):
    return list(ops.NDTest("GetStrain", tag))


def tangent(tag):
    return list(ops.NDTest("GetTangentStiffness", tag))


def init_tangent(tag):
    return list(ops.NDTest("GetInitialTangent", tag))


def response(tag, name):
    return list(ops.NDTest("GetResponse", tag, name))


def driver(tag, x):
    ops.NDTest("UpdateDoubleParameter", tag, 4002, x)


def snapshot(tag):
    return (stress(tag), tangent(tag), strain(tag), response(tag, "properties"),
            response(tag, "projectedStrain"), response(tag, "projectedStress"))


def props_at(x):
    """Law evaluated exactly as the native class does (E, K, A interpolated)."""
    import bisect
    i = bisect.bisect_left(LAW_T, x)
    if i == 0:
        E, K, A = LAW_E[0], LAW_K[0], LAW_A[0]
    elif i < len(LAW_T):
        a = (x - LAW_T[i - 1]) / (LAW_T[i] - LAW_T[i - 1])
        E = (1 - a) * LAW_E[i - 1] + a * LAW_E[i]
        K = (1 - a) * LAW_K[i - 1] + a * LAW_K[i]
        A = (1 - a) * LAW_A[i - 1] + a * LAW_A[i]
    else:
        E, K, A = LAW_E[-1], LAW_K[-1], LAW_A[-1]
    G = 3 * K * E / (9 * K - E)
    return K, G, A


# a loading path with properties changing at held and moving strain
PATH = [(x, [5.0 * v for v in e]) for x, e in [
    (0.0, [1e-4, -2e-5, 0, 3e-4, 0, 0]),
    (0.0, [2e-4, -4e-5, 0, 6e-4, 1e-4, 0]),     # plastic
    (0.5, [2e-4, -4e-5, 0, 6e-4, 1e-4, 0]),     # property change at held strain
    (1.0, [2e-4, -4e-5, 0, 6e-4, 1e-4, 0]),
    (1.5, [3e-4, 0, -1e-4, 9e-4, 1e-4, -2e-4]),
    (2.0, [3e-4, 0, -1e-4, 9e-4, 1e-4, -2e-4]),
    (2.5, [1e-4, 0, -1e-4, 2e-4, 0, -2e-4]),    # unloading
    (3.0, [-2e-4, 1e-4, 0, -5e-4, 0, 0]),       # reversal
]]


def q_of(s):
    p = (s[0] + s[1] + s[2]) / 3.0
    d = [s[0] - p, s[1] - p, s[2] - p]
    return math.sqrt(1.5 * (d[0] ** 2 + d[1] ** 2 + d[2] ** 2
                            + 2 * (s[3] ** 2 + s[4] ** 2 + s[5] ** 2)))


def run_path(tag, path):
    out = []
    for x, e in path:
        driver(tag, x)
        assert trial(tag, e) == 0
        assert commit(tag) == 0
        out.append(stress(tag))
    return out


# ---------------------------------------------------------------- conventions

def test_pure_shear_ordering():
    """Each engineering shear strain lands on its own tensor shear stress."""
    tag = build("elastic")
    for comp in (3, 4, 5):
        ops.NDTest("RevertToStart", tag)
        e = [0.0] * 6
        e[comp] = 2e-4
        assert trial(tag, e) == 0
        s = stress(tag)
        for k in range(6):
            expect = G0 * 2e-4 if k == comp else 0.0
            assert abs(s[k] - expect) <= 1e-9 * G0 * 2e-4, (comp, k, s)


# ------------------------------------------------------------------- physics

def test_elastic_closed_form_along_path():
    """Elastic reference: sigma = C(xi) eps_m at every committed step."""
    tag = build("elastic")
    for x, e in PATH:
        driver(tag, x)
        assert trial(tag, e) == 0 and commit(tag) == 0
        K, G, A = props_at(x)
        ref = matvec(C_iso(K, G), e)
        s = stress(tag)
        scale = max(abs(v) for v in ref)
        assert max(abs(a - b) for a, b in zip(s, ref)) <= 1e-12 * scale, (x, s, ref)


def test_sigma_ref_invariant_elastic():
    """sigma_ref = A sigma from a compatible (virgin) origin."""
    tag = build("elastic")
    for x, e in PATH:
        driver(tag, x)
        assert trial(tag, e) == 0 and commit(tag) == 0
        s, sp = stress(tag), response(tag, "projectedStress")
        A = props_at(x)[2]
        scale = max(abs(v) for v in sp)
        assert max(abs(a - A * b) for a, b in zip(sp, s)) <= 1e-12 * scale


def test_j2_path_yields_and_respects_current_strength():
    """The J2 path is plastic, and every committed state keeps q at or below
    the current strength fy / A(xi): a strength drop at held strain leaves no
    stress above it (the documented TimeVarying limitation, 153.7 MPa in A23).

    Exact checks: a bare J2 copy driven with the wrapper's projected strains
    gives the wrapper's projected stresses bitwise (the wrapper adds nothing
    to the return), and sigma_ref = A sigma (1e-12). Coarse physical check:
    q <= fy/A (1 + 1e-6). J2Plasticity shrinks the converged plastic
    multiplier by (1 - 1e-8), so it ends outside its surface by ~1e-8 of the
    plastic increment 2 G gamma, which here is up to tens of times q; a
    fixed 1e-8 relative bound on q is therefore not a property of J2."""
    tag = build("J2")
    bare = 90
    ops.NDTest("Copy", REF_J2, bare, "ThreeDimensional")
    yielded = 0
    for x, e in PATH:
        driver(tag, x)
        assert trial(tag, e) == 0 and commit(tag) == 0
        A = props_at(x)[2]
        sig, sref = stress(tag), response(tag, "projectedStress")
        assert ops.NDTest("TrySetStrain", bare, *response(tag, "projectedStrain")) == 0
        assert commit(bare) == 0
        assert stress(bare) == sref, x
        scale = max(abs(v) for v in sref)
        assert max(abs(a - A * b) for a, b in zip(sref, sig)) <= 1e-12 * scale, x
        q, qref = q_of(sig), q_of(sref)
        assert q <= FY / A * (1 + 1e-6), (x, q, FY / A)
        if qref >= FY * (1 - 1e-9):
            yielded += 1
    assert yielded >= 3


# --------------------------------------------------------------- transactions

def test_repeated_trial_is_identical():
    tag = build("J2")
    run_path(tag, PATH[:4])
    driver(tag, 2.0)
    e = [5e-4, 0, 0, 2e-3, 0, 0]          # plastic trial
    assert trial(tag, e) == 0
    first = snapshot(tag)
    assert trial(tag, [1e-4] * 6) == 0     # another trial in between
    assert trial(tag, e) == 0
    assert snapshot(tag) == first


def test_query_order_does_not_change_state():
    tag = build("J2")
    run_path(tag, PATH[:3])
    driver(tag, 1.5)
    assert trial(tag, [5e-4, 0, 0, 2e-3, 0, 0]) == 0
    a = tangent(tag); b = stress(tag); c = tangent(tag); d = stress(tag)
    assert a == c and b == d
    commit(tag)
    s1 = stress(tag)
    tag2 = build("J2")
    run_path(tag2, PATH[:3])
    driver(tag2, 1.5)
    assert trial(tag2, [5e-4, 0, 0, 2e-3, 0, 0]) == 0
    commit(tag2)
    assert stress(tag2) == s1


def test_invalid_input_leaves_state_untouched():
    tag = build("J2")
    run_path(tag, PATH[:4])
    driver(tag, 1.5)
    assert trial(tag, [3e-4, 0, 0, 1e-3, 0, 0]) == 0
    before = snapshot(tag)
    assert trial(tag, [float("nan"), 0, 0, 0, 0, 0]) != 0
    assert snapshot(tag) == before
    assert trial(tag, [1e300, 1e300, 0, 0, 0, 0]) != 0  # non-finite reference response
    assert snapshot(tag) == before
    # the kept trial still commits to the same state as an uninterrupted run
    assert commit(tag) == 0
    s = stress(tag)
    tag2 = build("J2")
    run_path(tag2, PATH[:4])
    driver(tag2, 1.5)
    assert trial(tag2, [3e-4, 0, 0, 1e-3, 0, 0]) == 0 and commit(tag2) == 0
    assert stress(tag2) == s


def test_revert_is_immediate():
    """Immediate queries after a revert, with NO trial in between."""
    tag = build("J2")
    run_path(tag, PATH[:5])
    committed = snapshot(tag)
    driver(tag, 2.5)
    assert trial(tag, [9e-4, 0, 0, 4e-3, 0, 0]) == 0   # plastic trial at other properties
    assert snapshot(tag) != committed
    assert ops.NDTest("RevertToLastCommit", tag) == 0
    assert snapshot(tag) == committed


def test_reference_j2_revert_is_not_immediate():
    """Documents why the wrapper never asks the reference after a revert."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("J2Plasticity", REF_J2, K0, G0, FY, FY, 0.0, 0.0, 0.0)
    m = 77
    ops.NDTest("Copy", REF_J2, m, "ThreeDimensional")
    assert ops.NDTest("TrySetStrain", m, 1e-4, 0, 0, 0, 0, 0) == 0
    assert commit(m) == 0
    s_n = stress(m)
    assert ops.NDTest("TrySetStrain", m, 5e-3, 0, 0, 0, 0, 0) == 0
    ops.NDTest("RevertToLastCommit", m)
    assert stress(m) != s_n          # stale trial stress: J2 revert is a no-op


def test_commit_after_revert_without_trial_continues_exactly():
    tag = build("J2")
    run_path(tag, PATH[:5])
    driver(tag, 2.5)
    assert trial(tag, [9e-4, 0, 0, 4e-3, 0, 0]) == 0
    ops.NDTest("RevertToLastCommit", tag)
    assert commit(tag) == 0              # nothing new: must not commit the stale reference trial
    rest = run_path(tag, PATH[5:])
    tag2 = build("J2")
    run_path(tag2, PATH[:5])
    assert run_path(tag2, PATH[5:]) == rest


def test_revert_restores_driver_and_eigenstrain():
    tag = build("J2")
    run_path(tag, PATH[:3])
    committed = snapshot(tag)
    driver(tag, 2.9)
    ops.NDTest("UpdateDoubleParameter", tag, 4001, 1e-4)
    assert trial(tag, PATH[2][1]) == 0
    ops.NDTest("RevertToLastCommit", tag)
    assert snapshot(tag) == committed
    # a trial with no new inputs maps to the committed projected strain
    # exactly; the reference J2 then re-returns a point on its yield surface,
    # which moves it within its own tolerance (see test_reference_j2_retrial)
    assert trial(tag, PATH[2][1]) == 0
    assert response(tag, "projectedStrain") == committed[4]
    scale = max(abs(v) for v in committed[0])
    assert max(abs(a - b) for a, b in zip(stress(tag), committed[0])) <= J2_TOL * scale


def test_revert_to_start_and_initial_tangent():
    tag = build("J2")
    driver(tag, 2.0)                       # non-identity map before any trial
    K, G, _ = props_at(2.0)
    expect = [v for row in C_iso(K, G) for v in row]
    got = init_tangent(tag)
    assert max(abs(a - b) for a, b in zip(got, expect)) <= 1e-12 * K
    # the driver update already evaluated a (zero-strain) trial, so the
    # tangent is Cref S / A = C up to roundoff
    assert max(abs(a - b) for a, b in zip(tangent(tag), got)) <= 1e-14 * K
    run_path(tag, PATH)
    assert ops.NDTest("RevertToStart", tag) == 0
    assert stress(tag) == [0.0] * 6 and strain(tag) == [0.0] * 6
    # after reset the domain-time driver (0) is back in charge: C(xi0)
    K, G, _ = props_at(0.0)
    expect = [v for row in C_iso(K, G) for v in row]
    got = init_tangent(tag)
    assert max(abs(a - b) for a, b in zip(got, expect)) <= 1e-12 * K
    # and the full path reproduces a fresh material bitwise
    first = run_path(tag, PATH)
    tag2 = build("J2")
    assert run_path(tag2, PATH) == first


def test_copy_is_deep_and_diverges_without_contamination():
    tag = build("J2")
    run_path(tag, PATH[:4])
    driver(tag, 1.5)
    assert trial(tag, PATH[4][1]) == 0      # pending trial is copied too
    ops.NDTest("Copy", tag, 500)
    assert snapshot(500) == snapshot(tag)
    # the copy follows a different path; the original is not affected
    commit(500)
    other = [(3.0, [-8e-4, 0, 0, -3e-3, 0, 0])]
    run_path(500, other)
    commit(tag)
    rest = run_path(tag, PATH[5:])
    tag2 = build("J2")
    run_path(tag2, PATH[:5])
    assert run_path(tag2, PATH[5:]) == rest


@pytest.mark.parametrize("pending", [False, True])
def test_send_recv_round_trip(tmp_path, pending):
    tag = build("J2")
    run_path(tag, PATH[:4])
    if pending:
        driver(tag, 1.5)
        assert trial(tag, PATH[4][1]) == 0
    before = snapshot(tag)
    assert ops.NDTest("SendRecv", tag, 600, str(tmp_path / "tvec")) == 0
    assert snapshot(600) == before             # immediate query, no re-evaluation
    if pending:
        assert commit(600) == 0 and commit(tag) == 0
    a = run_path(600, PATH[5:])
    b = run_path(tag, PATH[5:])
    assert a == b


def test_reference_j2_retrial():
    """Bare J2: a trial at the committed strain of a plastic state is not
    bit-identical to the committed stress (re-return within 1e-8 sigma_0)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("J2Plasticity", REF_J2, K0, G0, FY, FY, 0.0, 0.0, 0.0)
    m = 82
    ops.NDTest("Copy", REF_J2, m, "ThreeDimensional")
    e = [1e-3, -2e-4, 0, 3e-3, 5e-4, 0]
    assert ops.NDTest("TrySetStrain", m, *e) == 0 and commit(m) == 0
    s1 = stress(m)
    ops.NDTest("TrySetStrain", m, *e)
    s2 = stress(m)
    d = max(abs(a - b) for a, b in zip(s1, s2)) / max(abs(v) for v in s1)
    assert 0 < d <= J2_TOL


def test_reference_j2_contract():
    """Minimum contract of the reference: trials from committed, deep copy."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("J2Plasticity", REF_J2, K0, G0, FY, FY, 0.0, 0.0, 0.0)
    m = 80
    ops.NDTest("Copy", REF_J2, m, "ThreeDimensional")
    assert ops.NDTest("TrySetStrain", m, 3e-3, 0, 0, 0, 0, 0) == 0 and commit(m) == 0
    e = [4e-3, 0, 0, 1e-3, 0, 0]
    ops.NDTest("TrySetStrain", m, *e); s1 = stress(m)
    ops.NDTest("TrySetStrain", m, 9e-3, 0, 0, 0, 0, 0)
    ops.NDTest("TrySetStrain", m, *e)
    assert stress(m) == s1                   # trial integrates from committed state
    ops.NDTest("Copy", m, 81)
    ops.NDTest("TrySetStrain", 81, -5e-3, 0, 0, 0, 0, 0); commit(81)
    ops.NDTest("TrySetStrain", m, *e)
    assert stress(m) == s1                   # copy is independent


def test_strain_semantics_with_eigenstrain():
    """getStrain() is the total input strain; mechanical = total - eigen."""
    tag = build("elastic")
    ops.NDTest("UpdateDoubleParameter", tag, 4001, 2e-4)
    e = [5e-4, 1e-4, -1e-4, 2e-4, 0, 1e-4]
    assert trial(tag, e) == 0 and commit(tag) == 0
    assert strain(tag) == e
    assert response(tag, "eigenStrain") == [2e-4, 2e-4, 2e-4, 0, 0, 0]
    em = response(tag, "mechanicalStrain")
    assert em == [a - b for a, b in zip(e, [2e-4, 2e-4, 2e-4, 0, 0, 0])]
    K, G, _ = props_at(0.0)
    ref = matvec(C_iso(K, G), em)
    assert max(abs(a - b) for a, b in zip(stress(tag), ref)) <= 1e-12 * max(map(abs, ref))
