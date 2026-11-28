"""State-operation tests of nDMaterial TimeVarying (fix/timevarying-state).

They need the NDTest state commands (TrySetStrain, RevertToLastCommit, Copy,
SendRecv, GetInitialTangent), which ship with TimeVaryingEC; without them the
module is skipped. Reference material: perfect J2 (K, G in MPa).
"""
import math

import pytest

ops = pytest.importorskip("opensees")

K0, G0, FY = 175000.0, 80769.230769230766, 355.0
E0 = 9 * K0 * G0 / (3 * K0 + G0)
T = [0.0, 1.0, 2.0]
E = [E0, 0.7 * E0, 0.5 * E0]
K = [K0, 0.8 * K0, 0.6 * K0]
A = [1.0, 1.5, 2.5]
REF, TV = 1, 2


def _has_state_commands():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("J2Plasticity", 99, K0, G0, FY, FY, 0.0, 0.0, 0.0)
    try:
        # a no-op on the base class; an unknown NDTest command raises
        return ops.NDTest("RevertToLastCommit", 99) == 0
    except Exception:
        return False


pytestmark = pytest.mark.skipif(not _has_state_commands(),
                                reason="NDTest state commands (TimeVaryingEC branch) not available")


def build():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("J2Plasticity", REF, K0, G0, FY, FY, 0.0, 0.0, 0.0)
    ops.nDMaterial("TimeVarying", TV, REF, len(T), *T, *E, *K, *A)
    return TV


def drv(x):
    ops.NDTest("UpdateDoubleParameter", TV, 4002, x)


def trial(e):
    return ops.NDTest("TrySetStrain", TV, *e)


def commit(tag=TV):
    return ops.NDTest("TryCommit", tag)


def stress(tag=TV):
    return list(ops.NDTest("GetStress", tag))


def tangent(tag=TV):
    return list(ops.NDTest("GetTangentStiffness", tag))


PATH = [(0.0, [5e-4, -1e-4, 0, 1.5e-3, 0, 0]),
        (0.0, [1e-3, -2e-4, 0, 3e-3, 5e-4, 0]),      # plastic
        (0.5, [1e-3, -2e-4, 0, 3e-3, 5e-4, 0]),
        (1.0, [1.5e-3, 0, -5e-4, 4.5e-3, 5e-4, -1e-3]),
        (1.5, [5e-4, 0, -5e-4, 1e-3, 0, -1e-3]),
        (2.0, [-1e-3, 5e-4, 0, -2.5e-3, 0, 0])]


def run(path, tag=TV):
    out = []
    for x, e in path:
        ops.NDTest("UpdateDoubleParameter", tag, 4002, x)
        assert ops.NDTest("TrySetStrain", tag, *e) == 0 and commit(tag) == 0
        out.append(stress(tag))
    return out


def test_initial_tangent_before_first_trial():
    build()
    drv(1.0)
    C = ops.NDTest("GetInitialTangent", TV)
    assert all(math.isfinite(v) for v in C)
    G = 3 * K[1] * E[1] / (9 * K[1] - E[1])
    assert abs(C[21] - G) <= 1e-9 * G          # C(3,3) = G(xi)


def test_revert_is_immediate():
    build()
    run(PATH[:3])
    committed = (stress(), tangent())
    drv(1.5)
    assert trial([3e-3, 0, 0, 9e-3, 0, 0]) == 0
    ops.NDTest("RevertToLastCommit", TV)
    assert (stress(), tangent()) == committed


def test_commit_after_revert_without_trial():
    build()
    run(PATH[:3])
    drv(1.5)
    assert trial([3e-3, 0, 0, 9e-3, 0, 0]) == 0
    ops.NDTest("RevertToLastCommit", TV)
    assert commit() == 0
    rest = run(PATH[3:])
    build()
    run(PATH[:3])
    assert run(PATH[3:]) == rest


def test_rejected_trial_leaves_state():
    build()
    run(PATH[:2])
    drv(0.5)
    assert trial(PATH[2][1]) == 0
    before = (stress(), list(ops.NDTest("GetStrain", TV)))
    assert trial([float("nan"), 0, 0, 0, 0, 0]) != 0
    assert (stress(), list(ops.NDTest("GetStrain", TV))) == before


def test_send_recv_round_trip(tmp_path):
    build()
    run(PATH[:3])
    before = (stress(), tangent())
    assert ops.NDTest("SendRecv", TV, 600, str(tmp_path / "tv")) == 0
    assert (stress(600), tangent(600)) == before
    assert run(PATH[3:], 600) == run(PATH[3:], TV)


def test_copy_independent():
    build()
    run(PATH[:3])
    ops.NDTest("Copy", TV, 500)
    run([(2.0, [-3e-3, 0, 0, -9e-3, 0, 0])], 500)
    rest = run(PATH[3:])
    build()
    run(PATH[:3])
    assert run(PATH[3:]) == rest
