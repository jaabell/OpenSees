"""Corotational tangent of the ASD shell elements: geometric stiffness.

With the polar-aligned local frame the rotation gradient G was left at zero, so
the corotational tangent had no geometric stiffness while the residual, computed
on the exact frame, did. The tangent now carries both EICR geometric terms.

Model: 30 m cantilever strip with a tip mass and an axial load P held constant
after a static step. Euler load Pcr = pi^2 E I / (4 H^2). References:

  softening        out-of-plane sway, P = 0.5 Pcr: period ratio
                   T(P) / T(0) = 1 / sqrt(1 - P / Pcr) = sqrt(2).
  in-plane         sway in the plane of the strip, out-of-plane held, P = 0.3 Pcr:
                   period ratio 1 / sqrt(0.7). Exercises the in-plane spin row of G.
  newton           lateral step on the preloaded strip: a consistent tangent
                   converges in at most 8 iterations (33 without G), with the
                   second-order amplification 1 / (1 - P / Pcr) = 2. The
                   amplification is a control, it held before the change.
  zero load        control: without axial load the corotational period equals the
                   linear one.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

E, NU, T, B, H, M = 25.0e6, 0.2, 0.3, 1.0, 30.0, 100.0
PCR_OUT = math.pi ** 2 * E * (B * T ** 3 / 12.0) / (4.0 * H ** 2)
PCR_IN = math.pi ** 2 * E * (T * B ** 3 / 12.0) / (4.0 * H ** 2)

ELEMENTS = [("ASDShellQ4", False), ("ASDShellT3", True)]


def _strip(element, options, axial, tri, lateral_dof=2, nx=1, nz=10):
    """Period and lateral step of the strip after the axial load is held.

    The strip lies in the global XZ plane, so dof 2 is out of plane and dof 1 is
    in plane. Returns (period, Newton iterations of the step, lateral tip
    displacement).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    def tag(i, k):
        return 1 + i + k * (nx + 1)

    in_plane = lateral_dof == 1
    for k in range(nz + 1):
        for i in range(nx + 1):
            ops.node(tag(i, k), i * B / nx, 0.0, k * H / nz)
            if k == 0:
                ops.fix(tag(i, k), 1, 1, 1, 1, 1, 1)
            elif in_plane:
                ops.fix(tag(i, k), 0, 1, 0, 1, 0, 1)
    ops.section("ElasticMembranePlateSection", 1, E, NU, T, 0.0)
    e = 1
    for k in range(nz):
        for i in range(nx):
            a, b, c, d = tag(i, k), tag(i + 1, k), tag(i + 1, k + 1), tag(i, k + 1)
            if tri:
                ops.element(element, e, a, b, c, 1, *options)
                ops.element(element, e + 1, a, c, d, 1, *options)
                e += 2
            else:
                ops.element(element, e, a, b, c, d, 1, *options)
                e += 1
    top = [tag(i, nz) for i in range(nx + 1)]
    for n in top:
        m = M / len(top)
        ops.mass(n, m, 0.0 if in_plane else m, 0.0, 0.0, 0.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in top:
        ops.load(n, 0.0, 0.0, -axial / len(top), 0.0, 0.0, 0.0)
    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-10, 100)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)
    period = 2.0 * math.pi / math.sqrt(ops.eigen("-fullGenLapack", 1)[0])
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    force = [0.0] * 6
    force[lateral_dof - 1] = 1.0 / len(top)
    for n in top:
        ops.load(n, *force)
    assert ops.analyze(1) == 0
    return period, ops.testIter(), ops.nodeDisp(top[0], lateral_dof)


@pytest.mark.parametrize("element,tri", ELEMENTS)
def test_out_of_plane_softening_matches_euler(element, tri):
    try:
        t0, _, _ = _strip(element, ["-corotational"], 0.0, tri)
        tp, _, _ = _strip(element, ["-corotational"], 0.5 * PCR_OUT, tri)
    finally:
        ops.wipe()
    assert tp / t0 == pytest.approx(math.sqrt(2.0), rel=1.0e-2)


@pytest.mark.parametrize("element,tri", ELEMENTS)
def test_in_plane_softening_matches_euler(element, tri):
    try:
        t0, _, _ = _strip(element, ["-corotational"], 0.0, tri, lateral_dof=1, nx=4)
        tp, _, _ = _strip(element, ["-corotational"], 0.3 * PCR_IN, tri, lateral_dof=1, nx=4)
    finally:
        ops.wipe()
    assert tp / t0 == pytest.approx(1.0 / math.sqrt(0.7), rel=1.0e-2)


@pytest.mark.parametrize("element,tri", ELEMENTS)
def test_newton_on_a_p_delta_step(element, tri):
    try:
        _, _, u0 = _strip(element, ["-corotational"], 0.0, tri)
        _, iters, up = _strip(element, ["-corotational"], 0.5 * PCR_OUT, tri)
    finally:
        ops.wipe()
    assert up / u0 == pytest.approx(2.0, rel=2.0e-2)
    assert iters <= 8


@pytest.mark.parametrize("element,tri", ELEMENTS)
def test_zero_axial_load_control(element, tri):
    try:
        linear, _, _ = _strip(element, [], 0.0, tri)
        corotational, _, _ = _strip(element, ["-corotational"], 0.0, tri)
    finally:
        ops.wipe()
    assert corotational == pytest.approx(linear, rel=1.0e-7)
