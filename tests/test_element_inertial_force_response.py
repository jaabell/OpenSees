"""Element response "inertialForce" (Element::getResponse id 444444).

The response is the element force including inertia minus the Rayleigh
damping force minus the static resisting force, i.e. M * a. It is reached
through Element::setResponse by elements that fall back to the base class,
here dispBeamColumn (DispBeamColumn2d) with a lumped mass.

Model: 2D cantilever, node 1 fixed, node 2 free at L = 2.0, elastic section
E = 200, A = 1, I = 1, mass per unit length rho = 3.0. The lumped mass is
m = rho * L / 2 = 3.0 on each translational DOF and zero on the rotation.
One Newmark step (gamma = 0.5, beta = 0.25, dt = 0.1) under a linearly
ramped load (10, 5, 0) at node 2 gives a nonzero acceleration a2 at node 2.

Reference: inertialForce = [0, 0, 0, m * a2x, m * a2y, 0], with a2 read from
nodeAccel. The undamped case and a mass-proportional Rayleigh case
(alphaM = 0.4, so the damping force alphaM * M * v is nonzero and must be
removed) are both checked.

Before the fix the three accessors were evaluated in one expression with an
unspecified order, while getResistingForce and getResistingForceIncInertia
return the same element vector; with GCC's order the response was 0 (undamped)
instead of M * a.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


L = 2.0
RHO = 3.0
M_LUMP = 0.5 * RHO * L
REL_TOL = 1.0e-10


def _cantilever_step(alpha_m):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, 200.0, 1.0, 1.0)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Legendre", 1, 1, 3)
    ops.element("dispBeamColumn", 1, 1, 2, 1, 1, "-mass", RHO)
    if alpha_m != 0.0:
        ops.rayleigh(alpha_m, 0.0, 0.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0, 5.0, 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    assert ops.analyze(1, 0.1) == 0

    a2 = ops.nodeAccel(2)
    v2 = ops.nodeVel(2)
    got = list(ops.eleResponse(1, "inertialForce"))
    ops.wipe()
    return got, a2, v2


@pytest.mark.parametrize("alpha_m", [0.0, 0.4])
def test_inertial_force_equals_lumped_mass_times_accel(alpha_m):
    got, a2, v2 = _cantilever_step(alpha_m)
    assert abs(a2[0]) > 1.0e-3 and abs(a2[1]) > 1.0e-3, a2
    if alpha_m != 0.0:
        # the damping force removed from the response is not negligible
        assert abs(alpha_m * M_LUMP * v2[0]) > 1.0e-3, v2

    expected = [0.0, 0.0, 0.0, M_LUMP * a2[0], M_LUMP * a2[1], 0.0]
    assert len(got) == 6, got
    scale = max(abs(x) for x in expected)
    for g, e in zip(got, expected):
        assert abs(g - e) <= REL_TOL * scale, (got, expected)
