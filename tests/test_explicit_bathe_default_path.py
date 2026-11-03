"""ExplicitBathe: Python registration and the default path.

The reference values below were produced by the unmodified ExplicitBathe of
the base branch through the classic Tcl interpreter (OpenSees.exe), on the
same model: a two-DOF zeroLength chain (k = 100, masses 1.0 and 0.5), a
linear-in-time load of 1.0 on the free end, an initial velocity of 0.5 on the
free end, rayleigh 0.1 0 0 0, system Diagonal, algorithm Linear,
`integrator ExplicitBathe 0.54`, dt = 0.01, after steps 1, 10, 100 and 500.
The classic Tcl nodeDisp/nodeVel/nodeAccel return the value printed with 20
decimals ("%35.20f"), so a reference is the double rounded to 1e-20, not the
double itself: values near 1e-6 lose their last bits.

Checks:
  * the integrator is available from Python;
  * the default path matches the references to the precision they carry,
    |x - x_ref| <= 1e-20 + 2 ulp(x_ref);
  * enabling the critical time step estimate (row-sum and HRZ lumping) does
    not change the trajectory:
    the 500-step histories are equal as floats to the default run.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


# step: (u2, u3, v2, v3, a2, a3)
REFERENCE = {
    1: (6.2099999999999998e-06, 0.0049866784000000006,
        0.0024925175979999999, 0.49476372570423705,
        0.49722874, -1.0256708640000001),
    10: (0.006775963362791417, 0.035797711026313156,
         0.17543809399042704, 0.11800198502015731,
         2.2070112383049945, -5.616137199837028),
    100: (0.033296496397423492, 0.043778782558142097,
          -0.11517615164846932, 0.28192607179937773,
          -2.2698495582995615, -0.12473120051519704),
    500: (0.066129876884802485, 0.10274938998310312,
          0.15340941713605091, 0.12126149998048337,
          -2.9663979572977408, 2.6639815916116447),
}

DEFAULT = ("ExplicitBathe", 0.54)
VARIANTS = [
    ("ExplicitBathe", 0.54, 1),
    ("ExplicitBathe", 0.54, 1, "-lump", "hrz"),
]


def _run(integrator_args):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.node(3, 0.0)
    ops.fix(1, 1)
    ops.mass(2, 1.0)
    ops.mass(3, 0.5)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.element("zeroLength", 2, 2, 3, "-mat", 1, "-dir", 1)
    ops.setNodeVel(3, 1, 0.5, "-commit")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0)
    ops.rayleigh(0.1, 0.0, 0.0, 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    history = []
    for i in range(1, 501):
        assert ops.analyze(1, 0.01) == 0, "step %d failed" % i
        history.append((ops.nodeDisp(2, 1), ops.nodeDisp(3, 1),
                        ops.nodeVel(2, 1), ops.nodeVel(3, 1),
                        ops.nodeAccel(2, 1), ops.nodeAccel(3, 1)))
    return history


def test_default_path_matches_reference():
    history = _run(DEFAULT)
    for step, ref in REFERENCE.items():
        for got, exp in zip(history[step - 1], ref):
            assert abs(got - exp) <= 1.0e-20 + 2.0 * math.ulp(exp), (
                "step %d: %r != %r" % (step, got, exp))


@pytest.mark.parametrize("integrator_args", VARIANTS,
                         ids=[" ".join(str(a) for a in v[2:]) for v in VARIANTS])
def test_options_do_not_change_the_trajectory(integrator_args):
    assert _run(integrator_args) == _run(DEFAULT)
