"""ExplicitBathe critical time step through the shared per-element estimate.

One 2-D Truss, E = 100, A = 1, L = 1, rho = 1 (le/c = 0.1). Checks:
  * with the critical time step flag, criticalTimeStep returns le/c = 0.1
    (the estimate of the base read the mass after the stiffness and, since
    Truss returns one shared matrix from both, reported inf);
  * without the flag the value is <= 0;
  * -lump on a consistent-mass Truss: rowsum and hrz give 0.1, diagonal gives
    2/sqrt(600) = 0.08165.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _truss(integrator_args, cmass=0):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 1.0, 0.0)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0, "-cMass", cmass)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    assert ops.analyze(1, 1.0e-3) == 0
    return ops.criticalTimeStep()


def test_truss_critical_time_step():
    assert _truss(("ExplicitBathe", 0.54, 1)) == pytest.approx(0.1, rel=1e-9)


def test_no_value_without_flag():
    assert _truss(("ExplicitBathe", 0.54)) <= 0.0


@pytest.mark.parametrize("lump,expected", [
    ("rowsum", 0.1),
    ("hrz", 0.1),
    ("diagonal", 2.0 / math.sqrt(600.0)),
])
def test_lumping_on_consistent_truss(lump, expected):
    dt = _truss(("ExplicitBathe", 0.54, 1, "-lump", lump), cmass=1)
    assert dt == pytest.approx(expected, rel=1e-9)
