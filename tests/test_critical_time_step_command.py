"""criticalTimeStep command.

Checks that the command exists, that it returns a non-positive value for a
transient integrator that does not compute a critical time step (Newmark),
and that it fails when no transient integrator is set.
"""

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _sdof():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    ops.mass(2, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")


def test_non_explicit_integrator_reports_no_value():
    _sdof()
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    assert ops.analyze(1, 0.01) == 0
    assert ops.criticalTimeStep() <= 0.0


def test_without_transient_integrator_fails():
    _sdof()
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    with pytest.raises(ops.OpenSeesError):
        ops.criticalTimeStep()
