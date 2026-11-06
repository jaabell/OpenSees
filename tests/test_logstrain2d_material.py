"""nDMaterial LogStrain2D: command, inner validation and the small-strain guard.

LogStrain2D is driven through setTrialF(F) with a 2x2 deformation gradient by
a plane finite-strain element. No element in this source tree calls
setTrialF, so this file checks what is reachable from the interpreter:

  * the command builds over a 3D inner (ElasticIsotropic, J2Plasticity) in
    plane strain (default and -planeStrain) and plane stress (-planeStress,
    -PlaneStress);
  * a missing inner, or an inner without a 3D (order 6) view (ElasticPlaneStress),
    is refused without terminating the interpreter;
  * a small-strain quad that drives the material through setTrialStrain gets
    -1 with a message, so the analysis fails (the same quad with the inner
    material converges);
  * the "Cauchy" material response is answered by LogStrain2D with the three
    in-plane components.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _model(ndm=3, ndf=3):
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndf)


def _refused(*args):
    try:
        ops.nDMaterial(*args)
    except Exception:
        return True
    return False


def test_builds_both_plane_modes():
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("J2Plasticity", 2, 1.6e5, 7.7e4, 250.0, 250.0, 0.0, 1.0e3)
    ops.nDMaterial("LogStrain2D", 10, 1)
    ops.nDMaterial("LogStrain2D", 11, 1, "-planeStrain")
    ops.nDMaterial("LogStrain2D", 12, 1, "-planeStress")
    ops.nDMaterial("LogStrain2D", 13, 2, "-PlaneStress")


def test_invalid_inner_is_refused():
    _model()
    ops.nDMaterial("ElasticPlaneStress", 3, 2.0e5, 0.3, 0.0)
    assert _refused("LogStrain2D", 20, 999)
    assert _refused("LogStrain2D", 21, 3, "-planeStress")
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("LogStrain2D", 22, 1, "-planeStress")


def _quad_with(mat_tag, plane="PlaneStrain"):
    for i, xy in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(i, *xy)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(4, 1, 0)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, plane, mat_tag)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0.0)
    ops.load(3, 1.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.1)
    ops.test("NormDispIncr", 1.0e-8, 5)
    ops.algorithm("Newton")
    ops.analysis("Static")


@pytest.mark.parametrize("flag,plane", [("-planeStrain", "PlaneStrain"),
                                        ("-planeStress", "PlaneStress")])
def test_small_strain_element_is_refused(capfd, flag, plane):
    _model(2, 2)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1, flag)
    _quad_with(2, plane)
    ok = ops.analyze(1)
    err = capfd.readouterr().err
    assert ok < 0
    assert "setTrialStrain() is invalid" in err


def test_control_inner_alone_converges():
    _model(2, 2)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    _quad_with(1)
    assert ops.analyze(1) == 0


def test_cauchy_response_is_owned_by_wrapper():
    _model(2, 2)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    _quad_with(2)
    stress = ops.eleResponse(1, "material", 1, "Cauchy")
    assert len(stress) == 3
    assert all(s == 0.0 for s in stress)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
