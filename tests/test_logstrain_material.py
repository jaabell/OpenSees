"""nDMaterial LogStrain: command, inner validation and the small-strain guard.

LogStrain is driven through setTrialF(F) by a finite-strain element. No
element in this source tree calls setTrialF, so the stress update itself is
checked by tests/test_logstrain_kernel.cpp (the numerical core) and this file
checks what is reachable from the interpreter:

  * the command builds over 3D inner materials (ElasticIsotropic, J2Plasticity)
    with both spellings, LogStrain and LogStrainNDMaterial;
  * a missing inner, or an inner without a 3D (order 6) view (ElasticPlaneStress),
    is refused without terminating the interpreter;
  * a small-strain element (stdBrick) that drives the material through
    setTrialStrain gets -1 with a message, so the analysis fails instead of
    running with a strain vector fed into the finite-strain kinematics;
  * the "Cauchy" material response is answered by the wrapper (6 components),
    not by the inner material.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)


def _refused(*args):
    """True when the nDMaterial command fails (exception or no material)."""
    try:
        ops.nDMaterial(*args)
    except Exception:
        return True
    return False


def test_builds_over_3d_inner_materials():
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("J2Plasticity", 2, 1.6e5, 7.7e4, 250.0, 250.0, 0.0, 1.0e3)
    ops.nDMaterial("LogStrain", 10, 1)
    ops.nDMaterial("LogStrainNDMaterial", 11, 2)


def test_invalid_inner_is_refused():
    _model()
    ops.nDMaterial("ElasticPlaneStress", 3, 2.0e5, 0.3, 0.0)
    assert _refused("LogStrain", 20, 999)   # inner not defined
    assert _refused("LogStrain", 21, 3)     # no 3D view
    # the interpreter is still usable
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("LogStrain", 22, 1)


def _brick_with(mat_tag):
    coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    for i, xyz in enumerate(coords, 1):
        ops.node(i, *xyz)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, mat_tag)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.0, 0.0, -1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.1)
    ops.test("NormDispIncr", 1.0e-8, 5)
    ops.algorithm("Newton")
    ops.analysis("Static")


def test_small_strain_element_is_refused(capfd):
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("LogStrain", 2, 1)
    _brick_with(2)
    ok = ops.analyze(1)
    err = capfd.readouterr().err
    assert ok < 0
    assert "setTrialStrain() is invalid" in err


def test_control_inner_alone_converges():
    """Same brick with the inner material directly: the analysis converges."""
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    _brick_with(1)
    assert ops.analyze(1) == 0


def test_cauchy_response_is_owned_by_wrapper():
    _model()
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("LogStrain", 2, 1)
    _brick_with(2)
    stress = ops.eleResponse(1, "material", 1, "Cauchy")
    assert len(stress) == 6
    assert all(s == 0.0 for s in stress)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
