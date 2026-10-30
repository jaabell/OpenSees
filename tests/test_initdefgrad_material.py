"""nDMaterial InitDefGrad: command, options, inner validation and the F0 response.

InitDefGrad wraps a finite-strain material and forwards F F0^-1 through
setTrialF. No element in this source tree calls setTrialF, so the
re-referencing itself is not driven here. This file checks:

  * the command and its aliases (InitDefGradNDMaterial, StagedDefGrad,
    StagedDefGradNDMaterial) over a LogStrain inner, with -noInitF and -F0;
  * refusal, without terminating the interpreter, of a missing inner, a
    small-strain inner (ElasticIsotropic), an -F0 with fewer than 9
    components and an unknown option;
  * a singular -F0 is refused by the command; the constructor would exit
    the process, so this case runs in a child interpreter;
  * the "F0" material response returns the explicit -F0 row by row, and the
    identity before any capture;
  * a small-strain element (stdBrick) driving the material through
    setTrialStrain makes the analysis fail with a message.
"""
import subprocess
import sys
import textwrap

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

F0 = [1.02, 0.01, 0.0,
      0.0, 0.99, 0.03,
      0.0, 0.0, 1.01]


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.nDMaterial("LogStrain", 2, 1)


def _refused(*args):
    try:
        ops.nDMaterial(*args)
    except Exception:
        return True
    return False


def test_command_options_and_aliases():
    _model()
    ops.nDMaterial("InitDefGrad", 10, 2)
    ops.nDMaterial("InitDefGradNDMaterial", 11, 2, "-noInitF")
    ops.nDMaterial("StagedDefGrad", 12, 2, "-F0", *F0)
    ops.nDMaterial("StagedDefGradNDMaterial", 13, 2)


def test_invalid_input_is_refused():
    _model()
    assert _refused("InitDefGrad", 20, 999)              # inner not defined
    assert _refused("InitDefGrad", 21, 1)                # small-strain inner
    assert _refused("InitDefGrad", 22, 2, "-F0", 1.0, 0.0, 0.0)
    assert _refused("InitDefGrad", 24, 2, "-bogus")
    ops.nDMaterial("InitDefGrad", 25, 2)                 # interpreter still usable


def test_singular_f0_is_refused():
    script = textwrap.dedent("""
        try:
            import opensees as ops
        except ModuleNotFoundError:
            import openseespy.opensees as ops
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
        ops.nDMaterial("LogStrain", 2, 1)
        try:
            ops.nDMaterial("InitDefGrad", 3, 2, "-F0", 1, 0, 0, 0, 1, 0, 0, 0, 0)
            print("ACCEPTED")
        except Exception:
            print("REFUSED")
        ops.nDMaterial("InitDefGrad", 4, 2)
        print("ALIVE")
    """)
    run = subprocess.run([sys.executable, "-c", script], capture_output=True,
                         text=True, timeout=120)
    assert "REFUSED" in run.stdout, run.stdout + run.stderr
    assert "ALIVE" in run.stdout, run.stdout + run.stderr


def _brick_with(mat_tag):
    coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    for i, xyz in enumerate(coords, 1):
        ops.node(i, *xyz)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, mat_tag)


@pytest.mark.parametrize("opts,expected", [
    (("-F0", *F0), F0),
    ((), [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]),
])
def test_f0_response(opts, expected):
    _model()
    ops.nDMaterial("InitDefGrad", 3, 2, *opts)
    _brick_with(3)
    got = ops.eleResponse(1, "material", 1, "F0")
    assert got == pytest.approx(expected, abs=1e-15)


def test_small_strain_element_is_refused(capfd):
    _model()
    ops.nDMaterial("InitDefGrad", 3, 2)
    _brick_with(3)
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
    ok = ops.analyze(1)
    err = capfd.readouterr().err
    assert ok < 0
    assert "setTrialStrain() is invalid" in err


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
