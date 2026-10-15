"""Symmetric-storage systems warn when an element matrix is unsymmetric.

The systems listed in SYSTEMS store only the upper triangle and drop the
mirrored entry of every element matrix, so an unsymmetric tangent is solved
as the symmetric matrix built from its upper triangle. These cases check
that the discard is now reported, once, and that the results are
untouched.

  plastic   one plane-strain quad on a non-associated DruckerPrager
            (rho = 0.398, rho_bar = 0.1), three nodes fixed, the fourth
            loaded (1, -3) in four steps of 0.25. First yield is near load
            factor 0.55; after it the free 2x2 tangent is unsymmetric
            (about 9.3e3 against 7.6e3 off the diagonal). The warning must
            appear exactly once with each symmetric-storage system, and never
            with BandGeneral, which stores both triangles.
  elastic   two 1D trusses in series, E = 4, A = 1, L = 2 each (k = EA/L = 2),
            end load P = 1. Closed form u2 = P/k = 0.5, u3 = 2P/k = 1.0, both
            exact in binary. No warning, and the displacements equal the
            closed form exactly.

Each case runs in a child interpreter so that the warning can be read from
its stderr.
"""
import json
import subprocess
import sys

import pytest

WARNING = "unsymmetric element matrix"
SYSTEMS = ["ProfileSPD", "BandSPD", "SProfileSPD"]

CHILD = r'''
import json, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

case, system = sys.argv[1], sys.argv[2]
ops.wipe()
if case == "plastic":
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for n, (x, y) in enumerate([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], 1):
        ops.node(n, x, y)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.fix(4, 1, 1)
    ops.nDMaterial("DruckerPrager", 1, 27777.78, 9259.26, 5.0, 0.398, 0.1,
                   0.0, 0.0, 0.0, 0.0, 1000.0, 1.0, 0.0)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0, -3.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", 1.0e-10, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.25)
    ops.analysis("Static")
    for _ in range(4):
        ops.analyze(1)
    disp = ops.nodeDisp(3)
else:
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 2.0)
    ops.node(3, 4.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, 4.0)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.element("truss", 2, 2, 3, 1.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system(system)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.analyze(1)
    disp = [ops.nodeDisp(2, 1), ops.nodeDisp(3, 1)]
ops.wipe()
print("RESULT " + json.dumps(disp))
'''


def _run(case, system):
    proc = subprocess.run(
        [sys.executable, "-c", CHILD, case, system],
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    assert proc.returncode == 0, (
        f"{case}/{system}: child exited with {proc.returncode}\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, f"{case}/{system}: no result\n{proc.stdout}\n{proc.stderr}"
    return json.loads(lines[-1][len("RESULT "):]), proc.stderr


@pytest.mark.parametrize("system", SYSTEMS)
def test_unsymmetric_tangent_warns_once(system):
    _, err = _run("plastic", system)
    assert err.count(WARNING) == 1, err


def test_unsymmetric_tangent_no_warning_on_general_storage():
    _, err = _run("plastic", "BandGeneral")
    assert WARNING not in err, err


@pytest.mark.parametrize("system", SYSTEMS)
def test_symmetric_model_no_warning_and_exact_result(system):
    disp, err = _run("elastic", system)
    assert WARNING not in err, err
    assert disp == [0.5, 1.0]
