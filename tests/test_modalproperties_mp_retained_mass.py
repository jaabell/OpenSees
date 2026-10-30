"""modalProperties keeps the mass of multi-point constrained DOFs.

DomainModalProperties treated every DOF whose DOF_Group ID is -1 as fixed.
The Transformation handler also gives -1 to the constrained DOFs of
multi-point constraints, so mass placed on such a node dropped out of the
free mass, the participation factors and the effective masses, and
responseSpectrumAnalysis returned zero.

Model (ndm 2, ndf 3): a truss of stiffness k = E*A/L = 3000*10/10 = 3000
from the fixed node 1 to node 2; node 3 is tied to node 2 by equalDOF 2 3 1
and carries all the mass, m = 8. Translation Y and rotation are fixed at
nodes 2 and 3 by SP constraints. Closed form of the single free mode:

    lambda = k / m = 375,  total mass X = 8,  free mass X = 8,
    free mass Y = 0 (SP-fixed),  effective mass X = 8,
    participation factor X = sqrt(8)  (mass-normalised mode).

responseSpectrumAnalysis scales its modal response by these participation
factors, so a zero factor gives a zero response.

Before the fix the Transformation case gives free mass X = 0 and effective
mass X = 0. Penalty is a control that is right before and after the fix.
Each case runs in a child interpreter.
"""
import json
import subprocess
import sys

import pytest

try:
    import opensees as ops  # noqa: F401
except ModuleNotFoundError:
    import openseespy.opensees as ops  # noqa: F401

CHILD = r'''
import json, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

handler = sys.argv[1]
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)
ops.node(1, 0.0, 0.0)
ops.node(2, 10.0, 0.0)
ops.node(3, 10.0, 0.0)
ops.fix(1, 1, 1, 1)
ops.fix(2, 0, 1, 1)
ops.fix(3, 0, 1, 1)
ops.mass(3, 8.0, 8.0, 0.0)
ops.uniaxialMaterial("Elastic", 1, 3000.0)
ops.element("Truss", 1, 1, 2, 10.0, 1)
ops.equalDOF(2, 3, 1)
if handler == "Penalty":
    ops.constraints("Penalty", 1.0e10, 1.0e10)
else:
    ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("BandGeneral")
ops.integrator("LoadControl", 1.0)
ops.analysis("Static")
lam = ops.eigen("fullGenLapack", 1)
props = ops.modalProperties("-return")
print("RESULT " + json.dumps({"lam": lam, "props": props}), flush=True)
'''


def _run(handler):
    proc = subprocess.run([sys.executable, "-c", CHILD, handler],
                          stdin=subprocess.DEVNULL, capture_output=True,
                          text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr[-2000:]
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, proc.stdout + proc.stderr[-2000:]
    return json.loads(lines[-1][len("RESULT "):])


@pytest.mark.parametrize("handler,tol", [("Transformation", 1e-9),
                                         ("Penalty", 1e-4)])
def test_mass_on_constrained_node_is_free_mass(handler, tol):
    r = _run(handler)
    p = r["props"]
    assert r["lam"][0] == pytest.approx(375.0, rel=tol)
    assert p["totalMass"][0] == pytest.approx(8.0, rel=tol)
    assert p["totalFreeMass"][0] == pytest.approx(8.0, rel=tol)
    # The penalty handler numbers SP-constrained DOFs, so only Transformation
    # reports the fixed Y direction as having no free mass.
    if handler == "Transformation":
        assert p["totalFreeMass"][1] == pytest.approx(0.0, abs=tol)
    assert p["partiMassMX"][0] == pytest.approx(8.0, rel=tol)
    assert abs(p["partiFactorMX"][0]) == pytest.approx(8.0 ** 0.5, rel=tol)
