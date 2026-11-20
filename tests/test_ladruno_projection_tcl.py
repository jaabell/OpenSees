"""constraints LadrunoProjection from the classic Tcl interpreter.

Runs OpenSees(.exe) on three scripts with two masses (2 and 3) tied by
equalDOF, system Diagonal, algorithm Linear, CentralDifferenceLadruno:
- a constant force of 10 for 100 steps of 1e-3: analyze returns 0, the
  tie holds to 1e-12 and the acceleration is 10 / 5 = 2;
- a constraint chain: analyze returns nonzero (the error of the handler's
  handle() reaches DirectIntegrationAnalysis);
- initial displacements off the constraint manifold: analyze returns
  nonzero (the error of the integrator's domainChanged() reaches
  DirectIntegrationAnalysis).

Skipped when no Tcl executable is found next to the Python module or in
OPENSEES_EXE.
"""

import os
import subprocess
from pathlib import Path

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _tcl_executable():
    env = os.environ.get("OPENSEES_EXE")
    if env and Path(env).is_file():
        return Path(env)
    here = Path(ops.__file__).resolve().parent
    for name in ("OpenSees.exe", "OpenSees"):
        exe = here / name
        if exe.is_file():
            return exe
    return None


pytestmark = pytest.mark.skipif(_tcl_executable() is None,
                                reason="no OpenSees Tcl executable (set OPENSEES_EXE)")

_ANALYSIS = """
constraints LadrunoProjection
numberer Plain
system Diagonal
test NormDispIncr 1e-12 10
algorithm Linear
integrator CentralDifferenceLadruno
analysis Transient
"""

TIED = """
model basic -ndm 1 -ndf 1
node 1 0.0
mass 1 2.0
node 2 1.0
mass 2 3.0
equalDOF 1 2 1
timeSeries Constant 1
pattern Plain 1 1 { load 1 10.0 }
""" + _ANALYSIS + """
set ok [analyze 100 1.0e-3]
puts "RESULT $ok [nodeDisp 1 1] [nodeDisp 2 1] [nodeAccel 1 1]"
"""

CHAIN = """
model basic -ndm 1 -ndf 1
node 1 0.0
mass 1 1.0
node 2 1.0
mass 2 1.0
node 3 2.0
mass 3 1.0
equalDOF 1 2 1
equalDOF 2 3 1
""" + _ANALYSIS + """
set ok [analyze 1 1.0e-3]
puts "RESULT $ok"
"""

OFF_MANIFOLD = """
model basic -ndm 1 -ndf 1
node 1 0.0
mass 1 2.0
node 2 1.0
mass 2 3.0
equalDOF 1 2 1
setNodeDisp 1 1 1.0 -commit
setNodeDisp 2 1 2.0 -commit
""" + _ANALYSIS + """
set ok [analyze 1 1.0e-3]
puts "RESULT $ok"
"""


def _run(tmp_path, text):
    script = tmp_path / "projection.tcl"
    script.write_text(text)
    res = subprocess.run([str(_tcl_executable()), str(script)],
                         capture_output=True, text=True, cwd=tmp_path,
                         timeout=300, stdin=subprocess.DEVNULL,
                         env=dict(os.environ))
    out = res.stdout + res.stderr
    lines = [x for x in out.splitlines() if x.startswith("RESULT")]
    assert lines, out
    return [float(v) for v in lines[-1].split()[1:]], out


def test_tied_pair(tmp_path):
    vals, out = _run(tmp_path, TIED)
    ok, u1, u2, a1 = vals
    assert ok == 0.0, out
    assert abs(u1 - u2) < 1e-12
    assert abs(a1 - 2.0) < 1e-10


def test_chain_fails(tmp_path):
    vals, out = _run(tmp_path, CHAIN)
    assert vals[0] != 0.0, out


def test_initial_conditions_off_manifold_fail(tmp_path):
    vals, out = _run(tmp_path, OFF_MANIFOLD)
    assert vals[0] != 0.0, out
