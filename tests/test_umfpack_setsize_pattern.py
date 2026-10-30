"""UmfPack re-analysis after the degree-of-freedom graph changes.

UmfpackGenLinSOE::setSize reserved and appended to the column pointer and row
index arrays without clearing them, so a second setSize appended the new
sparsity pattern behind the old one. The assembly and the symbolic
factorization kept reading the old pattern (entries of a new element were not
assembled, a shrunk pattern left the value array too short).

Model: one-dimensional chain of unit-length trusses, nodes 0..3, node 0 fixed,
EA = 100 (spring stiffness k = 100), force 10 x lambda at node 3, LoadControl
step 1, Newton, constraints Transformation. analyze(1) reaches lambda = 1, the
graph changes, analyze(1) reaches lambda = 2. All responses are linear, so the
displacements are exact closed forms:

  chain        u = (0.1, 0.2, 0.3) x lambda            (three springs in series)
  chain + brace between nodes 1 and 3 (EA = 100, L = 2, k = 50):
               lambda = 1: u = (0.1, 0.15, 0.2)
               lambda = 2: u = (0.2, 0.3, 0.4)
  sp_add       sp on node 3 of 0.002 x t added after the first analyze, t = 2:
               u = (0.004 / 3, 0.008 / 3, 0.004)

  elem_add     brace added after the first analyze (pattern grows). Truss
               records the current elongation 0.2 as initial displacement and
               adds it to its length, L = 2.2, so k = 100 / 2.2 and the brace
               starts unstressed:
               lambda = 1: (0.1, 0.2, 0.3)
               lambda = 2: 50 x + (100 / 2.2) (x - 0.2) = 20 gives
                           x = u3 - u1 = 32 / 105, u = (0.2, 0.2 + x / 2, 0.2 + x)
  elem_remove  brace removed after the first analyze (pattern shrinks)
               lambda = 1: (0.1, 0.15, 0.2); lambda = 2: (0.2, 0.4, 0.6)
  sp_add       free equations drop from 3 to 2

Each case runs in a child interpreter so that a crash in a failing case is
reported as a test failure. Every case fails on the unmodified solver.
"""
import json
import os
import subprocess
import sys

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


ENGINE = os.path.abspath(ops.__file__)
REL_TOL = 1.0e-9

CHILD = r'''
import json, os, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

CASE = sys.argv[1]
out = {"engine": os.path.abspath(ops.__file__)}

ops.wipe()
ops.model("basic", "-ndm", 1, "-ndf", 1)
for i in range(4):
    ops.node(i, float(i))
ops.fix(0, 1)
ops.uniaxialMaterial("Elastic", 1, 100.0)
for i in range(1, 4):
    ops.element("Truss", i, i - 1, i, 1.0, 1)
if CASE == "elem_remove":
    ops.element("Truss", 7, 1, 3, 1.0, 1)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(3, 10.0)
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormUnbalance", 1.0e-10, 10)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0)
ops.analysis("Static")

def disp():
    return [ops.nodeDisp(i, 1) for i in (1, 2, 3)]

out["rc1"] = ops.analyze(1)
out["u1"] = disp()
if CASE == "sp_add":
    ops.timeSeries("Linear", 2, "-factor", 0.002)
    ops.pattern("Plain", 2, 2)
    ops.sp(3, 1, 1.0)
elif CASE == "elem_add":
    ops.element("Truss", 7, 1, 3, 1.0, 1)
elif CASE == "elem_remove":
    ops.remove("ele", 7)
out["rc2"] = ops.analyze(1)
out["u2"] = disp()
print("RESULT " + json.dumps(out))
'''

EXPECTED = {
    "sp_add": ([0.1, 0.2, 0.3], [0.004 / 3, 0.008 / 3, 0.004]),
    "elem_add": ([0.1, 0.2, 0.3], [0.2, 0.2 + 16.0 / 105.0, 0.2 + 32.0 / 105.0]),
    "elem_remove": ([0.1, 0.15, 0.2], [0.2, 0.4, 0.6]),
}


def _run(case):
    proc = subprocess.run(
        [sys.executable, "-c", CHILD, case],
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=120,
    )
    assert proc.returncode == 0, (
        f"{case}: child exited with {proc.returncode}\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, f"{case}: no result\n{proc.stdout}\n{proc.stderr}"
    out = json.loads(lines[-1][len("RESULT "):])
    assert os.path.normcase(out["engine"]) == os.path.normcase(ENGINE), out
    return out


@pytest.mark.parametrize("case", sorted(EXPECTED))
def test_reanalysis_after_graph_change(case):
    out = _run(case)
    first, second = EXPECTED[case]
    assert out["rc1"] == 0, out
    assert out["u1"] == pytest.approx(first, rel=REL_TOL), out
    assert out["rc2"] == 0, out
    assert out["u2"] == pytest.approx(second, rel=REL_TOL), out
