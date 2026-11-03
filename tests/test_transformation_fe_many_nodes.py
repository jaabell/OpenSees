"""TransformationFE with an element that has more than 64 nodes.

TransformationFE keeps one DOF count per element node in a shared buffer.
The buffer held 64 ints, so an element with more than 64 nodes wrote past it
on every tangent and residual assembly under constraints Transformation.
The damage is heap corruption, which usually ends the process at a later
allocation or at wipe.

Model: a KinematicCoupling element with a reference node and 150 slave nodes
(151 nodes, 87 entries past the old buffer) on a circle of radius 1. The
reference node is prescribed a translation c = (0.003, -0.002, 0.001) with
no rotation, so every slave must move by exactly c. The model is built,
solved and wiped three times in a child interpreter; the child must exit
with status 0 and report the exact slave displacements.

The overrun is undefined behaviour, so the unfixed build is not guaranteed
to fail on every platform or allocator; it fails when the overwritten heap
block is used or freed later, which the repeated wipe makes likely.
"""
import json
import subprocess
import sys

import pytest

try:
    import opensees as ops  # noqa: F401
except ModuleNotFoundError:
    import openseespy.opensees as ops  # noqa: F401

NSLAVE = 150
C = (0.003, -0.002, 0.001)

CHILD = r'''
import json, math, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

nslave = int(sys.argv[1])
c = json.loads(sys.argv[2])
worst = 0.0
for rep in range(3):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0, "-ndf", 6)
    slaves = []
    for i in range(nslave):
        a = 2.0 * math.pi * i / nslave
        ops.node(2 + i, math.cos(a), math.sin(a), 0.0)
        slaves.append(2 + i)
    ops.element("KinematicCoupling", 1, 1, nslave, *slaves, "-k", 1.0e7)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for d, v in enumerate(list(c) + [0.0, 0.0, 0.0]):
        ops.sp(1, d + 1, v)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        print(json.dumps({"ok": False}))
        sys.exit(0)
    for n in slaves:
        u = ops.nodeDisp(n)
        worst = max(worst, max(abs(u[k] - c[k]) for k in range(3)))
ops.wipe()
print(json.dumps({"ok": True, "worst": worst}))
'''


def test_transformation_fe_element_with_more_than_64_nodes():
    p = subprocess.run(
        [sys.executable, "-c", CHILD, str(NSLAVE), json.dumps(C)],
        capture_output=True, text=True, timeout=60)
    assert p.returncode == 0, p.stderr[-500:]
    out = json.loads(p.stdout.strip().splitlines()[-1])
    assert out["ok"]
    assert out["worst"] < 1.0e-9
