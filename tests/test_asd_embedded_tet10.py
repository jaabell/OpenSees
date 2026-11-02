"""ASDEmbeddedNodeElementTet10: a 10-node (quadratic) tetrahedron host for
ASDEmbeddedNodeElement.

The element's host dispatch is keyed on retained-node count (3 => triangle,
4 => tet4); a quadratic host fed through the existing command would either be
misread as one of those two (the exact bug a hex host already hit this way:
8 nodes silently read as "some triangle/tet subset") or need yet another
overload disambiguated by argument count. ASDEmbeddedNodeElementTet10 is a
separate command instead, so the two paths can never be confused.

The natural coordinates of the embedded node still come from the 4 corner
nodes alone (the element's existing linear machinery, unchanged): exact for
any straight-sided element, since a quadratic map with mid-edge nodes at the
true edge midpoint collapses algebraically to the same affine map as the
linear tet4. Only the DISPLACEMENT interpolation is quadratic.

  reproduces_quadratic_field   prescribing the exact nodal values of a known
                                quadratic field on all 10 host nodes, an
                                embedded node at a known interior point picks
                                up the field's exact value there (to machine
                                precision) -- confirms the quadratic shape
                                functions are the ones actually evaluated,
                                not a silent fallback to linear/corner-only
                                interpolation (which would give a visibly
                                different, wrong answer at the same point).
  tet4_path_unchanged          the original 3/4-node ASDEmbeddedNodeElement
                                command, same style of check: regression
                                guard that this change didn't touch it.
"""
import json
import os
import subprocess
import sys

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


ENGINE = os.path.abspath(ops.__file__)

# gmsh tet10 node ordering: 4 corners, then mid-edge nodes
# 5=(1,2) 6=(2,3) 7=(1,3) 8=(1,4) 9=(3,4) 10=(2,4) -- matching
# TenNodeTetrahedron's default, non-"Std", numbering.
CORNERS = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0),
           3: (0.0, 1.0, 0.0), 4: (0.0, 0.0, 1.0)}
EDGES = [(1, 2), (2, 3), (1, 3), (1, 4), (3, 4), (2, 4)]

# the embedded node, strictly interior (barycentric L1=0.55 L2=0.2 L3=0.1 L4=0.15)
C = (0.2, 0.1, 0.15)

CHILD = r'''
import json, os, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

SCEN = sys.argv[1]
out = {"engine": os.path.abspath(ops.__file__)}


def quad_field(x, y, z):
    # an arbitrary genuine quadratic: linear-only interpolation over the
    # corners would NOT reproduce this at an interior point
    return x * x + 2.0 * y - z


if SCEN == "tet10":
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    corners = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0),
               3: (0.0, 1.0, 0.0), 4: (0.0, 0.0, 1.0)}
    for tag, xyz in corners.items():
        ops.node(tag, *xyz)
    edges = [(1, 2), (2, 3), (1, 3), (1, 4), (3, 4), (2, 4)]
    for k, (a, b) in enumerate(edges):
        mid = tuple((corners[a][d] + corners[b][d]) / 2.0 for d in range(3))
        ops.node(5 + k, *mid)
    c = (0.2, 0.1, 0.15)
    ops.node(11, *c)

    node_xyz = dict(corners)
    for k, (a, b) in enumerate(edges):
        node_xyz[5 + k] = tuple((corners[a][d] + corners[b][d]) / 2.0 for d in range(3))

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for tag in range(1, 11):
        xx, yy, zz = node_xyz[tag]
        ops.sp(tag, 1, quad_field(xx, yy, zz))
        ops.sp(tag, 2, 0.0)
        ops.sp(tag, 3, 0.0)

    ops.element("ASDEmbeddedNodeElementTet10", 1, 11, *range(1, 11), "-K", 1.0e14)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormUnbalance", 1.0e-8, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    out["rc"] = ops.analyze(1)
    out["ux"] = ops.nodeDisp(11, 1)
    out["uy"] = ops.nodeDisp(11, 2)
    out["uz"] = ops.nodeDisp(11, 3)
    out["exact"] = quad_field(*c)

elif SCEN == "tet4":
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    corners = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0),
               3: (0.0, 1.0, 0.0), 4: (0.0, 0.0, 1.0)}
    for tag, xyz in corners.items():
        ops.node(tag, *xyz)
    c = (0.2, 0.1, 0.15)
    ops.node(5, *c)

    linear_field = {1: 0.0, 2: 1.0, 3: 2.0, 4: -1.0}  # corner ux values
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for tag in range(1, 5):
        ops.sp(tag, 1, linear_field[tag])
        ops.sp(tag, 2, 0.0)
        ops.sp(tag, 3, 0.0)

    ops.element("ASDEmbeddedNodeElement", 1, 5, 1, 2, 3, 4, "-K", 1.0e14)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormUnbalance", 1.0e-8, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    out["rc"] = ops.analyze(1)
    out["ux"] = ops.nodeDisp(5, 1)
    # exact LINEAR interpolation at L1=0.55 L2=0.2 L3=0.1 L4=0.15
    out["exact"] = 0.55 * 0.0 + 0.2 * 1.0 + 0.1 * 2.0 + 0.15 * (-1.0)

print("RESULT " + json.dumps(out), flush=True)
'''


def _run(scenario):
    proc = subprocess.run(
        [sys.executable, "-c", CHILD, scenario],
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=60,
    )
    assert proc.returncode == 0, (
        f"{scenario}: child exited with {proc.returncode}\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, f"{scenario}: no result\n{proc.stdout}\n{proc.stderr}"
    out = json.loads(lines[-1][len("RESULT "):])
    assert os.path.normcase(out["engine"]) == os.path.normcase(ENGINE), out
    return out


def test_tet10_reproduces_quadratic_field():
    out = _run("tet10")
    assert out["rc"] == 0, out["rc"]
    assert abs(out["ux"] - out["exact"]) < 1.0e-8, out
    assert abs(out["uy"]) < 1.0e-8, out
    assert abs(out["uz"]) < 1.0e-8, out


def test_tet4_path_unchanged():
    out = _run("tet4")
    assert out["rc"] == 0, out["rc"]
    assert abs(out["ux"] - out["exact"]) < 1.0e-8, out
