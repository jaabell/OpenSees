"""ASDAbsorbingBoundary3D: node-order handedness and degenerate bricks.

The free-field column of a vertical boundary element (addKff, addRff, and the
free-field damping built from addKff) is integrated over the 8-node brick taken
in the order the nodes are written, while the free-field DOF mapping and the
mass use the nodes sorted by coordinates. The volume measure was the signed
det J, so the same unit cube written in a mirrored order (det J = -1/8 at every
gauss point instead of +1/8) assembled a negative free-field stiffness and
residual. The reference is the right-handed order of the same geometry: the
mirrored orders must give the same element force and the same displacement
history.

A brick whose det J changes sign between gauss points (folded) or vanishes
(collapsed) has no valid volume measure in either order; setDomain now rejects
it on vertical boundaries with a message naming the element and its nodes.
Before, the folded brick was accepted and the collapsed one stopped later with
an unrelated distortion message.

One-element model: unit cube, G = 1000, v = 0.3, rho = 2, bottom nodes fixed,
one vertical boundary type per case (L, R, F, K). The element is switched
to the absorbing stage from a zero state and runs 20 linear Newmark steps
under a short top-node load pulse. Compared quantities: the top-node
displacements at every step and the final element force, indexed by node id. With the signed det J the mirrored orders diverge from the reference
(top displacements of order 1e+2 against 1e-4 for the L case).

Each case runs in a child interpreter, because an element that refuses its
input stops the process.
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

NODES = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
         (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
RIGHT = [1, 2, 3, 4, 5, 6, 7, 8]
MIRROR = [1, 4, 3, 2, 5, 8, 7, 6]      # mirrored about the x = y plane
TOP_FIRST = [5, 6, 7, 8, 1, 2, 3, 4]   # mirrored about z
FOLDED = [1, 2, 3, 4, 5, 7, 6, 8]      # top face written as a bow tie

CHILD = r'''
import json, os, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

case = json.loads(sys.argv[1])
out = {"engine": os.path.abspath(ops.__file__)}
order = case["order"]
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
for i, c in enumerate(case["coords"], 1):
    ops.node(i, *c)
ops.element("ASDAbsorbingBoundary3D", 1, *order, 1000.0, 0.3, 2.0, case["btype"])
for i in range(1, 5):
    ops.fix(i, 1, 1, 1)

# switch to the absorbing stage directly: zero stage-0 state
ops.setParameter("-val", 1, "-ele", 1, "stage")

# transient under a load pulse on the top nodes
ops.timeSeries("Path", 2, "-dt", 0.01, "-values", 0.0, 1.0, 0.0, 0.0)
ops.pattern("Plain", 2, 2)
ops.load(5, 1.0, 0.5, -2.0)
ops.load(6, 0.0, 0.3, -1.0)
ops.load(7, 0.5, 0.0, -1.0)
ops.load(8, 0.2, 0.1, -3.0)
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("BandGeneral")
ops.algorithm("Linear")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
hist = []
for k in range(20):
    ok = ops.analyze(1, 0.01)
    hist.append([ok] + [ops.nodeDisp(n, d) for n in (5, 6, 7, 8) for d in (1, 2, 3)])
out["hist"] = hist
force = ops.eleResponse(1, "force")
out["force"] = {str(order[k]): force[3 * k:3 * k + 3] for k in range(8)}
print("RESULT " + json.dumps(out))
'''


def run_case(order, btype, coords=NODES):
    case = json.dumps({"order": order, "btype": btype, "coords": coords})
    r = subprocess.run([sys.executable, "-c", CHILD, case],
                       stdin=subprocess.DEVNULL, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=120)
    text = r.stdout + r.stderr
    res = None
    for line in r.stdout.splitlines():
        if line.startswith("RESULT "):
            res = json.loads(line[len("RESULT "):])
            assert os.path.normcase(res["engine"]) == os.path.normcase(ENGINE), res
    return r.returncode, text, res


@pytest.mark.parametrize("btype", ["L", "R", "F", "K"])
def test_mirrored_order_same_response(btype):
    rc, text, ref = run_case(RIGHT, btype)
    assert rc == 0 and ref is not None, text
    assert all(h[0] == 0 for h in ref["hist"])
    assert max(abs(x) for h in ref["hist"] for x in h[1:]) > 0.0
    for order in (MIRROR, TOP_FIRST):
        rc, text, res = run_case(order, btype)
        assert rc == 0 and res is not None, text
        assert all(h[0] == 0 for h in res["hist"]), res["hist"]
        for n in map(str, range(1, 9)):
            assert res["force"][n] == pytest.approx(ref["force"][n], rel=1e-10, abs=1e-6)
        for h, h0 in zip(res["hist"], ref["hist"]):
            assert h[1:] == pytest.approx(h0[1:], rel=1e-8, abs=1e-18)


@pytest.mark.parametrize("btype", ["L", "R", "F", "K"])
def test_folded_brick_rejected(btype):
    # det J ranges from -0.0305 to +0.1138 over the gauss points
    rc, text, res = run_case(FOLDED, btype)
    assert rc != 0 and res is None, text
    assert "element 1 is a degenerate brick" in text
    assert "1 2 3 4 5 7 6 8" in text


def test_collapsed_brick_rejected():
    # zero height: det J = 0 at every gauss point
    coords = NODES[:4] + [(x, y, 0.0) for (x, y, _) in NODES[4:]]
    rc, text, res = run_case(RIGHT, "L", coords)
    assert rc != 0 and res is None, text
    assert "element 1 is a degenerate brick" in text


def test_bottom_boundary_not_checked():
    # a pure bottom element never integrates over the brick
    rc, text, res = run_case(FOLDED, "B")
    assert rc == 0 and res is not None, text
    assert "degenerate" not in text
