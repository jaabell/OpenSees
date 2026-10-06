"""Linear SOE accessors on an unsized system and on a zero-equation model.

Each case runs in a child interpreter: a regression ends the child with
exit(-1) or an access violation, which fails that case instead of the run.

  zero_fresh   every DOF fixed or sp-constrained from the start (0 equations)
  zero_shrink  3 equations, then the same SOE is resized to 0
  unsized      printA/printB before any analyze(), so setSize() never ran
  ok           a solved model still returns the real A and B
  shrink       FullGeneral resized 12 -> 6 equations reports the live size
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
EPS = (1.0e-3, -2.5e-4, 5.0e-4)   # prescribed ux, uy, uz on the x=1, y=1, z=1 faces

SYSTEMS_ZERO = ["FullGeneral", "BandGeneral", "BandSPD", "ProfileSPD",
                "SProfileSPD", "UmfPack"]
SYSTEMS_UNSIZED = ["FullGeneral", "Diagonal", "SparseSYM", "BandGeneral",
                   "BandSPD", "ProfileSPD", "SProfileSPD", "SparseGeneral",
                   "Pardiso", "UmfPack"]

CHILD = r'''
import json, os, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

SYS, SCEN = sys.argv[1], sys.argv[2]
EPS = %(EPS)r
out = {"engine": os.path.abspath(ops.__file__)}

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
try:
    ops.system(SYS)
except Exception:
    print("RESULT " + json.dumps({"unsupported": True, "engine": out["engine"]}))
    sys.exit(0)


def solver():
    ops.numberer("Plain")
    ops.system(SYS)
    ops.test("NormDispIncr", 1.0e-10, 20, 0)
    ops.algorithm("Newton")


def cube(free_corner):
    # unit stdBrick, 1/8-symmetry fixes, remaining DOFs driven by sp
    nodes = {1: ((0, 0, 0), (1, 1, 1)), 2: ((1, 0, 0), (0, 1, 1)),
             3: ((1, 1, 0), (0, 0, 1)), 4: ((0, 1, 0), (1, 0, 1)),
             5: ((0, 0, 1), (1, 1, 0)), 6: ((1, 0, 1), (0, 1, 0)),
             7: ((1, 1, 1), (0, 0, 0)), 8: ((0, 1, 1), (1, 0, 0))}
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (xyz, fx) in nodes.items():
        ops.node(tag, *map(float, xyz))
        if any(fx):
            ops.fix(tag, *fx)
    ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.25)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Path", 1, "-time", 0.0, 1.0, 2.0, "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    for tag, (xyz, fx) in nodes.items():
        if free_corner and tag == 7:
            continue
        for d in range(3):
            if xyz[d] == 1 and fx[d] == 0:
                ops.sp(tag, d + 1, EPS[d])
    ops.constraints("Transformation")
    solver()
    ops.integrator("LoadControl", 0.5)
    ops.analysis("Static")


def quads(nx, ny):
    # nx x ny plane-strain quads, bottom row fixed, one vertical load
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    tag = lambda i, j: 1 + i + j * (nx + 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(tag(i, j), float(i), float(j))
    for i in range(nx + 1):
        ops.fix(tag(i, 0), 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    e = 1
    for j in range(ny):
        for i in range(nx):
            ops.element("quad", e, tag(i, j), tag(i + 1, j), tag(i + 1, j + 1),
                        tag(i, j + 1), 1.0, "PlaneStrain", 1)
            e += 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tag(nx, ny), 0.0, -1.0)
    ops.constraints("Plain")
    solver()
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return tag


if SCEN == "zero_fresh":
    cube(free_corner=False)
    out["rc"] = ops.analyze(2)
    out["u7"] = ops.nodeDisp(7)
elif SCEN == "zero_shrink":
    cube(free_corner=True)
    out["rc1"] = ops.analyze(2)
    ops.timeSeries("Path", 2, "-time", 1.0, 2.0, "-values", 0.0, 1.0)
    ops.pattern("Plain", 2, 2)
    for d in range(3):
        ops.sp(7, d + 1, EPS[d])
    out["rc"] = ops.analyze(2)
    out["u7"] = ops.nodeDisp(7)
elif SCEN == "unsized":
    quads(1, 1)
    out["A"] = list(ops.printA("-ret"))
    out["B"] = list(ops.printB("-ret"))
elif SCEN == "ok":
    quads(1, 1)
    out["rc"] = ops.analyze(1)
    out["n"] = ops.systemSize()
    out["A"] = list(ops.printA("-ret"))
    out["B"] = list(ops.printB("-ret"))
elif SCEN == "shrink":
    tag = quads(2, 2)
    out["rc1"] = ops.analyze(1)
    out["n1"] = ops.systemSize()
    for i in range(3):
        ops.fix(tag(i, 2), 1, 1)
    out["rc"] = ops.analyze(1)
    out["n2"] = ops.systemSize()
    out["A2"] = len(ops.printA("-ret"))
    out["B2"] = len(ops.printB("-ret"))

print("RESULT " + json.dumps(out))
'''


def _run(system, scenario):
    proc = subprocess.run(
        [sys.executable, "-c", CHILD % {"EPS": EPS}, system, scenario],
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    assert proc.returncode == 0, (
        f"{system}/{scenario}: child exited with {proc.returncode}\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, f"{system}/{scenario}: no result\n{proc.stdout}\n{proc.stderr}"
    out = json.loads(lines[-1][len("RESULT "):])
    assert os.path.normcase(out["engine"]) == os.path.normcase(ENGINE), out
    if out.get("unsupported"):
        pytest.skip(f"system {system} is not available in this build")
    return out


@pytest.mark.parametrize("scenario,system",
                         [("zero_fresh", s) for s in SYSTEMS_ZERO + ["Diagonal"]]
                         + [("zero_shrink", s) for s in SYSTEMS_ZERO])
def test_zero_free_equations(scenario, system):
    out = _run(system, scenario)
    assert out["rc"] == 0, out
    assert all(abs(u - e) < 1e-12 for u, e in zip(out["u7"], EPS)), out


@pytest.mark.parametrize("system", SYSTEMS_UNSIZED)
def test_unsized_soe_print_returns_empty(system):
    out = _run(system, "unsized")
    assert out["A"] == [] and out["B"] == [], out


@pytest.mark.parametrize("system", [s for s in SYSTEMS_UNSIZED if s != "Diagonal"])
def test_sized_soe_print_returns_system(system):
    out = _run(system, "ok")
    assert out["rc"] == 0 and out["n"] == 4, out
    assert len(out["B"]) == 4, out
    if system == "FullGeneral":
        assert len(out["A"]) == 16 and any(v != 0.0 for v in out["A"]), out


def test_fullgeneral_shrink_reports_live_size():
    out = _run("FullGeneral", "shrink")
    assert out["rc1"] == 0 and out["rc"] == 0, out
    assert out["n2"] < out["n1"], out
    assert out["A2"] == out["n2"] ** 2 and out["B2"] == out["n2"], out
