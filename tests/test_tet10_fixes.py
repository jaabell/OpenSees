"""TenNodeTetrahedron output and error reporting.

Each case runs in a child interpreter so that stdout written by the element
(C++ std::cout) can be captured and a crash fails only that case.

  quiet    a one-element analysis that forms the initial stiffness writes
           nothing to stdout
  fail     a material that refuses setTrialStrain makes analyze() fail
  healthy  the same model with a working material converges (control)
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

CHILD = r'''
import json, os, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

SCEN = sys.argv[1]
out = {"engine": os.path.abspath(ops.__file__)}


def tet10(mats, e2=None):
    # unit corner tetrahedron, base face fixed, vertical load at the apex;
    # one TenNodeTetrahedron per material tag, all on the same ten nodes
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    xyz = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0),
           3: (0.0, 1.0, 0.0), 4: (0.0, 0.0, 1.0)}
    for tag, x in xyz.items():
        ops.node(tag, *x)
    for k, (a, b) in enumerate([(1, 2), (2, 3), (1, 3), (1, 4), (3, 4), (2, 4)]):
        ops.node(5 + k, *[(xyz[a][d] + xyz[b][d]) / 2.0 for d in range(3)])
    for tag in (1, 2, 3, 5, 6, 7):
        ops.fix(tag, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    if e2 is not None:
        ops.nDMaterial("ElasticIsotropic", 2, e2, 0.25)
        ops.nDMaterial("TimeVarying", 3, 2, 1, 0.0, 1000.0, 2000.0, 1.0)
    for e, m in enumerate(mats, start=1):
        ops.element("TenNodeTetrahedron", e, *range(1, 11), m)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(4, 0.0, 0.0, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 10, 0)


if SCEN == "quiet":
    tet10([1])
    ops.algorithm("Newton", "-initial")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    print("BEGIN", flush=True)
    out["rc"] = ops.analyze(1)
    out["u4"] = ops.nodeDisp(4)
elif SCEN in ("fail", "healthy"):
    # TimeVarying refuses setTrialStrain when the wrapped material has a
    # singular initial tangent (E = 0); a second element on the same nodes
    # keeps the system nonsingular, so only the refusal can stop analyze()
    tet10([1, 3], e2=0.0 if SCEN == "fail" else 1000.0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    out["rc"] = ops.analyze(1)

sys.stdout.flush()
print("RESULT " + json.dumps(out), flush=True)
'''


def _run(scenario):
    proc = subprocess.run(
        [sys.executable, "-c", CHILD, scenario],
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    assert proc.returncode == 0, (
        f"{scenario}: child exited with {proc.returncode}\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
    lines = proc.stdout.splitlines()
    results = [l for l in lines if l.startswith("RESULT ")]
    assert results, f"{scenario}: no result\n{proc.stdout}\n{proc.stderr}"
    out = json.loads(results[-1][len("RESULT "):])
    assert os.path.normcase(out["engine"]) == os.path.normcase(ENGINE), out
    out["stdout"] = lines
    return out


def test_tet10_writes_nothing_to_stdout():
    out = _run("quiet")
    assert out["rc"] == 0, out["rc"]
    assert out["u4"][2] > 0.0, out["u4"]
    lines = out["stdout"]
    begin = lines.index("BEGIN")
    stray = [l for l in lines[begin + 1:] if not l.startswith("RESULT ")]
    assert stray == [], f"{len(stray)} stray stdout lines, first: {stray[:5]}"


def test_tet10_reports_material_failure():
    assert _run("healthy")["rc"] == 0
    out = _run("fail")
    assert out["rc"] < 0, f"analyze() returned {out['rc']} with a failed material"
