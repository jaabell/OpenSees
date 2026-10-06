"""Small interpreter fixes: repeated eigen and printA -sparse option order.

  eigen_*        repeated openseespy eigen calls, with and without an analysis
  sparse_flags   openseespy printA with -sparse before or after -ret
  tcl_sparse     the same printA orderings in the classic Tcl interpreter

Python cases run in a child interpreter so that one failing case cannot take
down the run. The Tcl case needs the OpenSees executable named by the
OPENSEES_EXE environment variable and is skipped when it is unset or missing.
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
REL_TOL = 1.0e-8

CHILD = r'''
import json, os, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

SCEN = sys.argv[1]
out = {"engine": os.path.abspath(ops.__file__)}


def chain():
    # truss chain with lateral springs of distinct stiffness: 20 free DOFs,
    # no repeated eigenvalues
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 3000.0)
    ops.node(0, 0.0, 0.0)
    ops.fix(0, 1, 1)
    for i in range(1, 11):
        ops.node(i, 10.0 * i, 0.0)
        ops.mass(i, 5.0, 5.0)
        ops.element("Truss", i, i - 1, i, 10.0, 1)
        ops.node(100 + i, 10.0 * i, 10.0)
        ops.fix(100 + i, 1, 1)
        ops.element("Truss", 100 + i, 100 + i, i, 4.0 + 0.5 * i, 1)


def transient():
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1.0e-8, 10)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


def eig(*args):
    try:
        return list(ops.eigen(*args))
    except Exception as ex:
        return "raised: " + repr(ex)


if SCEN == "eigen_same_type":
    chain()
    out["runs"] = [eig(3), eig(3), eig(3)]
elif SCEN == "eigen_explicit_type":
    chain()
    out["runs"] = [eig("-genBandArpack", 3), eig("-genBandArpack", 3)]
elif SCEN == "eigen_switched_type":
    chain()
    out["runs"] = [eig("-genBandArpack", 3), eig("-fullGenLapack", 3),
                   eig("-genBandArpack", 3)]
elif SCEN == "eigen_with_analysis":
    chain()
    transient()
    out["runs"] = [eig("-fullGenLapack", 3), eig("-fullGenLapack", 3)]
elif SCEN == "eigen_then_analysis":
    chain()
    first = eig(3)
    transient()
    out["runs"] = [first, eig(3), eig(3)]
elif SCEN == "sparse_flags":
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, x, y in ((1, 0.0, 0.0), (2, 1.0, 0.0), (3, 1.0, 1.0), (4, 0.0, 1.0)):
        ops.node(tag, x, y)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, -1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    out["rc"] = ops.analyze(1)
    forms = {
        "ret_sparse0": ("-ret", "-sparse", 0),
        "sparse_ret": ("-sparse", "-ret"),
        "ret_sparse1": ("-ret", "-sparse", 1),
    }
    got = {}
    for name, args in forms.items():
        try:
            r = ops.printA(*args)
            rows = list(r["rowIndices"])
            got[name] = ["ok", len(rows), min(rows) if rows else None]
        except Exception as ex:
            got[name] = ["raised", repr(ex), None]
    out["forms"] = got

print("RESULT " + json.dumps(out))
'''

TCL_DECK = r'''
set nfail 0
proc check {name ok detail} {
    global nfail
    if {$ok} {
        puts "PASS $name"
    } else {
        incr nfail
        puts "FAIL $name -- $detail"
    }
}

model basic -ndm 2 -ndf 2
node 1 0.0 0.0
node 2 1.0 0.0
node 3 1.0 1.0
node 4 0.0 1.0
fix 1 1 1
fix 2 1 1
nDMaterial ElasticIsotropic 1 1000.0 0.25
element quad 1 1 2 3 4 1.0 PlaneStrain 1
timeSeries Linear 1
pattern Plain 1 1 { load 3 0.0 -1.0 }
constraints Plain
numberer Plain
system FullGeneral
test NormDispIncr 1.0e-10 20 0
algorithm Newton
integrator LoadControl 1.0
analysis Static
analyze 1

set rc [catch {printA -ret -sparse 0} msg]
check ret_sparse0 [expr {$rc == 0}] $msg
set rc [catch {printA -sparse -ret} msg]
check sparse_ret [expr {$rc == 0}] $msg
set rc [catch {printA -ret -sparse 1} msg]
check ret_sparse1 [expr {$rc == 0}] $msg
set rc [catch {printA -sparse} msg]
check bare_sparse [expr {$rc == 0}] $msg

puts "SELF-TEST: $nfail failure(s)"
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
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, f"{scenario}: no result\n{proc.stdout}\n{proc.stderr}"
    out = json.loads(lines[-1][len("RESULT "):])
    assert os.path.normcase(out["engine"]) == os.path.normcase(ENGINE), out
    return out


def _assert_same_modes(runs, tol=REL_TOL):
    for r in runs:
        assert isinstance(r, list) and len(r) == 3, runs
        assert all(x > 0.0 for x in r), runs
    for r in runs[1:]:
        for x, y in zip(runs[0], r):
            assert abs(x - y) <= tol * abs(x), runs


@pytest.mark.parametrize("scenario", ["eigen_same_type", "eigen_explicit_type",
                                      "eigen_with_analysis", "eigen_then_analysis"])
def test_repeated_eigen(scenario):
    _assert_same_modes(_run(scenario)["runs"])


def test_repeated_eigen_switched_type():
    runs = _run("eigen_switched_type")["runs"]
    _assert_same_modes(runs, tol=1.0e-6)
    _assert_same_modes([runs[0], runs[2]])


def test_printa_sparse_option_order_python():
    out = _run("sparse_flags")
    assert out["rc"] == 0, out
    f = out["forms"]
    for name in ("ret_sparse0", "sparse_ret", "ret_sparse1"):
        assert f[name][0] == "ok", out
    assert f["ret_sparse0"][1] == f["sparse_ret"][1] == f["ret_sparse1"][1] > 0, out
    assert f["ret_sparse0"][2] == 0 and f["sparse_ret"][2] == 0, out
    assert f["ret_sparse1"][2] == 1, out


def _tcl_exe():
    exe = os.environ.get("OPENSEES_EXE")
    return exe if exe and os.path.isfile(exe) else None


@pytest.mark.skipif(_tcl_exe() is None,
                    reason="OPENSEES_EXE is unset or does not name a file")
def test_printa_sparse_option_order_tcl(tmp_path):
    deck = tmp_path / "printa_sparse.tcl"
    deck.write_text(TCL_DECK)
    proc = subprocess.run(
        [_tcl_exe(), str(deck)], cwd=str(tmp_path),
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    out = proc.stdout + proc.stderr
    assert "SELF-TEST:" in out, out
    assert out.count("PASS ") + out.count("FAIL ") == 4, out
    assert "SELF-TEST: 0 failure(s)" in out, out
