"""Interpreter fixes: Tcl string arguments and Python module re-import.

  tcl_string_from_all   OPS_GetStringFromAll in the classic Tcl interpreter
                        fills the buffer the caller passes in
  reimport              the Python module removed from sys.modules and
                        imported again, then a normal interpreter exit

KinematicCoupling reads its -k value and its -dof components through
OPS_GetStringFromAll into a local buffer. Reference values, 2D model with
the reference node fixed and two free slaves tied in x and y:

  -k 1000.0        eleResponse kt = 1000.0; a unit load in x on a slave,
                   whose only stiffness is the tie, gives u = P/Kt = 1.0e-3
  -dof 1 2         slaves with 3 DOFs, rotation left free: 2 x 2 = 4 tied
                   DOFs (6 with the default component list)

Before the fix the buffer is never written under Tcl, so the element parser
reads uninitialised memory and rejects the option; the element command fails.

The Python module declares a nonzero module state size, so every import
after removal from sys.modules runs its init function again. Before the fix
each init registers the exit handler with Py_AtExit, and the handler
dereferences and deletes the static module object without clearing it: the
second call at interpreter exit wipes and deletes a freed object, and the
child process exits with an access violation or heap corruption status
instead of 0. The child re-imports 40 times (more than the 32 Py_AtExit
slots), solves an elastic truss, u = P L / (E A) = 10 * 2 / (1000 * 4) =
5.0e-3, and must exit with code 0. A single import is checked as a control.

The Tcl case needs the OpenSees executable named by the OPENSEES_EXE
environment variable and is skipped when it is unset or missing.
"""
import importlib.machinery
import json
import os
import subprocess
import sys

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _extension_module_name():
    """Name of the compiled OpenSees module in sys.modules, or None."""
    for name, mod in list(sys.modules.items()):
        f = getattr(mod, "__file__", None) or ""
        if (hasattr(mod, "wipe") and hasattr(mod, "analyze")
                and f.endswith(tuple(importlib.machinery.EXTENSION_SUFFIXES))):
            return name
    return None


EXT_NAME = _extension_module_name()

CHILD = r'''
import importlib, json, sys
name = sys.argv[1]
count = int(sys.argv[2])
ops = importlib.import_module(name)
for _ in range(count):
    del sys.modules[name]
    ops = importlib.import_module(name)
ops.wipe()
ops.model("basic", "-ndm", 1, "-ndf", 1)
ops.node(1, 0.0)
ops.node(2, 2.0)
ops.fix(1, 1)
ops.uniaxialMaterial("Elastic", 1, 1000.0)
ops.element("Truss", 1, 1, 2, 4.0, 1)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, 10.0)
ops.constraints("Plain")
ops.numberer("Plain")
ops.system("FullGeneral")
ops.algorithm("Linear")
ops.integrator("LoadControl", 1.0)
ops.analysis("Static")
rc = ops.analyze(1)
print("RESULT " + json.dumps({"rc": rc, "u": ops.nodeDisp(2, 1)}))
sys.stdout.flush()
'''


@pytest.mark.skipif(EXT_NAME is None,
                    reason="the OpenSees module is not a compiled extension")
@pytest.mark.parametrize("count", [0, 40])
def test_python_reimport_clean_exit(count):
    proc = subprocess.run(
        [sys.executable, "-c", CHILD, EXT_NAME, str(count)],
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    detail = (f"{count} re-imports: child exited with {proc.returncode}\n"
              f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
    assert lines, detail
    out = json.loads(lines[-1][len("RESULT "):])
    assert out["rc"] == 0, out
    assert out["u"] == pytest.approx(5.0e-3, rel=1.0e-12), out
    assert proc.returncode == 0, detail


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

# -k value read through OPS_GetStringFromAll
model basic -ndm 2 -ndf 2
node 1 0.0 0.0 -ndf 3
node 2 1.0 0.0
node 3 0.0 1.0
fix 1 1 1 1
set rc [catch {element KinematicCoupling 1 1 2 2 3 -k 1000.0} msg]
check k_element [expr {$rc == 0}] $msg
if {$rc == 0} {
    set kt [eleResponse 1 kt]
    check k_value [expr {abs($kt - 1000.0) <= 1.0e-9}] "kt = $kt"
    timeSeries Linear 1
    pattern Plain 1 1 { load 2 1.0 0.0 }
    constraints Plain
    numberer Plain
    system FullGeneral
    test NormDispIncr 1.0e-12 10 0
    algorithm Newton
    integrator LoadControl 1.0
    analysis Static
    set arc [analyze 1]
    set ux [nodeDisp 2 1]
    check k_disp [expr {$arc == 0 && abs($ux - 1.0e-3) <= 1.0e-12}] "analyze $arc, ux = $ux"
} else {
    check k_value 0 "no element"
    check k_disp 0 "no element"
}

# -dof components read through OPS_GetStringFromAll
wipe
model basic -ndm 2 -ndf 3
node 1 0.0 0.0
node 2 1.0 0.0
node 3 0.0 1.0
fix 1 1 1 1
set rc [catch {element KinematicCoupling 1 1 2 2 3 -dof 1 2 -k 1000.0} msg]
check dof_element [expr {$rc == 0}] $msg
if {$rc == 0} {
    set nt [eleResponse 1 tiedDOFs]
    check dof_count [expr {abs($nt - 4.0) <= 1.0e-12}] "tiedDOFs = $nt"
} else {
    check dof_count 0 "no element"
}

puts "SELF-TEST: $nfail failure(s)"
'''


def _tcl_exe():
    exe = os.environ.get("OPENSEES_EXE")
    return exe if exe and os.path.isfile(exe) else None


@pytest.mark.skipif(_tcl_exe() is None,
                    reason="OPENSEES_EXE is unset or does not name a file")
def test_tcl_string_from_all_fills_buffer(tmp_path):
    deck = tmp_path / "string_from_all.tcl"
    deck.write_text(TCL_DECK)
    proc = subprocess.run(
        [_tcl_exe(), str(deck)], cwd=str(tmp_path),
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    out = proc.stdout + proc.stderr
    assert "SELF-TEST:" in out, out
    assert out.count("PASS ") + out.count("FAIL ") == 5, out
    assert "SELF-TEST: 0 failure(s)" in out, out
