"""ASDConcrete3D and ASDConcrete1D laws do not survive wipe under a reused tag.

Each material keeps a process-wide registry of the original hardening laws,
keyed by material tag and filled only when the slot is empty. With
-autoRegularization and a characteristic length different from lch_ref,
regularize() copies the stored law back over the instance. Before the registry
was cleared on wipe, the first material defined with tag N owned that tag for
the life of the process, so a later model reusing the tag ran on the earlier
backbones, elastic modulus included.

  3D  one 1 x 1 x 2 stdBrick, rollers, top load P = 7200
  1D  one truss of length L = 2 and area 1, axial compression P = 7200

Model 1 defines the tag with E = 3e6 and is wiped. Model 2 defines the same tag
with E = 30e6. Closed form on model 2: stress = -P/A = -7200,
uz = -P L / (E A) = -4.8e-4, support reaction = P. Both stay inside the elastic
branch of the backbones. Without the fix model 2 runs with E = 3e6.

Python cases run in a child interpreter. The Tcl cases need the OpenSees
executable named by the OPENSEES_EXE environment variable and are skipped when
it is unset or missing.
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


P, L = 7200.0, 2.0
E2 = 30.0e6

CHILD_3D = r'''
import json
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops
P, E, NU = 7200.0, 30.0e6, 0.2
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
ops.nDMaterial("ASDConcrete3D", 1, 3.0e6, 0.2, "-fc", 3.0e3, "-ft", 3.0e2)
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
xy = [(0, 0), (1, 0), (1, 1), (0, 1)]
for i, (x, y) in enumerate(xy):
    ops.node(1 + i, float(x), float(y), 0.0)
    ops.node(5 + i, float(x), float(y), 2.0)
ops.nDMaterial("ASDConcrete3D", 1, E, NU,
    "-Te", 0.0, 9e-05, 0.00015, 0.00507, 0.0250501, 0.250501,
    "-Ts", 0.0, 2700.0, 3000.0, 600.0, 3.0, 3.0,
    "-Td", 0.0, 0.0, 0.0, 0.960552268244576, 0.9999800399998403, 0.9999995660869531,
    "-Ce", 0.0, 0.0005, 0.001, 0.0015, 0.002, 0.18327272727272728, 0.18377272727272728,
    "-Cs", 0.0, 15000.0, 24852.8137423857, 28923.048454132637, 30000.0, 3000.0, 3000.0,
    "-Cd", 0.0, 0.0, 0.0, 0.0, 0.0, 0.9981618744961699, 0.9981786141748574,
    "-autoRegularization", 0.03112529146138255)
ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
ops.fix(1, 1, 1, 1)
ops.fix(2, 0, 1, 1)
ops.fix(3, 0, 0, 1)
ops.fix(4, 1, 0, 1)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
for n in (5, 6, 7, 8):
    ops.load(n, 0.0, 0.0, -P / 4.0)
ops.constraints("Plain")
ops.numberer("RCM")
ops.system("BandGeneral")
ops.test("NormDispIncr", 1e-10, 20, 0)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0)
ops.analysis("Static")
ok = ops.analyze(1)
res = {"ok": ok}
if ok == 0:
    res["szz"] = [ops.eleResponse(1, "material", gp, "stress")[2] for gp in range(1, 9)]
    res["uz"] = [ops.nodeDisp(n, 3) for n in (5, 6, 7, 8)]
    ops.reactions()
    res["rz"] = sum(ops.nodeReaction(n, 3) for n in (1, 2, 3, 4))
print("RESULT " + json.dumps(res))
'''


def _run_child(code):
    proc = subprocess.run(
        [sys.executable, "-c", code], stdin=subprocess.DEVNULL,
        capture_output=True, text=True, encoding="utf-8", errors="replace",
        timeout=120,
    )
    out = proc.stdout + proc.stderr
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, out
    return json.loads(lines[-1][len("RESULT "):]), out


def test_asdconcrete3d_same_tag_after_wipe():
    res, out = _run_child(CHILD_3D)
    assert res["ok"] == 0, out
    for gp, s in enumerate(res["szz"], start=1):
        assert s == pytest.approx(-P, rel=1e-5), f"gauss point {gp}: {s}"
    for uz in res["uz"]:
        assert uz == pytest.approx(-P * L / E2, rel=1e-5)
    assert res["rz"] == pytest.approx(P, rel=1e-5)


TCL_3D = """
set nFail 0
proc check {name got want} {
    global nFail
    if {abs($got - $want) <= 1e-5 * abs($want)} {
        puts "PASS $name"
    } else {
        puts "FAIL $name got=$got want=$want"
        incr nFail
    }
}
wipe
model BasicBuilder -ndm 3 -ndf 3
nDMaterial ASDConcrete3D 1 3.0e6 0.2 -fc 3.0e3 -ft 3.0e2
wipe
model BasicBuilder -ndm 3 -ndf 3
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 2.0
node 6 1.0 0.0 2.0
node 7 1.0 1.0 2.0
node 8 0.0 1.0 2.0
nDMaterial ASDConcrete3D 1 30.0e6 0.2 \\
    -Te 0.0 9e-05 0.00015 0.00507 0.0250501 0.250501 \\
    -Ts 0.0 2700.0 3000.0 600.0 3.0 3.0 \\
    -Td 0.0 0.0 0.0 0.960552268244576 0.9999800399998403 0.9999995660869531 \\
    -Ce 0.0 0.0005 0.001 0.0015 0.002 0.18327272727272728 0.18377272727272728 \\
    -Cs 0.0 15000.0 24852.8137423857 28923.048454132637 30000.0 3000.0 3000.0 \\
    -Cd 0.0 0.0 0.0 0.0 0.0 0.9981618744961699 0.9981786141748574 \\
    -autoRegularization 0.03112529146138255
element stdBrick 1 1 2 3 4 5 6 7 8 1
fix 1 1 1 1
fix 2 0 1 1
fix 3 0 0 1
fix 4 1 0 1
timeSeries Linear 1
pattern Plain 1 1 {
    foreach n {5 6 7 8} { load $n 0.0 0.0 -1800.0 }
}
constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1e-10 20 0
algorithm Newton
integrator LoadControl 1.0
analysis Static
set ok [analyze 1]
if {$ok == 0} { puts "PASS converge" } else { puts "FAIL converge ok=$ok"; incr nFail }
if {$ok == 0} {
    for {set gp 1} {$gp <= 8} {incr gp} {
        check "szz$gp" [lindex [eleResponse 1 material $gp stress] 2] -7200.0
    }
    check "uz5" [nodeDisp 5 3] [expr {-7200.0 * 2.0 / 30.0e6}]
}
puts "WIPE_3D: $nFail failure(s)"
"""


def _tcl_exe():
    exe = os.environ.get("OPENSEES_EXE")
    return exe if exe and os.path.isfile(exe) else None


def _run_tcl(tmp_path, name, text):
    deck = tmp_path / name
    deck.write_text(text)
    proc = subprocess.run(
        [_tcl_exe(), str(deck)], cwd=str(tmp_path),
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    return proc.stdout + proc.stderr


@pytest.mark.skipif(_tcl_exe() is None,
                    reason="OPENSEES_EXE is unset or does not name a file")
def test_asdconcrete3d_same_tag_after_wipe_tcl(tmp_path):
    out = _run_tcl(tmp_path, "wipe_3d.tcl", TCL_3D)
    assert "WIPE_3D: 0 failure(s)" in out, out
    assert out.count("PASS ") == 10, out
