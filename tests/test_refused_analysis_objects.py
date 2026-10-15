"""A refused system, numberer, constraints, test or integrator command is an
error and leaves the previous object in force.

Each case defines a complete set of analysis objects, issues one command that
creates nothing (an unknown integrator name, or a known name whose factory
rejects its options), and then rebuilds the analysis from the stored objects.

Model: one truss, E = 1, A = 1, L = 1, load P = 1 at node 2, so k = EA/L = 1.
The previous objects are LoadControl 0.1 and FixedNumIter 3. One step must
give t = 0.1, u2 = 0.1 * P / k = 0.1 and testIter = 3. A stored object that was
nulled is replaced by a default when the analysis is rebuilt: the analysis
command prints "no LinearSOE specified" (or the numberer, handler, integrator
equivalent), the integrator default LoadControl 1.0 gives t = 1.0, and the
test default NormUnbalance 1e-6 25 converges this linear model in 1 iteration.

Before the fix the refused commands returned 0, so OpenSeesPy raised nothing,
and system, numberer, constraints and test installed a null object.

The Tcl cases need the OpenSees executable named by the OPENSEES_EXE
environment variable and are skipped when it is unset or missing. They use
the same model; the transient cases use Newmark 0.5 0.25 with a unit mass at
node 2, and three steps of 0.01 give t = 0.03.
"""
import os
import subprocess

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


DT = 0.1
DEFAULT_WARNINGS = (
    "no LinearSOE specified",
    "no Numberer specified",
    "no ConstraintHandler yet specified",
    "no Integrator specified",
)


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("FixedNumIter", 3)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", DT)


# (command, arguments for which no object is created)
REFUSED = [
    pytest.param("system", ("PFEM", "-noSuchOption"), id="system-PFEM-unknown-option"),
    pytest.param("numberer", ("ParallelRCM",), id="numberer-ParallelRCM-serial"),
    pytest.param("constraints", ("Penalty",), id="constraints-Penalty-no-alpha"),
    pytest.param("test", ("NormUnbalance",), id="test-NormUnbalance-no-tol"),
    pytest.param("integrator", ("LoadControl",), id="integrator-LoadControl-no-dlam"),
    pytest.param("integrator", ("NoSuchIntegrator", 1.0), id="integrator-unknown-name"),
]


@pytest.mark.parametrize("cmd,refused", REFUSED)
@pytest.mark.parametrize("analysis_first", [False, True],
                         ids=["no-analysis", "with-analysis"])
def test_refused_object_raises_and_keeps_previous(cmd, refused, analysis_first, capfd):
    try:
        _model()
        if analysis_first:
            ops.analysis("Static")
        capfd.readouterr()
        with pytest.raises(Exception):
            getattr(ops, cmd)(*refused)
        assert "left unchanged" in capfd.readouterr().err

        ops.analysis("Static")
        err = capfd.readouterr().err
        for warning in DEFAULT_WARNINGS:
            assert warning not in err, err
        assert ops.analyze(1) == 0
        assert ops.getTime() == pytest.approx(DT, rel=1e-12)
        assert ops.nodeDisp(2, 1) == pytest.approx(DT, rel=1e-12)
        assert ops.testIter() == 3
    finally:
        ops.wipe()


# --- Tcl interpreter ----------------------------------------------------------

TCL_MODEL = """
model basic -ndm 1 -ndf 1
node 1 0.0
node 2 1.0
fix 1 1
mass 2 1.0
uniaxialMaterial Elastic 1 1.0
element Truss 1 1 2 1.0 1
timeSeries Linear 1
pattern Plain 1 1 { load 2 1.0 }
constraints Transformation
numberer Plain
system FullGeneral
test FixedNumIter 3
algorithm Newton
"""

STATIC = ("integrator LoadControl 0.1\n", "analysis Static", "analyze 1", 0.1)
TRANSIENT = ("integrator Newmark 0.5 0.25\n", "analysis Transient", "analyze 3 0.01", 0.03)

# name -> (analysis setup, refused command, analysis defined before it)
TCL_CASES = {
    "system-unknown-with-analysis": (STATIC, "system NoSuchSystem", True),
    "system-PFEM-unknown-option": (STATIC, "system PFEM -noSuchOption", False),
    "constraints-Auto-missing-value": (STATIC, "constraints Auto -autoPenalty", False),
    "integrator-Newmark-no-args": (TRANSIENT, "integrator Newmark", False),
    "integrator-Newmark-no-args-with-analysis": (TRANSIENT, "integrator Newmark", True),
}


def _tcl_exe():
    exe = os.environ.get("OPENSEES_EXE")
    return exe if exe and os.path.isfile(exe) else None


@pytest.mark.skipif(_tcl_exe() is None,
                    reason="OPENSEES_EXE is unset or does not name a file")
@pytest.mark.parametrize("case", list(TCL_CASES))
def test_refused_object_tcl(case, tmp_path):
    """Before the fix system NoSuchSystem with an analysis present returned
    TCL_OK and deleted the system the analysis still used (the next analyze
    crashed), and integrator Newmark without arguments with an analysis
    present dereferenced a null integrator (segmentation fault)."""
    (integrator, analysis, analyze, t_expected), refused, analysis_first = TCL_CASES[case]
    deck = tmp_path / "refused.tcl"
    deck.write_text(
        TCL_MODEL + integrator
        + (analysis + "\n" if analysis_first else "")
        + "puts \"RC [catch {" + refused + "} msg]\"\n"
        + analysis + "\n" + analyze + "\n"
        + "puts \"TIME [getTime]\"\n"
        + "puts \"ITER [testIter]\"\n"
        + "puts DECK-DONE\n")
    proc = subprocess.run(
        [_tcl_exe(), str(deck)], cwd=str(tmp_path),
        stdin=subprocess.DEVNULL, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=300,
    )
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, f"exit code {proc.returncode}\n{out}"
    assert "DECK-DONE" in out, out
    assert "RC 1" in out, out
    for warning in DEFAULT_WARNINGS:
        assert warning not in out, out
    lines = {l.split()[0]: l.split()[1] for l in out.splitlines()
             if l.split() and l.split()[0] in ("TIME", "ITER")}
    assert float(lines["TIME"]) == pytest.approx(t_expected, rel=1e-12), out
    assert int(lines["ITER"]) == 3, out
