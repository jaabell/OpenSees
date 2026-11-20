"""constraints LadrunoProjection: conditioning of L^T M L and the
small-rotation note.

A two-retained-DOF equationConstraint (u3 = 0.5 u1 + 0.5 u2) puts both
retained DOFs in one group, so L^T M L is 2x2 with a condition number close
to the mass ratio of the retained DOFs.

Checks:
- equal masses: the analysis runs;
- mass ratio 1e13 (condition number above 1e12): refused, analyze() fails;
- mass ratio 1e10 (condition number between 1e8 and 1e12): a note is
  printed and the analysis runs;
- rigidLink -beam with the retained node rotated past 0.1 rad: the note on
  the frozen small-rotation lever arm is printed; at a rotation well below
  0.1 rad it is not.
"""

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


PROJ = "LadrunoProjection"
CDL = "CentralDifferenceLadruno"


def _analysis():
    ops.constraints(PROJ)
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")


def _run_eq_group(m1, m2, m3, nsteps=3, dt=1.0e-3):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.node(3, 2.0); ops.mass(3, m3)
    ops.equationConstraint(3, 1, 1.0, 1, 1, -0.5, 2, 1, -0.5)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 1.0)
    _analysis()
    rc = 0
    for _ in range(nsteps):
        rc = ops.analyze(1, dt)
        if rc != 0:
            break
    return rc


def test_well_conditioned_group_runs():
    assert _run_eq_group(1.0, 1.0, 1.0) == 0


def test_near_singular_group_refused():
    assert _run_eq_group(1.0e13, 1.0, 1.0) != 0


def test_ill_conditioned_group_runs():
    assert _run_eq_group(1.0e10, 1.0, 1.0) == 0


# Output is captured with capfd (file descriptor level): ops.logFile would
# redirect opserr for the rest of the process.
def _run_rigidlink_beam(moment, nsteps, dt):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.mass(2, 1.0, 1.0, 1.0)
    ops.rigidLink("beam", 1, 2)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, moment)
    _analysis()
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
    return abs(ops.nodeDisp(1, 3))


def test_note_on_large_rotation(capfd):
    rz = _run_rigidlink_beam(50.0, 400, 5.0e-3)
    out = capfd.readouterr()
    txt = (out.err + out.out).lower()
    assert rz > 0.1, rz
    assert "small rotation" in txt and "has rotated" in txt, txt[-800:]


def test_no_note_on_small_rotation(capfd):
    rz = _run_rigidlink_beam(1.0e-4, 20, 1.0e-3)
    out = capfd.readouterr()
    txt = (out.err + out.out).lower()
    assert rz < 0.1, rz
    assert "has rotated" not in txt, txt
