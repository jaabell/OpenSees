"""Prescribed displacements under constraints Auto, sp -subtractInit, and Path -useLast.

  auto_sp        two-bar chain, node 1 fixed, node 3 driven to u3 = 3.0 by an sp,
                 Newton with NormDispIncr. Equal bars give the closed form
                 u2 = u3 / 2 = 1.5 for every constraint handler.
  subtract_init  same chain, a first stage with P = 40 at node 3 held constant
                 (k = E*A/L = 100, u2 = 0.4, u3 = 0.8), then sp 3 1 1.0
                 -subtractInit. The prescribed value is incremental, so
                 u3 = 0.8 + 1.0 = 1.8 and u2 = u3 / 2 = 0.9. Without the
                 flag taking effect u3 = 1.0 and u2 = 0.5.
  path_use_last  Path series with points (0, 0) and (1, 2) and -useLast,
                 given as -time/-values, as -fileTime/-filePath, and as
                 -dt/-values. At t = 2, past the last point, the load factor
                 is the last value 2.0, so P = 10 at node 3 gives
                 u3 = 2*P*L/(E*A)*2 = 0.4. Without -useLast reaching the
                 series the factor is 0.0 and u3 = 0.0. The -dt form already
                 forwarded the flag and is a control.

Each case builds its own model after ops.wipe().
"""
import os

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E = 100.0
A = 1.0
L = 1.0


def bar_chain():
    # nodes 1-2-3 on the x axis, two equal trusses, y DOFs fixed
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.node(3, 2.0 * L, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(3, 0, 1)
    ops.element("Truss", 1, 1, 2, A, 1)
    ops.element("Truss", 2, 2, 3, A, 1)


def static_analysis(handler, test_type="NormDispIncr"):
    if handler == "Penalty":
        ops.constraints("Penalty", 1.0e12, 1.0e12)
    else:
        ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test(test_type, 1.0e-10, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


@pytest.mark.parametrize("handler", ["Auto", "Transformation", "Penalty"])
@pytest.mark.parametrize("test_type", ["NormDispIncr", "EnergyIncr"])
def test_auto_sp_prescribed_displacement(handler, test_type):
    """u3 = 3.0 prescribed; u2 = 1.5 (equal bars). Auto returned u2 = 0.0."""
    bar_chain()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(3, 1, 3.0)
    static_analysis(handler, test_type)
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(3, 1) == pytest.approx(3.0, rel=1.0e-6)
    assert ops.nodeDisp(2, 1) == pytest.approx(1.5, rel=1.0e-6)


@pytest.mark.parametrize("handler", ["Transformation", "Penalty"])
def test_sp_subtract_init_is_incremental(handler):
    """u3 = 0.8 + 1.0 = 1.8 and u2 = 0.9. The flag was ignored: u3 = 1.0."""
    bar_chain()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 40.0, 0.0)
    static_analysis(handler)
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(3, 1) == pytest.approx(0.8, rel=1.0e-6)
    ops.loadConst("-time", 0.0)

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.sp(3, 1, 1.0, "-subtractInit")
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(3, 1) == pytest.approx(1.8, rel=1.0e-6)
    assert ops.nodeDisp(2, 1) == pytest.approx(0.9, rel=1.0e-6)


def path_series(form, tmp_path):
    if form == "time_values":
        ops.timeSeries("Path", 1, "-time", 0.0, 1.0,
                       "-values", 0.0, 2.0, "-useLast")
    elif form == "files":
        ftime = os.path.join(str(tmp_path), "time.txt")
        fpath = os.path.join(str(tmp_path), "path.txt")
        with open(ftime, "w") as f:
            f.write("0.0\n1.0\n")
        with open(fpath, "w") as f:
            f.write("0.0\n2.0\n")
        ops.timeSeries("Path", 1, "-fileTime", ftime, "-filePath", fpath,
                       "-useLast")
    else:
        ops.timeSeries("Path", 1, "-dt", 1.0, "-values", 0.0, 2.0,
                       "-useLast")


@pytest.mark.parametrize("form", ["time_values", "files", "dt"])
def test_path_series_use_last(form, tmp_path):
    """At t = 2 the factor holds the last value 2.0 and u3 = 0.4 (was 0.0)."""
    bar_chain()
    path_series(form, tmp_path)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 10.0, 0.0)
    static_analysis("Transformation")
    assert ops.analyze(2) == 0
    assert ops.getTime() == pytest.approx(2.0)
    assert ops.getLoadFactor(1) == pytest.approx(2.0, rel=1.0e-12)
    assert ops.nodeDisp(3, 1) == pytest.approx(0.4, rel=1.0e-6)
