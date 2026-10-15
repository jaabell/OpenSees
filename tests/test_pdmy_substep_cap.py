"""PressureDependMultiYield: bounded substep count for huge trial strains.

setSubStrainRate() sizes the substep loop of getStress() as the octahedral
shear strain increment over 1e-4 (PressureDependMultiYield) or 1e-5
(PressureDependMultiYield02 and 03), or the volumetric strain increment over
1e-5, converted to int with no upper bound.

Model: one 1 x 1 plane-strain quad, rollers on the bottom and left edges,
top nodes tied in y, isotropic confinement P0 = 100 by nodal loads in stage
0 (elastic), then stage 1 at the same loads. A second pattern prescribes a
vertical displacement increment -D on the top (penalty constraints, so the
increment enters through the first Newton iterate). That iterate is a
uniform strain increment eps_yy = -D, eps_xx = gamma_xy = 0, so

    octahedralShear(1) = 2/sqrt(3) * sqrt(2/3) * D = 0.9428 * D.

huge   D = 2.0e4: octahedral shear 18856, i.e. 1.886e8 substeps per Gauss
       point and getStress() call for PressureDependMultiYield and 1.886e9
       (below INT_MAX) for PressureDependMultiYield02/03. Before the fix the
       formUnbalance() after that update runs this loop at 4 Gauss points,
       about a minute even at 0.1 microseconds per substep, so the child is
       still running at the 20 s timeout. With the fix setTrialStrain()
       returns -1 (limit 1e5 substeps), the quad propagates it, analyze()
       returns a negative code in well under a second and the domain is
       reverted: the top node is back at its committed position.
small  D = 5.0e-4: 4 (PressureDependMultiYield) or 47 substeps from the
       strain criteria, far below the limit, on a code path the fix does not
       change. The step converges with the top node at u_y - D (penalty
       error below 1e-9), the mean lateral stress equal to -P0 (lateral
       equilibrium) and the mean vertical stress equal to the value of the
       unmodified code (relative tolerance 1e-4).

Each case runs in a child interpreter under a wall-clock timeout so that the
unbounded loop fails the test instead of hanging the session.
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
TIMEOUT = 20.0
P0 = 100.0

# mean vertical stress after the small step, unmodified code
SYY_SMALL = {
    "PressureDependMultiYield": -239.23405793931204,
    "PressureDependMultiYield02": -241.29693501184508,
    "PressureDependMultiYield03": -241.28999251703226,
}

CHILD = r'''
import json, os, sys, time
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

MAT, D, P0 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
PARAMS = {
    # nd rho G B phi gammaPeak refP d PTAng contrac dilat1 dilat2
    # liquefac1 liquefac2 liquefac3 NYS e cs1 cs2 cs3 pa c
    "PressureDependMultiYield": [2, 2.0, 1.3e5, 2.6e5, 40.0, 0.1, 101.0, 0.5,
                                 26.0, 0.013, 0.3, 3.0, 5.0, 0.0, 1.0,
                                 20, 0.6, 0.9, 0.02, 0.7, 101.0, 0.1],
    # nd rho G B phi gammaPeak refP d PTAng contrac1 contrac3 dilat1 dilat3
    # NYS contrac2 dilat2 liquefac1 liquefac2 e cs1 cs2 cs3 pa c
    "PressureDependMultiYield02": [2, 2.0, 1.3e5, 2.6e5, 40.0, 0.1, 101.0,
                                   0.5, 26.0, 0.013, 0.0, 0.3, 0.0,
                                   20, 5.0, 3.0, 1.0, 0.0, 0.6, 0.9, 0.02,
                                   0.7, 101.0, 0.1],
    # nd rho G B phi gammaPeak refP d PTAng mType ca cb cc cd ce da db dc
    # NYS liquefac1 liquefac2 pa c
    "PressureDependMultiYield03": [2, 2.0, 1.3e5, 2.6e5, 40.0, 0.1, 101.0,
                                   0.5, 26.0, 0, 0.013, 0.0, 0.3, 0.0, 0.0,
                                   0.3, 3.0, 0.0, 20, 1.0, 0.0, 101.0, 0.1],
}

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.nDMaterial(MAT, 1, *PARAMS[MAT])
ops.node(1, 0.0, 0.0)
ops.node(2, 1.0, 0.0)
ops.node(3, 1.0, 1.0)
ops.node(4, 0.0, 1.0)
ops.fix(1, 1, 1)
ops.fix(2, 0, 1)
ops.fix(4, 1, 0)
ops.equalDOF(3, 4, 2)
ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, -0.5 * P0, 0.0)
ops.load(3, -0.5 * P0, -0.5 * P0)
ops.load(4, 0.0, -0.5 * P0)
ops.constraints("Penalty", 1.0e12, 1.0e12)
ops.numberer("Plain")
ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-9, 100)
ops.algorithm("KrylovNewton")
ops.integrator("LoadControl", 0.1)
ops.analysis("Static")
rc0 = ops.analyze(10)
ops.loadConst("-time", 0.0)
ops.integrator("LoadControl", 0.0)
ops.updateMaterialStage("-material", 1, "-stage", 1)
rc1 = ops.analyze(1)
ops.loadConst("-time", 0.0)

uy = ops.nodeDisp(3, 2)
ops.timeSeries("Constant", 2)
ops.pattern("Plain", 2, 2)
ops.sp(3, 2, uy - D)
ops.integrator("LoadControl", 1.0)
ops.analysis("Static")

t0 = time.perf_counter()
rc = ops.analyze(1)
dt = time.perf_counter() - t0

out = {"engine": os.path.abspath(ops.__file__),
       "rc0": rc0, "rc1": rc1, "rc": rc, "dt": dt, "uy": uy,
       "uy_after": ops.nodeDisp(3, 2),
       "stress": list(ops.eleResponse(1, "stresses"))}
ops.wipe()
print("RESULT " + json.dumps(out), flush=True)
'''


def _run(mat, d):
    try:
        proc = subprocess.run(
            [sys.executable, "-c", CHILD, mat, repr(d), repr(P0)],
            stdin=subprocess.DEVNULL, capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        pytest.fail(f"{mat}: analyze() still running after {TIMEOUT:.0f} s "
                    f"for a trial strain increment of {d}")
    assert proc.returncode == 0, (
        f"{mat}: child exited with {proc.returncode}\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, f"{mat}: no result\n{proc.stdout}\n{proc.stderr}"
    out = json.loads(lines[-1][len("RESULT "):])
    assert os.path.normcase(out["engine"]) == os.path.normcase(ENGINE), out
    assert out["rc0"] == 0 and out["rc1"] == 0, out
    return out, proc.stdout + proc.stderr


@pytest.mark.parametrize("mat", sorted(SYY_SMALL))
def test_huge_trial_strain_fails_fast(mat):
    out, log = _run(mat, 2.0e4)
    assert out["rc"] < 0, out
    assert out["dt"] < 1.0, out
    assert out["uy_after"] == out["uy"], out
    assert "substeps; increment rejected" in log, log


@pytest.mark.parametrize("mat", sorted(SYY_SMALL))
def test_small_increment_unchanged(mat):
    d = 5.0e-4
    out, log = _run(mat, d)
    assert out["rc"] == 0, out
    assert out["uy_after"] == pytest.approx(out["uy"] - d, rel=0.0, abs=1.0e-8)
    s = out["stress"]
    assert len(s) == 12, out
    assert sum(s[0::3]) / 4.0 == pytest.approx(-P0, rel=1.0e-6), out
    assert sum(s[1::3]) / 4.0 == pytest.approx(SYY_SMALL[mat], rel=1.0e-4), out
    assert "increment rejected" not in log, log
