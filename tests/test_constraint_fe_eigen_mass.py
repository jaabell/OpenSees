"""Constraint FE_Elements must not contribute mass to eigen or modalDamping.

Model (1D, ndf 1): a truss of axial stiffness k = E*A/L = 3000*10/10 = 3000
from a fixed node 0 to node 1 (mass 5). Node 2 (mass 3) sits at node 1 and is
tied to it by equalDOF 1 2 1. Twelve independent oscillators (k = 3000 to a
fixed node, mass 3/j for j = 1..12, eigenvalue 1000*j) keep the eigen problem
well above the ARPACK subspace size nev + 8 = 10 under every handler (13
equations under Transformation, more under the others). Closed-form
eigenvalues of the tied system:

    k / (5 + 3) = 375,   1000,   2000,   ...,   12000

A static Newton step with a load F on node 2 stretches the penalty spring of
the tie (Penalty, Auto) or loads the multiplier (Lagrange), so the constraint
FE_Element holds a residual of magnitude F when eigen is called. The eigen at
the initial state (zero residual) is a control that is right before and after
the fix.

Each case runs in a child interpreter, so global state and a crash stay
isolated. Penalty and Auto use a penalty of 1e11: the penalty approximation
error of the first eigenvalue is k*m2/(alpha*(m1+m2)) = 1.1e-8, inside the
1e-6 tolerance. Lagrange is exact. Lagrange is not combined with
modalDamping: IncrementalIntegrator::setupModal asks every DOF_Group for
eigenvectors, and the multiplier DOF_Group has no node.
"""

import json
import subprocess
import sys

import pytest

try:
    import opensees as ops  # noqa: F401
except ModuleNotFoundError:
    import openseespy.opensees as ops  # noqa: F401


REL_TOL = 1.0e-6
EIGEN_TIED = [375.0, 1000.0]
EIGEN_UNTIED = [3000.0 / 5.0, 1000.0]

CHILD = r'''
import json, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

case, handler, solver = sys.argv[1], sys.argv[2], sys.argv[3]
loads = [float(x) for x in sys.argv[4].split(",")]

ALPHA = 1.0e11
HANDLERS = {
    "Transformation": ("Transformation",),
    "Penalty": ("Penalty", ALPHA, ALPHA),
    "Auto": ("Auto", "-userPenalty", ALPHA),
    "Lagrange": ("Lagrange",),
}


def model(tie):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.uniaxialMaterial("Elastic", 1, 3000.0)
    ops.node(0, 0.0)
    ops.fix(0, 1)
    ops.node(1, 10.0)
    ops.mass(1, 5.0)
    ops.element("Truss", 1, 0, 1, 10.0, 1)
    if tie:
        ops.node(2, 10.0)
        ops.mass(2, 3.0)
        ops.equalDOF(1, 2, 1)
    for j in range(1, 13):
        i, m = 10 + j, 3.0 / j
        ops.node(100 + i, 0.0)
        ops.fix(100 + i, 1)
        ops.node(i, 10.0)
        ops.mass(i, m)
        ops.element("Truss", i, 100 + i, i, 10.0, 1)


def static(tie, load):
    ops.constraints(*HANDLERS[handler])
    ops.numberer("RCM")
    ops.system("BandGeneral")
    # displacement 0.33e-3*F; the tolerance scales with it
    ops.test("NormDispIncr", 1.0e-8 * load, 10)
    ops.algorithm("Newton")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2 if tie else 1, load)
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def eig(nev):
    try:
        if solver == "default":
            return list(ops.eigen(nev))
        return list(ops.eigen(solver, nev))
    except Exception as ex:
        return "raised: " + repr(ex)


out = {"runs": []}
if case in ("eigen", "eigen_sp"):
    tie = case == "eigen"
    for load in loads:
        model(tie)
        static(tie, load)
        run = {"load": load, "k0": eig(2)}
        run["ok"] = ops.analyze(1)
        run["eig"] = eig(2)
        out["runs"].append(run)
elif case in ("modal", "undamped"):
    model(True)
    static(True, loads[0])
    out["k0"] = eig(2)
    out["ok"] = ops.analyze(1)
    out["eig"] = eig(2)
    ops.loadConst("-time", 0.0)
    u0 = ops.nodeDisp(1, 1)
    if case == "modal":
        ops.modalDamping(0.05)
    ops.timeSeries("Trig", 2, 0.0, 1.0, 0.5, "-factor", 200.0)
    ops.pattern("UniformExcitation", 3, 1, "-accel", 2)
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    hist = []
    for _ in range(200):
        if ops.analyze(1, 0.005) != 0:
            break
        hist.append(ops.nodeDisp(1, 1) - u0)
    out["hist"] = hist
print("RESULT " + json.dumps(out), flush=True)
'''


def _run(case, handler, loads, solver="default"):
    arg = ",".join(repr(float(x)) for x in loads)
    proc = subprocess.run(
        [sys.executable, "-c", CHILD, case, handler, solver, arg],
        stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode == 0, (
        f"{case} {handler} {solver}: child exited with {proc.returncode}\n"
        f"{proc.stderr[-2000:]}")
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, f"no result\n{proc.stdout}\n{proc.stderr[-2000:]}"
    return json.loads(lines[-1][len("RESULT "):])


def _close(got, ref, tol=REL_TOL):
    return (isinstance(got, list) and len(got) == len(ref)
            and all(abs(g - r) <= tol * abs(r) for g, r in zip(got, ref)))


def _check_eigen_runs(r, tol=REL_TOL, ref=EIGEN_TIED):
    assert r["runs"], r
    for run in r["runs"]:
        assert _close(run["k0"], ref, tol), run
        assert run["ok"] == 0, run
        assert _close(run["eig"], ref, tol), run


# ARPACK reads M*x through getM_Force whenever the assembled M is not
# diagonal, which the constraint FE makes it. Before the fix the residual of
# the last step entered M*x as a constant vector of size F. The result of a
# polluted Lanczos run varies irregularly with F, hence several load levels.
# Under Lagrange the pollution is small, hence the large load and the tight
# tolerance (the fixed result is exact).
ARPACK_CASES = {
    "Penalty": ([1.0e3, 1.0e4, 1.0e6], REL_TOL),
    "Auto": ([1.0e3, 1.0e4, 1.0e6], REL_TOL),
    "Lagrange": ([1.0e8], 1.0e-9),
}


@pytest.mark.parametrize("handler", list(ARPACK_CASES))
def test_arpack_eigen_after_newton_step(handler):
    """ARPACK eigenvalues after a Newton step equal the closed form 375, 1000."""
    loads, tol = ARPACK_CASES[handler]
    _check_eigen_runs(_run("eigen", handler, loads), tol)


@pytest.mark.parametrize("handler", ["Penalty", "Auto"])
def test_modal_damping_after_newton_step(handler):
    """A modalDamping transient built on the eigen after a Newton preload
    (F = 1e3; eigen at the initial state first, as in the eigen test)
    matches the Transformation handler history of node 1 (200 steps of
    0.005 s, measured from the preloaded position). The undamped run checks
    that the damping is active (it differs by more than 1 % of the peak)."""
    ref = _run("modal", "Transformation", [1.0e3])
    und = _run("undamped", "Transformation", [1.0e3])
    got = _run("modal", handler, [1.0e3])
    for r in (ref, und, got):
        assert r["ok"] == 0, r
        assert len(r["hist"]) == 200, r
    assert _close(ref["eig"], EIGEN_TIED), ref
    peak = max(abs(x) for x in ref["hist"])
    assert max(abs(a - b) for a, b in zip(und["hist"], ref["hist"])) > 0.01 * peak
    err = max(abs(a - b) for a, b in zip(got["hist"], ref["hist"]))
    assert err <= 1.0e-5 * peak, (handler, err, peak, got["eig"])
    assert _close(got["eig"], EIGEN_TIED), got


# The eigen mass loops added the tangent of every FE_Element to M. A
# constraint FE returns its constraint stiffness whatever it is asked for,
# so M received alpha*C^T*C (Penalty, Auto), alpha at every Penalty SP dof,
# and the [0 C^T; C 0] coupling block (Lagrange). The explicit-M solver then
# returns spurious modes near +1 or -1 already at the initial state.
@pytest.mark.parametrize("handler", ["Penalty", "Auto", "Lagrange"])
def test_full_lapack_eigen_constraint_has_no_mass(handler):
    """-fullGenLapack eigenvalues at the initial state and after a Newton
    step (F = 1e3) equal the closed form 375, 1000."""
    _check_eigen_runs(_run("eigen", handler, [1.0e3], "-fullGenLapack"))


@pytest.mark.parametrize("solver", ["default", "-fullGenLapack"])
def test_penalty_sp_only_eigen(solver):
    """Without the tie (fix only, Penalty handler) the eigenvalues are
    3000/5 = 600 and 1000. Before the fix PenaltySP_FE put alpha on the
    diagonal of M at every fixed dof, which also keeps M diagonal for ARPACK:
    each fixed dof became a spurious low mode (eigenvalue alpha/alpha = 1)."""
    _check_eigen_runs(_run("eigen_sp", "Penalty", [1.0e3], solver),
                      ref=EIGEN_UNTIED)
