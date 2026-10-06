"""Arpack eigen solvers: unconverged modes and a model that grows.

Each case runs in a child interpreter, so a crash fails that case only.
The reference eigenvalues come from -fullGenLapack on the same model.
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


REL_TOL = 1.0e-8

CHILD = r'''
import json, sys
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

case, solver, nev = sys.argv[1], sys.argv[2], int(sys.argv[3])


def add_bays(first, last):
    # truss chain along x; each node also has a lateral spring in y.
    # Equal springs and masses give one eigenvalue of multiplicity
    # (last - first) at k/m = 240.
    for i in range(first, last):
        ops.node(i, 10.0 * i, 0.0)
        ops.mass(i, 5.0, 5.0)
        ops.element("Truss", i, i - 1, i, 10.0, 1)
        ops.node(100000 + i, 10.0 * i, 10.0)
        ops.fix(100000 + i, 1, 1)
        ops.element("Truss", 100000 + i, 100000 + i, i, 4.0, 1)


def model(nbays):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 3000.0)
    ops.node(0, 0.0, 0.0)
    ops.fix(0, 1, 1)
    add_bays(1, nbays + 1)


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


out = {}
if case == "cluster":
    model(10)
    transient()
    out["arpack"] = eig(*([solver] if solver != "default" else []), nev)
    out["lapack"] = eig("-fullGenLapack", nev)
elif case == "grow":
    model(10)
    transient()
    out["before"] = eig(*([solver] if solver != "default" else []), nev)
    add_bays(11, 200)
    out["arpack"] = eig(*([solver] if solver != "default" else []), nev)
    out["lapack"] = eig("-fullGenLapack", nev)
print("RESULT " + json.dumps(out), flush=True)
'''


def _run(case, solver, nev):
    proc = subprocess.run(
        [sys.executable, "-c", CHILD, case, solver, str(nev)],
        stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=600,
    )
    assert proc.returncode == 0, (
        f"{case} {solver} {nev}: child exited with {proc.returncode}\n"
        f"{proc.stderr[-2000:]}")
    lines = [l for l in proc.stdout.splitlines() if l.startswith("RESULT ")]
    assert lines, f"no result\n{proc.stdout}\n{proc.stderr[-2000:]}"
    return json.loads(lines[-1][len("RESULT "):])


def _same(a, b):
    return (isinstance(a, list) and isinstance(b, list) and len(a) == len(b)
            and all(abs(x - y) <= REL_TOL * abs(y) for x, y in zip(a, b)))


@pytest.mark.parametrize("solver", ["default", "-genBandArpack"])
def test_unconverged_modes_are_not_returned(solver):
    """Arpack either returns the correct eigenvalues or reports failure;
    it never returns values it has not computed."""
    r = _run("cluster", solver, 3)
    assert isinstance(r["lapack"], list)
    if isinstance(r["arpack"], list):
        assert _same(r["arpack"], r["lapack"]), r


@pytest.mark.parametrize("solver", ["default", "-genBandArpack"])
def test_eigen_after_model_grows(solver):
    """eigen after nodes and elements are added without wipe."""
    r = _run("grow", solver, 4)
    assert isinstance(r["lapack"], list)
    assert _same(r["arpack"], r["lapack"]), r


@pytest.mark.parametrize("nev", [3, 6])
def test_repeated_eigenvalue_copies(nev):
    """The third to tenth eigenvalues are equal (240). Arpack must return
    every requested copy, not skip to the next distinct values."""
    r = _run("cluster", "default", nev)
    assert _same(r["arpack"], r["lapack"]), r
