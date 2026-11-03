"""ElasticTimoshenkoBeam3d must use its geometric transformation.

The element formed its own rotation-only matrix from the axes of the
transformation, and its geometric flag came only from -geomNonlinear. A PDelta
transformation gave no chord stiffness, a Corotational transformation was
linear, and -jntOffset rigid arms were ignored. recvSelf reset the flag from the
transformation type, so a restored Linear + -geomNonlinear element lost P-Delta.

All cases are a vertical cantilever of length L along global Z, base at node 1,
vecxz = (1, 0, 0): global X is local z (bending about Iy, shear area Avz) and
global Y is local y (bending about Iz, shear area Avy). Reference values:

  P-Delta      f = L^3/(3 E I) + L/(G Av) is the tip compliance. With compression
               P and tip shear H the chord stiffness of a PDelta transformation
               gives u = H f / (1 - P f / L). The unmodified element gives H f.
  Corotational a tip moment M on one element gives equal and opposite end
               moments, so the shear term cancels: theta = M L / (E I), the chord
               keeps length L and turns by theta/2. Tip displacement
               u = L sin(theta/2), w = L (cos(theta/2) - 1). The unmodified
               element gives w = 0. At theta = 1e-6 one solve with the initial
               stiffness gives u = L theta / 2.
  Rigid arms   -jntOffset is the same system as an element between the arm ends,
               tied to the nodes with rigidLink beam.
  Element      with Av -> infinity the element is ElasticBeam3d, which uses the
  loads        same transformation calls.
  Restore      a Linear + -geomNonlinear element gives the same tip displacement
               after a database save and restore.
"""
import math
import os
import tempfile

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E, G = 30.0e6, 12.0e6
A, JX, IY, IZ = 0.12, 2.0e-3, 9.0e-4, 1.6e-3
AVY, AVZ = 0.05, 0.03
RIGID = 1.0e9      # shear area that makes the element shear rigid
L = 4.0
OFFSETS = (0.2, -0.1, 0.3, -0.15, 0.25, -0.4)
LOAD = (3.0, -2.0, 5.0, 0.7, -0.4, 0.3)


@pytest.fixture(autouse=True)
def clean_domain():
    ops.wipe()
    yield
    ops.wipe()


def _cantilever(ele, transf, offsets=None, avy=AVY, avz=AVZ, extra=()):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, L)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    args = [transf, 1, 1.0, 0.0, 0.0]
    if offsets is not None:
        args += ["-jntOffset", *offsets]
    ops.geomTransf(*args)
    if ele == "timoshenko":
        ops.element("ElasticTimoshenkoBeam", 1, 1, 2, E, G, A, JX, IY, IZ, avy, avz, 1, *extra)
    else:
        ops.element("elasticBeamColumn", 1, 1, 2, A, E, G, JX, IY, IZ, 1, *extra)


def _solve(load_fn, nsteps=1, algorithm=("Newton",), iterations=50):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    load_fn()
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-11, iterations)
    ops.algorithm(*algorithm)
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0


def _tip():
    return [ops.nodeDisp(2, i) for i in range(1, 7)]


def _compliance(inertia, shear_area):
    return L**3 / (3.0 * E * inertia) + L / (G * shear_area)


@pytest.mark.parametrize("dof,inertia,shear_area", [(1, IY, AVZ), (2, IZ, AVY)],
                         ids=["bending-Iy", "bending-Iz"])
def test_pdelta_transf_gives_chord_softening(dof, inertia, shear_area):
    f = _compliance(inertia, shear_area)
    H, P = 1.0, 0.3 * L / f
    _cantilever("timoshenko", "PDelta")
    tip = [H if i == dof else 0.0 for i in (1, 2)]
    _solve(lambda: ops.load(2, tip[0], tip[1], -P, 0.0, 0.0, 0.0))
    assert ops.nodeDisp(2, dof) == pytest.approx(H * f / (1.0 - 0.3), rel=1e-8)


def test_corotational_tip_moment_closed_form():
    theta = 0.5
    M = theta * E * IY / L
    _cantilever("timoshenko", "Corotational")
    _solve(lambda: ops.load(2, 0.0, 0.0, 0.0, 0.0, M, 0.0), nsteps=5)
    assert ops.nodeDisp(2, 1) == pytest.approx(L * math.sin(theta / 2.0), rel=1e-6)
    assert ops.nodeDisp(2, 3) == pytest.approx(L * (math.cos(theta / 2.0) - 1.0), rel=1e-6)
    assert ops.nodeDisp(2, 5) == pytest.approx(theta, rel=1e-6)

    # one solve with the initial stiffness at a small rotation, u = L theta/2
    theta = 1.0e-6
    M = theta * E * IY / L
    _cantilever("timoshenko", "Corotational")
    _solve(lambda: ops.load(2, 0.0, 0.0, 0.0, 0.0, M, 0.0), algorithm=("Linear", "-initial"))
    assert ops.nodeDisp(2, 1) == pytest.approx(L * theta / 2.0, rel=1e-5)


def _explicit_rigid_links(extra=()):
    """The flexible element between the arm ends, rigid links to the nodes."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, L)
    ops.node(11, *OFFSETS[:3])
    ops.node(12, OFFSETS[3], OFFSETS[4], L + OFFSETS[5])
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.rigidLink("beam", 1, 11)
    ops.rigidLink("beam", 2, 12)
    ops.geomTransf("Linear", 1, 1.0, 0.0, 0.0)
    ops.element("ElasticTimoshenkoBeam", 1, 11, 12, E, G, A, JX, IY, IZ, AVY, AVZ, 1, *extra)


def test_joint_offsets_match_explicit_rigid_links():
    _cantilever("timoshenko", "Linear", offsets=OFFSETS)
    _solve(lambda: ops.load(2, *LOAD))
    with_offsets = _tip()
    _explicit_rigid_links()
    _solve(lambda: ops.load(2, *LOAD))
    assert with_offsets == pytest.approx(_tip(), rel=1e-7)


def test_geomnonlinear_pdelta_offsets_match_legacy_path_with_rigid_links():
    """PDelta chord + basic P-delta matrix = the full local geometric matrix."""
    load = (1.0, 0.5, -0.3 * L / _compliance(IY, AVZ), 0.0, 0.0, 0.0)
    _cantilever("timoshenko", "PDelta", offsets=OFFSETS, extra=("-geomNonlinear",))
    _solve(lambda: ops.load(2, *load))
    with_offsets = _tip()
    _explicit_rigid_links(extra=("-geomNonlinear",))
    _solve(lambda: ops.load(2, *load))
    assert with_offsets == pytest.approx(_tip(), rel=1e-7)


@pytest.mark.parametrize("transf,offsets", [("Linear", OFFSETS), ("PDelta", OFFSETS),
                                            ("Corotational", None)],
                         ids=["Linear-offsets", "PDelta-offsets", "Corotational"])
def test_uniform_load_matches_elastic_beam(transf, offsets):
    """CorotCrdTransf3d has no rigid joint zones, so that case has no offsets."""
    def run(ele):
        _cantilever(ele, transf, offsets=offsets, avy=RIGID, avz=RIGID)
        _solve(lambda: ops.eleLoad("-ele", 1, "-type", "-beamUniform", 2000.0, -1500.0, 800.0))
        ops.reactions()
        return _tip() + [ops.nodeReaction(1, i) for i in range(1, 7)]

    assert run("timoshenko") == pytest.approx(run("elastic"), rel=1e-6, abs=1e-7)


def test_local_force_with_offsets_and_load_matches_elastic_beam():
    def run(ele):
        _cantilever(ele, "Linear", offsets=OFFSETS, avy=RIGID, avz=RIGID)
        def load():
            ops.load(2, *LOAD)
            ops.eleLoad("-ele", 1, "-type", "-beamUniform", 2.0, -1.5, 0.8)
        _solve(load)
        return ops.eleResponse(1, "localForce")

    assert run("timoshenko") == pytest.approx(run("elastic"), rel=1e-6, abs=1e-8)


def test_restore_keeps_geomnonlinear_flag():
    H, P = 1.0, 0.3 * L / _compliance(IY, AVZ)

    def model():
        _cantilever("timoshenko", "Linear", extra=("-geomNonlinear",))

    def run():
        _solve(lambda: ops.load(2, H, 0.0, -P, 0.0, 0.0, 0.0))
        return ops.nodeDisp(2, 1)

    model()
    reference = run()
    assert reference > 1.2 * H * _compliance(IY, AVZ)

    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "timoshenko_db")
        model()
        ops.database("File", path)
        ops.save(1)
        ops.wipe()
        model()
        ops.database("File", path)
        ops.restore(1)
        restored = run()
        ops.wipe()
    assert restored == pytest.approx(reference, rel=1e-9)
