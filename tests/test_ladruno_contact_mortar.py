"""Mortar segment-to-segment contact with augmented Lagrangian (Uzawa) normal enforcement.

Checks:

- Penalty: a slave facet pressed on a fixed coincident master settles at the penalty
  penetration P/(epsN A) (1e-7), with static Newton converging; also for a non-matching
  2x2 slave on one master quad and a coarse slave quad over a 4x4 master (full coverage).
- A mortar pair whose facets do not overlap leaves a Newmark run bitwise unchanged.
- Augmentation: a held-load loop (ladrunoBeginAugment, analyze at a zero load increment,
  ladrunoEndAugment) drives ladrunoMortarPenetration below 1e-9 for epsN = 1e4 and 1e6,
  and the converged slave position is the same for both (1e-7); the penalty penetration
  is O(1/epsN). The same on a deformable stdBrick block.
- The equation count does not change across augmentations (4 free DOFs).
- Release: after an augmented contact is pulled off, the penetration query is 0 and the
  slave rides its springs alone (1e-6).
- Recorders: the augmentation sweep adds no recorder sample and leaves the recorded time
  sequence unchanged; without the begin/end bracket the same loop adds three samples per
  step; an empty begin/end pair leaves the output byte-identical.
- The held re-solves do not advance the load factor or re-apply the load (reaction equal
  to the load to 1e-6).
- Command surface: ladrunoContactInfo returns 4 values with the mortar count last;
  mortar-only options without -mortar, node-to-segment options with -mortar and a -slave
  node set as the mortar slave are refused.

Every test fails on the base with an unknown command.
"""
import os
import tempfile

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


EPS = 1.0e6
PTOT = 400.0


def _analyze_augmented(maxAug=20, augTol=1.0e-9, query=None):
    """Held-load Uzawa loop: re-solve and commit at a zero load increment until the
    measure drops below augTol. Returns (status, augmentations, measure); status 1 when
    maxAug is reached, < 0 for a failed re-solve."""
    if query is None:
        query = ops.ladrunoMortarPenetration
    ops.integrator("LoadControl", 0.0)
    try:
        ops.ladrunoBeginAugment()
        pen = query()
        for k in range(1, int(maxAug) + 1):
            rc = ops.analyze(1)
            if rc != 0:
                return rc, k, pen
            pen = query()
            if pen < augTol:
                return 0, k, pen
        return 1, int(maxAug), pen
    finally:
        ops.ladrunoEndAugment()


def _unit_square(tags, z):
    for t, x, y in zip(tags, (0.0, 1.0, 1.0, 0.0), (0.0, 0.0, 1.0, 1.0)):
        ops.node(t, x, y, z)
    return tags


def _static_setup(tol=1.0e-11, maxit=50):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, maxit, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _facet_on_master(epsN, series="Linear"):
    """One slave quad seeded 1e-4 into a coincident fixed master quad, x and y fixed,
    pressed by PTOT. Returns the slave tags."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    _unit_square([5, 6, 7, 8], -1.0e-4)
    for t in (5, 6, 7, 8):
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries(series, 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -PTOT / 4.0)
    _static_setup()
    return [5, 6, 7, 8]


# ------------------------------------------------------------------ penalty
def test_matched_penalty_penetration():
    _facet_on_master(EPS)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    pos = -1.0e-4 + ops.nodeDisp(5, 3)
    assert abs(pos - (-PTOT / EPS)) < 1.0e-7
    for t in (6, 7, 8):
        assert abs(-1.0e-4 + ops.nodeDisp(t, 3) - pos) < 1.0e-9


def test_nonmatched_penalty_penetration():
    """2x2 slave mesh on one master quad, consistent nodal loads: uniform depth P/epsN."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    z0 = -1.0e-4
    nid = {}
    k = 11
    for j in range(3):
        for i in range(3):
            nid[(i, j)] = k
            ops.node(k, i * 0.5, j * 0.5, z0)
            ops.fix(k, 1, 1, 0)
            k += 1
    sflat = []
    for j in range(2):
        for i in range(2):
            sflat += [nid[(i, j)], nid[(i + 1, j)], nid[(i + 1, j + 1)], nid[(i, j + 1)]]
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, *sflat)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for (i, j), n in nid.items():
        ne = (1 if i in (0, 2) else 0) + (1 if j in (0, 2) else 0)
        w = {2: 1.0 / 16.0, 1: 1.0 / 8.0, 0: 1.0 / 4.0}[ne]   # tributary areas
        ops.load(n, 0.0, 0.0, -PTOT * w)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    for n in nid.values():
        assert abs(z0 + ops.nodeDisp(n, 3) - (-PTOT / EPS)) < 1.0e-6


def test_coarse_slave_fine_master_full_coverage():
    """One slave quad over a 4x4 master grid: every master facet is paired."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nm = 4
    mid = {}
    k = 1
    for j in range(nm + 1):
        for i in range(nm + 1):
            mid[(i, j)] = k
            ops.node(k, i / nm, j / nm, 0.0)
            ops.fix(k, 1, 1, 1)
            k += 1
    mflat = []
    for j in range(nm):
        for i in range(nm):
            mflat += [mid[(i, j)], mid[(i + 1, j)], mid[(i + 1, j + 1)], mid[(i, j + 1)]]
    _unit_square([101, 102, 103, 104], -1.0e-4)
    for t in (101, 102, 103, 104):
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, *mflat)
    ops.contactSurface(2, "-slave-segments", 4, 101, 102, 103, 104)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (101, 102, 103, 104):
        ops.load(t, 0.0, 0.0, -PTOT / 4.0)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    for t in (101, 102, 103, 104):
        assert abs(-1.0e-4 + ops.nodeDisp(t, 3) - (-PTOT / EPS)) < 1.0e-7


def test_disjoint_mortar_pair_is_inert():
    def run(define):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _unit_square([1, 2, 3, 4], 0.0)
        for t in (1, 2, 3, 4):
            ops.fix(t, 1, 1, 1)
        for t, x, y in [(5, 5.0, 0.0), (6, 6.0, 0.0), (7, 6.0, 1.0), (8, 5.0, 1.0)]:
            ops.node(t, x, y, 0.10)
            ops.fix(t, 1, 1, 0)
            ops.mass(t, 1.0, 1.0, 1.0)
            ops.setNodeVel(t, 3, 0.2, "-commit")
        if define:
            ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
            ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
            ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0)
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25, 0)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        out = []
        for _ in range(8):
            assert ops.analyze(1, 1.0e-3) == 0
            out.append(tuple(ops.nodeDisp(t, 3) for t in (5, 6, 7, 8)))
        return out
    assert run(True) == run(False)


# ------------------------------------------------------------- augmentation
def _augmented_facet(epsN, augTol=1.0e-9, maxAug=40):
    _facet_on_master(epsN)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    n_before = ops.systemSize()
    pen_penalty = ops.ladrunoMortarPenetration()
    status, n_aug, pen = _analyze_augmented(maxAug=maxAug, augTol=augTol)
    assert status == 0, f"augmentation did not converge (status {status}, pen {pen:.2e})"
    return pen_penalty, pen, -1.0e-4 + ops.nodeDisp(5, 3), n_before, ops.systemSize()


def test_penetration_is_epsN_independent_after_augmentation():
    augTol = 1.0e-9
    lo = _augmented_facet(1.0e4, augTol)
    hi = _augmented_facet(1.0e6, augTol)
    assert lo[0] > 50.0 * hi[0]                  # penalty: O(1/epsN)
    assert lo[1] < augTol and hi[1] < augTol
    assert abs(lo[2] - hi[2]) < 1.0e-7
    assert abs(lo[2]) < 1.0e-7


def test_equation_count_constant_across_augmentations():
    res = _augmented_facet(1.0e5)
    assert res[3] == res[4] == 4


def test_release_after_augmentation():
    epsN, ks = 1.0e5, 1.0e3
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    _unit_square([5, 6, 7, 8], -1.0e-4)
    ops.uniaxialMaterial("Elastic", 1, ks)
    for i, t in enumerate((5, 6, 7, 8)):
        ops.node(100 + t, *[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)][i])
        ops.fix(100 + t, 1, 1, 1)
        ops.fix(t, 1, 1, 0)
        ops.element("zeroLength", 200 + t, 100 + t, t, "-mat", 1, "-dir", 3)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -PTOT / 4.0)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    status, _, _ = _analyze_augmented(maxAug=40, augTol=1.0e-9)
    assert status == 0
    assert ops.ladrunoMortarPenetration() < 1.0e-9
    ops.pattern("Plain", 2, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, +2.0 * PTOT / 4.0)   # net load now pulls the slave away
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    ops.integrator("LoadControl", 0.0)
    for _ in range(5):
        assert ops.analyze(1) == 0
    assert ops.ladrunoMortarPenetration() == 0.0
    assert abs(ops.nodeDisp(5, 3) - (PTOT / 4.0) / ks) < 1.0e-6


def _elastic_block_alm(epsN, augTol=1.0e-9):
    """Unit stdBrick block, bottom face the mortar slave on a fixed master quad, top
    face loaded; uniaxial (x, y fixed, nu = 0)."""
    L, E = 1.0, 2.0e4
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    z0 = -1.0e-4
    pts = [(0, 0, z0), (L, 0, z0), (L, L, z0), (0, L, z0),
           (0, 0, z0 + L), (L, 0, z0 + L), (L, L, z0 + L), (0, L, z0 + L)]
    tags = list(range(5, 13))
    for t, (x, y, z) in zip(tags, pts):
        ops.node(t, float(x), float(y), float(z))
        ops.fix(t, 1, 1, 0)
    ops.nDMaterial("ElasticIsotropic", 1, E, 0.0)
    ops.element("stdBrick", 1, *tags, 1)
    bottom, top = tags[:4], tags[4:]
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, *bottom)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    P = 1.0e-3 * E * L * L
    for t in top:
        ops.load(t, 0.0, 0.0, -P / 4.0)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    penalty_bottom = z0 + ops.nodeDisp(bottom[0], 3)
    status, _, _ = _analyze_augmented(maxAug=60, augTol=augTol)
    assert status == 0
    return penalty_bottom, z0 + ops.nodeDisp(bottom[0], 3)


def test_augmentation_on_a_deformable_block():
    lo = _elastic_block_alm(1.0e4)
    hi = _elastic_block_alm(1.0e6)
    assert lo[0] < -1.0e-6
    assert abs(lo[0]) > 10.0 * abs(hi[0])
    assert abs(lo[1]) < 1.0e-7 and abs(hi[1]) < 1.0e-7
    assert abs(lo[1] - hi[1]) < 1.0e-7


# ------------------------------------------------------------------ recorders
def _record(with_aug, bracket=True, noop=False, nsteps=4, epsN=1.0e4):
    fd, path = tempfile.mkstemp(suffix=".out", prefix="mortar_rec_")
    os.close(fd)
    try:
        _facet_on_master(epsN)
        ops.recorder("Node", "-file", path, "-time", "-node", 5, "-dof", 3, "disp")
        for _ in range(nsteps):
            ops.integrator("LoadControl", 1.0)
            if noop:
                ops.ladrunoBeginAugment()
                ops.ladrunoEndAugment()
            assert ops.analyze(1) == 0
            if with_aug and bracket:
                _analyze_augmented(maxAug=40, augTol=1.0e-10)
            elif with_aug:
                ops.integrator("LoadControl", 0.0)
                for _ in range(3):
                    assert ops.analyze(1) == 0
        ops.wipe()
        with open(path) as f:
            return [line.strip() for line in f if line.strip()]
    finally:
        if os.path.exists(path):
            os.remove(path)


def test_augmentation_adds_no_recorder_samples():
    base = [float(s.split()[0]) for s in _record(False)]
    aug = [float(s.split()[0]) for s in _record(True)]
    assert len(base) == 4
    assert aug == base


def test_unbracketed_loop_adds_recorder_samples():
    assert len(_record(True, bracket=False)) == 4 + 4 * 3


def test_empty_bracket_is_byte_identical():
    assert _record(False, noop=True) == _record(False)


def test_held_resolves_do_not_reapply_the_load():
    _facet_on_master(1.0e4)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    assert ops.getTime() == 1.0
    status, n_aug, _ = _analyze_augmented(maxAug=60, augTol=1.0e-10)
    assert status == 0 and n_aug >= 1
    assert ops.getTime() == 1.0
    ops.reactions()
    r = sum(ops.nodeReaction(t, 3) for t in (5, 6, 7, 8))
    assert abs(r - PTOT) < 1.0e-6 * PTOT
    ops.ladrunoBeginAugment()
    try:
        ops.integrator("LoadControl", 0.0)
        for _ in range(10):
            assert ops.analyze(1) == 0
            assert ops.getTime() == 1.0
    finally:
        ops.ladrunoEndAugment()
    ops.reactions()
    r = sum(ops.nodeReaction(t, 3) for t in (5, 6, 7, 8))
    assert abs(r - PTOT) < 1.0e-6 * PTOT


# ------------------------------------------------------------- command surface
def _two_facets():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(1 + i, float(x), float(y), 0.0)
        ops.node(5 + i, float(x), float(y), 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contactSurface(3, "-slave", 5, 6, 7, 8)


def test_contact_info_counts_mortar_contacts():
    _two_facets()
    ops.contact(1, 1, 3, 1.0e5)
    ops.contact(2, 1, 2, "-mortar", "-epsN", 1.0e5, "-ngp", 4, "-augTol", 1e-9, "-maxAug", 30)
    info = ops.ladrunoContactInfo()
    assert len(info) == 4
    assert info[0] == 1 and info[3] == 1


@pytest.mark.parametrize("args", [
    (1, 1, 3, 1.0e5, "-epsN", 1.0e5),                         # mortar option without -mortar
    (1, 1, 3, 1.0e5, "-ngp", 4),
    (1, 1, 2, "-mortar", "-epsN", 1.0e5, "-geomtan"),         # node-to-segment options
    (1, 1, 2, "-mortar", "-epsN", 1.0e5, "-cell", 2.0),
    (1, 1, 2, "-mortar", "-epsN", 1.0e5, "-visc", 1.0),
    (1, 1, 3, "-mortar", "-epsN", 1.0e5),                     # node-set slave
    (1, 1, 2, "-mortar", "-epsN", -1.0),                      # negative penalty
])
def test_refusals(args):
    _two_facets()
    with pytest.raises(ops.OpenSeesError):
        ops.contact(*args)


def test_mortar_without_penalty_fails_at_handle():
    _two_facets()
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6, 7, 8):
        ops.fix(t, 1, 1, 0)
    ops.contact(1, 1, 2, "-mortar")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, 0.0, -1.0)
    _static_setup()
    ops.integrator("LoadControl", 1.0)
    try:
        rc = ops.analyze(1)
    except Exception:
        rc = -1
    assert rc < 0
