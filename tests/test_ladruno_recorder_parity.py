"""recorder ladruno: values against recorder mpco and against analytic data.

Each model runs `recorder mpco` and `recorder ladruno` with identical requests
in one process; the two files are flattened by information, not by layout
(ladruno_format.iter_step_slices reads both the chunked .ladruno DATA and the
per-step .mpco DATA/STEP_<k>), and every common value must agree to 1e-12
(absolute) + 1e-9 (relative):

  * nodal: displacement and reactionForce of a 2-D truss frame;
  * multi-stage: the same frame with a node and an element added between two
    analyze() calls (two MODEL_STAGEs);
  * element: FourNodeQuad `stress`, dispBeamColumn `section.fiber.stress`,
    Truss `axialForce`/`basicForces` and zeroLength `force`, ShellMITC4 `forces`;
  * eigen: `-N modesOfVibration` eigenvectors of a portal frame (stage KIND
    "eigen", DATA/STEP_<k>/MODE_<k> layout).

Without an mpco oracle:

  * envelope: `-envelope` MIN/MAX/ABSMAX/ARG_STEP equal the reduction of the
    time series of the same run;
  * shell layers: `material.fiber.stress` and `section.fiber.stress` on a
    LayeredShell both expand to nGP x nLayers x 5 columns with equal values; a
    bare `section` / `material` verb records no bucket and does not crash;
  * precision: `-precision f32` stores float32 DATA (STORED_PRECISION "f32")
    within 6e-7 relative of the f64 run;
  * quadrature: GP_PARAM of quad, tri31, stdBrick, quad9n, TenNodeTetrahedron
    and 20NodeBrick equals the element's own Gauss rule in its own order, and
    GLOBAL_GP_COORDS equals x(GP_PARAM) rebuilt from the nodes;
  * local axes: beam LOCAL_AXES frames are unit quaternions whose local x is
    the element axis (skipped when no element answers the `localAxes`
    response).

Every written file also passes ladruno_format.validate().
"""
import os

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

h5py = pytest.importorskip("h5py")

h5py = pytest.importorskip("h5py")  # the reader needs h5py
import ladruno_basis  # noqa: E402
import ladruno_format as lf  # noqa: E402

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

ATOL = 1e-12
RTOL = 1e-9


def _norm(name):
    return "".join(c for c in name.lower() if c.isalnum())


def _nodal(path):
    """{(node, result:comp, step): value} over ON_NODES, .mpco or .ladruno."""
    out = {}
    with h5py.File(path, "r") as f:
        for stage in [k for k in f if k.startswith("MODEL_STAGE")]:
            base = f.get(f"{stage}/RESULTS/ON_NODES")
            if base is None:
                continue
            for rname in base:
                grp = base[rname]
                if "ID" not in grp or "DATA" not in grp or "modesofvibration" in _norm(rname):
                    continue
                ids = np.asarray(grp["ID"][...]).reshape(-1)
                comps = [c for c in lf._attr(grp, "COMPONENTS").split(",") if c]
                for step, arr in lf.iter_step_slices(grp):
                    for r, tag in enumerate(ids):
                        for c, cn in enumerate(comps):
                            if c < arr.shape[1]:
                                out[(int(tag), f"{_norm(rname)}:{cn}", step)] = float(arr[r, c])
    return out


def _element(path):
    """{(element, result, column, step): value} over ON_ELEMENTS buckets."""
    out = {}
    with h5py.File(path, "r") as f:
        for stage in [k for k in f if k.startswith("MODEL_STAGE")]:
            base = f.get(f"{stage}/RESULTS/ON_ELEMENTS")
            if base is None:
                continue
            for rname in base:
                rgrp = base[rname]
                if not isinstance(rgrp, h5py.Group):
                    continue
                for bname in rgrp:
                    bgrp = rgrp[bname]
                    if not isinstance(bgrp, h5py.Group) or "ID" not in bgrp or "DATA" not in bgrp:
                        continue
                    ids = np.asarray(bgrp["ID"][...]).reshape(-1)
                    for step, arr in lf.iter_step_slices(bgrp):
                        for r, tag in enumerate(ids):
                            for j in range(arr.shape[1]):
                                out[(int(tag), _norm(rname), j, step)] = float(arr[r, j])
    return out


def _modes(path):
    """{(node, MODE_k, comp): value} over the modesOfVibration groups."""
    out = {}
    with h5py.File(path, "r") as f:
        for s in [k for k in f if k.startswith("MODEL_STAGE")]:
            base = f.get(f"{s}/RESULTS/ON_NODES")
            if base is None:
                continue
            for rname in base:
                if "modesofvibration" not in _norm(rname):
                    continue
                grp = base[rname]
                ids = np.asarray(grp["ID"][...]).reshape(-1)
                data = grp["DATA"]
                for step in data:
                    for mk in data[step]:
                        arr = np.atleast_2d(np.asarray(data[step][mk][...]))
                        for r, tag in enumerate(ids):
                            for c in range(arr.shape[1]):
                                out[(int(tag), mk, c)] = float(arr[r, c])
    return out


def _assert_parity(ref, new, what):
    assert ref, f"no {what} values read from the .mpco reference"
    assert new, f"no {what} values read from the .ladruno file"
    common = set(ref) & set(new)
    assert common, f"{what}: no common keys (ref e.g. {sorted(ref)[:2]}, new e.g. {sorted(new)[:2]})"
    bad = [(k, ref[k], new[k]) for k in sorted(common)
           if abs(ref[k] - new[k]) > ATOL + RTOL * abs(ref[k])]
    assert not bad, f"{what}: {len(bad)} mismatches, e.g. {bad[:5]}"
    assert any(abs(ref[k]) > 0.0 for k in common), f"{what}: every compared value is zero"
    return len(common)


def _static(steps, dlam, system="BandSPD"):
    ops.system(system)
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", dlam)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(steps) == 0


def _truss_frame():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 2.0, 0.0)
    ops.node(3, 4.0, 0.0)
    ops.node(4, 2.0, 2.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, 2.0e5)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.element("truss", 2, 2, 3, 1.0, 1)
    ops.element("truss", 3, 1, 4, 1.0, 1)
    ops.element("truss", 4, 2, 4, 1.0, 1)
    ops.element("truss", 5, 3, 4, 1.0, 1)


def _both(tmp_path, stem, *req):
    ref = str(tmp_path / f"{stem}.mpco")
    new = str(tmp_path / f"{stem}.ladruno")
    ops.recorder("mpco", ref, *req, "-T", "dt", 0.0)
    ops.recorder("ladruno", new, *req, "-T", "dt", 0.0)
    return ref, new


def _valid(path):
    problems = lf.validate(path)
    assert problems == [], problems


# ---------------------------------------------------------------------------
# parity against recorder mpco
# ---------------------------------------------------------------------------

def test_nodal_parity(tmp_path):
    _truss_frame()
    ref, new = _both(tmp_path, "nodal", "-N", "displacement", "reactionForce")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(4, 0.0, -1.0e3)
    _static(5, 0.2)
    ops.wipe()
    _valid(new)
    n = _assert_parity(_nodal(ref), _nodal(new), "nodal")
    assert n == 4 * 4 * 5          # 4 nodes x (2 disp + 2 reaction) x 5 steps


def test_multistage_parity(tmp_path):
    _truss_frame()
    ref, new = _both(tmp_path, "stages", "-N", "displacement", "reactionForce")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(4, 0.0, -1.0e3)
    _static(3, 0.2)
    ops.node(5, 2.0, 4.0)
    ops.fix(5, 1, 1)
    ops.element("truss", 6, 4, 5, 1.0, 1)
    assert ops.analyze(3) == 0
    ops.wipe()
    _valid(new)
    with h5py.File(new, "r") as f:
        assert len([k for k in f if k.startswith("MODEL_STAGE")]) == 2
    _assert_parity(_nodal(ref), _nodal(new), "multi-stage nodal")


def test_quad_and_fiber_section_parity(tmp_path):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, x, y in ((1, 0, 0), (2, 1, 0), (3, 2, 0), (4, 0, 1), (5, 1, 1), (6, 2, 1)):
        ops.node(tag, float(x), float(y))
    ops.fix(1, 1, 1)
    ops.fix(4, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    ops.element("quad", 10, 1, 2, 5, 4, 1.0, "PlaneStress", 1)
    ops.element("quad", 11, 2, 3, 6, 5, 1.0, "PlaneStress", 1)
    ref, new = _both(tmp_path, "quad", "-E", "stress")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0e2, 0.0)
    ops.load(6, 1.0e2, 0.0)
    _static(4, 0.25)
    ops.wipe()
    _valid(new)
    _assert_parity(_element(ref), _element(new), "FourNodeQuad stress")

    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 2.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, 2.0e5)
    ops.section("Fiber", 1)
    ops.patch("rect", 1, 2, 3, -0.5, -0.25, 0.5, 0.25)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 3)
    ops.element("dispBeamColumn", 10, 1, 2, 1, 1)
    ops.element("dispBeamColumn", 11, 2, 3, 1, 1)
    bref, bnew = _both(tmp_path, "beam", "-E", "section.fiber.stress")
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(3, 0.0, -1.0e1, 0.0)
    _static(4, 0.25)
    ops.wipe()
    _valid(bnew)
    _assert_parity(_element(bref), _element(bnew), "dispBeamColumn section.fiber.stress")


def test_truss_and_zero_length_force_parity(tmp_path):
    _truss_frame()
    ops.node(5, 0.0, 4.0)
    ops.node(6, 0.0, 4.0)
    ops.fix(6, 1, 1)
    ops.uniaxialMaterial("Elastic", 2, 5.0e3)
    ops.element("zeroLength", 20, 6, 5, "-mat", 2, 2, "-dir", 1, 2)
    ref, new = _both(tmp_path, "force", "-E", "axialForce", "basicForces", "force",
                     "-N", "displacement")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(4, 0.0, -1.0e3)
    ops.load(5, 1.0e2, -2.0e2)
    _static(5, 0.2)
    ops.wipe()
    _valid(new)
    _assert_parity(_element(ref), _element(new), "Truss/zeroLength forces")


def _shell_mesh(section_cmd):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nid = 0
    for j in range(3):
        for i in range(3):
            nid += 1
            ops.node(nid, float(i), float(j), 0.0)
    for n in (1, 2, 3):
        ops.fix(n, 1, 1, 1, 1, 1, 1)
    section_cmd()

    def n(i, j):
        return j * 3 + i + 1
    etag = 100
    for j in range(2):
        for i in range(2):
            etag += 1
            ops.element("ShellMITC4", etag, n(i, j), n(i + 1, j), n(i + 1, j + 1), n(i, j + 1), 1)


def _shell_load(steps, dlam):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for nn in (7, 8, 9):
        ops.load(nn, 0.0, 0.0, -1.0e1, 0.0, 0.0, 0.0)
    _static(steps, dlam, system="BandGeneral")
    ops.wipe()


def _layered_section():
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.25)
    ops.nDMaterial("PlateFiber", 2, 1)
    layers = []
    for _ in range(3):
        layers += [2, 0.1 / 3]
    ops.section("LayeredShell", 1, 3, *layers)


def test_shell_forces_parity(tmp_path):
    _shell_mesh(lambda: ops.section("ElasticMembranePlateSection", 1, 2.0e5, 0.25, 0.1, 0.0))
    ref, new = _both(tmp_path, "shell", "-E", "forces", "-N", "displacement", "rotation")
    _shell_load(5, 0.2)
    _valid(new)
    _assert_parity(_element(ref), _element(new), "ShellMITC4 forces")
    _assert_parity(_nodal(ref), _nodal(new), "shell nodal")


def test_eigen_modes_parity(tmp_path):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 4.0, 0.0)
    ops.node(3, 0.0, 3.0, "-mass", 2.0, 2.0, 1.0e-4)
    ops.node(4, 4.0, 3.0, "-mass", 2.0, 2.0, 1.0e-4)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    ops.element("elasticBeamColumn", 1, 1, 3, 1.0e2, 2.0e5, 1.0e3, 1)
    ops.element("elasticBeamColumn", 2, 2, 4, 1.0e2, 2.0e5, 1.0e3, 1)
    ops.element("elasticBeamColumn", 3, 3, 4, 1.0e2, 2.0e5, 1.0e3, 1)
    ref, new = _both(tmp_path, "eigen", "-N", "modesOfVibration")
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.eigen(3)
    ops.record()
    ops.wipe()
    _valid(new)
    with h5py.File(new, "r") as f:
        kinds = [str(lf._attr(f[s], "KIND")) for s in f if s.startswith("MODEL_STAGE")]
    assert "eigen" in kinds, kinds
    n = _assert_parity(_modes(ref), _modes(new), "modesOfVibration")
    assert n >= 2 * 3 * 2          # free nodes x modes x translational components


# ---------------------------------------------------------------------------
# checks without an mpco oracle
# ---------------------------------------------------------------------------

def test_envelope_equals_time_series_reduction(tmp_path):
    ts_path = str(tmp_path / "ts.ladruno")
    env_path = str(tmp_path / "env.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 100.0, 0.0)
    ops.node(3, 200.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(3, 0, 1)
    ops.mass(2, 1.0, 1.0)
    ops.mass(3, 1.0, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.element("Truss", 2, 2, 3, 1.0, 1)
    ops.recorder("ladruno", ts_path, "-N", "displacement", "-T", "dt", 0.0)
    ops.recorder("ladruno", env_path, "-N", "displacement", "-T", "dt", 0.0, "-envelope")
    ops.timeSeries("Trig", 1, 0.0, 100.0, 0.5)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 5.0, 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    for _ in range(40):
        assert ops.analyze(1, 0.02) == 0
    ops.wipe()
    _valid(env_path)

    with lf.LadrunoReader(ts_path) as r:
        grp = r.f[f"{r.stages()[0]}/RESULTS/ON_NODES/DISPLACEMENT"]
        comps = [c for c in lf._attr(grp, "COMPONENTS").split(",") if c]
        mn = mx = am = arg = None
        for k, arr in lf.iter_step_slices(grp):
            a = np.atleast_2d(np.asarray(arr, dtype=float))
            if mn is None:
                mn, mx, am = a.copy(), a.copy(), np.abs(a)
                arg = np.full(a.shape, k, dtype=int)
            else:
                mn, mx = np.minimum(mn, a), np.maximum(mx, a)
                upd = np.abs(a) > am                       # first occurrence wins
                am = np.where(upd, np.abs(a), am)
                arg = np.where(upd, k, arg)
    with lf.LadrunoReader(env_path) as r:
        env = r.envelopes(r.stages()[0])["ON_NODES"]["DISPLACEMENT"]
    assert np.any(am > 0.0)
    for key, want in (("min", mn), ("max", mx), ("absmax", am), ("arg_step", arg)):
        got = np.atleast_2d(np.asarray(env[key]))
        assert got.shape == want.shape and np.allclose(got, want, atol=1e-9), key
    assert env.get("components") == comps


def test_shell_layer_verbs(tmp_path):
    f_mat = str(tmp_path / "material.ladruno")
    f_sec = str(tmp_path / "section.ladruno")
    _shell_mesh(_layered_section)
    ops.recorder("ladruno", f_mat, "-E", "material.fiber.stress", "-N", "displacement")
    ops.recorder("ladruno", f_sec, "-E", "section.fiber.stress", "-N", "displacement")
    _shell_load(2, 0.5)
    ncols = []
    for path in (f_mat, f_sec):
        _valid(path)
        with lf.LadrunoReader(path) as r:
            base = r.f[f"{r.stages()[0]}/RESULTS/ON_ELEMENTS"]
            n, ids = 0, 0
            for rname in base:
                for bname in base[rname]:
                    b = base[rname][bname]
                    n += len(lf.decode_columns(b["COLUMN_MAP"]))
                    ids = max(ids, int(b["ID"].shape[0]))
            ncols.append((n, ids))
    assert ncols == [(4 * 3 * 5, 4)] * 2, ncols   # nGP x nLayers x nComp, 4 shells
    def strip(path):
        # (element, gp, section, fiber, component, step), result name dropped
        with lf.LadrunoReader(path) as r:
            nm = lf.normalize_element(r, r.stages()[0])
        return {(t, gp, sec, fib, rc.split(":", 1)[1], step): v
                for (t, gp, sec, fib, rc, step), v in nm.items()}
    vm, vs = strip(f_mat), strip(f_sec)
    assert vm.keys() == vs.keys() and any(abs(v) > 1e-12 for v in vm.values())
    assert max(abs(vm[k] - vs[k]) for k in vm) <= 1e-12


def test_shell_bare_verbs_record_nothing(tmp_path):
    paths = [str(tmp_path / "bare_section.ladruno"), str(tmp_path / "bare_material.ladruno")]
    _shell_mesh(_layered_section)
    ops.recorder("ladruno", paths[0], "-E", "section", "-N", "displacement")
    ops.recorder("ladruno", paths[1], "-E", "material", "-N", "displacement")
    _shell_load(2, 0.5)
    for path in paths:
        _valid(path)
        with lf.LadrunoReader(path) as r:
            stage = r.stages()[0]
            base = r.f[f"{stage}/RESULTS/ON_ELEMENTS"]
            assert sum(len(base[k]) for k in base) == 0
            assert lf.normalize_nodal(r, stage)


def test_precision_f32(tmp_path):
    f64, f32 = str(tmp_path / "f64.ladruno"), str(tmp_path / "f32.ladruno")
    nx, ny = 20, 6

    def nid(i, j):
        return 1 + i + j * (nx + 1)
    for path, prec in ((f64, "f64"), (f32, "f32")):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j), i * 6.0 / nx, j * 2.0 / ny)
        for j in range(ny + 1):
            ops.fix(nid(0, j), 1, 1)
        ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.25, 2.0)
        etag = 1
        for j in range(ny):
            for i in range(nx):
                ops.element("quad", etag, nid(i, j), nid(i + 1, j), nid(i + 1, j + 1),
                            nid(i, j + 1), 0.1, "PlaneStress", 1)
                etag += 1
        ops.recorder("ladruno", path, "-N", "displacement", "-E", "stress",
                     "-precision", prec, "-T", "dt", 0.0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(nid(nx, ny), 0.0, -40.0)
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("BandGeneral")
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        assert ops.analyze(120, 1.0e-3) == 0
        ops.wipe()
        _valid(path)
    with lf.LadrunoReader(f32) as r:
        assert r.stored_precision() == "f32"
    with lf.LadrunoReader(f64) as r:
        assert r.stored_precision() == "f64"
    with h5py.File(f32, "r") as f:
        st = f[[k for k in f if k.startswith("MODEL_STAGE")][0]]
        assert st["RESULTS/ON_NODES/DISPLACEMENT/DATA"].dtype == np.float32
    def by_step_index(d):
        # the commit tag keeps counting across wipe in one interpreter, so the
        # two runs are aligned by the position of each step, not by its tag
        order = {s: i for i, s in enumerate(sorted({k[-1] for k in d}))}
        return {k[:-1] + (order[k[-1]],): v for k, v in d.items()}
    for read in (_nodal, _element):
        a, b = by_step_index(read(f64)), by_step_index(read(f32))
        assert a.keys() == b.keys() and a
        bad = [k for k in a if abs(a[k] - b[k]) > 1e-30 + 6.0e-7 * abs(a[k])]
        assert not bad, bad[:5]


_A = 1.0 / np.sqrt(3.0)
_B = np.sqrt(0.6)
_AL = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0
_BE = (5.0 - np.sqrt(5.0)) / 20.0
_GP = {
    "quad": [[-_A, -_A], [_A, -_A], [_A, _A], [-_A, _A]],
    "tri31": [[1.0 / 3.0, 1.0 / 3.0]],
    "stdBrick": [[-_A, -_A, -_A], [-_A, -_A, _A], [-_A, _A, -_A], [-_A, _A, _A],
                 [_A, -_A, -_A], [_A, -_A, _A], [_A, _A, -_A], [_A, _A, _A]],
    "quad9n": [[-_B, -_B], [_B, -_B], [_B, _B], [-_B, _B],
               [0, -_B], [_B, 0], [0, _B], [-_B, 0], [0, 0]],
    "TenNodeTetrahedron": [[_AL, _BE, _BE], [_BE, _AL, _BE], [_BE, _BE, _AL], [_BE, _BE, _BE]],
    "20NodeBrick": [
        [-_B, -_B, -_B], [_B, -_B, -_B], [_B, _B, -_B], [-_B, _B, -_B],
        [-_B, -_B, _B], [_B, -_B, _B], [_B, _B, _B], [-_B, _B, _B],
        [0, -_B, -_B], [_B, 0, -_B], [0, _B, -_B], [-_B, 0, -_B],
        [0, -_B, _B], [_B, 0, _B], [0, _B, _B], [-_B, 0, _B],
        [-_B, -_B, 0], [_B, -_B, 0], [_B, _B, 0], [-_B, _B, 0],
        [_B, 0, 0], [0, _B, 0], [0, 0, _B], [-_B, 0, 0],
        [0, -_B, 0], [0, 0, -_B], [0, 0, 0]],
}
_NODES = {
    "quad": (2, [(0, 0), (1, 0), (1, 1), (0, 1)]),
    "tri31": (2, [(0, 0), (1, 0), (0, 1)]),
    "stdBrick": (3, [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                     (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]),
    "quad9n": (2, [(0, 0), (1, 0), (1, 1), (0, 1), (.5, 0), (1, .5), (.5, 1), (0, .5), (.5, .5)]),
    "TenNodeTetrahedron": (3, [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (.5, 0, 0),
                               (.5, .5, 0), (0, .5, 0), (0, 0, .5), (0, .5, .5), (.5, 0, .5)]),
    "20NodeBrick": (3, [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1),
                        (1, 1, 1), (0, 1, 1), (.5, 0, 0), (1, .5, 0), (.5, 1, 0), (0, .5, 0),
                        (.5, 0, 1), (1, .5, 1), (.5, 1, 1), (0, .5, 1), (0, 0, .5),
                        (1, 0, .5), (1, 1, .5), (0, 1, .5)]),
}


@pytest.mark.parametrize("ele", list(_GP))
def test_standard_quadrature_and_global_gp_coords(tmp_path, ele):
    ndm, xyz = _NODES[ele]
    path = str(tmp_path / f"{ele}.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndm)
    for k, c in enumerate(xyz, start=1):
        ops.node(k, *[float(v) for v in c])
        if c[ndm - 1] == 0:                  # fix the y = 0 (2-D) or z = 0 (3-D) face
            ops.fix(k, *([1] * ndm))
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    tags = list(range(1, len(xyz) + 1))
    if ndm == 2:
        ops.element(ele, 1, *tags, 1.0, "PlaneStress", 1)
    else:
        ops.element(ele, 1, *tags, 1)
    ops.recorder("ladruno", path, "-N", "displacement")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for k, c in enumerate(xyz, start=1):
        if c[ndm - 1] == 1:
            ops.load(k, *([0.0] * (ndm - 1)), 1.0e2)
    _static(2, 0.5, system="BandGeneral")
    ops.wipe()
    _valid(path)
    assert lf.round_trip_oracle(path) == []
    with lf.LadrunoReader(path) as r:
        groups = r.element_groups(r.stages()[0])
    assert len(groups) == 1
    g = next(iter(groups.values()))
    want = np.array(_GP[ele], dtype=float)
    assert g["gp_param"] is not None and g["gp_param"].shape == want.shape
    assert np.allclose(g["gp_param"], want, atol=1e-9)
    assert g["ndir"] == want.shape[1] and g["num_gp"] == want.shape[0]
    assert "global_gp_coords" in g
    X = np.array(xyz, dtype=float)
    ggp = g["global_gp_coords"].reshape(1, want.shape[0], ndm)[0]
    for k in range(want.shape[0]):
        N = ladruno_basis.shape_functions(g["topology"], len(xyz), want[k])
        assert np.max(np.abs(N @ X - ggp[k])) <= 1e-12


def _quat_rows(q):
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def test_beam_local_axes(tmp_path):
    path = str(tmp_path / "axes.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for tag, c in ((1, (0, 0, 0)), (2, (1, 1, 0)), (3, (0, 0, 1)), (4, (1, 1, 1)),
                   (5, (0, 0, 2)), (6, (1, 1, 2))):
        ops.node(tag, *[float(v) for v in c])
    for n in (1, 3, 5):
        ops.fix(n, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 1, 1, 2, 1.0, 1000.0, 400.0, 1.0, 1.0, 1.0, 1)
    ops.section("Elastic", 1, 1000.0, 1.0, 1.0, 1.0, 400.0, 1.0)
    ops.beamIntegration("Legendre", 1, 1, 3)
    ops.element("dispBeamColumn", 2, 3, 4, 1, 1)
    ops.element("forceBeamColumn", 3, 5, 6, 1, 1)
    ops.recorder("ladruno", path, "-N", "displacement", "-kind", "transient")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 4, 6):
        ops.load(n, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)
    _static(1, 1.0, system="BandGeneral")
    ops.wipe()
    axis = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    frames = 0
    with h5py.File(path, "r") as f:
        stage = [k for k in f if k.startswith("MODEL_STAGE")][0]
        assert str(lf._attr(f[stage], "KIND")) == "transient"
        la = f.get(f"{stage}/MODEL/LOCAL_AXES")
        if la is None or len(la) == 0:
            pytest.skip("no element answers the localAxes response in this build")
        for gname in la:
            ids = la[gname]["ID"][...].ravel()
            frame = la[gname]["FRAME"][...]
            assert frame.shape == (ids.shape[0], 4)
            for q in frame:
                assert abs(np.linalg.norm(q) - 1.0) <= 1e-9
                R = _quat_rows(q)
                assert np.max(np.abs(R @ R.T - np.eye(3))) <= 1e-9
                assert np.max(np.abs(R[0, :] - axis)) <= 1e-9   # local x = element axis
                frames += 1
    assert frames > 0
