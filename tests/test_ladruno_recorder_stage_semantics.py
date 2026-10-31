"""recorder ladruno: when a MODEL_STAGE opens, and what it is labelled.

  * a moved domain-change stamp opens a new MODEL_STAGE only when the
    node/element/pressure-constraint SET changed: a pattern holding SPs does not
    copy the model again; a new node does.
  * the energy integrals are recorder state integrated on every commit: they
    continue across a stage change (no rate x t jump) and do not depend on the
    -T sampling.
  * INFO carries RUN_ID and RUN_ID_SCOPE (LADRUNO_RUN_ID -> scope "user").
  * a user `domainChange` (OPS_domainChange in SRC/interpreter/
    OpenSeesMiscCommands.cpp) opens a stage on an unchanged model; a
    static/transient flip opens one without it; a barrier before the first
    record, or a repeated kind, opens none; `-kind` only sets the label; an
    `eigen` + `record` snapshot opens no stage; and a free-field deck shape
    (gravity, eigen, loadConst, a new node and element, domainChange,
    rayleigh, UniformExcitation) gives [static, transient].

The classic Tcl twin of the stage cases is test_ladruno_recorder_stage_kind_tcl.py.
"""
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

h5py = pytest.importorskip("h5py")

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


def _truss():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 1.0e3)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 1.0, 0.0); ops.fix(2, 0, 1); ops.mass(2, 1.0, 0.0)
    ops.node(3, 2.0, 0.0); ops.fix(3, 0, 1); ops.mass(3, 1.0, 0.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.element("Truss", 2, 2, 3, 1.0, 1)


def _stages(path):
    with h5py.File(path, "r") as f:
        return sorted(k for k in f if k.startswith("MODEL_STAGE"))


def _static(dlam):
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
    ops.integrator("LoadControl", dlam); ops.analysis("Static")


def test_sp_pattern_does_not_start_a_new_stage(tmp_path):
    """A new pattern with an imposed displacement bumps the stamp (Domain::addSP_
    Constraint) but changes no node/element: one MODEL_STAGE, rows continue."""
    out = str(tmp_path / "sp.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(2) == 0
    ops.pattern("Plain", 2, 1)
    ops.sp(3, 1, 0.01)
    assert ops.analyze(2) == 0
    ops.wipe()
    stages = _stages(out)
    assert len(stages) == 1, stages
    with h5py.File(out, "r") as f:
        assert f[f"{stages[0]}/RESULTS/ON_NODES/DISPLACEMENT/DATA"].shape[0] == 4


def test_topology_change_starts_a_new_stage(tmp_path):
    out = str(tmp_path / "topo.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(2) == 0
    ops.node(99, 5.0, 5.0); ops.fix(99, 1, 1)          # a new node: topology changed
    assert ops.analyze(2) == 0
    ops.wipe()
    assert len(_stages(out)) == 2, _stages(out)


def _dynamic_energy(out, extra, steps_a=40, steps_b=40, add_node=False):
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement", "-G", "energy", *extra)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
    assert ops.analyze(steps_a, 0.01) == 0
    if add_node:
        ops.node(99, 5.0, 5.0); ops.fix(99, 1, 1)
    assert ops.analyze(steps_b, 0.01) == 0
    ops.wipe()
    rows = []
    with h5py.File(out, "r") as f:
        for st in sorted((k for k in f if k.startswith("MODEL_STAGE")),
                         key=lambda k: int(k[len("MODEL_STAGE["):-1])):
            g = f[f"{st}/RESULTS/ON_DOMAIN/energyBalance"]
            rows.append((g["TIME"][...], g["DATA"][:, 0, :]))   # [T x 6] KE,IE,DW,ULW,RES,ERR
    return rows


def test_energy_continues_across_a_stage_change(tmp_path):
    rows = _dynamic_energy(str(tmp_path / "energy.ladruno"), [], add_node=True)
    assert len(rows) == 2, len(rows)
    (t1, e1), (t2, e2) = rows
    ie_last, ie_next = e1[-1, 1], e2[0, 1]
    ulw_last, ulw_next = e1[-1, 3], e2[0, 3]
    # the next commit adds ONE step of work, not a restart near rate x t_total
    step_ulw = abs(e1[-1, 3] - e1[-2, 3])
    assert abs(ulw_next - ulw_last) < 5.0 * step_ulw + 1e-12, (ulw_last, ulw_next, step_ulw)
    assert ie_next > 0.5 * ie_last, (ie_last, ie_next)
    # and the closure still holds: |RES| small next to the external work
    assert abs(e2[-1, 4]) < 0.05 * abs(e2[-1, 3]) + 1e-12, e2[-1]


def test_energy_does_not_depend_on_T_sampling(tmp_path):
    every = _dynamic_energy(str(tmp_path / "e1.ladruno"), ["-T", "nsteps", 1])
    sparse = _dynamic_energy(str(tmp_path / "e5.ladruno"), ["-T", "nsteps", 5])
    (t_a, e_a), (t_b, e_b) = every[0], sparse[0]
    # compare at the common sample times
    common = np.intersect1d(np.round(t_a, 9), np.round(t_b, 9))
    assert len(common) >= 10, common
    ia = np.isin(np.round(t_a, 9), common)
    ib = np.isin(np.round(t_b, 9), common)
    np.testing.assert_allclose(e_b[ib][:, 1:4], e_a[ia][:, 1:4], rtol=1e-9, atol=1e-12)


def test_info_carries_a_run_id(tmp_path):
    out = str(tmp_path / "runid.ladruno")
    mod = ops.__name__
    here = os.path.dirname(os.path.abspath(ops.__file__))
    moddir = os.path.dirname(here) if "." in mod else here
    script = "\n".join([
        "import os, sys",
        f"sys.path.insert(0, {moddir!r})",
        f"os.add_dll_directory({here!r})" if hasattr(os, "add_dll_directory") else "",
        f"import {mod} as ops" if "." not in mod
        else f"import {mod.rsplit('.', 1)[0]}.opensees as ops",
        textwrap.dedent(f"""
            ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
            ops.node(1, 0.0, 0.0); ops.node(2, 1.0, 0.0); ops.fix(1, 1, 1); ops.fix(2, 0, 1)
            ops.uniaxialMaterial("Elastic", 1, 1.0e3); ops.element("Truss", 1, 1, 2, 1.0, 1)
            ops.recorder("ladruno", {out!r}, "-N", "displacement")
            ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(2, 1.0, 0.0)
            ops.system("FullGeneral"); ops.integrator("LoadControl", 1.0)
            ops.algorithm("Linear"); ops.analysis("Static"); ops.analyze(1); ops.wipe()
        """),
    ])
    launcher = ("SLURM_JOB_ID", "SLURM_STEP_ID", "OMPI_MCA_ess_base_jobid",
                "PMI_SIZE", "PMI_RANK", "OMPI_COMM_WORLD_SIZE", "OMPI_COMM_WORLD_RANK",
                "SLURM_NTASKS", "SLURM_PROCID")
    env = {k: v for k, v in os.environ.items() if k not in launcher}
    env.update(LADRUNO_RUN_ID="run-abc", HDF5_USE_FILE_LOCKING="FALSE")
    r = subprocess.run([sys.executable, "-c", script], env=env,
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr

    def s(v):
        v = v.flat[0] if hasattr(v, "flat") else v
        return v.decode() if isinstance(v, bytes) else str(v)

    with h5py.File(out, "r") as f:
        assert s(f["INFO"].attrs["RUN_ID"]) == "run-abc"
        assert s(f["INFO"].attrs["RUN_ID_SCOPE"]) == "user"


# --- user stage barriers and per-stage KIND -----------------------------------

def _kinds(path):
    """[(KIND, recorded rows)] per MODEL_STAGE in id order, and INFO's minor."""
    def s(v):
        v = np.atleast_1d(v)[0]
        return v.decode() if isinstance(v, bytes) else str(v)
    with h5py.File(path, "r") as f:
        ids = sorted(int(k[len("MODEL_STAGE["):-1]) for k in f if k.startswith("MODEL_STAGE["))
        out = [(s(f[f"MODEL_STAGE[{i}]"].attrs["KIND"]),
                f[f"MODEL_STAGE[{i}]/RESULTS/ON_NODES/DISPLACEMENT/DATA"].shape[0])
               for i in ids]
        minor = int(np.atleast_1d(f["INFO"].attrs["FORMAT_MINOR_VERSION"])[0])
    return out, minor


def _transient():
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")


def _gravity_then_transient(out, barrier):
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(4) == 0
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    if barrier:
        ops.domainChange()
    ops.timeSeries("Constant", 2); ops.pattern("Plain", 2, 2); ops.load(2, 1.0, 0.0)
    _transient()
    assert ops.analyze(5, 0.01) == 0
    ops.wipe()
    return _kinds(out)


def test_user_domainchange_splits_static_then_transient(tmp_path):
    """a user `domainChange` barrier opens a stage on an unchanged model, and
    each stage carries the kind of the analysis that ran in it."""
    got, minor = _gravity_then_transient(str(tmp_path / "b.ladruno"), barrier=True)
    assert got == [("static", 4), ("transient", 5)], got
    assert minor == 1


def test_kind_flip_splits_without_domainchange(tmp_path):
    """No stage mixes kinds: a static -> transient flip opens a stage by itself."""
    got, _ = _gravity_then_transient(str(tmp_path / "a.ladruno"), barrier=False)
    assert got == [("static", 4), ("transient", 5)], got


def test_domainchange_splits_two_static_stages(tmp_path):
    out = str(tmp_path / "ss.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(4) == 0
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 3); ops.pattern("Plain", 3, 3); ops.load(2, 0.5, 0.0)
    ops.domainChange()
    assert ops.analyze(3) == 0
    ops.wipe()
    got, _ = _kinds(out)
    assert got == [("static", 4), ("static", 3)], got


def test_repeated_kind_and_first_record_open_no_stage(tmp_path):
    """A barrier before the first record, and static -> wipeAnalysis -> static,
    open nothing; a transient-only run reads transient."""
    out = str(tmp_path / "rep.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.domainChange()
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(2) == 0
    ops.wipeAnalysis()
    _static(0.25)
    assert ops.analyze(2) == 0
    ops.wipe()
    got, _ = _kinds(out)
    assert got == [("static", 4)], got

    out = str(tmp_path / "tr.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _transient()
    assert ops.analyze(3, 0.01) == 0
    ops.wipe()
    got, _ = _kinds(out)
    assert got == [("transient", 3)], got


def test_kind_override_is_a_label_the_flip_still_splits(tmp_path):
    """`-kind static` sets the label; static -> transient still opens a stage."""
    out = str(tmp_path / "ov.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement", "-kind", "static")
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(4) == 0
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    _transient()
    assert ops.analyze(5, 0.01) == 0
    ops.wipe()
    got, _ = _kinds(out)
    assert got == [("static", 4), ("static", 5)], got


def test_eigen_snapshot_opens_no_stage(tmp_path):
    """`eigen; record` then a static analysis: ONE static stage (the classic-Tcl
    twin is the eigen_snapshot case of test_ladruno_recorder_stage_kind_tcl.py)."""
    out = str(tmp_path / "eig.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.eigen("-fullGenLapack", 1)
    ops.record()
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(4) == 0
    ops.wipe()
    got, _ = _kinds(out)
    assert got == [("static", 5)], got


def test_freefield_shape_transient_stage_is_kind_transient(tmp_path):
    """free-field deck shape: the recorder is created before both stages and
    samples on ``-T dt``; gravity runs Static, `eigen` feeds rayleigh while
    the static analysis is alive, the stage close is loadConst + wipeAnalysis,
    and the earthquake stage adds a node and an element (stamp and
    fingerprint move), fires the `domainChange` barrier, binds rayleigh,
    drives a UniformExcitation and runs Transient. The transient stage is
    KIND transient. The classic-Tcl twin is the freefield_shape case of
    test_ladruno_recorder_stage_kind_tcl.py."""
    out = str(tmp_path / "ff.ladruno")
    _truss()
    ops.recorder("ladruno", out, "-N", "displacement", "-T", "dt", 0.02)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
    _static(0.25)
    assert ops.analyze(4) == 0
    lam = ops.eigen("-fullGenLapack", 1)
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    ops.node(4, 3.0, 0.0); ops.fix(4, 0, 1); ops.mass(4, 1.0, 0.0)
    ops.element("Truss", 3, 3, 4, 1.0, 1)
    ops.domainChange()
    w = lam[0] ** 0.5
    ops.rayleigh(0.05 * w, 0.0, 0.05 / w, 0.0)
    ops.timeSeries("Path", 4, "-dt", 0.01, "-values",
                   0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0, -0.5, 0.0, 0.5, 1.0, 0.5)
    ops.pattern("UniformExcitation", 4, 1, "-accel", 4)
    _transient()
    assert ops.analyze(10, 0.01) == 0
    ops.wipe()
    got, minor = _kinds(out)
    # 4 static steps 0.25 apart all sample; 10 transient steps of 0.01 on the
    # 0.02 grid sample at 0.01, 0.03, 0.05, 0.07, 0.09.
    assert got == [("static", 4), ("transient", 5)], got
    assert minor == 1
