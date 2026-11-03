"""recorder Monitor: the live SWMR HDF5 sink.

  * smoke: a 200-step linear SDOF transient streamed with `recorder Monitor`
    gives FORMAT "ladruno-monitor", COLUMNS ["node2.disp.dof1"], 200 frames,
    a strictly increasing STEP and finite FRAMES;
  * observation only: the streamed channel equals a Node recorder on the same
    node/dof to 1e-9 (the Node recorder writes 16 digits);
  * -every K keeps one commit in K;
  * the file reopens in SWMR read mode and refreshes (tests/ladruno_monitor_reader.py);
  * -dof 0 is refused (DOFs are 1-based);
  * openseesmp rank 1 of 2 (PMI_SIZE=2, PMI_RANK=1, set before import in a
    subprocess) writes mon.part-1.h5, not a shared mon.h5;
  * the reader alone, on a synthetic sink: metadata, incremental cursor,
    rejection of a foreign HDF5 file.
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

h5py = pytest.importorskip("h5py")  # the reader needs h5py
from ladruno_monitor_reader import MonitorReader  # noqa: E402

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

NSTEPS = 200
DT = 0.01


def _sdof(sink=None, node_file=None, every=1):
    """A linear SDOF under a ramp load. Node 2 is the free mass."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.mass(2, 2.0)
    ops.uniaxialMaterial("Elastic", 1, 50.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0)
    if node_file is not None:
        ops.recorder("Node", "-file", node_file, "-time", "-precision", 16,
                     "-node", 2, "-dof", 1, "disp")
    if sink is not None:
        ops.recorder("Monitor", "-node", 2, "-dof", 1, "-resp", "disp",
                     "-sink", sink, "-every", every)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    for _ in range(NSTEPS):
        assert ops.analyze(1, DT) == 0
    ops.wipe()


def _s(v):
    return v.decode() if isinstance(v, bytes) else v


def test_monitor_smoke_and_node_recorder_parity(tmp_path):
    sink = str(tmp_path / "mon.h5")
    ref = str(tmp_path / "node.out")
    _sdof(sink=sink, node_file=ref)
    with h5py.File(sink, "r") as f:
        assert _s(f.attrs["FORMAT"]) == "ladruno-monitor"
        cols = [_s(c) for c in f["COLUMNS"][:]]
        step = f["STEP"][:]
        frames = f["FRAMES"][:]
    assert cols == ["node2.disp.dof1"]
    assert len(step) == NSTEPS and frames.shape == (NSTEPS, 1)
    assert np.all(np.diff(step) > 0)
    assert np.all(np.isfinite(frames)) and np.any(frames != 0.0)
    disp = np.loadtxt(ref)[:, 1]
    n = min(len(disp), NSTEPS)
    assert np.max(np.abs(disp[:n] - frames[:n, 0])) < 1e-9

    m = MonitorReader(sink).meta()
    assert m["columns"] == cols and m["nframes"] == NSTEPS
    with h5py.File(sink, "r", swmr=True) as f:
        ds = f["STEP"]
        ds.id.refresh()
        assert ds.shape[0] == NSTEPS


def test_monitor_every(tmp_path):
    sink = str(tmp_path / "every.h5")
    _sdof(sink=sink, every=10)
    with h5py.File(sink, "r") as f:
        step = f["STEP"][:]
    assert len(step) == NSTEPS // 10
    assert np.all(np.diff(step) == 10)


def test_monitor_rejects_dof_zero(tmp_path, capfd):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    try:
        ops.recorder("Monitor", "-node", 1, "-dof", 0, "-sink", str(tmp_path / "m.h5"))
    except Exception:
        pass                         # openseespy raises on a refused command
    log = "".join(capfd.readouterr())
    ops.wipe()
    assert "1-based" in log, log


def test_monitor_writes_a_per_rank_sink(tmp_path):
    sink = str(tmp_path / "mon.h5")
    mod = ops.__name__
    here = os.path.dirname(os.path.abspath(ops.__file__))
    moddir = os.path.dirname(here) if "." in mod else here
    body = textwrap.dedent(f"""
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.uniaxialMaterial("Elastic", 1, 1.0e6)
        ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
        ops.node(2, 1.0, 0.0); ops.fix(2, 0, 1)
        ops.element("Truss", 1, 1, 2, 1.0, 1)
        ops.recorder("Monitor", "-node", 2, "-dof", 1, "-sink", {sink!r})
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(2, 1.0, 0.0)
        ops.system("FullGeneral"); ops.integrator("LoadControl", 0.5)
        ops.algorithm("Linear"); ops.analysis("Static"); ops.analyze(2)
        ops.wipe()
        print("MONITOR_CHILD_DONE", flush=True)
    """)
    script = "\n".join([
        "import os, sys",
        f"sys.path.insert(0, {moddir!r})",
        f"os.add_dll_directory({here!r})" if hasattr(os, "add_dll_directory") else "",
        f"import {mod} as ops" if "." not in mod
        else f"import {mod.rsplit('.', 1)[0]}.opensees as ops",
        body,
    ])
    launcher = ("PMI_SIZE", "PMI_RANK", "OMPI_COMM_WORLD_SIZE", "OMPI_COMM_WORLD_RANK",
                "SLURM_NTASKS", "SLURM_PROCID", "SLURM_STEP_ID")
    env = {k: v for k, v in os.environ.items() if k not in launcher}
    env.update(PMI_SIZE="2", PMI_RANK="1", HDF5_USE_FILE_LOCKING="FALSE")
    r = subprocess.run([sys.executable, "-c", script], env=env, capture_output=True,
                       text=True, timeout=180, cwd=str(tmp_path))
    assert "MONITOR_CHILD_DONE" in r.stdout, r.stdout + r.stderr
    assert os.path.exists(str(tmp_path / "mon.part-1.h5")), os.listdir(tmp_path)
    assert not os.path.exists(sink), os.listdir(tmp_path)


# --- the reader alone, on a synthetic sink -----------------------------------

def _write_sink(path, columns, nframes, units=""):
    with h5py.File(path, "w", libver="latest") as f:
        f.attrs["FORMAT"] = "ladruno-monitor"
        f.attrs["FORMAT_VERSION"] = 1
        f.attrs["GENERATOR"] = "test"
        if units:
            f.attrs["units"] = units
        vstr = h5py.special_dtype(vlen=str)
        f.create_dataset("COLUMNS", data=np.array(columns, dtype=object), dtype=vstr)
        step = np.arange(1, nframes + 1, dtype="i4")
        t = np.arange(1, nframes + 1, dtype="f8") * 0.01
        frames = np.column_stack([np.sin(t * (j + 1)) for j in range(len(columns))])
        f.create_dataset("STEP", data=step, maxshape=(None,), chunks=(256,))
        f.create_dataset("TIME", data=t, maxshape=(None,), chunks=(256,))
        f.create_dataset("FRAMES", data=frames.astype("f8"),
                         maxshape=(None, len(columns)), chunks=(256, len(columns)))


def test_reader_meta_and_incremental(tmp_path):
    p = str(tmp_path / "s.h5")
    cols = ["node2.disp.dof1", "node5.disp.dof1"]
    _write_sink(p, cols, 300, units="m")
    r = MonitorReader(p)
    m = r.meta()
    assert m["columns"] == cols
    assert m["nframes"] == 300 and m["nColumns"] == 2
    assert m["format"] == "ladruno-monitor" and m["units"] == "m"
    d0 = r.frames_since(0)
    assert d0["n"] == 300 and len(d0["step"]) == 300
    assert len(d0["rows"]) == 300 and len(d0["rows"][0]) == 2
    d1 = r.frames_since(290)
    assert d1["since"] == 290 and len(d1["step"]) == 10
    assert d1["step"][0] == 291
    assert r.frames_since(300)["step"] == []
    assert r.frames_since(999)["step"] == []


def test_reader_rejects_foreign_h5(tmp_path):
    p = str(tmp_path / "x.h5")
    with h5py.File(p, "w") as f:
        f.create_dataset("nope", data=[1, 2, 3])
    with pytest.raises((ValueError, KeyError)):
        MonitorReader(p).meta()


def test_reader_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        MonitorReader(str(tmp_path / "absent.h5")).meta()
