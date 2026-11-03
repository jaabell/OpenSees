"""recorder ladruno: robustness of the write path.

  * modes after wipe: `eigen` -> `wipe` -> new model with `-N modesOfVibration`
    -> `analyze` must not exit the process. `wipe` does not reset the
    interpreter's numEigen, so the recorder gates its modal path on the domain
    spectrum (Domain::getNumEigenvalues), not on numEigen.
  * emptied region: removing every node of the `-R` region opens an EMPTY
    MODEL_STAGE (EMPTY_PARTITION=1, zero-length NODES/ID) and the run goes on;
    no Response of a deleted element is touched.
  * duplicate request (`-N displacement displacement`): recorded once, with one
    notice, and DATA keeps T rows.
  * `-G -T nsteps 2`: -G hands the option token back, so -T is honoured.
  * `-T dt 0.1` with LoadControl 0.1: accumulated time (0.30000000000000004)
    does not skip samples.
  * launcher environment: SLURM_NTASKS/SLURM_PROCID without an srun step is a
    sequential run (no part file); PMI_SIZE > 1 without PMI_RANK is refused.

The first two would exit or corrupt the interpreter on a defective build, so
they run in a SUBPROCESS and assert a marker printed after the last command.
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

DONE = "LADRUNO_CHILD_DONE"


def _child_script(body):
    """`body` (uses `ops`) wrapped to import THIS build of opensees."""
    mod = ops.__name__
    here = os.path.dirname(os.path.abspath(ops.__file__))
    moddir = os.path.dirname(here) if "." in mod else here
    lines = ["import os, sys", f"sys.path.insert(0, {moddir!r})"]
    if hasattr(os, "add_dll_directory"):
        lines.append(f"os.add_dll_directory({here!r})")
    lines += [f"import {mod} as ops" if "." not in mod
              else f"import {mod.rsplit('.', 1)[0]}.opensees as ops",
              textwrap.dedent(body), f"print({DONE!r}, flush=True)"]
    return "\n".join(lines)


def _run_child(body, tmp_path, env_extra=None):
    """Run `body` in a fresh interpreter on THIS build of opensees."""
    launcher = ("PMI_SIZE", "PMI_RANK", "OMPI_COMM_WORLD_SIZE", "OMPI_COMM_WORLD_RANK",
                "SLURM_NTASKS", "SLURM_PROCID", "SLURM_STEP_ID")
    env = {k: v for k, v in os.environ.items() if k not in launcher}
    env.update(env_extra or {}, HDF5_USE_FILE_LOCKING="FALSE")
    return subprocess.run([sys.executable, "-c", _child_script(body)], env=env,
                          capture_output=True, text=True, timeout=180,
                          cwd=str(tmp_path))


# Two x-direction trusses with lumped x-masses: 2 free DOFs, a well-posed
# generalized eigenproblem for -fullGenLapack.
TRUSS_MODEL = """
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.uniaxialMaterial("Elastic", 1, 1.0e6)
ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
ops.node(2, 1.0, 0.0); ops.fix(2, 0, 1); ops.mass(2, 1.0, 0.0)
ops.node(3, 2.0, 0.0); ops.fix(3, 0, 1); ops.mass(3, 1.0, 0.0)
ops.element("Truss", 1, 1, 2, 1.0, 1)
ops.element("Truss", 2, 2, 3, 1.0, 1)
"""

STATIC = """
ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(3, 1.0, 0.0)
ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0); ops.analysis("Static")
"""


def test_modes_after_wipe_do_not_exit_the_process(tmp_path):
    out = str(tmp_path / "modes.ladruno")
    body = TRUSS_MODEL + """
ops.eigen("-fullGenLapack", 1)       # leaves the interpreter's numEigen = 1
""" + TRUSS_MODEL + f"""
ops.recorder("ladruno", {out!r}, "-N", "displacement", "modesOfVibration")
""" + STATIC + """
ops.analyze(1)                       # fresh domain: no spectrum -> no modal write
ops.wipeAnalysis()
ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
ops.eigen("-fullGenLapack", 1)       # a real spectrum now exists
ops.record()                         # -> modes are written
ops.wipe()
"""
    r = _run_child(body, tmp_path)
    log = r.stdout + r.stderr
    assert DONE in r.stdout, f"process died (rc={r.returncode}):\n{log}"
    assert "Eigenvalues were never set" not in log, log
    with h5py.File(out, "r") as f:
        modes = []
        f.visit(lambda name: modes.append(name)
                if name.rsplit("/", 1)[-1].startswith("MODE_") else None)
    assert modes, "modesOfVibration wrote no MODE_<k> dataset after the real eigen"


def test_emptied_region_opens_an_empty_stage(tmp_path):
    out = str(tmp_path / "region.ladruno")
    body = f"""
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.uniaxialMaterial("Elastic", 1, 1.0e6)
# region 1: two trusses (nodes 1-2-3); a separate truss 10-11 keeps the model alive
ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
ops.node(2, 1.0, 0.0); ops.fix(2, 0, 1)
ops.node(3, 2.0, 0.0); ops.fix(3, 0, 1)
ops.node(10, 0.0, 5.0); ops.fix(10, 1, 1)
ops.node(11, 1.0, 5.0); ops.fix(11, 0, 1)
ops.element("Truss", 1, 1, 2, 1.0, 1)
ops.element("Truss", 2, 2, 3, 1.0, 1)
ops.element("Truss", 3, 10, 11, 1.0, 1)
ops.region(1, "-ele", 1, 2)
ops.recorder("ladruno", {out!r}, "-R", 1, "-N", "displacement", "-E", "force")
ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(11, 1.0, 0.0)
ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
ops.integrator("LoadControl", 0.1); ops.analysis("Static")
ops.analyze(2)
# remove the whole region: its elements, its nodes' SPs, its nodes
ops.remove("element", 1); ops.remove("element", 2)
ops.remove("sp", 1, 1); ops.remove("sp", 1, 2)
ops.remove("sp", 2, 2); ops.remove("sp", 3, 2)
ops.remove("node", 1); ops.remove("node", 2); ops.remove("node", 3)
ops.analyze(3)                       # rebuild fails once, then 2 more commits
ops.wipe()
"""
    r = _run_child(body, tmp_path)
    log = r.stdout + r.stderr
    assert DONE in r.stdout, f"process died (rc={r.returncode}):\n{log}"
    # a region with no nodes left is a valid EMPTY stage, not a failed rebuild
    assert "written as an empty partition" in log, log
    assert "could not be written" not in log, log
    with h5py.File(out, "r") as f:
        stages = sorted((k for k in f if k.startswith("MODEL_STAGE")),
                        key=lambda k: int(k[len("MODEL_STAGE["):-1]))
        assert len(stages) == 2, stages
        first, second = f[stages[0]], f[stages[1]]
        # the first stage is intact: 2 rows of displacement
        assert first["RESULTS/ON_NODES/DISPLACEMENT/DATA"].shape[0] == 2
        assert int(np.asarray(second.attrs["EMPTY_PARTITION"]).flat[0]) == 1
        assert second["MODEL/NODES/ID"].shape[0] == 0


def _stage(f):
    return f[[k for k in f if k.startswith("MODEL_STAGE")][0]]


@pytest.mark.parametrize("request_args", [
    ("-N", "displacement", "displacement"),
    ("-N", "displacement", "-E", "force", "force"),
])
def test_duplicate_request_is_recorded_once(tmp_path, capfd, request_args):
    """repeated -N/-E tokens are dropped at parse time with a notice (the sink
    also refuses a pre-existing group)."""
    out = str(tmp_path / "dup.ladruno")
    exec(TRUSS_MODEL, {"ops": ops})
    ops.recorder("ladruno", out, *request_args)
    exec(STATIC, {"ops": ops})
    ops.analyze(3)
    ops.wipe()
    log = "".join(capfd.readouterr())
    assert log.count("requested twice") == 1, log
    with h5py.File(out, "r") as f:
        g = _stage(f)["RESULTS/ON_NODES/DISPLACEMENT"]
        assert g["DATA"].shape[0] == 3, g["DATA"].shape          # 6 with two sinks
        assert g["TIME"].shape[0] == 3 and g["STEP"].shape[0] == 3
        if "force" in request_args:
            buckets = list(_stage(f)["RESULTS/ON_ELEMENTS/force"].values())
            assert buckets and all(b["DATA"].shape[0] == 3 for b in buckets)


def test_G_does_not_swallow_the_next_option(tmp_path, capfd):
    """`-G -T nsteps 2`: -T is honoured."""
    out = str(tmp_path / "g.ladruno")
    exec(TRUSS_MODEL, {"ops": ops})
    ops.recorder("ladruno", out, "-N", "displacement", "-G", "-T", "nsteps", 2)
    exec(STATIC, {"ops": ops})
    ops.analyze(4)                  # commits 1..4: records at 1 (first) and 3
    ops.wipe()
    with h5py.File(out, "r") as f:
        n = _stage(f)["RESULTS/ON_NODES/DISPLACEMENT/DATA"].shape[0]
    assert n == 2, f"-T nsteps 2 lost: {n} rows for 4 steps"


def test_T_dt_tolerates_accumulated_time(tmp_path):
    """LoadControl 0.1 x 10 with -T dt 0.1: time sums to 0.30000000000000004,
    0.4 - that = 0.0999...98 < 0.1, which an exact comparison would skip."""
    out = str(tmp_path / "dt.ladruno")
    exec(TRUSS_MODEL, {"ops": ops})
    ops.recorder("ladruno", out, "-N", "displacement", "-T", "dt", 0.1)
    exec(STATIC.replace('ops.integrator("LoadControl", 1.0)',
                        'ops.integrator("LoadControl", 0.1)'), {"ops": ops})
    ops.analyze(10)
    ops.wipe()
    with h5py.File(out, "r") as f:
        t = _stage(f)["RESULTS/ON_NODES/DISPLACEMENT/TIME"][...]
    assert len(t) == 10, t


# --- launcher environment ---------------------------------------------------

def test_sequential_run_inside_sbatch_is_not_partitioned(tmp_path):
    """`sbatch --ntasks=4` exports SLURM_NTASKS/SLURM_PROCID into the batch shell;
    without an srun step the run is sequential and keeps its filename."""
    out = str(tmp_path / "sbatch.ladruno")
    body = TRUSS_MODEL + f"""
ops.recorder("ladruno", {out!r}, "-N", "displacement")
""" + STATIC + "\nops.analyze(1)\nops.wipe()\n"
    r = _run_child(body, tmp_path, {"SLURM_NTASKS": "4", "SLURM_PROCID": "0"})
    assert DONE in r.stdout, r.stdout + r.stderr
    assert os.path.exists(out), os.listdir(tmp_path)
    with h5py.File(out, "r") as f:
        assert int(f["INFO"].attrs["PARTITIONED"].flat[0]) == 0


def test_size_without_rank_is_refused(tmp_path):
    """PMI_SIZE=4 with no PMI_RANK: refused loudly, no file (every rank would
    otherwise write part-0)."""
    out = str(tmp_path / "rank.ladruno")
    body = TRUSS_MODEL + f"""
ops.recorder("ladruno", {out!r}, "-N", "displacement")
""" + STATIC + "\nops.analyze(1)\nops.wipe()\n"
    r = _run_child(body, tmp_path, {"PMI_SIZE": "4"})
    log = r.stdout + r.stderr
    assert DONE in r.stdout, log
    assert "missing or not in" in log, log
    assert not [p for p in os.listdir(tmp_path) if p.endswith(".ladruno")], os.listdir(tmp_path)
