"""recorder ladruno: MODEL_STAGE boundaries and KIND under classic Tcl.

Runs tests/ladruno_stage_kind_model.tcl (one case per run) with the classic Tcl
interpreter and checks the (KIND, recorded rows) list of every MODEL_STAGE, the
strictly increasing stage ids, INFO/FORMAT_MINOR_VERSION == 1 and the format
validator. The cases cover the hooks of SRC/tcl/commands.cpp: a user
`domainChange` opens a stage on an unchanged model, a static/transient flip
opens one, internal re-stamps (SP pattern, eleLoad) do not, and the transient
analysis that `eigen` builds by itself is not a running analysis until
`analysis` or `analyze`.

The interpreter is taken from the OPENSEES_TCL environment variable or found
as `OpenSees` on PATH; the test is skipped when none is available.
"""
import os
import shutil
import subprocess

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

h5py = pytest.importorskip("h5py")  # the reader needs h5py
import ladruno_format as lf  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DECK = os.path.join(HERE, "ladruno_stage_kind_model.tcl")
MARKER = "STAGE_KIND_MODEL OK"

EXPECT = {
    "dc_static_transient": [("static", 4), ("transient", 5)],
    "flip_static_transient": [("static", 4), ("transient", 5)],
    "dc_static_static": [("static", 4), ("static", 3)],
    "repeat_static": [("static", 7)],
    "internal_restamp": [("static", 6)],
    "dc_first": [("static", 4)],
    "transient_only": [("transient", 5)],
    "record_before_analysis": [("transient", 6)],
    "override_kind": [("static", 4), ("static", 5)],
    "override_kind_dc": [("transient", 4), ("transient", 3)],
    "eigen_snapshot": [("static", 5)],
    "eigen_snapshot_transient": [("transient", 6)],
    "flip_then_topology": [("static", 4), ("transient", 5), ("transient", 2)],
    "freefield_shape": [("static", 4), ("transient", 5)],
}


def _interpreter():
    exe = os.environ.get("OPENSEES_TCL") or shutil.which("OpenSees")
    if not exe:
        pytest.skip("no classic Tcl OpenSees interpreter (set OPENSEES_TCL)")
    return exe


def _s(v):
    v = np.atleast_1d(v)[0]
    return v.decode() if isinstance(v, bytes) else str(v)


@pytest.mark.parametrize("case", list(EXPECT))
def test_stage_kind_tcl(tmp_path, case):
    exe = _interpreter()
    out = str(tmp_path / f"{case}.ladruno").replace("\\", "/")
    env = dict(os.environ, HDF5_USE_FILE_LOCKING="FALSE")
    r = subprocess.run([exe, DECK, out, case], capture_output=True, text=True,
                       env=env, cwd=str(tmp_path), timeout=300)
    log = r.stdout + r.stderr
    assert f"{MARKER} {case}" in log, log
    with h5py.File(out, "r") as f:
        ids = sorted(int(k[len("MODEL_STAGE["):-1]) for k in f if k.startswith("MODEL_STAGE["))
        got = []
        for i in ids:
            g = f[f"MODEL_STAGE[{i}]"]
            d = g.get("RESULTS/ON_NODES/DISPLACEMENT/DATA")
            got.append((_s(g.attrs["KIND"]), 0 if d is None else int(d.shape[0])))
        minor = int(np.atleast_1d(f["INFO"].attrs["FORMAT_MINOR_VERSION"])[0])
    assert got == EXPECT[case], (got, ids)
    assert ids == sorted(set(ids))
    assert minor == 1
    assert lf.validate(out) == []
