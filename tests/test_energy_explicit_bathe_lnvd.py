"""EnergyBalance -v2: local non-viscous damping work of ExplicitBathe in E_lnvd.

ExplicitBathe -lnvd alpha removes alpha |a| sign(v) from each solved sub-step
acceleration (Cundall 1987 local damping). That force sits outside every
element and node damping term, so the six-column RES drifts by its
dissipation. ExplicitBathe now accumulates alpha |F| |v| over each sub-step
and publishes it to the LNVD_WORK channel on commit.

Checked on a free-vibration truss bar (10 elements, p = 0.54, alpha = 0.4,
lumped mass, Diagonal system): E_lnvd is positive and monotone, equals the
six-column RES drift within 30 %, and the -v2 RES drift is below 20 % of the
six-column drift. A run without -lnvd shows no E_lnvd column.

The integrator is created through the classic Tcl interpreter, where
ExplicitBathe is registered on this base, so the test needs the OpenSees
executable (OPENSEES_EXE, or OpenSees.exe next to the Python module).
"""
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _tcl_executable():
    env = os.environ.get("OPENSEES_EXE")
    if env and Path(env).is_file():
        return Path(env)
    here = Path(ops.__file__).resolve().parent
    for name in ("OpenSees.exe", "OpenSees"):
        exe = here / name
        if exe.is_file():
            return exe
    return None


EXE = _tcl_executable()
pytestmark = pytest.mark.skipif(EXE is None,
                                reason="no OpenSees Tcl executable (set OPENSEES_EXE)")

N, L, E, A, RHO, V0 = 10, 1.0, 100.0, 1.0, 1.0, 1.0


def _script(efile, v2, alpha):
    le = L / N
    dt = 0.2 * le / (E / RHO) ** 0.5
    lines = ["model basic -ndm 2 -ndf 2"]
    for i in range(N + 1):
        lines.append(f"node {i + 1} {i * le} 0.0")
    lines.append("fix 1 1 1")
    for i in range(1, N + 1):
        lines.append(f"fix {i + 1} 0 1")
    lines.append(f"uniaxialMaterial Elastic 1 {E}")
    for i in range(N):
        lines.append(f"element Truss {i + 1} {i + 1} {i + 2} {A} 1 -rho {RHO}")
    for i in range(1, N + 1):
        lines.append(f"setNodeVel {i + 1} 1 {V0} -commit")
    lines += [
        "constraints Transformation",
        "numberer Plain",
        "system Diagonal",
        "test NormDispIncr 1.0e-12 1",
        "algorithm Linear",
        # the integrator before the recorder: -v2 fixes its columns at the
        # first record from the channels declared by then
        "integrator ExplicitBathe 0.54" + (f" 0 -lnvd {alpha}" if alpha > 0 else ""),
        f"recorder EnergyBalance -file {{{efile}}} -time" + (" -v2" if v2 else ""),
        "analysis Transient",
        f"analyze 400 {dt}",
        "wipe",
    ]
    return "\n".join(lines) + "\n"


def _run(tmp_path, tag, v2, alpha):
    efile = (tmp_path / f"{tag}.txt").as_posix()
    tcl = tmp_path / f"{tag}.tcl"
    tcl.write_text(_script(efile, v2, alpha))
    r = subprocess.run([str(EXE), str(tcl)], capture_output=True, text=True,
                       timeout=300, cwd=str(tmp_path))
    data = np.atleast_2d(np.loadtxt(efile))
    m = re.findall(r"columns =( time)? \[model: ([^\]]+)\]", r.stdout + r.stderr)
    cols = None
    if m:
        has_time, names = m[-1]
        cols = (["time"] if has_time else []) + names.split()
    return data, cols


def test_lnvd_work_closes_the_balance(tmp_path):
    alpha = 0.4
    d_six, _ = _run(tmp_path, "six", v2=False, alpha=alpha)
    d_v2, cols = _run(tmp_path, "v2", v2=True, alpha=alpha)
    assert cols is not None, "no EnergyBalance column echo captured"

    # six-column layout: time KE IE DW ULW RES ERR% -> RES is col 5
    drift_six = d_six[-1, 5] - d_six[0, 5]
    e_lnvd = d_v2[:, cols.index("E_lnvd")]
    res_v2 = d_v2[:, cols.index("RES")]
    drift_v2 = res_v2[-1] - res_v2[0]

    assert e_lnvd[-1] > 0.0
    assert np.all(np.diff(e_lnvd) > -1e-14), "E_lnvd must be monotone"
    assert drift_six > 0.0
    assert 0.7 * drift_six < e_lnvd[-1] < 1.3 * drift_six, (e_lnvd[-1], drift_six)
    assert abs(drift_v2) < 0.2 * abs(drift_six), (drift_v2, drift_six)


def test_no_lnvd_column_without_local_damping(tmp_path):
    _, cols = _run(tmp_path, "plain", v2=True, alpha=0.0)
    assert cols is not None, "no EnergyBalance column echo captured"
    assert "E_lnvd" not in cols
