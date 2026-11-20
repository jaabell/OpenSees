"""EnergyBalance -v2: modal-damping dissipation in E_modal.

Modal damping (modalDamping) is applied by the integrator, outside any
element or node damping matrix, so the six-column RES drifts by the modal
dissipation. IncrementalIntegrator::commit() publishes v^T C_modal v,
integrated over each step, to the MODAL_WORK channel. Checked on a
free-vibration truss bar (20 elements, zeta = 0.05 on mode 1, Newmark average
acceleration): E_modal is positive and monotone, matches the six-column RES
drift within 30 %, and the -v2 RES drift is below 20 % of the six-column
drift.
"""
import math
import re

import numpy as np

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def v2_cols(capfd):
    """Column names from the recorder's echoed layout (the -v2 output is
    parsed by name). Reads the output captured since the last readouterr()."""
    cap = capfd.readouterr()
    text = cap.out + cap.err
    m = re.findall(r"columns =( time)? \[model: ([^\]]+)\]", text)
    assert m, "no EnergyBalanceRecorder column echo captured"
    has_time, names = m[-1]
    return (["time"] if has_time else []) + names.split()


# ---------------------------------------------------------------------------
# Modal damping: E_modal closes what the six-column layout leaves in RES
# ---------------------------------------------------------------------------
def run_modal_bar(tmp_path, fname, v2, zeta):
    N, L, E, A, rho, v0 = 20, 10.0, 100.0, 1.0, 1.0, 2.0
    efile = str(tmp_path / fname)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    le = L / N
    for i in range(N + 1):
        ops.node(i + 1, i * le, 0.0)
    ops.fix(1, 1, 1)
    for i in range(1, N + 1):
        ops.fix(i + 1, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(N):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", rho)
    for i in range(1, N + 1):
        ops.setNodeVel(i + 1, 1, v0, "-commit")
    ops.eigen(1)
    ops.modalDamping(zeta)
    rec = ["EnergyBalance", "-file", efile, "-time"]
    if v2:
        rec.append("-v2")
    ops.recorder(*rec)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    assert ops.analyze(400, 0.2 * le / math.sqrt(E / rho)) == 0
    ops.wipe()
    return np.atleast_2d(np.loadtxt(efile))


def test_modal_damping_closure(tmp_path, capfd):
    zeta = 0.05
    d_six = run_modal_bar(tmp_path, "modal_six.txt", v2=False, zeta=zeta)
    capfd.readouterr()
    d_v2 = run_modal_bar(tmp_path, "modal_v2.txt", v2=True, zeta=zeta)
    cols = v2_cols(capfd)

    # six-column layout: time KE IE DW ULW RES ERR% -> RES is col 5
    res_six = d_six[:, 5]
    drift_six = res_six[-1] - res_six[0]
    assert drift_six > 0.0, drift_six

    assert "E_modal" in cols, cols
    e_modal = d_v2[:, cols.index("E_modal")]
    res_v2 = d_v2[:, cols.index("RES")]
    drift_v2 = res_v2[-1] - res_v2[0]

    assert e_modal[-1] > 0.0
    assert np.all(np.diff(e_modal) > -1e-9)
    assert 0.7 * drift_six < e_modal[-1] < 1.3 * drift_six, (e_modal[-1], drift_six)
    assert abs(drift_v2) < 0.2 * abs(drift_six), (drift_v2, drift_six)
