"""EnergyBalance -v2: absorbing-boundary injection work in E_inject.

-v2 model block:

    time  KE_ele KE_nod IE DW_ele DW_nod ULW  [E_inject] ...  RES ERR%

Absorbing boundaries apply their incident-wave forcing inside the element
resisting force, so the six-column IE integral contains that injection work
with resisting-force sign. ASDAbsorbingBoundary2D/3D and LysmerTriangle
publish it at commit; -v2 reports IE without it and the injected energy as
E_inject. Checked, for an ASDAbsorbingBoundary2D compliant base (bottom,
Ricker velocity) and an ASDAbsorbingBoundary2D lateral free-field column
(initial velocity), and a LysmerTriangle base driven by
eleLoad -type -lysmerVelocityLoader (Ricker velocity), Newmark average
acceleration:
  * the six-column IE carries a leak of at least 5 % of the energy scale;
  * the -v2 IE is below 5 % of the energy scale once the motion has
    settled (ASD cases), and RES is the same in both layouts (Lysmer case);
  * E_inject equals IE_v2 - IE_six_column (the same leak through two paths);
  * a zero-amplitude control gives E_inject ~ 0.
"""
import math

import numpy as np

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def run_transient(efile, v2, nstep, dt, system="BandGeneral"):
    rec = ["EnergyBalance", "-file", efile, "-time"]
    if v2:
        rec.append("-v2")
    ops.recorder(*rec)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(system)
    ops.test("NormDispIncr", 1.0e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0
    ops.wipe()
    return np.atleast_2d(np.loadtxt(efile))


def check_injection(d_six, d_v2, d_ctl, leak_fraction=0.05, ie_settles=True):
    # E_inject is the first channel column (col 7) whenever an absorbing
    # element is present: the recorder writes E_inject before any other
    # channel.
    assert d_v2.shape[1] >= 10
    ie_six = d_six[-1, 2]
    ie_v2 = d_v2[-1, 3]
    e_inject = d_v2[-1, 7]
    eref = max(np.abs(d_v2[:, 1:8]).max(), 1e-30)
    assert abs(ie_six) > leak_fraction * eref, (ie_six, eref)
    if ie_settles:
        assert abs(ie_v2) < 0.05 * eref, (ie_v2, eref)
    diff = ie_v2 - ie_six
    assert abs(e_inject) > 0.0
    assert abs(e_inject - diff) < 0.1 * max(abs(diff), 1e-30), (e_inject, diff)
    assert abs(d_ctl[-1, 7]) < 1e-9, d_ctl[-1, 7]


# ---------------------------------------------------------------------------
# ASDAbsorbingBoundary2D, bottom: compliant-base incident velocity
# ---------------------------------------------------------------------------
def build_asd_column(amp_scale, n_layers=3):
    E, nu, rho = 2.0 * 2000.0 * 200.0 ** 2 * 1.25, 0.25, 2000.0
    G = 2000.0 * 200.0 ** 2
    dt = 1.0e-3
    nstep = 60
    t = np.arange(nstep + 1) * dt
    f0, t0, amp = 10.0, 0.02, 0.05
    a = (np.pi * f0 * (t - t0)) ** 2
    vals = (amp * amp_scale * (1.0 - 2.0 * a) * np.exp(-a)).tolist()

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1001, 0.0, -1.0)
    ops.node(1002, 1.0, -1.0)
    ops.fix(1001, 1, 1)
    ops.fix(1002, 1, 1)
    for row in range(n_layers + 1):
        ops.node(1 + 2 * row, 0.0, float(row))
        ops.node(2 + 2 * row, 1.0, float(row))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for row in range(n_layers):
        n1, n2 = 1 + 2 * row, 2 + 2 * row
        n3, n4 = 2 + 2 * (row + 1), 1 + 2 * (row + 1)
        ops.element("quad", 200 + row, n1, n2, n3, n4, 1.0, "PlaneStrain", 1)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *vals)
    ops.element("ASDAbsorbingBoundary2D", 101, 1001, 1002, 2, 1,
                G, nu, rho, 1.0, "B", "-fx", 1)
    ops.setParameter("-val", 1, "-ele", 101, "stage")
    return dt, nstep


def run_asd_column(tmp_path, fname, v2, amp_scale):
    dt, nstep = build_asd_column(amp_scale)
    return run_transient(str(tmp_path / fname), v2, nstep, dt)


def test_asd_absorbing_injection_closure(tmp_path):
    d_six = run_asd_column(tmp_path, "asd_six.txt", v2=False, amp_scale=1.0)
    d_v2 = run_asd_column(tmp_path, "asd_v2.txt", v2=True, amp_scale=1.0)
    d_ctl = run_asd_column(tmp_path, "asd_ctl.txt", v2=True, amp_scale=0.0)
    check_injection(d_six, d_v2, d_ctl)


# ---------------------------------------------------------------------------
# ASDAbsorbingBoundary2D, lateral: free-field column transfer
#
# Soil quad column at x in [0, 1] and a free-field node column at x = -1 with
# one "L" boundary element per layer. Driven by an initial velocity on both
# DOFs of every free node (both axes must move for the cross-axis transfer
# to do first-order work). Mass-proportional Rayleigh damping lets the
# free-field column decay as well.
# ---------------------------------------------------------------------------
def build_asd_lateral_column(amp_scale, n_layers=3):
    E, nu, rho = 2.0 * 2000.0 * 200.0 ** 2 * 1.25, 0.25, 2000.0
    G = 2000.0 * 200.0 ** 2
    dt = 1.0e-3
    nstep = 300
    v0 = 0.5 * amp_scale

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for row in range(n_layers + 1):
        ops.node(1 + 2 * row, 0.0, float(row))
        ops.node(2 + 2 * row, 1.0, float(row))
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    for row in range(n_layers + 1):
        ops.node(1000 + row, -1.0, float(row))
    ops.fix(1000, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for row in range(n_layers):
        n1, n2 = 1 + 2 * row, 2 + 2 * row
        n3, n4 = 2 + 2 * (row + 1), 1 + 2 * (row + 1)
        ops.element("quad", 200 + row, n1, n2, n3, n4, 1.0, "PlaneStrain", 1)
    for row in range(n_layers):
        ops.element("ASDAbsorbingBoundary2D", 100 + row,
                    1000 + row, 1000 + row + 1,
                    1 + 2 * row, 1 + 2 * (row + 1),
                    G, nu, rho, 1.0, "L")
    for row in range(n_layers):
        ops.setParameter("-val", 1, "-ele", 100 + row, "stage")
    for row in range(1, n_layers + 1):
        for nd in (1 + 2 * row, 2 + 2 * row, 1000 + row):
            ops.setNodeVel(nd, 1, v0, "-commit")
            ops.setNodeVel(nd, 2, v0, "-commit")
    ops.rayleigh(20.0, 0.0, 0.0, 0.0)
    return dt, nstep


def run_asd_lateral_column(tmp_path, fname, v2, amp_scale):
    dt, nstep = build_asd_lateral_column(amp_scale)
    return run_transient(str(tmp_path / fname), v2, nstep, dt)


def test_asd_lateral_freefield_closure(tmp_path):
    d_six = run_asd_lateral_column(tmp_path, "lat_six.txt", v2=False, amp_scale=1.0)
    d_v2 = run_asd_lateral_column(tmp_path, "lat_v2.txt", v2=True, amp_scale=1.0)
    d_ctl = run_asd_lateral_column(tmp_path, "lat_ctl.txt", v2=True, amp_scale=0.0)
    check_injection(d_six, d_v2, d_ctl)


# ---------------------------------------------------------------------------
# LysmerTriangle base driven by eleLoad -type -lysmerVelocityLoader
#
# 1 x 1 x NL soil column of stdBrick, x-polarized shear (y and z fixed), base
# face covered by two LysmerTriangle (stage 0, pure dashpot), incident
# x-velocity given by a Ricker pulse through the velocity loader.
# ---------------------------------------------------------------------------
def build_lysmer_column(amp_scale, nl=5):
    rho, vs, nu = 2000.0, 200.0, 0.25
    G = rho * vs * vs
    E = 2.0 * G * (1.0 + nu)
    vp = vs * math.sqrt((2.0 - 2.0 * nu) / (1.0 - 2.0 * nu))
    dt = 1.0e-3
    nstep = 200
    f0, t0 = 10.0, 0.1
    vals = []
    for i in range(nstep + 1):
        a = (math.pi * f0 * (i * dt - t0)) ** 2
        vals.append(amp_scale * 0.1 * (1.0 - 2.0 * a) * math.exp(-a))

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for k in range(nl + 1):
        b = 4 * k
        ops.node(b + 1, 0.0, 0.0, float(k))
        ops.node(b + 2, 1.0, 0.0, float(k))
        ops.node(b + 3, 1.0, 1.0, float(k))
        ops.node(b + 4, 0.0, 1.0, float(k))
    for n in range(1, 4 * (nl + 1) + 1):
        ops.fix(n, 0, 1, 1)
    for k in range(nl):
        b, t = 4 * k, 4 * (k + 1)
        ops.element("stdBrick", k + 1, b + 1, b + 2, b + 3, b + 4,
                    t + 1, t + 2, t + 3, t + 4, 1)
    ops.element("LysmerTriangle", 101, 1, 2, 3, rho, vp, vs)
    ops.element("LysmerTriangle", 102, 1, 3, 4, rho, vp, vs)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *vals)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 101, 102, "-type", "-lysmerVelocityLoader", 1)
    return dt, nstep


def run_lysmer_column(tmp_path, fname, v2, amp_scale):
    dt, nstep = build_lysmer_column(amp_scale)
    return run_transient(str(tmp_path / fname), v2, nstep, dt)


def test_lysmer_velocity_loader_injection_closure(tmp_path):
    d_six = run_lysmer_column(tmp_path, "lys_six.txt", v2=False, amp_scale=1.0)
    d_v2 = run_lysmer_column(tmp_path, "lys_v2.txt", v2=True, amp_scale=1.0)
    d_ctl = run_lysmer_column(tmp_path, "lys_ctl.txt", v2=True, amp_scale=0.0)
    # the column still holds strain energy at the end of the run, so IE_v2
    # is not expected to vanish here
    check_injection(d_six, d_v2, d_ctl, ie_settles=False)
    # the incident wave puts energy into the column
    assert d_v2[:, 7].max() > 0.0
    # moving L between IE and E_inject leaves RES unchanged
    scale = np.abs(d_six[:, 1:5]).max()
    np.testing.assert_allclose(d_v2[:, -2], d_six[:, 5], rtol=0.0, atol=1e-6 * scale)
