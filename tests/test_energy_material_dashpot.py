"""EnergyBalance recorder: a material dashpot is counted once (in IE).

A rate-dependent material (Viscous, Parallel(..., Viscous)) puts its viscous
force in getStress(), hence in the element resisting force F, so its power is
already in IE = int F.v dt. The element getDamp() also returns the material
damping tangent; booking v^T C v of it in DW as well would count the dashpot
twice and give RES = -(dashpot dissipation). DW is the power of the damping
force the element adds outside F (g = F_incInertia - F - M a); Rayleigh-only
elements keep v^T C v.

SDOF, m = 1, k = 100, c = 2 on zeroLength or truss, Newmark average
acceleration, a load pulse followed by free vibration (4000 steps of 1 ms):

  T1  Elastic + Viscous (one zeroLength, two materials), Parallel(Elastic,
      Viscous) -factors, and a separate dashpot element: DW = 0,
      IE = trapz(k u v) + trapz(c v^2), RES ~ 0, and the same RES and IE+DW
      as the Rayleigh twin (betaK k = c).
  T2  Rayleigh damping only (region betaK on the spring):
      DW = trapz(betaK k v^2), RES ~ 0.
  T3  Truss with Parallel(Elastic, Viscous) and Rayleigh betaK on the same
      element: DW holds only the Rayleigh part, IE the dashpot.
  T4  Nonlinear dashpot c|v|^a, a = 0.5: DW = 0, RES ~ 0.
  T5  -v2: DW_ele follows the six-column DW; RES identical.
"""
import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


M, K, C, DT, N = 1.0, 100.0, 2.0, 0.001, 4000


def run(tmp, name, model, ndm=1):
    d = tmp / name
    d.mkdir()
    e1 = str(d / "e.txt")
    e2 = str(d / "e2.txt")
    ops.wipe()
    if ndm == 1:
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0)
        ops.node(2, 0.0, "-mass", M)
        ops.fix(1, 1)
    else:
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 1.0, 0.0, "-mass", M, M)
        ops.fix(1, 1, 1)
        ops.fix(2, 0, 1)
    model()
    ops.timeSeries("Path", 1, "-time", 0.0, 0.05, 0.1, 100.0,
                   "-values", 0.0, 10.0, 0.0, 0.0)
    ops.pattern("Plain", 1, 1)
    if ndm == 1:
        ops.load(2, 1.0)
    else:
        ops.load(2, 1.0, 0.0)
    ops.recorder("EnergyBalance", "-time", "-precision", 17, "-file", e1)
    ops.recorder("EnergyBalance", "-time", "-precision", 17, "-v2", "-file", e2)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-14, 30)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    u = [0.0]
    v = [0.0]
    for _ in range(N):
        assert ops.analyze(1, DT) == 0
        u.append(ops.nodeDisp(2, 1))
        v.append(ops.nodeVel(2, 1))
    ops.wipe()
    return np.array(u), np.array(v), np.loadtxt(e1), np.loadtxt(e2)


def trapz(f):
    return float(np.sum(0.5 * (f[:-1] + f[1:]) * DT))


def check_closure(e, d_visc):
    """RES is the trapezoidal-quadrature residual of the Newmark run, about
    1e-7 of the energy scale here. Counting the dashpot twice would make
    RES = -D, about 0.5 of the scale."""
    ref = np.max(np.abs(e[:, 1:5]).sum(axis=1))
    assert abs(e[-1, 5]) < 1e-6 * ref, (e[-1], d_visc)


def rayleigh_twin():
    ops.uniaxialMaterial("Elastic", 1, K)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1, "-doRayleigh", 1)
    ops.region(1, "-ele", 1, "-rayleigh", 0.0, C / K, 0.0, 0.0)


def two_materials():
    ops.uniaxialMaterial("Elastic", 1, K)
    ops.uniaxialMaterial("Viscous", 2, C, 1.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, 2, "-dir", 1, 1, "-doRayleigh", 0)


def parallel_factors():
    ops.uniaxialMaterial("Elastic", 11, 1.0)
    ops.uniaxialMaterial("Viscous", 12, 1.0, 1.0)
    ops.uniaxialMaterial("Parallel", 3, 11, 12, "-factors", K, C)
    ops.element("zeroLength", 1, 1, 2, "-mat", 3, "-dir", 1)


def separate_dashpot():
    ops.uniaxialMaterial("Elastic", 1, K)
    ops.uniaxialMaterial("Viscous", 2, C, 1.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.element("zeroLength", 2, 1, 2, "-mat", 2, "-dir", 1)


DASHPOT_VARIANTS = {
    "two_materials": two_materials,
    "parallel_factors": parallel_factors,
    "separate_dashpot": separate_dashpot,
}


@pytest.mark.parametrize("variant", sorted(DASHPOT_VARIANTS))
def test_T1_linear_dashpot_booked_once(tmp_path, variant):
    u, v, e, e2 = run(tmp_path, variant, DASHPOT_VARIANTS[variant])
    d_visc = trapz(C * v * v)
    ie_expect = trapz(K * u * v) + d_visc
    assert d_visc > 1e-3
    assert abs(e[-1, 3]) < 1e-14 * d_visc, ("DW must be 0", e[-1])
    assert abs(e[-1, 2] - ie_expect) < 1e-9 * d_visc, (e[-1, 2], ie_expect)
    check_closure(e, d_visc)
    # same dynamics with the dashpot as Rayleigh damping (betaK k = c): same
    # RES (quadrature residual), and DW of the twin equals the dissipation
    ur, vr, er, _ = run(tmp_path, variant + "_rayleigh_twin", rayleigh_twin)
    assert np.max(np.abs(vr - v)) < 1e-12 * np.max(np.abs(v))
    assert abs(e[-1, 5] - er[-1, 5]) < 1e-12 * d_visc, (e[-1, 5], er[-1, 5])
    assert abs((e[-1, 2] + e[-1, 3]) - (er[-1, 2] + er[-1, 3])) < 1e-12 * d_visc


def test_T2_rayleigh_spring_unchanged_path(tmp_path):
    beta = 0.02  # c_eq = beta k = 2.0, the same as the dashpot

    def model():
        ops.uniaxialMaterial("Elastic", 1, K)
        ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1, "-doRayleigh", 1)
        ops.region(1, "-ele", 1, "-rayleigh", 0.0, beta, 0.0, 0.0)

    u, v, e, e2 = run(tmp_path, "rayleigh", model)
    dw_expect = trapz(beta * K * v * v)
    assert abs(e[-1, 3] - dw_expect) < 1e-10 * dw_expect, (e[-1, 3], dw_expect)
    assert abs(e[-1, 2] - trapz(K * u * v)) < 1e-10 * dw_expect
    check_closure(e, dw_expect)


def test_T3_truss_dashpot_plus_rayleigh(tmp_path):
    beta = 0.01

    def model():
        ops.uniaxialMaterial("Elastic", 11, 1.0)
        ops.uniaxialMaterial("Viscous", 12, 1.0, 1.0)
        ops.uniaxialMaterial("Parallel", 3, 11, 12, "-factors", K, C)
        ops.element("truss", 1, 1, 2, 1.0, 3, "-doRayleigh", 1)
        ops.region(1, "-ele", 1, "-rayleigh", 0.0, beta, 0.0, 0.0)

    u, v, e, e2 = run(tmp_path, "truss_mixed", model, ndm=2)
    d_visc = trapz(C * v * v)
    dw_expect = trapz(beta * K * v * v)
    assert abs(e[-1, 3] - dw_expect) < 1e-8 * dw_expect, (e[-1, 3], dw_expect, d_visc)
    assert abs(e[-1, 2] - (trapz(K * u * v) + d_visc)) < 1e-8 * d_visc
    check_closure(e, d_visc)


def test_T4_nonlinear_dashpot(tmp_path):
    a = 0.5

    def model():
        ops.uniaxialMaterial("Elastic", 1, K)
        ops.uniaxialMaterial("Viscous", 2, C, a)
        ops.element("zeroLength", 1, 1, 2, "-mat", 1, 2, "-dir", 1, 1)

    u, v, e, e2 = run(tmp_path, "nonlinear", model)
    d_visc = trapz(C * np.abs(v) ** a * np.abs(v))
    assert abs(e[-1, 3]) < 1e-14 * d_visc, e[-1]
    assert abs(e[-1, 2] - (trapz(K * u * v) + d_visc)) < 1e-8 * d_visc
    check_closure(e, d_visc)


def test_T5_v2_columns_follow(tmp_path):
    u, v, e, e2 = run(tmp_path, "v2", parallel_factors)
    # -v2: time KE_ele KE_nod IE DW_ele DW_nod ULW [channels] RES ERR
    assert abs(e2[-1, 4]) == 0.0 and abs(e2[-1, 5]) == 0.0, e2[-1]
    assert abs(e2[-1, 3] - e[-1, 2]) < 1e-12 * abs(e[-1, 2])
    assert abs(e2[-1, -2] - e[-1, 5]) < 1e-12 * abs(e[-1, 2])
