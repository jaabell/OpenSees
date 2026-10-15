"""ZeroLengthND assembles the full, possibly unsymmetric, material tangent.

One zeroLengthND element in 2D (ndf = 3, default orientation, so the
transformation is the identity and the element strain is u2 - u1 =
(eps_xx, eps_yy, gamma_xy)). Node 1 is fixed, node 2 keeps ux and the
rotation free, so the assembled system is the 2x2 block of the plane strain
material tangent on (eps_xx, gamma_xy).

The material is DruckerPrager in plane strain, perfectly plastic
(Kinf = Ko = H = 0, theta = 1), non-associated (rho = 0.4, rho_bar = 0.1),
with K = 30000, G = 10000, sigma_y = 5. A pure shear trial strain
gamma = 1.0e-3 (yield at sigma_y / (sqrt(3) G) = 2.887e-4) gives
I1 = 0, n = (0, 0, 0, 1/sqrt(2), 0, 0), a single active cone and the
consistent tangent

    Cep = Ce - (2G n + 3K rho_bar I)(2G n + 3K rho I)^T / d
             - (4 G^2 dgamma / |s_tr|) (IIdev - n n^T),
    d = 2G + 9 K rho rho_bar = 30800,

so that
    C(xx,xx) = K + 4G/3 - 9 K^2 rho rho_bar / d - (2/3) 4 G^2 dgamma / |s_tr|
             = 26655.196
    C(xx,xy) = -3 sqrt(2) G K rho_bar / d = -4132.442
    C(xy,xx) = -3 sqrt(2) G K rho     / d = -16529.769
    C(xy,xy) = G - 2 G^2 / d              =  3506.494

The element used to assemble only the lower triangle and mirror it, so both
off-diagonal terms came out as C(xy,xx) = -16529.769.
"""
import math

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


K_BULK = 30000.0
G_SHEAR = 10000.0
SIG_Y = 5.0
RHO = 0.4
RHO_BAR = 0.1


def _closed_form(gamma):
    d = 2.0 * G_SHEAR + 9.0 * K_BULK * RHO * RHO_BAR
    s_tr = math.sqrt(2.0) * G_SHEAR * gamma
    dgamma = (s_tr - math.sqrt(2.0 / 3.0) * SIG_Y) / d
    soft = 4.0 * G_SHEAR**2 * dgamma / s_tr
    c00 = (K_BULK + 4.0 * G_SHEAR / 3.0
           - 9.0 * K_BULK**2 * RHO * RHO_BAR / d - 2.0 / 3.0 * soft)
    c02 = -3.0 * math.sqrt(2.0) * G_SHEAR * K_BULK * RHO_BAR / d
    c20 = -3.0 * math.sqrt(2.0) * G_SHEAR * K_BULK * RHO / d
    c22 = G_SHEAR - 2.0 * G_SHEAR**2 / d
    return np.array([[c00, c02], [c20, c22]])


def _tangent_at_shear(gamma):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)
    ops.nDMaterial("DruckerPrager", 1, K_BULK, G_SHEAR, SIG_Y, RHO, RHO_BAR,
                   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    ops.element("zeroLengthND", 1, 1, 2, 1)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0           # sizes the system at zero strain

    # trial state: pure shear gamma_xy = gamma, eps_xx = eps_yy = 0
    ops.setNodeDisp(2, 3, gamma)
    flat = np.asarray(ops.printA("-ret"), dtype=float)
    assert flat.size == 4
    return flat.reshape(2, 2, order="F")  # Matrix storage is column major


def test_plastic_shear_tangent_is_not_symmetrized():
    try:
        A = _tangent_at_shear(1.0e-3)
        ref = _closed_form(1.0e-3)
        # equation 0 is ux of node 2 (eps_xx), equation 1 its rotation
        # (gamma_xy)
        assert A[1, 0] == pytest.approx(ref[1, 0], rel=1.0e-6)
        assert A[0, 1] == pytest.approx(ref[0, 1], rel=1.0e-6)
        assert A[0, 0] == pytest.approx(ref[0, 0], rel=1.0e-6)
        assert A[1, 1] == pytest.approx(ref[1, 1], rel=1.0e-6)
        assert abs(A[0, 1] - A[1, 0]) > 1.0e4
    finally:
        ops.wipe()


def test_elastic_shear_tangent_unchanged():
    try:
        A = _tangent_at_shear(1.0e-4)   # below yield
        ref = np.array([[K_BULK + 4.0 * G_SHEAR / 3.0, 0.0],
                        [0.0, G_SHEAR]])
        assert np.allclose(A, ref, rtol=1.0e-12, atol=1.0e-9)
    finally:
        ops.wipe()
