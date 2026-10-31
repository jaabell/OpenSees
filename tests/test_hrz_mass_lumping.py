"""HRZ (Hinton-Rock-Zienkiewicz) lumping of the critical-time-step pencil.

Checks through `integrator CentralDifferenceLadruno -cfl -lump hrz`:
  * one 2-D Truss with consistent mass (rho L/6 [2 1; 1 2] per direction):
    HRZ conserves the mass (m/2 per node), so dt_cr = le/c = 0.1 equals the
    row sum, while the consistent diagonal (m/3 per node) gives
    2/sqrt(600) = 0.08165;
  * one 2-D Truss with lumped mass: HRZ, row sum and diagonal agree;
  * one 3-D elasticBeamColumn with consistent mass (rotational DOFs): HRZ is
    finite and positive and differs from the consistent diagonal (rel > 1e-3),
    so it did not fall back to the diagonal.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

CDL = "CentralDifferenceLadruno"


def _analysis(lump):
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl", "-lump", lump)
    ops.analysis("Transient")


def _truss_dtcr(lump, cmass):
    # E = 100, A = 1, L = 1, rho = 1: k = 100, m = 1, c = 10, le/c = 0.1
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 1.0, 0.0)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0, "-cMass", cmass)
    _analysis(lump)
    assert ops.analyze(1, 1.0e-3) == 0
    return ops.criticalTimeStep()


def _beam_dtcr(lump):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.node(2, 1.0, 0.0, 0.0)
    E, G, A, Jx, Iy, Iz, rho = 200.0, 80.0, 1.0, 1.0, 1.0, 1.0, 1.0
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
                "-mass", rho, "-cMass")
    _analysis(lump)
    assert ops.analyze(1, 1.0e-6) == 0
    return ops.criticalTimeStep()


def test_hrz_conserves_mass_on_consistent_truss():
    h = _truss_dtcr("hrz", 1)
    r = _truss_dtcr("rowsum", 1)
    d = _truss_dtcr("diagonal", 1)
    assert h == pytest.approx(0.1, rel=1e-9)
    assert r == pytest.approx(0.1, rel=1e-9)
    assert d == pytest.approx(2.0 / math.sqrt(600.0), rel=1e-9)


def test_hrz_equals_other_lumps_on_lumped_truss():
    values = [_truss_dtcr(lump, 0) for lump in ("hrz", "rowsum", "diagonal")]
    for v in values:
        assert v == pytest.approx(0.1, rel=1e-9)


def test_hrz_on_consistent_beam():
    h = _beam_dtcr("hrz")
    d = _beam_dtcr("diagonal")
    assert math.isfinite(h) and h > 0.0, h
    assert math.isfinite(d) and d > 0.0, d
    assert abs(h - d) / d > 1.0e-3, "HRZ %.12g equals the diagonal %.12g" % (h, d)
