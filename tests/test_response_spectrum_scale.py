"""responseSpectrumAnalysis -scale multiplies every spectral ordinate.

The option was parsed into m_scale and never read, so the modal displacements
were identical with and without it. Modal displacements are linear in Sa, so
the reference is the same run without -scale: with -scale s every nodal
displacement equals s times the unscaled one. Unscaled and scaled runs use
the same model, spectrum and mode, so the ratio is exactly s up to rounding.

Cases: single mode (-mode k), the full mode loop (last mode left committed,
-scale given before the lists), the timeSeries form, and -scale 1.0 (unchanged).
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

TN = [0.02, 0.1, 0.3, 0.6, 1.0, 3.0]
SA = [0.40, 1.00, 1.20, 0.90, 0.60, 0.20]
NODES = (2, 3, 4)


def _setup():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.fix(1, 1)
    for i, m in enumerate((2.0, 1.5, 1.0), start=2):
        ops.node(i, float(i - 1))
        ops.mass(i, m)
    for e, k in enumerate((1200.0, 800.0, 500.0), start=1):
        ops.uniaxialMaterial("Elastic", e, k)
        ops.element("Truss", e, e, e + 1, 1.0, e)
    ops.eigen("-fullGenLapack", 3)
    ops.modalProperties()


def _field():
    return [ops.nodeDisp(n, 1) for n in NODES]


def _assert_scaled(base, scaled, s):
    assert any(abs(b) > 1e-12 for b in base), "reference field is zero"
    for b, v in zip(base, scaled):
        assert v == pytest.approx(s * b, rel=1e-10, abs=1e-15), (base, scaled, s)


@pytest.mark.parametrize("mode", [1, 2, 3])
def test_scale_single_mode(mode):
    _setup()
    ops.responseSpectrumAnalysis(1, "-Tn", *TN, "-Sa", *SA, "-mode", mode)
    base = _field()
    ops.responseSpectrumAnalysis(1, "-Tn", *TN, "-Sa", *SA, "-mode", mode,
                                 "-scale", 9.81)
    _assert_scaled(base, _field(), 9.81)


def test_scale_full_loop_scale_first():
    _setup()
    ops.responseSpectrumAnalysis(1, "-Tn", *TN, "-Sa", *SA)
    base = _field()
    ops.responseSpectrumAnalysis(1, "-scale", 0.25, "-Tn", *TN, "-Sa", *SA)
    _assert_scaled(base, _field(), 0.25)


def test_scale_timeseries_form():
    _setup()
    ops.timeSeries("Path", 7, "-time", *TN, "-values", *SA)
    ops.responseSpectrumAnalysis(7, 1, "-mode", 1)
    base = _field()
    ops.responseSpectrumAnalysis(7, 1, "-mode", 1, "-scale", 3.0)
    _assert_scaled(base, _field(), 3.0)


def test_scale_one_is_unchanged():
    _setup()
    ops.responseSpectrumAnalysis(1, "-Tn", *TN, "-Sa", *SA, "-mode", 2)
    base = _field()
    ops.responseSpectrumAnalysis(1, "-Tn", *TN, "-Sa", *SA, "-mode", 2,
                                 "-scale", 1.0)
    assert _field() == base
