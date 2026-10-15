"""Arc-length integrators report a failed domain update.

Model: one truss in 1D, L = 1, A = 1, node 1 fixed, reference load -1.0 at
node 2, uniaxialMaterial ReinforcingSteel (fy = 60, fu = 90, Es = 29000,
Esh = 1000, esh = 0.01, eult = 0.1). ReinforcingSteel::setTrialStrain
returns -1 for any trial strain below -0.95, so Truss::update and
AnalysisModel::updateDomain() return a negative value for such a strain.

With alpha = 0 (ArcLength, ArcLength1) or psi_u = 1, psi_f = 0, u_ref = 1
(HSConstraint) and a single free DOF, the predictor of newStep() is
dlambda = s / |dUhat| and dU = dUhat * dlambda = -s exactly, independent of
the tangent. An arc length s = 2.0 gives a trial strain of -2.0, which the
material refuses.

Before the fix newStep() and update() discard the return value of
updateDomain(); with algorithm Linear the step is reported as converged and
analyze(1) returns 0. With the fix newStep() returns -1 and analyze(1)
returns -2 (StaticAnalysis: integrator failed in newStep).

Control: s = 1e-4 stays elastic (yield strain 60/29000 = 2.07e-3) and
analyze(1) returns 0. For ArcLength and ArcLength1 the corrector keeps
|dU| = s, so the displacement of node 2 is -s. The HSConstraint corrector
does not hold that constraint exactly, so only a compressive displacement
inside the elastic range is checked for it.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


INTEGRATORS = {
    "ArcLength": lambda s: ("ArcLength", s, 0.0),
    "ArcLength1": lambda s: ("ArcLength1", s, 0.0),
    "HSConstraint": lambda s: ("HSConstraint", s, 1.0, 0.0, 1.0),
}


def _build(integrator, s):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("ReinforcingSteel", 1,
                         60.0, 90.0, 29000.0, 1000.0, 0.01, 0.1)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator(*INTEGRATORS[integrator](s))
    ops.algorithm("Linear")
    ops.analysis("Static")


@pytest.mark.parametrize("integrator", sorted(INTEGRATORS))
def test_refused_trial_strain_fails_analyze(integrator):
    try:
        _build(integrator, 2.0)
        assert ops.analyze(1) < 0
    finally:
        ops.wipe()


@pytest.mark.parametrize("integrator", sorted(INTEGRATORS))
def test_elastic_control_still_solves(integrator):
    try:
        s = 1.0e-4
        _build(integrator, s)
        assert ops.analyze(1) == 0
        u = ops.nodeDisp(2, 1)
        if integrator == "HSConstraint":
            assert -60.0 / 29000.0 < u < 0.0
        else:
            assert u == pytest.approx(-s, rel=1.0e-6)
    finally:
        ops.wipe()
