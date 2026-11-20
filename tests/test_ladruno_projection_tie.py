"""constraints LadrunoProjection with LadrunoTie mesh ties (explicit).

LadrunoTie writes one EQ_Constraint per slave node and translational DOF
(u_s = sum_i N_i u_m,i). Under CentralDifferenceLadruno the projection
handler enforces them on the acceleration, so the tie adds no stiffness.

Mesh: one SSPbrick (lumped mass) master block, x in [0, 1], and two SSPbrick
slave blocks, x in [1, 2], stacked in z. The slave nodes at z = 0.5 lie at
the middle of the master facet edges (weights 0.5 / 0.5).

Checks:
- the critical time step with the tie equals the one without it (1e-9);
- under an axial pull on the slave end the tie holds to 1e-10 for 200
  steps and the slave end moves;
- a massless slave node is refused at analyze() when the constraints
  command comes after LadrunoTie.

Skipped when the LadrunoTie command is not available.
"""

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


pytestmark = pytest.mark.skipif(not hasattr(ops, "LadrunoTie"),
                                reason="LadrunoTie command not available")

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"

MASTER_FACET = [2, 3, 7, 6]
SLAVE_NODES = [11, 14, 15, 18, 19, 22]
TIE_RELATIONS = {
    11: [(2, 1.0)], 14: [(3, 1.0)], 19: [(6, 1.0)], 22: [(7, 1.0)],
    15: [(2, 0.5), (6, 0.5)], 18: [(3, 0.5), (7, 0.5)],
}
_MASTER = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
           5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_SLAVE = {11: (1, 0, 0.0), 12: (2, 0, 0.0), 13: (2, 1, 0.0), 14: (1, 1, 0.0),
          15: (1, 0, 0.5), 16: (2, 0, 0.5), 17: (2, 1, 0.5), 18: (1, 1, 0.5),
          19: (1, 0, 1.0), 20: (2, 0, 1.0), 21: (2, 1, 1.0), 22: (1, 1, 1.0)}
_ALL = {**_MASTER, **_SLAVE}


def _analysis(handler):
    ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")


def _build(tie, with_load):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.0, 1.0)
    ops.element("SSPbrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("SSPbrick", 2, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    ops.element("SSPbrick", 3, 15, 16, 17, 18, 19, 20, 21, 22, 1)
    if tie:
        ops.LadrunoTie("-slaveNodes", len(SLAVE_NODES), *SLAVE_NODES,
                       "-masterFacets", 4, 1, *MASTER_FACET)
    for t, (x, y, z) in _ALL.items():
        if abs(x) < 1e-12:
            ops.fix(t, 1, 1, 1)
        elif tie and t in SLAVE_NODES:
            continue                  # tied DOFs must not be fixed
        else:
            ops.fix(t, 0, 1, 1)
    if with_load:
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for t in (20, 21):
            ops.load(t, 5.0, 0.0, 0.0)
    _analysis(PROJ if tie else "Plain")


def test_tie_does_not_change_critical_time_step():
    _build(tie=True, with_load=False)
    assert ops.analyze(1, 1e-4) == 0
    dt_tied = ops.criticalTimeStep()
    _build(tie=False, with_load=False)
    assert ops.analyze(1, 1e-4) == 0
    dt_free = ops.criticalTimeStep()
    assert dt_tied > 0 and dt_free > 0
    assert abs(dt_tied / dt_free - 1.0) < 1e-9


def test_tie_holds_every_step():
    _build(tie=True, with_load=True)
    assert ops.analyze(1, 1e-4) == 0
    dt = 0.4 * ops.criticalTimeStep()
    tie_err = 0.0
    for _ in range(200):
        assert ops.analyze(1, dt) == 0
        for s, rel in TIE_RELATIONS.items():
            us = ops.nodeDisp(s, 1)
            um = sum(w * ops.nodeDisp(m, 1) for (m, w) in rel)
            tie_err = max(tie_err, abs(us - um))
    assert tie_err < 1e-10, tie_err
    assert abs(ops.nodeDisp(20, 1)) > 1e-9


def test_massless_slave_refused_constraints_after_tie():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _MASTER.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.0, 1.0)
    ops.element("SSPbrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for t, (x, y, z) in _MASTER.items():
        if abs(x) < 1e-12:
            ops.fix(t, 1, 1, 1)
    ops.node(99, 1.0, 0.5, 0.5)       # on the facet, no mass, no element
    ops.LadrunoTie("-slaveNodes", 1, 99, "-masterFacets", 4, 1, *MASTER_FACET)
    try:
        _analysis(PROJ)
        ran = ops.analyze(1, 1e-4) == 0
    except Exception:
        ran = False
    assert not ran
