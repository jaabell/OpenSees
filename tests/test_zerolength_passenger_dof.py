"""zeroLength between 3D nodes with different numbers of DOFs (passenger DOFs).

Node 1 is fixed; node 2 is coincident, free along the spring direction, joined by a zeroLength
spring (K = 1e4) and loaded along it. On a node with ndf > 3 the fourth DOF is imposed by `sp`,
so that the element is seen not to read it (the answer is that of a (3,3) pair) and not to write
it (the reaction on it is exactly zero).

Checks:

- A no-tension (ENT) spring across an ndf-3 / ndf-4 pair carries compression: K u balances the
  load (1e-12 relative), the passenger DOF keeps its imposed value with zero reaction, and the
  element `force` response has 3 + 4 slots with the translational slots +-P and zero elsewhere.
- Under tension the ENT spring opens (zero element force) while a soft parallel spring carries the
  load (1e-9).
- The pairs (3,3), (4,4), (4,3) and (3,6) give the same displacement (1e-12) and an element-sized
  `force` response with zeros on every non-translational slot.
- A rotational direction on a passenger pair is refused (element disabled); a (2,3) pair still
  gets the stock "differing dof at ends" refusal.

The (4,4), (4,3), (3,4) and (3,6) cases fail on the base: the element is refused there and
left without its transformation matrix, and the analysis stops with an access violation.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


K_SPRING = 1.0e4
P_LOAD = 250.0
PASSENGER_VAL = 7.0     # imposed on DOF 4 of an ndf-4 node: never read, never written


def _spring_pair(ndf_a, ndf_b, direction=3, load=+P_LOAD, ent=True):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', ndf_a)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, *([1] * ndf_a))
    ops.model('basic', '-ndm', 3, '-ndf', ndf_b)
    ops.node(2, 0.0, 0.0, 0.0)
    fx = [1] * ndf_b
    fx[direction - 1] = 0
    if ndf_b > 3:
        fx[3] = 0          # DOF 4 is imposed by sp below, not fixed
    ops.fix(2, *fx)
    if ent:
        ops.uniaxialMaterial('ENT', 1, K_SPRING)
    else:
        ops.uniaxialMaterial('Elastic', 1, K_SPRING)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', direction)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ld = [0.0] * ndf_b
    ld[direction - 1] = load
    ops.load(2, *ld)
    if ndf_b > 3:
        ops.sp(2, 4, PASSENGER_VAL)
    ops.constraints('Transformation')   # Plain would make the DOF-4 sp homogeneous
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.25)
    ops.analysis('Static')
    rc = 0
    for _ in range(4):
        rc = ops.analyze(1)
        if rc != 0:
            break
    u = ops.nodeDisp(2, direction)
    p4 = r4 = None
    if ndf_b > 3:
        p4 = ops.nodeDisp(2, 4)
        ops.reactions()
        r4 = ops.nodeReaction(2, 4)
    force = list(ops.eleResponse(1, 'force'))
    return rc, u, p4, r4, force


def test_ent_spring_3_4_carries_compression():
    rc, u, p4, r4, force = _spring_pair(3, 4, load=-P_LOAD)
    assert rc == 0
    assert abs(u + P_LOAD / K_SPRING) <= 1.0e-12 * P_LOAD / K_SPRING, u
    assert p4 == PASSENGER_VAL, p4
    assert r4 == 0.0, r4
    assert len(force) == 7, force
    assert abs(force[2] - P_LOAD) <= 1.0e-9 * P_LOAD, force
    assert abs(force[5] + P_LOAD) <= 1.0e-9 * P_LOAD, force
    assert force[6] == 0.0, force


def test_ent_spring_3_4_opens_under_tension():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.model('basic', '-ndm', 3, '-ndf', 4)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 0, 0)          # DOF 4 imposed by sp below
    # a soft parallel spring keeps the tension DOF non-singular, so the opened
    # ENT spring is measured rather than inferred from a failed solve
    ops.uniaxialMaterial('ENT', 1, K_SPRING)
    ops.uniaxialMaterial('Elastic', 2, 1.0e-3 * K_SPRING)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 3)
    ops.element('zeroLength', 2, 1, 2, '-mat', 2, '-dir', 3)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 0.0, 0.0, +P_LOAD, 0.0)
    ops.sp(2, 4, PASSENGER_VAL)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 30, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.25)
    ops.analysis('Static')
    for _ in range(4):
        assert ops.analyze(1) == 0
    u = ops.nodeDisp(2, 3)
    assert abs(u - P_LOAD / (1.0e-3 * K_SPRING)) <= 1.0e-9 * u, u
    ent = list(ops.eleResponse(1, 'force'))
    assert max(abs(v) for v in ent) == 0.0, ent
    assert ops.nodeDisp(2, 4) == PASSENGER_VAL
    ops.reactions()
    assert ops.nodeReaction(2, 4) == 0.0


@pytest.mark.parametrize('pair', [(3, 3), (4, 4), (4, 3), (3, 6)])
def test_same_answer_on_every_pair(pair):
    rc, u, p4, r4, force = _spring_pair(*pair, load=-P_LOAD)
    assert rc == 0, pair
    assert abs(u + P_LOAD / K_SPRING) <= 1.0e-12 * P_LOAD / K_SPRING, (pair, u)
    if p4 is not None:
        assert p4 == PASSENGER_VAL and r4 == 0.0, (pair, p4, r4)
    na, nb = pair
    assert len(force) == na + nb, (pair, force)
    assert abs(force[2] - P_LOAD) <= 1.0e-9 * P_LOAD, (pair, force)
    assert abs(force[na + 2] + P_LOAD) <= 1.0e-9 * P_LOAD, (pair, force)
    extras = [force[i] for i in range(na + nb) if i not in (0, 1, 2, na, na + 1, na + 2)]
    assert all(v == 0.0 for v in extras), (pair, force)


def test_rotational_dir_is_refused_on_a_passenger_pair(capfd):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.model('basic', '-ndm', 3, '-ndf', 4)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 1, 0)
    ops.uniaxialMaterial('Elastic', 1, K_SPRING)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 4)
    out = capfd.readouterr()
    assert 'element disabled' in out.err + out.out


def test_ndf_2_pair_is_still_refused(capfd):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 2)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(2, 1, 1, 0)
    ops.uniaxialMaterial('Elastic', 1, K_SPRING)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 3)
    assert 'differing dof at ends' in ''.join(capfd.readouterr())
