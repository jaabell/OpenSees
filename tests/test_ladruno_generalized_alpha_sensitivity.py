"""DDM response sensitivities of the LadrunoGeneralizedAlpha integrator.

LadrunoGeneralizedAlpha is the generalized-alpha integrator (Chung and
Hulbert, 1993) with Direct Differentiation Method sensitivities (Conte,
Vijalapura and Meghella, 2003). The checks use a single-DOF truss under a sine
pulse with alphaM = 0.9 and alphaF = 0.8, so both weightings of the
sensitivity residual are active:

1. Without sensitivities the displacement history is identical, bit for bit,
   to that of GeneralizedAlpha with the same parameters.
2. d(u)/dE from the DDM matches a central finite difference of the response,
   undamped and with mass-proportional Rayleigh damping. Relative error below
   1e-5. The undamped case exercises the alphaM-weighted mass term, which is
   the derivative of the inertia evaluated at t + alphaM*dt by
   GeneralizedAlpha::update().
3. d(u)/dfy for a yielding Steel01 truss matches a central finite difference
   with Newton and with ModifiedNewton; the latter checks that the
   sensitivity solves use a tangent re-formed at the converged state.
4. At alphaM = alphaF = 1 the DDM equals the Newmark (gamma = 0.5,
   beta = 0.25) DDM.

On a build without LadrunoGeneralizedAlpha the integrator command is refused
and every test fails.
"""
try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
import pytest

# single-DOF truss: fixed node 1, truss, node 2 with mass, free in x only
E0, A, L, M = 2.0e3, 1.5, 1.0, 1.0   # k = EA/L = 3000, omega ~ 54.8 rad/s
DT, NSTEPS = 5.0e-4, 80
ALPHA_M, ALPHA_F = 0.9, 0.8
FY0, B0 = 20.0, 0.02                 # Steel01 case: yield force FY0*A = 30


def run(intArgs, E=E0, sens=False, alphaRay=0.0, steel=None, algo='Newton'):
   """Run the transient analysis.

   Returns the history of u2x and, with sens=True, the final DDM
   sensitivity d(u2x)/d(parameter 1). steel=(fy, b) uses a Steel01 truss
   with parameter fy instead of an elastic one with parameter E.
   """
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 2)
   ops.node(1, 0.0, 0.0)
   ops.node(2, L, 0.0)
   ops.fix(1, 1, 1)
   ops.fix(2, 0, 1)
   ops.mass(2, M, M)

   if steel is None:
      ops.uniaxialMaterial('Elastic', 1, E)
      factor = 100.0
   else:
      ops.uniaxialMaterial('Steel01', 1, steel[0], E0, steel[1])
      factor = 2.0*FY0*A   # twice the nominal yield force, independent of fy
   ops.element('Truss', 1, 1, 2, A, 1)

   ops.timeSeries('Trig', 1, 0.0, DT*NSTEPS, DT*NSTEPS, '-factor', factor)
   ops.pattern('Plain', 1, 1)
   ops.load(2, 1.0, 0.0)

   if sens:
      if steel is None:
         ops.parameter(1, 'element', 1, 'E')
      else:
         ops.parameter(1, 'element', 1, 'fy')

   if alphaRay > 0.0:
      ops.rayleigh(alphaRay, 0.0, 0.0, 0.0)   # C = alphaRay*M, independent of E

   ops.constraints('Plain')
   ops.numberer('Plain')
   ops.system('BandGeneral')
   ops.test('NormDispIncr', 1.0e-12, 100)
   ops.algorithm(algo)
   ops.integrator(*intArgs)
   ops.analysis('Transient')
   if sens:
      ops.sensitivityAlgorithm('-computeAtEachStep')

   hist = []
   for i in range(NSTEPS):
      assert ops.analyze(1, DT) == 0
      hist.append(ops.nodeDisp(2, 1))

   ddm = ops.sensNodeDisp(2, 1, 1) if sens else None
   ops.wipe()
   return hist, ddm


def test_without_sensitivity_identical_to_generalized_alpha():
   ref, _ = run(('GeneralizedAlpha', ALPHA_M, ALPHA_F))
   new, _ = run(('LadrunoGeneralizedAlpha', ALPHA_M, ALPHA_F))
   assert new == ref


@pytest.mark.parametrize('alphaRay', [0.0, 4.0], ids=['undamped', 'massDamped'])
def test_ddm_matches_finite_difference_E(alphaRay):
   args = ('LadrunoGeneralizedAlpha', ALPHA_M, ALPHA_F)
   _, ddm = run(args, sens=True, alphaRay=alphaRay)
   dE = 1.0e-6*E0
   up, _ = run(args, E=E0+dE, alphaRay=alphaRay)
   um, _ = run(args, E=E0-dE, alphaRay=alphaRay)
   fd = (up[-1] - um[-1])/(2.0*dE)
   assert abs(fd) > 0.0
   assert abs(ddm - fd) <= 1.0e-5*abs(fd)


@pytest.mark.parametrize('algo', ['Newton', 'ModifiedNewton'])
def test_ddm_matches_finite_difference_fy_steel01(algo):
   args = ('LadrunoGeneralizedAlpha', ALPHA_M, ALPHA_F)
   hist, ddm = run(args, sens=True, steel=(FY0, B0), algo=algo)
   assert max(abs(u) for u in hist) > FY0/E0*L   # the truss yields
   dfy = 1.0e-6*FY0
   up, _ = run(args, steel=(FY0+dfy, B0), algo=algo)
   um, _ = run(args, steel=(FY0-dfy, B0), algo=algo)
   fd = (up[-1] - um[-1])/(2.0*dfy)
   assert abs(fd) > 0.0
   assert abs(ddm - fd) <= 1.0e-6*abs(fd)


def test_alpha_one_equals_newmark_ddm():
   _, ddmGA = run(('LadrunoGeneralizedAlpha', 1.0, 1.0), sens=True)
   _, ddmNM = run(('Newmark', 0.5, 0.25), sens=True)
   assert abs(ddmNM) > 0.0
   assert abs(ddmGA - ddmNM) <= 1.0e-8*abs(ddmNM)
