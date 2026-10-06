try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops
from math import pi, sin, log
import pytest

# Undamped linear SDOF in free vibration from an initial velocity, integrated
# with the generalized-alpha method (Chung and Hulbert, 1993). With alphaM and
# alphaF taken from the spectral radius rho_inf and the default gamma and beta,
# the scheme is second-order accurate for every rho_inf, so the end-time
# displacement error must fall by about 4 each time dt is halved. If the
# inertia term is evaluated at t+dt instead of t+alphaM*dt the scheme drops to
# first order whenever alphaM != 1.

m = 1.0
omega = 2*pi
k = m*omega**2
v0 = 1.0
tEnd = 1.1   # generic phase, so both amplitude and period errors show up
dts = [0.01, 0.005, 0.0025]

def exact_disp(t):
   return v0/omega*sin(omega*t)

def end_disp(alphaM, alphaF, dt):
   ops.wipe()
   ops.model('basic','-ndm',1,'-ndf',1)

   ops.node(1,0.0); ops.fix(1,1)
   ops.node(2,0.0); ops.mass(2,m)

   ops.uniaxialMaterial('Elastic',1,k)
   ops.element('zeroLength',1,1,2,'-mat',1,'-dir',1)

   # u0 = 0 and a0 = 0 are consistent for the undamped free vibration
   ops.setNodeVel(2,1,v0,'-commit')

   ops.constraints('Plain')
   ops.numberer('Plain')
   ops.system('FullGeneral')
   ops.test('NormDispIncr',1.0e-12,50)
   ops.algorithm('Newton')
   ops.integrator('GeneralizedAlpha',alphaM,alphaF)
   ops.analysis('Transient','-noWarnings')

   nSteps = round(tEnd/dt)
   assert ops.analyze(nSteps,dt) == 0
   return ops.nodeDisp(2,1)

@pytest.mark.parametrize('rho', [1.0, 0.8, 0.5, 0.0])
def test_second_order_convergence(rho):
   alphaM = (2.0 - rho)/(1.0 + rho)
   alphaF = 1.0/(1.0 + rho)

   err = [abs(end_disp(alphaM,alphaF,dt) - exact_disp(tEnd)) for dt in dts]

   for i in range(len(dts) - 1):
      assert err[i+1] < err[i]
      order = log(err[i]/err[i+1])/log(dts[i]/dts[i+1])
      assert order >= 1.8
