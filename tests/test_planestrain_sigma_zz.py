"""Out-of-plane stress sigma_zz of plane-strain elements: "stressesPlaneStrain".

A plane-strain material computes sigma_zz under eps_zz = 0, but the legacy
"stresses" response returns only [sxx, syy, sxy] per integration point. The
"stressesPlaneStrain" response of FourNodeQuad and Tri31 returns
[sxx, syy, sxy, szz] per integration point, with szz read from
NDMaterial::getStressZZ().

Checks:
  1. elastic (ElasticIsotropic, DruckerPrager below yield): szz equals
     nu*(sxx + syy) to 1e-9 on quad and Tri31,
  2. the first three components equal the legacy "stresses" response bitwise,
  3. J2 plastic (perfect plasticity, oedometric push past yield): von Mises
     built with the reported szz equals sigma_y to 1e-8, while the elastic
     reconstruction nu*(sxx + syy) misses it by more than 5 percent,
  4. the PlaneStrain wrapper around the 3D J2 material reports the same
     szz as the dedicated J2PlaneStrain material,
  5. a plane-stress material, which has no sigma_zz override, reports NaN.

On a build without the response the element returns no response and every
test fails.
"""
import math

import pytest

try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops

E, NU = 30000.0, 0.3
DELTA = 1.0e-4  # prescribed top settlement for the elastic runs

QUAD_NODES = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
TRI_NODES = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.0, 1.0)}


def _mises(sxx, syy, sxy, szz):
   m = (sxx + syy + szz) / 3.0
   dxx, dyy, dzz = sxx - m, syy - m, szz - m
   return math.sqrt(1.5 * (dxx * dxx + dyy * dyy + dzz * dzz + 2.0 * sxy * sxy))


def _oedometer(nodes):
   """eps_xx = 0 (all x fixed), bottom y fixed, top y prescribed later."""
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 2)
   for tag, (x, y) in nodes.items():
      ops.node(tag, x, y)
      ops.fix(tag, 1, 1 if y == 0.0 else 0)


def _drive(top_nodes, delta, nsteps=10):
   ops.timeSeries('Linear', 1)
   ops.pattern('Plain', 1, 1)
   for n in top_nodes:
      ops.sp(n, 2, -delta)
   # every free DOF is prescribed, so Penalty keeps a non-empty system
   ops.constraints('Penalty', 1.0e15, 1.0e15)
   ops.numberer('Plain')
   ops.system('FullGeneral')
   ops.test('NormDispIncr', 1.0e-12, 50, 0)
   ops.algorithm('Newton')
   ops.integrator('LoadControl', 1.0 / nsteps)
   ops.analysis('Static')
   assert ops.analyze(nsteps) == 0


def _check_elastic(ele_tag, ngp, nu, rel=1.0e-9):
   s3 = ops.eleResponse(ele_tag, 'stresses')
   s4 = ops.eleResponse(ele_tag, 'stressesPlaneStrain')
   assert len(s3) == 3 * ngp
   assert len(s4) == 4 * ngp
   for g in range(ngp):
      sxx, syy, sxy, szz = s4[4 * g:4 * g + 4]
      assert (sxx, syy, sxy) == tuple(s3[3 * g:3 * g + 3])
      ref = nu * (sxx + syy)
      assert abs(ref) > 0.0
      assert abs(szz - ref) <= rel * abs(ref)


def test_quad_elastic_szz():
   _oedometer(QUAD_NODES)
   ops.nDMaterial('ElasticIsotropic', 1, E, NU)
   ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', 1)
   _drive([3, 4], DELTA)
   _check_elastic(1, 4, NU)


def test_tri31_elastic_szz():
   _oedometer(TRI_NODES)
   ops.nDMaterial('ElasticIsotropic', 1, E, NU)
   ops.element('Tri31', 1, 1, 2, 3, 1.0, 'PlaneStrain', 1)
   _drive([3], DELTA)
   _check_elastic(1, 1, NU)


def test_stress_alias():
   _oedometer(QUAD_NODES)
   ops.nDMaterial('ElasticIsotropic', 1, E, NU)
   ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', 1)
   _drive([3, 4], DELTA)
   s4 = ops.eleResponse(1, 'stressesPlaneStrain')
   assert len(s4) == 16
   assert ops.eleResponse(1, 'stressPlaneStrain') == s4


def test_drucker_prager_elastic_szz():
   # oedometric compression stays inside the cone for these constants:
   # sqrt(2/3)*2G < 3*K*rho, so the path is elastic at any strain
   K, G = 27777.78, 9259.26
   _oedometer(QUAD_NODES)
   # K G sigma_y rho rho_bar Kinf Ko delta1 delta2 H theta density
   ops.nDMaterial('DruckerPrager', 1, K, G, 5.0, 0.398, 0.1,
                  0.0, 0.0, 0.0, 0.0, 1000.0, 1.0, 0.0)
   ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', 1)
   _drive([3, 4], 1.0e-3)
   _check_elastic(1, 4, (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G)))


def test_j2_plastic_szz_on_yield_surface():
   K, G, SIG0 = 25000.0, 11500.0, 30.0
   _oedometer(QUAD_NODES)
   # perfect plasticity: sig0 == sigInf, H = 0
   ops.nDMaterial('J2Plasticity', 1, K, G, SIG0, SIG0, 0.0, 0.0)
   ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', 1)
   _drive([3, 4], 0.02, nsteps=40)

   nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
   s4 = ops.eleResponse(1, 'stressesPlaneStrain')
   assert len(s4) == 16
   for g in range(4):
      sxx, syy, sxy, szz = s4[4 * g:4 * g + 4]
      assert math.isfinite(szz)
      assert _mises(sxx, syy, sxy, szz) == pytest.approx(SIG0, rel=1.0e-8)
      elastic_guess = _mises(sxx, syy, sxy, nu * (sxx + syy))
      assert abs(elastic_guess - SIG0) > 0.05 * SIG0


def test_planestrain_wrapper_matches_dedicated_j2():
   K, G, SIG0 = 25000.0, 11500.0, 30.0

   def run(wrapped):
      _oedometer(QUAD_NODES)
      ops.nDMaterial('J2Plasticity', 1, K, G, SIG0, SIG0, 0.0, 0.0)
      mat = 1
      if wrapped:
         ops.nDMaterial('PlaneStrain', 2, 1)
         mat = 2
      ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', mat)
      _drive([3, 4], 0.02, nsteps=40)
      return ops.eleResponse(1, 'stressesPlaneStrain')

   direct, wrapped = run(False), run(True)
   assert len(direct) == 16 and len(wrapped) == 16
   assert all(math.isfinite(v) for v in wrapped)
   assert wrapped == pytest.approx(direct, rel=1.0e-10)


def test_plane_stress_reports_nan():
   _oedometer(QUAD_NODES)
   ops.nDMaterial('ElasticIsotropic', 1, E, NU)
   ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStress', 1)
   _drive([3, 4], DELTA)
   s4 = ops.eleResponse(1, 'stressesPlaneStrain')
   assert len(s4) == 16
   for g in range(4):
      assert math.isnan(s4[4 * g + 3])
      assert all(math.isfinite(v) for v in s4[4 * g:4 * g + 3])
