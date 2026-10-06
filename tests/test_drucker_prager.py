import json
import math
import os
import subprocess
import sys

import pytest

try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops

# The DruckerPrager nDMaterial (UW, Petek, Mackenzie-Helnwein and Arduino) has
# two yield surfaces: the cone f1 = ||eta|| + rho*I1 - sqrt(2/3)*K(alpha1) and
# the tension cutoff f2 = I1 - T(alpha2), with T(0) = sqrt(2/3)*sigma_y/rho.
#
# A single stdBrick on the unit cube is driven to a known homogeneous strain by
# prescribing u_i = eps_ij x_j on all 24 DOF (Penalty constraints). All eight
# Gauss points then carry the same state. A second stage starts from the state
# committed by the first.
#
# Before the return map was corrected, some of these paths never left the
# material's surface-selection loop, so every case runs in a subprocess with a
# time limit.

K = 1.0e4
G = 4.0e3
SY = 0.2
PHI = 20.0
RHO = 2.0*math.sin(math.radians(PHI))/(math.sqrt(3.0)*(3.0 - math.sin(math.radians(PHI))))
T_CUT = math.sqrt(2.0/3.0)*SY/RHO

# Hardening deck: with isotropic hardening the corner of the two surfaces is
# not the apex of the cone, so the consistent tangent there is not zero.
H_HARD = 6.0e4
THETA = 1.0

HEX = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
       (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
PENALTY = 1.0e12


def mode(eps):
   # nodal displacements of the homogeneous strain eps
   # (exx, eyy, ezz, gxy, gyz, gzx), engineering shear
   exx, eyy, ezz, gxy, gyz, gzx = eps
   u = []
   for (x, y, z) in HEX:
      u += [exx*x + 0.5*gxy*y + 0.5*gzx*z,
            0.5*gxy*x + eyy*y + 0.5*gyz*z,
            0.5*gzx*x + 0.5*gyz*y + ezz*z]
   return u


def apply_strain(eps, tag):
   ops.timeSeries('Linear', tag)
   ops.pattern('Plain', tag, tag)
   u = mode(eps)
   for i in range(8):
      for d in range(3):
         ops.sp(i + 1, d + 1, u[3*i + d])
   ops.constraints('Penalty', PENALTY, PENALTY)
   ops.numberer('Plain')
   ops.system('FullGeneral')
   ops.test('NormDispIncr', 1.0e-12, 50, 0)
   ops.algorithm('Newton')
   ops.integrator('LoadControl', 1.0)
   ops.analysis('Static', '-noWarnings')
   assert ops.analyze(1) == 0


def drive(stages, hard=0.0, theta=0.0):
   ops.wipe()
   ops.model('basic', '-ndm', 3, '-ndf', 3)
   for i, (x, y, z) in enumerate(HEX):
      ops.node(i + 1, x, y, z)
   # K G sigma_y rho rho_bar Kinf Ko delta1 delta2 H theta density
   ops.nDMaterial('DruckerPrager', 1, K, G, SY, RHO, 0.0,
                  0.0, 0.0, 0.0, 0.0, hard, theta, 0.0)
   ops.element('stdBrick', 1, *range(1, 9), 1)
   for k, eps in enumerate(stages, start=1):
      if k > 1:
         ops.remove('loadPattern', k - 1)
         ops.setTime(0.0)
         ops.wipeAnalysis()
      apply_strain(eps, k)


def stress(gp=1):
   return list(ops.eleResponse(1, 'material', gp, 'stress'))


def tangent():
   """Material tangent of the last step, read from the element stiffness.

   The stdBrick stiffness uses the tangent the material computed in the last
   Newton iteration of the step, i.e. the consistent tangent at the converged
   strain. For homogeneous strain modes u_j, u_k on the unit cube,
   u_j^T K_e u_k = C_jk. printA returns the matrix column by column, and the
   penalty constraints add PENALTY to its diagonal.
   """
   A = ops.printA('-ret')
   n = 24
   Ke = [[A[j*n + i] for j in range(n)] for i in range(n)]
   for i in range(n):
      Ke[i][i] -= PENALTY
   modes = []
   for k in range(6):
      e = [0.0]*6
      e[k] = 1.0
      modes.append(mode(e))
   C = [[0.0]*6 for _ in range(6)]
   for k in range(6):
      Ku = [sum(Ke[r][c]*modes[k][c] for c in range(n)) for r in range(n)]
      for j in range(6):
         C[j][k] = sum(modes[j][r]*Ku[r] for r in range(n))
   return C


# ---------------------------------------------------------------------------
# Tension cutoff: after a plastic step the mean stress is at the cutoff, not
# above it.
# ---------------------------------------------------------------------------

CUTOFF_CASES = {
   # hydrostatic tension, I1_trial = 9.0 > T: both surfaces active
   'hydrostatic': ([(1.0e-4, 1.0e-4, 1.0e-4, 0.0, 0.0, 0.0)], 0.0, 0.0),
   # tension with deviatoric strain: both surfaces active
   'corner': ([(3.0e-4, 1.0e-4, 1.0e-4, 0.0, 0.0, 0.0)], 0.0, 0.0),
   # hardened by shear, then shear removed and tension applied: the trial
   # state is inside the cone and above the cutoff (cutoff surface only)
   'cutoff_only': ([(0.0, 0.0, 0.0, 8.0e-4, 0.0, 0.0),
                    (5.0e-5, 5.0e-5, 5.0e-5, 0.0, 0.0, 0.0)], H_HARD, THETA),
}


def run_isolated(func, case, timeout=20):
   code = ('import sys, json; sys.path.insert(0, %r); '
           'import test_drucker_prager as t; '
           'print("RESULT" + json.dumps(t.%s(%r)))'
           % (os.path.dirname(os.path.abspath(__file__)), func, case))
   try:
      p = subprocess.run([sys.executable, '-c', code], capture_output=True,
                         text=True, timeout=timeout)
   except subprocess.TimeoutExpired:
      pytest.fail(f'{func}({case!r}) did not finish in {timeout} s')
   for line in p.stdout.splitlines():
      if line.startswith('RESULT'):
         return json.loads(line[6:])
   pytest.fail(f'{func}({case!r}) failed (exit {p.returncode}):\n'
               f'{p.stdout[-2000:]}\n{p.stderr[-2000:]}')


def cutoff_I1(case):
   stages, hard, theta = CUTOFF_CASES[case]
   drive(stages, hard, theta)
   return [sum(stress(gp)[:3]) for gp in range(1, 9)]


@pytest.mark.parametrize('case', list(CUTOFF_CASES))
def test_tension_cutoff_enforced(case):
   for gp, I1 in enumerate(run_isolated('cutoff_I1', case), start=1):
      assert I1 <= T_CUT*(1.0 + 1.0e-8), (
         f'GP {gp}: I1 = {I1} above the tension cutoff T = {T_CUT}')
      assert I1 == pytest.approx(T_CUT, rel=1.0e-8)


# ---------------------------------------------------------------------------
# Consistent tangent: central difference of the return map from the same
# committed state, compared column by column with the element-assembled
# tangent at the returned state.
# ---------------------------------------------------------------------------

FD_H = 1.0e-7
CASES = {
   # cone only: virgin material sheared past yield
   'cone': [(0.0, 0.0, 0.0, 8.0e-4, 0.0, 0.0)],
   # corner of cone and cutoff, away from the apex
   'corner': [(0.0, 0.0, 0.0, 8.0e-4, 0.0, 0.0),
              (1.0e-4, 1.0e-4, 1.0e-4, 8.0e-4, 0.0, 0.0)],
}


def tangent_fd(case):
   # The last stage is perturbed; the earlier stages are repeated exactly, so
   # every run starts the last step from the same committed state.
   stages = CASES[case]
   drive(stages, H_HARD, THETA)
   C = tangent()
   FD = [[0.0]*6 for _ in range(6)]
   for k in range(6):
      sig = []
      for sgn in (1.0, -1.0):
         last = list(stages[-1])
         last[k] += sgn*FD_H
         drive(stages[:-1] + [tuple(last)], H_HARD, THETA)
         sig.append(stress())
      for i in range(6):
         FD[i][k] = (sig[0][i] - sig[1][i])/(2.0*FD_H)
   return C, FD


@pytest.mark.parametrize('case', list(CASES))
def test_consistent_tangent_central_difference(case):
   C, FD = run_isolated('tangent_fd', case)
   for k in range(6):
      an = [C[i][k] for i in range(6)]
      fd = [FD[i][k] for i in range(6)]
      scale = max(max(abs(v) for v in an), max(abs(v) for v in fd))
      for i in range(6):
         assert abs(fd[i] - an[i]) <= 1.0e-3*scale, (
            f'{case}: C[{i}][{k}] = {an[i]}, central difference {fd[i]}')


# ---------------------------------------------------------------------------
# DruckerPragerPlaneStrain initial tangent: it is the elastic matrix whatever
# the loading history. A single plane-strain quad is pushed past yield, then
# one Newton step on the initial stiffness solves a small probe load, so the
# probe displacement is K_initial^-1 * probe. It must equal the probe
# displacement of a fresh, never loaded quad.
# ---------------------------------------------------------------------------

PS_K, PS_G = 27777.78, 9259.26
PS_RHO, PS_RHO_BAR = 0.398, 0.1
PS_SY, PS_H = 5.0, 1000.0


def plane_strain_probe(nsteps):
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 2)
   ops.nDMaterial('DruckerPrager', 1, PS_K, PS_G, PS_SY, PS_RHO, PS_RHO_BAR,
                  0.0, 0.0, 0.0, 0.0, PS_H, 1.0, 0.0)
   ops.node(1, 0.0, 0.0)
   ops.node(2, 1.0, 0.0)
   ops.node(3, 1.0, 1.0)
   ops.node(4, 0.0, 1.0)
   ops.element('quad', 1, 1, 2, 3, 4, 1.0, 'PlaneStrain', 1)
   ops.fix(1, 1, 1)
   ops.fix(2, 1, 1)
   ops.constraints('Plain')
   ops.numberer('Plain')
   ops.system('FullGeneral')
   ops.integrator('LoadControl', 1.0)
   if nsteps > 0:
      # shear and compression on the top edge, past first yield
      ops.timeSeries('Linear', 1)
      ops.pattern('Plain', 1, 1)
      ops.load(3, 1.5, -2.0)
      ops.load(4, 1.5, -2.0)
      ops.test('NormDispIncr', 1.0e-9, 60, 0)
      ops.algorithm('Newton')
      ops.integrator('LoadControl', 1.0/nsteps)
      ops.analysis('Static', '-noWarnings')
      for _ in range(nsteps):
         assert ops.analyze(1) == 0
      ops.loadConst('-time', 0.0)
   strain_before = list(ops.eleResponse(1, 'material', 1, 'strain'))

   ops.timeSeries('Linear', 2)
   ops.pattern('Plain', 2, 2)
   ops.load(3, 0.02, -0.02)
   ops.load(4, 0.02, -0.02)
   pre = [ops.nodeDisp(n, d) for n in (3, 4) for d in (1, 2)]
   # one solve on the initial stiffness, accepted without iterating
   ops.algorithm('Newton', '-initial')
   ops.test('NormDispIncr', 1.0e30, 1, 0)
   ops.integrator('LoadControl', 1.0)
   if nsteps == 0:
      ops.analysis('Static', '-noWarnings')
   assert ops.analyze(1) == 0
   post = [ops.nodeDisp(n, d) for n in (3, 4) for d in (1, 2)]
   return [b - a for a, b in zip(pre, post)], strain_before


def test_plane_strain_initial_tangent_is_elastic():
   du_fresh, _ = run_isolated('plane_strain_probe', 0)
   du_loaded, strain = run_isolated('plane_strain_probe', 30)
   assert max(abs(e) for e in strain) > 1.0e-4
   scale = max(abs(v) for v in du_fresh)
   assert scale > 1.0e-8
   for a, b in zip(du_fresh, du_loaded):
      assert abs(a - b) <= 1.0e-6*scale, (
         f'initial-stiffness probe: fresh {du_fresh}, after yield {du_loaded}')
