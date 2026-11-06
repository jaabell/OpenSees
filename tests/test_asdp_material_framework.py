"""ASDPlasticMaterial3D: integrators, Voigt convention and component consistency.

One stdBrick under a prescribed homogeneous strain (all nodes driven), so the
material state is the same at every Gauss point and is checked directly, plus
one brick with a fixed base for the global Newton iteration.

  radial_return_in_shear      a von Mises step with axial + shear strain returns
                              radially (exact Backward Euler answer) with every
                              integrator: the shear part of the plastic correction
                              uses the right Voigt convention
  full_backward_euler_options Full_Backward_Euler with both derivative modes gives
                              the same answer as Backward_Euler where the cutting
                              plane is exact (von Mises, linear hardening)
  mc_extension_corner         Backward_Euler converges with step refinement to the
                              RK45 answer on a path from the triaxial extension
                              meridian (the Mohr-Coulomb corner, now rounded)
  mc_tension_apex             a trial stress beyond the tensile apex (hydrostatic
                              extension, and extension plus shear) returns to
                              p = c cot(phi), q = 0 with Backward_Euler and
                              Full_Backward_Euler, for psi = 0 (no volumetric flow)
                              and psi > 0
  dp_hardening_consistency    Backward_Euler and RK45 agree for Drucker-Prager with
                              kinematic and cohesion hardening (consistent f and
                              hardening terms)
  algorithmic_tangent         with tangent_type Algorithmic the non-homogeneous brick
                              needs far fewer global iterations than Secant
  continuum_tangent_per_point the Continuum tangent converges in the non-homogeneous
                              brick (each Gauss point has its own tangent)
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E, NU = 50e6, 0.3


def vm_material(opts, Hiso=0.0, Hkin=0.0):
    ops.nDMaterial('ASDPlasticMaterial3D', 1, 'VonMises_YF', 'VonMises_PF', 'LinearIsotropic3D_EL',
                   'BackStress(TensorLinearHardeningFunction):YieldStress(ScalarLinearHardeningFunction):',
                   'Begin_Internal_Variables', 'BackStress', 0., 0., 0., 0., 0., 0., 'YieldStress', 20000.,
                   'End_Internal_Variables',
                   'Begin_Model_Parameters', 'InitialP0', 0., 'MassDensity', 0., 'PoissonsRatio', NU,
                   'ScalarLinearHardeningParameter', Hiso, 'TensorLinearHardeningParameter', Hkin,
                   'YoungsModulus', E, 'End_Model_Parameters',
                   'Begin_Integration_Options', 'f_absolute_tol', 1e-8, 'stress_absolute_tol', 1e-6,
                   'n_max_iterations', 60, 'return_to_yield_surface', 'Disabled', *opts,
                   'End_Integration_Options')


def mc_material(opts, psi=0.0):
    ops.nDMaterial('ASDPlasticMaterial3D', 1, 'MohrCoulomb_YF', 'MohrCoulomb_PF', 'LinearIsotropic3D_EL',
                   'BackStress(NullHardeningTensorFunction):',
                   'Begin_Model_Parameters', 'MC_phi', 30., 'MC_c', 10000., 'MC_ds', 1e-6, 'MC_psi', psi,
                   'YoungsModulus', E, 'PoissonsRatio', NU, 'MassDensity', 0., 'InitialP0', -1e-3,
                   'End_Model_Parameters',
                   'Begin_Integration_Options', 'f_absolute_tol', 1e-6, 'stress_absolute_tol', 1e-4,
                   'n_max_iterations', 60, 'return_to_yield_surface', 'Disabled', *opts,
                   'End_Integration_Options')


def dp_material(opts):
    import math
    sp, cp = math.sin(math.radians(30.)), math.cos(math.radians(30.))
    eta = 6 * sp / (math.sqrt(3) * (3 - sp))
    xi = 6 * 10000. * cp / (math.sqrt(3) * (3 - sp))
    ops.nDMaterial('ASDPlasticMaterial3D', 1, 'DruckerPrager_YF', 'DruckerPrager_PF', 'LinearIsotropic3D_EL',
                   'BackStress(TensorLinearHardeningFunction):DP_cohesion(ScalarLinearHardeningFunction):',
                   'Begin_Model_Parameters', 'DP_xi_c', xi, 'DP_eta', eta, 'DP_etabar', eta / 2,
                   'TensorLinearHardeningParameter', 3e6, 'ScalarLinearHardeningParameter', 1e6,
                   'YoungsModulus', E, 'PoissonsRatio', NU, 'MassDensity', 0., 'InitialP0', -1e-3,
                   'End_Model_Parameters',
                   'Begin_Integration_Options', 'f_absolute_tol', 1e-6, 'stress_absolute_tol', 1e-4,
                   'n_max_iterations', 60, 'return_to_yield_surface', 'Disabled', *opts,
                   'End_Integration_Options')


def strain_path(make_material, legs):
    """Drive one brick through piecewise-linear homogeneous strain legs.

    legs: list of (total strain [exx, eyy, ezz, gxy, gyz, gxz], steps). Returns the final stress.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    xyz = [(x, y, z) for z in (0, 1) for y in (0, 1) for x in (0, 1)]
    for i, c in enumerate(xyz):
        ops.node(i + 1, *c)
    # one free dof so the system of equations is never empty
    ops.node(9, 5., 5., 5.)
    ops.fix(9, 0, 1, 1)
    ops.uniaxialMaterial('Elastic', 99, 1.0)
    ops.element('zeroLength', 2, 9, 1, '-mat', 99, '-dir', 1)
    make_material()
    ops.element('stdBrick', 1, 1, 2, 4, 3, 5, 6, 8, 7, 1)

    def disp(c, e):
        x, y, z = c
        exx, eyy, ezz, gxy, gyz, gxz = e
        return (exx * x + gxy * y + gxz * z, eyy * y + gyz * z, ezz * z)

    times = [0.0] + [float(k + 1) for k in range(len(legs))] + [len(legs) + 1.0]
    tag = 1
    for i, c in enumerate(xyz):
        values = [(0., 0., 0.)] + [disp(c, e) for e, _ in legs]
        values.append(values[-1])
        for d in range(3):
            ops.timeSeries('Path', tag, '-time', *times, '-values', *[v[d] for v in values])
            ops.pattern('Plain', tag, tag)
            ops.sp(i + 1, d + 1, 1.0)
            tag += 1
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-10, 50)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    for _, n in legs:
        ops.integrator('LoadControl', 1.0 / n)
        for _ in range(n):
            assert ops.analyze(1) == 0
    return ops.eleResponse(1, 'material', 1, 'stress')


def rel_diff(a, b):
    sc = max(abs(x) for x in b)
    return max(abs(x - y) for x, y in zip(a, b)) / sc


BE = ('integration_method', 'Backward_Euler', 'tangent_type', 'Secant')
RK = ('integration_method', 'Runge_Kutta_45_Error_Control', 'tangent_type', 'Elastic')
FBE_N = ('integration_method', 'Full_Backward_Euler', 'derivatives', 'Numerical', 'tangent_type', 'Algorithmic')
FBE_CF = ('integration_method', 'Full_Backward_Euler', 'derivatives', 'Closed_Form', 'tangent_type', 'Algorithmic')


@pytest.mark.parametrize('opts', [BE, FBE_N, FBE_CF, RK,
                                  ('integration_method', 'Forward_Euler', 'tangent_type', 'Elastic')])
def test_radial_return_in_shear(opts):
    # eps_zz = 0.002, gamma_xz = 0.004 in one step, perfectly plastic von Mises:
    # the exact Backward Euler answer is the radial return of the trial deviator
    s = strain_path(lambda: vm_material(opts), [([0, 0, 0.002, 0, 0, 0.004], 1)])
    assert s == pytest.approx([80000., 80000., 90000., 0., 0., 10000.], rel=1e-8, abs=1e-6)


@pytest.mark.parametrize('opts', [FBE_N, FBE_CF])
def test_full_backward_euler_options(opts):
    legs = [([-5e-4] * 3 + [3e-3, 0, 0], 5), ([-5e-4, -5e-4, -1.5e-3, 3e-3, 3e-3, 0], 10)]
    ref = strain_path(lambda: vm_material(BE, Hiso=2e6, Hkin=5e6), legs)
    s = strain_path(lambda: vm_material(opts, Hiso=2e6, Hkin=5e6), legs)
    assert rel_diff(s, ref) < 1e-8


def test_mc_extension_corner():
    # confined, then axial extension (onto the extension meridian), then shear
    legs = lambda n: [([-5e-4] * 3 + [0, 0, 0], 10), ([-5e-4, -5e-4, 1e-3, 0, 0, 0], 10),
                      ([-5e-4, -5e-4, 1e-3, 3e-3, 0, 0], n)]
    ref = strain_path(lambda: mc_material(RK), legs(200))
    e = [rel_diff(strain_path(lambda: mc_material(BE), legs(n)), ref) for n in (10, 40, 160)]
    assert e[2] < e[0] and e[2] < 5e-3


@pytest.mark.parametrize('opts', [BE, FBE_N])
@pytest.mark.parametrize('psi', [0.0, 10.0])
@pytest.mark.parametrize('legs', [[([3e-4] * 3 + [0, 0, 0], 10), ([4e-4] * 3 + [0, 0, 0], 10)],
                                  [([0, 0, 0, 1e-3, 0, 0], 10), ([6e-4] * 3 + [1e-3, 0, 0], 10)]])
def test_mc_tension_apex(opts, psi, legs):
    apex = 10000. / 0.5773502691896257   # c cot(30 deg)
    s = strain_path(lambda: mc_material(opts, psi=psi), legs)
    assert s == pytest.approx([apex] * 3 + [0.] * 3, rel=1e-9, abs=1e-6)


def test_dp_hardening_consistency():
    legs = [([-5e-4] * 3 + [3e-3, 0, 0], 10), ([-5e-4, -5e-4, -1.5e-3, 3e-3, 3e-3, 0], 400)]
    a = strain_path(lambda: dp_material(BE), legs)
    b = strain_path(lambda: dp_material(RK), legs)
    assert rel_diff(a, b) < 1e-3


def brick_iterations(make_material):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    i = 1
    for z in (0, 1):
        for y in (0, 1):
            for x in (0, 1):
                ops.node(i, float(x), float(y), float(z))
                i += 1
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    make_material()
    ops.element('stdBrick', 1, 1, 2, 4, 3, 5, 6, 8, 7, 1)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.5, 0.0, -1.0)
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-10, 100)
    ops.algorithm('Newton')
    ops.integrator('DisplacementControl', 8, 3, -2e-4)
    ops.analysis('Static')
    total = 0
    for _ in range(20):
        assert ops.analyze(1) == 0
        total += ops.testIter()
    return total, ops.nodeDisp(8, 1)


def test_algorithmic_tangent():
    n_secant, u_secant = brick_iterations(lambda: vm_material(BE, Hiso=2e6, Hkin=5e6))
    n_alg, u_alg = brick_iterations(lambda: vm_material(FBE_N, Hiso=2e6, Hkin=5e6))
    assert u_alg == pytest.approx(u_secant, rel=1e-6)
    assert n_alg < n_secant / 3


def test_continuum_tangent_per_point():
    # the existing Continuum tangent converges in a non-homogeneous model only if
    # each Gauss point keeps its own tangent
    n_secant, u_secant = brick_iterations(lambda: vm_material(BE, Hiso=2e6, Hkin=5e6))
    cont = ('integration_method', 'Backward_Euler', 'tangent_type', 'Continuum')
    n_cont, u_cont = brick_iterations(lambda: vm_material(cont, Hiso=2e6, Hkin=5e6))
    assert u_cont == pytest.approx(u_secant, rel=1e-6)
    assert n_cont < n_secant
