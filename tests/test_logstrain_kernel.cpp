// Standalone checks of SRC/material/nD/LogStrainKernel.h, the numerical core
// of the LogStrain nDMaterial. The kernel has no OpenSees dependency, so this
// test compiles against the header alone. An isotropic linear elastic law in
// Hencky strain (K, G) stands in for the inner material.
//
// Checks:
//   1. Principal stretch F = diag(l1, l2, l3): the Cauchy stress equals the
//      closed form sigma = [K tr(e) I + 2G dev(e)] / J with e_i = ln(l_i).
//   2. Pure rotation: zero stress.
//   3. Objectivity: sigma(Q F) = Q sigma(F) Q^T.
//   4. Spatial modulus: c_ijkl = (1/J) d(tau_ij)/dh at F_h = (I + h e_k e_l) F,
//      by central differences, for all 81 components, at distinct, two-fold
//      and three-fold repeated principal stretches and at a nearly repeated
//      pair (the analytic-limit branches of the log derivative).
//   5. The 6x6 output of assemble_material equals the Voigt projection of the
//      full tensor.
//   6. be_from_hencky_voigt inverts hencky_voigt.

#include <LogStrainKernel.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace logstrain_kernel;

static const double K = 1500.0;
static const double G = 700.0;
static int failures = 0;

static void check(bool ok, const char *what, double got, double ref)
{
  if (!ok) {
    std::printf("FAIL: %s (got %.15g, expected %.15g)\n", what, got, ref);
    failures++;
  }
}

// tau = K tr(e) I + 2G dev(e); engineering-shear Voigt in and out.
static void elastic(const double e6[6], double tau6[6], double D6[36])
{
  double tr = e6[0] + e6[1] + e6[2];
  for (int i = 0; i < 3; i++) tau6[i] = K*tr + 2.0*G*(e6[i] - tr/3.0);
  for (int i = 3; i < 6; i++) tau6[i] = G*e6[i];
  for (int k = 0; k < 36; k++) D6[k] = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      D6[6*i + j] = K + 2.0*G*((i == j ? 1.0 : 0.0) - 1.0/3.0);
  D6[6*3+3] = G; D6[6*4+4] = G; D6[6*5+5] = G;
}

// Kirchhoff stress tau (Voigt) of the elastic Hencky law at total F.
static void kirchhoff(const double F[9], double tau6[6], double Betr[9], double D6[36])
{
  const double I9[9] = {1,0,0, 0,1,0, 0,0,1};
  trial_Be(F, I9, I9, Betr);
  double e6[6];
  hencky_voigt(Betr, e6);
  elastic(e6, tau6, D6);
}

static void cauchy(const double F[9], double sig6[6], double c6[36])
{
  double tau6[6], Betr[9], D6[36];
  kirchhoff(F, tau6, Betr, D6);
  assemble_material(Betr, D6, tau6, mat3_det(F), sig6, c6);
}

static double maxabs(const double *a, int n)
{
  double m = 0.0;
  for (int i = 0; i < n; i++) if (std::fabs(a[i]) > m) m = std::fabs(a[i]);
  return m;
}

static void test_principal_stretch()
{
  const double l[3] = {1.25, 0.95, 0.90};
  double F[9] = {l[0],0,0, 0,l[1],0, 0,0,l[2]};
  double sig6[6], c6[36];
  cauchy(F, sig6, c6);
  double e[3] = {std::log(l[0]), std::log(l[1]), std::log(l[2])};
  double tr = e[0] + e[1] + e[2], J = l[0]*l[1]*l[2];
  for (int i = 0; i < 3; i++) {
    double ref = (K*tr + 2.0*G*(e[i] - tr/3.0)) / J;
    check(std::fabs(sig6[i] - ref) <= 1.0e-12*std::fabs(ref), "principal stretch: normal stress", sig6[i], ref);
  }
  for (int i = 3; i < 6; i++)
    check(std::fabs(sig6[i]) <= 1.0e-12*K, "principal stretch: shear stress", sig6[i], 0.0);
}

static void rotation(double angle, double Q[9])
{
  // rotation about the axis (1,2,3)/sqrt(14) (Rodrigues formula)
  double n[3] = {1.0, 2.0, 3.0};
  double nn = std::sqrt(14.0);
  for (int i = 0; i < 3; i++) n[i] /= nn;
  double c = std::cos(angle), s = std::sin(angle);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double eps = 0.0;
      if (i != j) {
        int k = 3 - i - j;
        double sgn = ((i == 0 && j == 1) || (i == 1 && j == 2) || (i == 2 && j == 0)) ? 1.0 : -1.0;
        eps = -sgn*n[k];
      }
      Q[3*i+j] = c*(i == j ? 1.0 : 0.0) + (1.0 - c)*n[i]*n[j] + s*eps;
    }
}

static void voigt_to_tensor(const double v[6], double T[9])
{
  T[0] = v[0]; T[4] = v[1]; T[8] = v[2];
  T[1] = T[3] = v[3]; T[5] = T[7] = v[4]; T[2] = T[6] = v[5];
}

static void test_rotation_and_objectivity()
{
  double Q[9];
  rotation(0.7, Q);
  double sig6[6], c6[36];
  cauchy(Q, sig6, c6);
  check(maxabs(sig6, 6) <= 1.0e-10*K, "pure rotation: zero stress", maxabs(sig6, 6), 0.0);

  double F[9] = {1.20, 0.10, -0.05,  0.02, 0.90, 0.08,  0.03, -0.04, 1.10};
  double QF[9];
  mat3_mul(Q, F, QF);
  double s1[6], s2[6];
  cauchy(F, s1, c6);
  cauchy(QF, s2, c6);
  double S1[9], S2[9], QS[9], QT[9], QSQT[9];
  voigt_to_tensor(s1, S1);
  voigt_to_tensor(s2, S2);
  mat3_transpose(Q, QT);
  mat3_mul(Q, S1, QS);
  mat3_mul(QS, QT, QSQT);
  double scale = maxabs(S1, 9);
  for (int i = 0; i < 9; i++)
    check(std::fabs(S2[i] - QSQT[i]) <= 1.0e-10*scale, "objectivity: sigma(QF) = Q sigma Q^T", S2[i], QSQT[i]);
}

static void test_tangent_fd(const char *name, const double F[9])
{
  double tau6[6], Betr[9], D6[36];
  kirchhoff(F, tau6, Betr, D6);
  double J = mat3_det(F);
  double c[3][3][3][3];
  spatial_tangent_full(Betr, D6, J, c);

  double cmax = 0.0;
  for (int i = 0; i < 81; i++) {
    double v = std::fabs((&c[0][0][0][0])[i]);
    if (v > cmax) cmax = v;
  }

  const double h = 1.0e-6;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++) {
      double tp[6], tm[6], Bt[9], Dt[36];
      double Fp[9], Fm[9], Gp[9], Gm[9];
      for (int a = 0; a < 9; a++) { Gp[a] = (a % 4 == 0) ? 1.0 : 0.0; Gm[a] = Gp[a]; }
      Gp[3*k+l] += h;
      Gm[3*k+l] -= h;
      mat3_mul(Gp, F, Fp);
      mat3_mul(Gm, F, Fm);
      kirchhoff(Fp, tp, Bt, Dt);
      kirchhoff(Fm, tm, Bt, Dt);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          int I = VI(i, j);
          double fd = (tp[I] - tm[I]) / (2.0*h) / J;
          if (std::fabs(c[i][j][k][l] - fd) > 1.0e-6*cmax) {
            char what[96];
            std::snprintf(what, sizeof(what), "%s: c[%d][%d][%d][%d] vs central difference", name, i, j, k, l);
            check(false, what, c[i][j][k][l], fd);
          }
        }
    }

  // the 6x6 projection returned by assemble_material
  double sig6[6], c6[36];
  assemble_material(Betr, D6, tau6, J, sig6, c6);
  for (int I = 0; I < 6; I++)
    for (int Jn = 0; Jn < 6; Jn++) {
      double ref = c[REP_I[I]][REP_J[I]][REP_I[Jn]][REP_J[Jn]];
      if (std::fabs(c6[6*I+Jn] - ref) > 1.0e-12*cmax)
        check(false, "6x6 projection of c", c6[6*I+Jn], ref);
    }
}

static void test_exp_log_roundtrip()
{
  double F[9] = {1.20, 0.10, -0.05,  0.02, 0.90, 0.08,  0.03, -0.04, 1.10};
  const double I9[9] = {1,0,0, 0,1,0, 0,0,1};
  double Be[9], e6[6], Be2[9];
  trial_Be(F, I9, I9, Be);
  hencky_voigt(Be, e6);
  be_from_hencky_voigt(e6, Be2);
  for (int i = 0; i < 9; i++)
    check(std::fabs(Be2[i] - Be[i]) <= 1.0e-12, "exp(2 eps) recovers Be", Be2[i], Be[i]);
}

int main()
{
  test_principal_stretch();
  test_rotation_and_objectivity();

  const double generic[9] = {1.20, 0.10, -0.05,  0.02, 0.90, 0.08,  0.03, -0.04, 1.10};
  const double distinct[9] = {1.25, 0.0, 0.0,  0.0, 0.95, 0.0,  0.0, 0.0, 0.90};
  const double twofold[9] = {1.2, 0.0, 0.0,  0.0, 1.2, 0.0,  0.0, 0.0, 0.8};
  const double threefold[9] = {1.1, 0.0, 0.0,  0.0, 1.1, 0.0,  0.0, 0.0, 1.1};
  const double nearpair[9] = {1.2, 0.0, 0.0,  0.0, 1.2 + 1.0e-11, 0.0,  0.0, 0.0, 0.8};
  test_tangent_fd("generic", generic);
  test_tangent_fd("distinct", distinct);
  test_tangent_fd("two-fold", twofold);
  test_tangent_fd("three-fold", threefold);
  test_tangent_fd("nearly two-fold", nearpair);

  test_exp_log_roundtrip();

  if (failures != 0) {
    std::printf("%d check(s) failed\n", failures);
    return EXIT_FAILURE;
  }
  std::printf("all LogStrainKernel checks passed\n");
  return EXIT_SUCCESS;
}
