/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell (Ladruño project)
//
// References:
//   E.A. de Souza Neto, D. Peric, D.R.J. Owen (2008). Computational Methods
//     for Plasticity: Theory and Applications. Wiley, Chichester.
//     DOI 10.1002/9780470694626. Box 14.3, Sections 14.4 and 14.5.
//   J.C. Simo (1992). Algorithms for static and dynamic multiplicative
//     plasticity that preserve the classical return mapping schemes of the
//     infinitesimal theory. Computer Methods in Applied Mechanics and
//     Engineering 99(1), 61-112. DOI 10.1016/0045-7825(92)90123-2.
//   A.L. Eterovic, K.J. Bathe (1990). A hyperelastic-based large strain
//     elasto-plastic constitutive formulation with combined isotropic-
//     kinematic hardening using the logarithmic stress and strain measures.
//     International Journal for Numerical Methods in Engineering 30(6),
//     1099-1114. DOI 10.1002/nme.1620300602.
//   C. Miehe, N. Apel, M. Lambrecht (2002). Anisotropic additive plasticity
//     in the logarithmic strain space: modular kinematic formulation and
//     implementation based on incremental minimization principles for
//     standard materials. Computer Methods in Applied Mechanics and
//     Engineering 191(47-48), 5383-5425. DOI 10.1016/S0045-7825(02)00438-3.
//   C. Miehe (1998). Comparison of two algorithms for the computation of
//     fourth-order isotropic tensor functions. Computers & Structures 66(1),
//     37-43. DOI 10.1016/S0045-7949(97)00073-4.
//
// Description: LogStrainNDMaterial lifts an unchanged small-strain 3D
// NDMaterial to finite strain in logarithmic (Hencky) strain space, by the
// spatial multiplicative scheme of Box 14.3 of de Souza Neto et al. (2008):
//
//   Be_trial = F_inc Be_n F_inc^T                       Eq. (14.94)
//   eps_trial = 1/2 ln(Be_trial)                        Eq. (14.30)
//   inner small-strain law at eps_trial -> Kirchhoff tau, tangent D
//   sigma = tau / J                                     Box 14.3 (iv)
//   c = (1/2J) [D : L : B]                              Eq. (14.99)
//
// The inner return map is reused unchanged: with the Hencky elastic strain
// and the exponential-map plastic update, the finite-strain update has the
// format of the infinitesimal one (Simo 1992; Eterovic and Bathe 1990).
// Miehe, Apel and Lambrecht (2002) describe the same modular idea in the
// Lagrangian logarithmic strain space.
//
// The class implements the FiniteStrainNDMaterial contract: setTrialF(F)
// with the total F, getStress() = Cauchy stress, getSpatialTangentTensor()
// = full c_ijkl, getTangent() = its lossy 6x6 projection. The geometric
// term -sigma_il delta_jk of Eq. (14.99) is left to the element.
//
// State: the material owns, per integration point, the committed elastic
// left Cauchy-Green tensor Be_n, the committed F_n and the committed strain
// last fed to the inner material, and wraps the inner commit / revert /
// sendSelf / recvSelf.
//
// Plastic inner materials: a stock small-strain plastic NDMaterial treats
// the strain it receives as total strain and subtracts its own committed
// plastic strain. Since Be_n already carries the elastic history, the inner
// is fed eps_feed = eps_feed_n + (eps_trial - eps_e_n), so the inner's own
// subtraction returns exactly eps_trial (the invariant
// eps_p_inner = eps_feed_n - eps_e_n is maintained on commit). The updated
// elastic strain is recovered as eps_e_{n+1} = C_e : tau with C_e the
// inverse of the inner initial tangent, and Be_{n+1} = exp(2 eps_e_{n+1}).
//
// Limitations:
//   * The recovery eps_e = C_e : tau is exact only for an inner whose
//     elastic law is linear with constant stiffness (elastic and
//     rate-independent elastoplastic materials such as ElasticIsotropic,
//     J2Plasticity, Drucker-Prager). A stiffness-degrading (damage) inner,
//     tau = (1-d) D_e : eps_e, makes the recovered elastic strain shrink by
//     (1-d) and the committed Be drift; such materials need their own
//     FiniteStrainNDMaterial subclass.
//   * A kinematic-hardening back stress stored by the inner in a fixed
//     frame does not co-rotate, so combined hardening is not objective under
//     large rotation (de Souza Neto et al. 2008, Section 14.11). Isotropic
//     hardening is exact.
//   * The inner law is evaluated with Kirchhoff stress; pressure-sensitive
//     yield at large volume change sees tau, not sigma.
//   * The inner must provide a 3D (order 6) copy.
//
// Command: nDMaterial LogStrain $tag $innerTag

#ifndef LogStrainNDMaterial_h
#define LogStrainNDMaterial_h

#include <FiniteStrainNDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class LogStrainNDMaterial : public FiniteStrainNDMaterial
{
 public:
  LogStrainNDMaterial(int tag, NDMaterial &inner);
  LogStrainNDMaterial();
  ~LogStrainNDMaterial();

  const char *getClassType(void) const { return "LogStrainNDMaterial"; }

  // finite-strain entry point (FiniteStrainNDMaterial contract)
  int setTrialF(const Matrix &F);          // total deformation gradient
  const Vector &getStress(void);           // Cauchy stress (6)
  const Matrix &getTangent(void);          // 6x6 projection of c (lossy)
  int getSpatialTangentTensor(double c[3][3][3][3]);  // full c_ijkl
  const Matrix &getInitialTangent(void);
  double getJ(void) const { return Jdet; }

  const Vector &getStrain(void);           // trial Hencky strain (6)
  double getRho(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial *getCopy(void);
  NDMaterial *getCopy(const char *type);
  const char *getType(void) const { return "ThreeDimensional"; }
  int getOrder(void) const { return 6; }

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  // "stress", "strain" and "tangent" return the finite-strain measures of
  // this material (Cauchy stress, Hencky strain, spatial c), not the inner
  // Kirchhoff stress, fed strain and small-strain D, which differ by det F.
  // Any other response is delegated to the inner material.
  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &matInfo);

 private:
  NDMaterial *theMaterial;     // inner small-strain 3D material (own copy)

  // committed state
  double Fn[9];                // committed deformation gradient (row-major)
  double Be_n[9];              // committed elastic left Cauchy-Green tensor
  double epsFeed_n[6];         // committed strain last fed to the inner
                               // (engineering-shear Voigt)

  // trial state staged for commit
  double Ftrial9[9];           // trial F
  double BeTrial9[9];          // trial Be (used by the spatial tangent)
  double Be_trialUpd[9];       // Be to commit = exp(2 eps_e_{n+1})
  double epsFeedTrial[6];      // trial fed strain

  // cached trial output
  Vector sigmaCauchy;          // Cauchy stress (6)
  Vector henckyStrain;         // trial Hencky strain (6, engineering shear)
  Matrix aTangent;             // 6x6 projection of c
  double Jdet;                 // det F

  void setIdentity(double M[9]);
};

#endif
