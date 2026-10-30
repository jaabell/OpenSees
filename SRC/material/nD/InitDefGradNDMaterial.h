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
//   E.H. Lee (1969). Elastic-plastic deformation at finite strains. Journal
//     of Applied Mechanics 36(1), 1-6. DOI 10.1115/1.3564580.
//   E.A. de Souza Neto, D. Peric, D.R.J. Owen (2008). Computational Methods
//     for Plasticity: Theory and Applications. Wiley, Chichester.
//     DOI 10.1002/9780470694626. Chapter 14.
//
// Description: InitDefGradNDMaterial is the multiplicative, finite-strain
// counterpart of InitStrainNDMaterial / InitStressNDMaterial. It makes a
// finite-strain continuum element added to a deformed model (staged
// construction: a new member, a concrete lift, a fill layer) start
// stress-free in the current configuration.
//
// Displacements are measured from the original mesh coordinates X, so an
// element added when the nodes have moved by u0 computes F from X and
// starts with the whole accumulated deformation as stress. Subtracting u0
// additively is exact at small strain but not objective at finite strain:
// with F_naive = F - F0 + I, a later rigid rotation F = Q F0 gives
// Q F0 - F0 + I, which differs from a rotation for Q != I and produces
// stress.
//
// The wrapper instead stores the deformation gradient at activation,
//
//     F0 = F at the first setTrialF,
//
// and passes the inner finite-strain material the relative gradient
//
//     F_rel = F F0^-1,
//
// the multiplicative split F = F_rel F0 (Lee 1969). A later rotation
// F = Q F0 gives F_rel = Q, so C_rel = F_rel^T F_rel = I and the stress is
// zero. An element active from the start has F0 = I and F_rel = F.
//
// For a spatial inner material (setTrialF in, Cauchy stress out, as
// LogStrainNDMaterial) the wrapper needs no further transformation:
//   1. The step increment does not depend on F0:
//        F_rel F_rel,n^-1 = (F F0^-1)(F_n F0^-1)^-1 = F F_n^-1.
//      F0 changes only the stress-free reference state.
//   2. The inner returns sigma = tau_rel / J_rel, which is the Cauchy stress
//      in the current configuration; no pull-back with F0 is needed. getJ()
//      returns J_rel = det F / det F0; an element integrating on the current
//      configuration uses its own det F.
// The wrapper therefore captures F0 once, forwards F_rel and delegates
// stress, tangent and state to the inner material.
//
// The inner must be a FiniteStrainNDMaterial (for example LogStrain).
// getStrain() returns the inner strain relative to activation. The
// captured F0 is available through the "F0" material response.
//
// Command: nDMaterial InitDefGrad $tag $innerTag <-noInitF>
//            <-F0 F11 F12 F13 F21 F22 F23 F31 F32 F33>
//   -noInitF  pass-through (F0 = I, never captured)
//   -F0       explicit F0 (row-major); skips the capture at the first call
// Aliases: InitDefGradNDMaterial, StagedDefGrad, StagedDefGradNDMaterial.
#ifndef InitDefGradNDMaterial_h
#define InitDefGradNDMaterial_h

#include <FiniteStrainNDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class InitDefGradNDMaterial : public FiniteStrainNDMaterial
{
 public:
  // active=false  -> pass-through (F0 = I, never captured); identical to the
  //                  inner material (-noInitF).
  // F0given!=0    -> explicit birth gradient (9 row-major components); skips the
  //                  auto-capture. Otherwise F0 is captured on the FIRST setTrialF.
  InitDefGradNDMaterial(int tag, NDMaterial &inner, bool active = true,
                        const double *F0given = 0);
  InitDefGradNDMaterial();
  ~InitDefGradNDMaterial();

  const char *getClassType(void) const { return "InitDefGradNDMaterial"; }

  // finite-strain entry point (FiniteStrainNDMaterial contract)
  int setTrialF(const Matrix &F);                     // total F in; forwards F F0^-1
  const Vector &getStress(void);                      // Cauchy stress (delegated)
  const Matrix &getTangent(void);                     // 6x6 projection of c (delegated)
  int getSpatialTangentTensor(double c[3][3][3][3]);  // full c_ijkl (delegated)
  const Matrix &getInitialTangent(void);
  double getJ(void) const;                            // J_rel = det(F F0^-1)

  const Vector &getStrain(void);                      // relative Hencky (delegated)
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
  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &matInfo);

 private:
  FiniteStrainNDMaterial *theMaterial;  // inner finite-strain material (deep copy)

  bool   active;        // false: pass-through (F0 = I)
  bool   captured;      // F0 has been snapshot (true once birth is seen / given)
  bool   f0Explicit;    // F0 supplied at construction (revertToStart keeps it;
                        // auto-capture mode re-arms instead)
  double F0[9];         // birth deformation gradient (row-major), captured at birth
  double invF0[9];      // F0^-1 (rebuilt from F0 on construction and recvSelf)

  Matrix Frel;          // scratch 3x3 passed to the inner setTrialF
  Vector birthF0;       // setResponse("F0") channel (9, row-major)

  // 3x3 row-major helpers
  static void   setIdentity9(double M[9]);
  static double det9(const double M[9]);
  static int    inv9(const double M[9], double Mi[9]);  // returns -1 if singular
  static void   mul9(const double A[9], const double B[9], double C[9]); // C = A B
};

#endif
