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
//     DOI 10.1002/9780470694626. Box 14.3, Sections 14.5 and 14.7.
//   J.C. Simo (1992). Algorithms for static and dynamic multiplicative
//     plasticity that preserve the classical return mapping schemes of the
//     infinitesimal theory. Computer Methods in Applied Mechanics and
//     Engineering 99(1), 61-112. DOI 10.1016/0045-7825(92)90123-2.
//
// Description: LogStrain2D is the plane counterpart of LogStrainNDMaterial.
// It presents an order-3 PlaneStrain or PlaneStress material over an
// unchanged small-strain 3D NDMaterial by lifting the element's 2x2
// deformation gradient to 3x3 and driving an internal LogStrainNDMaterial,
// which owns the full 3D elastic left Cauchy-Green state:
//
//   Plane strain: F33 = 1, no out-of-plane shear. The 3D material is driven
//   directly and the in-plane block of sigma and c is returned. The full
//   3x3 Be is carried, so an out-of-plane plastic strain (eps_33 from
//   deviatoric J2 flow) is tracked.
//
//   Plane stress: the thickness stretch lambda = F33 is solved from
//   sigma_33(lambda) = 0 by a local Newton iteration (finite-difference
//   d sigma_33 / d lambda, started from the committed lambda), and the
//   out-of-plane direction is condensed out of the spatial tangent
//   (delta sigma_33 = 0), de Souza Neto et al. (2008), Section 14.7. For a
//   linear elastic inner the thickness strain is
//   eps_33 = -nu/(1-nu) (eps_11 + eps_22) in Hencky measure.
//
// FiniteStrainND2DMaterial contract: getStress() = in-plane Cauchy stress
// {11,22,12}; getSpatialTangentTensor2D(c) = full in-plane modulus (the
// element forms a = c - sigma_il delta_jk); getTangent() = lossy 3x3
// projection. The inner must provide a 3D (order 6) copy.
//
// Command: nDMaterial LogStrain2D $tag $innerTag <-planeStrain|-planeStress>
#ifndef LogStrain2D_h
#define LogStrain2D_h

#include <FiniteStrainND2DMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class LogStrainNDMaterial;

class LogStrain2D : public FiniteStrainND2DMaterial
{
 public:
  enum PlaneType { PlaneStrain = 1, PlaneStress = 2 };

  LogStrain2D(int tag, NDMaterial &inner3D, int planeType);
  LogStrain2D();
  ~LogStrain2D();

  const char *getClassType(void) const { return "LogStrain2D"; }

  // finite-strain entry point (FiniteStrainND2DMaterial contract)
  int setTrialF(const Matrix &F);              // total 2x2 deformation gradient
  const Vector &getStress(void);               // in-plane Cauchy stress (3)
  const Matrix &getTangent(void);              // 3x3 projection of c (lossy)
  int getSpatialTangentTensor2D(double c[2][2][2][2]);  // FULL in-plane c_ijkl
  const Matrix &getInitialTangent(void);
  double getJ(void) const { return Jarea; }

  const Vector &getStrain(void);               // in-plane Hencky strain (3-vector)
  double getRho(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial *getCopy(void);
  NDMaterial *getCopy(const char *type);
  const char *getType(void) const {
    return (planeType == PlaneStress) ? "PlaneStress" : "PlaneStrain";
  }
  int getOrder(void) const { return 3; }

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &matInfo);

 private:
  // internal 3D Hencky material (owns the 3D inner copy and Be_n)
  LogStrainNDMaterial *the3D;
  int planeType;                 // PlaneStrain (1) or PlaneStress (2)

  // committed state
  double lam_n;                  // committed thickness stretch F33 (plane stress)
  // trial state staged for commit
  double lamTrial;

  // cached trial outputs (in-plane, order 3)
  Vector sigma3;                 // in-plane Cauchy stress {11,22,12}
  Vector hencky3;                // in-plane Hencky strain {11,22,12}
  Matrix tangent3;               // 3x3 projection of c (lossy)
  double c2[2][2][2][2];         // FULL in-plane spatial modulus (element channel)
  double Jarea;                  // det of the 2x2 F

  // helpers
  int  driveLifted(const Matrix &F2, double F33);   // build F3D, call the3D->setTrialF
  int  buildOutputs(const Matrix &F2, bool condense); // fill sigma3/tangent3/c2; <0 = refuse
  static void voigt3FromTensor(const double c[2][2][2][2], Matrix &D);
};

#endif
