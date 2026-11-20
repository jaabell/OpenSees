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
//     DOI 10.1002/9780470694626. Chapter 14 (Box 14.3, Sections 14.5, 14.7).
//
// Description: abstract base class for two-dimensional (plane) finite-strain
// NDMaterials, the plane counterpart of FiniteStrainNDMaterial. It adds a
// 2x2 deformation-gradient entry point setTrialF(F):
//
//   * The element calls setTrialF(F) with the TOTAL in-plane deformation
//     gradient F = I + dU/dX (2x2). The material lifts F to 3x3 internally
//     (F33 = 1 for plane strain; F33 solved from sigma_33 = 0 for plane
//     stress). setTrialStrain() is disabled and returns -1.
//
//   * getStress() returns the in-plane Cauchy stress as a 3-vector in the
//     order {11, 22, 12}.
//
//   * getSpatialTangentTensor2D(c) returns the full in-plane spatial
//     constitutive modulus c_ijkl (i,j,k,l in {0,1}). The element forms
//     a_ijkl = c_ijkl - sigma_il delta_jk exactly as in 3D. For plane stress
//     c is statically condensed (delta sigma_33 = 0), so the same formula
//     yields the consistent plane-stress tangent. getTangent() (3x3) is the
//     lossy minor-symmetric projection, for output only.
//
//   * getJ() returns the in-plane Jacobian det(F) of the 2x2 gradient.
//
//   * getType() reports "PlaneStrain" or "PlaneStress"; getOrder() returns 3.

#ifndef FiniteStrainND2DMaterial_h
#define FiniteStrainND2DMaterial_h

#include <NDMaterial.h>
#include <OPS_Globals.h>

class Matrix;
class Vector;

class FiniteStrainND2DMaterial : public NDMaterial
{
 public:
  FiniteStrainND2DMaterial(int tag, int classTag) : NDMaterial(tag, classTag) {}
  FiniteStrainND2DMaterial() : NDMaterial() {}
  virtual ~FiniteStrainND2DMaterial() {}

  // Set the trial state from the total 2x2 deformation gradient F.
  // Returns 0 on success, < 0 on failure (det F <= 0, plane-stress
  // non-convergence).
  virtual int setTrialF(const Matrix &F) = 0;

  // In-plane Cauchy stress (3-vector, {11,22,12}).
  virtual const Vector &getStress(void) = 0;

  // Minor-symmetric 3x3 projection of the in-plane spatial modulus; lossy in
  // the (k,l) pair, for output only.
  virtual const Matrix &getTangent(void) = 0;

  // Full in-plane spatial constitutive modulus c_ijkl. Returns 0 on success;
  // the default -1 means "not supported".
  virtual int getSpatialTangentTensor2D(double c[2][2][2][2]) { return -1; }

  // In-plane Jacobian det(F) of the 2x2 gradient (> 0).
  virtual double getJ(void) const = 0;

  // Small-strain entry points are disabled (see FiniteStrainNDMaterial).
  virtual int setTrialStrain(const Vector & /*v*/) override {
    opserr << "FiniteStrainND2DMaterial: setTrialStrain() is invalid for a "
              "finite-strain material; the element must call setTrialF(F).\n";
    return -1;
  }
  virtual int setTrialStrain(const Vector & /*v*/, const Vector & /*r*/) override {
    opserr << "FiniteStrainND2DMaterial: setTrialStrain() is invalid for a "
              "finite-strain material; the element must call setTrialF(F).\n";
    return -1;
  }
};

#endif
