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
//     DOI 10.1002/9780470694626. Chapter 14 (Box 14.3, Section 14.5).
//   J.C. Simo (1992). Algorithms for static and dynamic multiplicative
//     plasticity that preserve the classical return mapping schemes of the
//     infinitesimal theory. Computer Methods in Applied Mechanics and
//     Engineering 99(1), 61-112. DOI 10.1016/0045-7825(92)90123-2.
//   C. Miehe (1998). Comparison of two algorithms for the computation of
//     fourth-order isotropic tensor functions. Computers & Structures 66(1),
//     37-43. DOI 10.1016/S0045-7949(97)00073-4.
//
// Description: abstract base class for three-dimensional finite-strain
// NDMaterials driven by the deformation gradient F instead of a small-strain
// vector. It adds one entry point to the NDMaterial interface, setTrialF(F),
// and fixes the stress and tangent measures exchanged with a finite-strain
// element:
//
//   * The element calls setTrialF(F) with the TOTAL deformation gradient F
//     (3x3) at the integration point. A material that needs the increment
//     F_inc = F F_n^-1 computes it from its own committed F_n. The element
//     never calls setTrialStrain() on a finite-strain material; both
//     setTrialStrain() overloads return -1 with a message.
//
//   * getStress() returns the Cauchy stress sigma as a 6-vector in the 3D
//     Voigt order {11, 22, 33, 12, 23, 31}. The element assembles the
//     spatial internal force, the integral of B^T sigma over the current
//     volume.
//
//   * The spatial constitutive modulus is c_ijkl = (1/2J) [D : L : B]_ijkl,
//     the material part of the spatial tangent of Eq. (14.99) of de Souza
//     Neto et al. (2008). c is not minor-symmetric in (k,l), so a 6x6 Voigt
//     matrix cannot hold it. getTangent() returns the minor-symmetric 6x6
//     projection, which is suitable for output only. An element building a
//     consistent tangent calls getSpatialTangentTensor(c), forms
//         a_ijkl = c_ijkl - sigma_il delta_jk
//     (the geometric term is element-owned) and assembles
//         K_(a,i)(b,k) = integral of dN_a/dx_j a_ijkl dN_b/dx_l dv.
//     The full non-symmetric c and the full geometric term must both be kept
//     for the assembled K to be symmetric.
//
//   * getJ() returns J = det F.
//
//   * getType() reports "ThreeDimensional" and getOrder() returns 6.
//
// State ownership, commit and revert of any wrapped material, and
// sendSelf/recvSelf are the responsibility of the concrete subclass.

#ifndef FiniteStrainNDMaterial_h
#define FiniteStrainNDMaterial_h

#include <NDMaterial.h>
#include <OPS_Globals.h>

class Matrix;
class Vector;

class FiniteStrainNDMaterial : public NDMaterial
{
 public:
  FiniteStrainNDMaterial(int tag, int classTag) : NDMaterial(tag, classTag) {}
  FiniteStrainNDMaterial() : NDMaterial() {}
  virtual ~FiniteStrainNDMaterial() {}

  // Set the trial state from the total deformation gradient F (3x3).
  // Returns 0 on success, < 0 on failure (det F <= 0, local non-convergence),
  // so the element can reject the step.
  virtual int setTrialF(const Matrix &F) = 0;

  // Optional increment form: F_inc maps configuration n onto n+1
  // (F_{n+1} = F_inc F_n). Not supported unless a subclass overrides it.
  virtual int setTrialFincr(const Matrix &Fincr) {
    opserr << "FiniteStrainNDMaterial::setTrialFincr() not supported by this "
              "material; pass the total F via setTrialF().\n";
    return -1;
  }

  // Cauchy stress (6-vector, Voigt {11,22,33,12,23,31}).
  virtual const Vector &getStress(void) = 0;

  // Minor-symmetric 6x6 projection of the spatial constitutive modulus c.
  // Lossy in the (k,l) pair; for output only.
  virtual const Matrix &getTangent(void) = 0;

  // Full fourth-order spatial constitutive modulus c_ijkl (81 components).
  // Returns 0 on success; the default -1 means "not supported".
  virtual int getSpatialTangentTensor(double c[3][3][3][3]) { return -1; }

  // J = det F (> 0). Kirchhoff stress is tau = J sigma.
  virtual double getJ(void) const = 0;

  // The small-strain entry points are disabled: a finite-strain material is
  // driven by F, and a strain vector fed by a small-strain element would be
  // interpreted silently with the wrong kinematics.
  virtual int setTrialStrain(const Vector & /*v*/) override {
    opserr << "FiniteStrainNDMaterial: setTrialStrain() is invalid for a "
              "finite-strain material; the element must call setTrialF(F).\n";
    return -1;
  }
  virtual int setTrialStrain(const Vector & /*v*/, const Vector & /*r*/) override {
    opserr << "FiniteStrainNDMaterial: setTrialStrain() is invalid for a "
              "finite-strain material; the element must call setTrialF(F).\n";
    return -1;
  }
};

#endif
