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
// Description: StagedStrainNDMaterial is a small-strain wrapper that makes a
// continuum element added to a deformed model (staged construction) start
// stress-free in the current configuration, for any material dimension. It
// is the small-strain counterpart of InitDefGrad.
//
// InitStrainNDMaterial adds a fixed prestrain known at construction. For
// staged activation it is limited in three ways: it is 3D only (getOrder()
// is 6 and getCopy() accepts only "ThreeDimensional"; a 2D inner exits the
// program), the strain is fixed at construction, so one material tag
// cannot serve elements activated at different deformed states, and it
// adds a prestrain rather than subtracting the activation strain.
// StagedStrain is dimension-general, captures the activation strain per
// material instance and subtracts it. InitStrain is unchanged.
//
// At activation the wrapper stores
//
//     eps0 = eps at the first setTrialStrain (the committed strain B u0)
//
// and passes the inner material the relative strain
//
//     eps_rel = eps - eps0.
//
// At activation eps_rel = 0, so a fresh inner returns zero stress and has
// no plastic history. The tangent is passed through unchanged
// (d sigma / d eps = d sigma / d eps_rel for constant eps0). The additive
// offset is exact at small strain; at finite strain it is not objective,
// and InitDefGrad applies instead.
//
// Dimension: eps0 is sized to the order of the inner copy (3 for
// PlaneStrain / PlaneStress, 4 for AxiSymmetric, 6 for ThreeDimensional),
// and getCopy(type) asks the inner for the requested view, so one tag
// serves 2D and 3D elements.
//
// Capture: an element creates one wrapper copy per integration point, and
// the first setTrialStrain of an element added to the domain is evaluated
// at the committed deformed state, so each copy captures its own eps0.
// An element active from the start captures eps0 = 0.
//
// Nesting with InitStrain: element -> StagedStrain( InitStrain(mat, eps_pre) )
// starts the element at sigma(eps_pre) in the deformed configuration. With
// the opposite nesting the prestrain is captured into eps0 and cancels, so
// StagedStrain must be the outer wrapper.
//
// Command: nDMaterial StagedStrain $tag $innerTag <-noInit> <-eps0 e1 e2 ...>
//   -noInit   pass-through (eps0 = 0, never captured)
//   -eps0     explicit activation strain (skips the capture); the number of
//             components must equal the order of the inner copy; the sign
//             convention is subtractive, opposite to InitStrain
// Alias: StagedStrainNDMaterial.
#ifndef StagedStrainNDMaterial_h
#define StagedStrainNDMaterial_h

#include <NDMaterial.h>
#include <Vector.h>

class StagedStrainNDMaterial : public NDMaterial
{
 public:
  // Takes ownership of innerOwned (a typed copy, e.g. from
  // inner.getCopy("ThreeDimensional")) without copying it again. A null
  // innerOwned leaves the object invalid (isValid() false) instead of exiting.
  //   active=false  -> pass-through (eps0 = 0; -noInit).
  //   eps0Given!=0  -> explicit activation strain (size must equal the inner
  //                    order); otherwise eps0 is captured at the first
  //                    setTrialStrain.
  StagedStrainNDMaterial(int tag, NDMaterial *innerOwned, bool active = true,
                         const Vector *eps0Given = 0);
  StagedStrainNDMaterial();

  bool isValid(void) const { return theMaterial != 0; }   // construction succeeded?
  ~StagedStrainNDMaterial();

  const char *getClassType(void) const { return "StagedStrainNDMaterial"; }

  // strain entry points (any order)
  int setTrialStrain(const Vector &v);
  int setTrialStrain(const Vector &v, const Vector &rate);
  int setTrialStrainIncr(const Vector &v);
  int setTrialStrainIncr(const Vector &v, const Vector &rate);

  const Vector &getStress(void);
  const Matrix &getTangent(void);
  const Matrix &getInitialTangent(void);
  const Vector &getStrain(void);      // relative strain eps - eps0
  double getRho(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial *getCopy(void);
  NDMaterial *getCopy(const char *type);   // adapts the inner to the requested view
  const char *getType(void) const;         // delegates to the inner
  int getOrder(void) const;                // delegates to the inner

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &matInfo);

  // sensitivity (DDM), as in InitStrainNDMaterial
  const Vector &getStressSensitivity(int gradIndex, bool conditional);
  int commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads);

 private:
  NDMaterial *theMaterial;   // inner small-strain material (deep copy)

  bool   active;             // false: pass-through (eps0 = 0)
  bool   captured;           // eps0 has been captured
  bool   eps0Explicit;       // eps0 given at construction (kept by revertToStart)
  Vector eps0;               // activation strain (sized to the inner order)
  Vector totalStrain;        // last total strain seen (for the totalStrain response)

  void sizeBuffers(int n);   // size eps0 and totalStrain to an order
  // capture eps0 at the first call, forward eps - eps0 (rate optional)
  int  forwardTrial(const Vector &v, const Vector *rate);
  // copy path: adopt a typed inner copy; keep the captured eps0 when the
  // order matches, otherwise the copy captures its own
  NDMaterial *adoptCopy(NDMaterial *innerCopy);
};

#endif
