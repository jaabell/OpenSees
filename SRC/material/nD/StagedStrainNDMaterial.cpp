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
// StagedStrainNDMaterial: small-strain, dimension-general staged-activation
// wrapper. See StagedStrainNDMaterial.h for the formulation.

#include <StagedStrainNDMaterial.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <Information.h>
#include <Response.h>
#include <MaterialResponse.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>

// =========================================================================== //
//  Factory:  nDMaterial StagedStrain $tag $innerTag <-noInit> <-eps0 e1 e2 ...>//
// =========================================================================== //
void *OPS_StagedStrainNDMaterial(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING invalid args: nDMaterial StagedStrain $tag $innerTag "
              "<-noInit> <-eps0 e1 e2 ...>\n";
    return 0;
  }
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid ints: nDMaterial StagedStrain $tag $innerTag\n";
    return 0;
  }

  NDMaterial *inner = OPS_getNDMaterial(iData[1]);
  if (inner == 0) {
    opserr << "WARNING nDMaterial StagedStrain " << iData[0]
           << " : inner nDMaterial " << iData[1] << " not found\n";
    return 0;
  }
  bool   active = true;
  bool   haveEps0 = false;
  Vector eps0in;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-noInit") == 0) {
      active = false;
    } else if (strcmp(opt, "-eps0") == 0) {
      int n = OPS_GetNumRemainingInputArgs();   // remaining args = the eps0 components
      if (n < 1) {
        opserr << "WARNING nDMaterial StagedStrain " << iData[0]
               << " : -eps0 needs at least one component\n";
        return 0;
      }
      double *buf = new double[n];
      if (OPS_GetDoubleInput(&n, buf) != 0) {
        opserr << "WARNING nDMaterial StagedStrain " << iData[0]
               << " : invalid -eps0 components\n";
        delete[] buf;
        return 0;
      }
      eps0in.resize(n);
      for (int i = 0; i < n; i++) eps0in(i) = buf[i];
      delete[] buf;
      haveEps0 = true;
    } else {
      opserr << "WARNING nDMaterial StagedStrain " << iData[0]
             << " : unknown option '" << opt << "'\n";
      return 0;
    }
  }

  // Take one "ThreeDimensional" copy of the inner, after option parsing so an
  // option error cannot leak it. The no-argument getCopy() of the generic
  // material base classes exits the program, so it is not used here; 2D and
  // axisymmetric views are derived from this 3D copy by getCopy(type). An
  // invalid inner makes the command return 0.
  NDMaterial *inner3D = inner->getCopy("ThreeDimensional");
  if (inner3D == 0) {
    opserr << "WARNING nDMaterial StagedStrain " << iData[0]
           << " : inner nDMaterial " << iData[1]
           << " has no ThreeDimensional view (StagedStrain needs a 3D-capable inner)\n";
    return 0;
  }
  // adopt inner3D (no second copy); the ctor takes ownership.
  StagedStrainNDMaterial *m =
    new StagedStrainNDMaterial(iData[0], inner3D, active, haveEps0 ? &eps0in : 0);
  if (m == 0 || !m->isValid()) { if (m) delete m; return 0; }   // graceful
  return m;
}

// =========================================================================== //
//  Construction (eps0 and totalStrain are sized at the first setTrialStrain,  //
//  to the order the element uses)                                            //
// =========================================================================== //
StagedStrainNDMaterial::StagedStrainNDMaterial(int tag, NDMaterial *innerOwned,
                                               bool active_, const Vector *eps0Given)
  : NDMaterial(tag, ND_TAG_StagedStrainNDMaterial),
    theMaterial(innerOwned), active(active_), captured(false), eps0Explicit(false),
    eps0(), totalStrain()
{
  // innerOwned is adopted without a second copy. A null innerOwned leaves the
  // object invalid (isValid() false) so the caller can refuse it.
  if (theMaterial != 0 && eps0Given != 0) {   // explicit eps0: captured now
    eps0 = *eps0Given;
    totalStrain.resize(eps0.Size()); totalStrain.Zero();
    captured     = true;
    eps0Explicit = true;
  }
}

StagedStrainNDMaterial::StagedStrainNDMaterial()
  : NDMaterial(0, ND_TAG_StagedStrainNDMaterial),
    theMaterial(0), active(true), captured(false), eps0Explicit(false),
    eps0(), totalStrain()
{
}

StagedStrainNDMaterial::~StagedStrainNDMaterial()
{
  if (theMaterial) delete theMaterial;
}

// size eps0 and totalStrain to order n (keep an explicit eps0 of that size)
void StagedStrainNDMaterial::sizeBuffers(int n)
{
  if (eps0.Size() != n) {
    Vector keep = eps0;
    eps0.resize(n); eps0.Zero();
    if (eps0Explicit && keep.Size() == n) eps0 = keep;   // keep a matching explicit eps0
  }
  if (totalStrain.Size() != n) { totalStrain.resize(n); totalStrain.Zero(); }
}

// =========================================================================== //
//  setTrialStrain: capture eps0 at the first call, forward eps - eps0         //
// =========================================================================== //
// Shared by both setTrialStrain overloads. eps0 is captured at the first
// call and the inner receives eps_rel = eps - eps0. The first evaluation of
// an element added to the domain is made at the committed deformed state
// (the predictor has not moved the displacements yet), so eps0 is the
// activation strain. A first trial strain that is later reverted would be
// kept as eps0; revertToStart clears the capture.
int StagedStrainNDMaterial::forwardTrial(const Vector &v, const Vector *rate)
{
  const int n = v.Size();
  // size eps0 and totalStrain so the subtraction below stays in bounds
  if (totalStrain.Size() != n || eps0.Size() != n) sizeBuffers(n);
  totalStrain = v;

  if (!active)                     // pass-through (-noInit)
    return rate ? theMaterial->setTrialStrain(v, *rate)
                : theMaterial->setTrialStrain(v);

  if (!captured) { eps0 = v; captured = true; }   // capture the activation strain

  Vector rel(v);                   // eps
  rel.addVector(1.0, eps0, -1.0);  // eps_rel = eps - eps0
  return rate ? theMaterial->setTrialStrain(rel, *rate)
              : theMaterial->setTrialStrain(rel);
}

int StagedStrainNDMaterial::setTrialStrain(const Vector &v)
{
  return forwardTrial(v, 0);
}

int StagedStrainNDMaterial::setTrialStrain(const Vector &v, const Vector &rate)
{
  return forwardTrial(v, &rate);
}

int StagedStrainNDMaterial::setTrialStrainIncr(const Vector &v)
{
  // Rebuild the total strain from the strain last passed to the inner:
  //   active and captured: the inner holds eps - eps0, so add eps0
  //   -noInit: the inner holds eps
  // then add the increment. Forwarding only the increment would drop the
  // accumulated strain.
  Vector total(theMaterial->getStrain());      // eps_rel (active) or eps (-noInit)
  if (active && captured) total.addVector(1.0, eps0, 1.0);   // + eps0
  total.addVector(1.0, v, 1.0);                 // + increment
  return setTrialStrain(total);
}

int StagedStrainNDMaterial::setTrialStrainIncr(const Vector &v, const Vector &/*rate*/)
{
  return setTrialStrainIncr(v);
}

// =========================================================================== //
//  Query: delegated to the inner, which works in eps_rel = eps - eps0          //
// =========================================================================== //
const Vector &StagedStrainNDMaterial::getStress(void)         { return theMaterial->getStress(); }
const Matrix &StagedStrainNDMaterial::getTangent(void)        { return theMaterial->getTangent(); }
const Matrix &StagedStrainNDMaterial::getInitialTangent(void) { return theMaterial->getInitialTangent(); }
const Vector &StagedStrainNDMaterial::getStrain(void)         { return theMaterial->getStrain(); }  // eps_rel
double        StagedStrainNDMaterial::getRho(void)            { return theMaterial->getRho(); }

const char *StagedStrainNDMaterial::getType(void) const  { return theMaterial->getType(); }
int         StagedStrainNDMaterial::getOrder(void) const { return theMaterial->getOrder(); }

// =========================================================================== //
//  State: eps0 is fixed after capture; the inner owns its own history          //
// =========================================================================== //
int StagedStrainNDMaterial::commitState(void)        { return theMaterial->commitState(); }
int StagedStrainNDMaterial::revertToLastCommit(void) { return theMaterial->revertToLastCommit(); }

int StagedStrainNDMaterial::revertToStart(void)
{
  // Back to the constructed state: a captured eps0 is cleared and captured
  // again at the next setTrialStrain; an explicit eps0 is kept.
  if (!eps0Explicit) {
    captured = false;
    eps0.Zero();
  }
  return theMaterial->revertToStart();
}

// =========================================================================== //
//  Copies (getCopy(type) adapts the inner to the requested dimensional view)   //
// =========================================================================== //
// Shared by getCopy() and getCopy(type). Adopts the typed innerCopy as is
// (converting it to 3D would give a 2D element an order-6 material). The
// captured eps0 is carried when the order matches, so a copy of a captured
// wrapper keeps its reference; on an order mismatch the copy captures its
// own (with a warning if an explicit eps0 is dropped).
NDMaterial *StagedStrainNDMaterial::adoptCopy(NDMaterial *innerCopy)
{
  if (innerCopy == 0) return 0;
  const int order = innerCopy->getOrder();
  // the constructor takes ownership; the copy starts uncaptured
  StagedStrainNDMaterial *c =
    new StagedStrainNDMaterial(this->getTag(), innerCopy, active, 0);
  if (c == 0 || !c->isValid()) { if (c) delete c; return 0; }

  if (captured && eps0.Size() == order) {           // keep the captured eps0
    c->captured     = true;
    c->eps0Explicit = eps0Explicit;
    c->eps0         = eps0;
    c->sizeBuffers(order);
  } else if (eps0Explicit && eps0.Size() != order) {  // explicit eps0 does not fit
    opserr << "StagedStrainNDMaterial::getCopy - explicit -eps0 size " << eps0.Size()
           << " != requested order " << order
           << "; ignoring it (the copy captures eps0 at its first setTrialStrain)\n";
  }
  return c;   // otherwise the copy captures eps0 at its first setTrialStrain
}

NDMaterial *StagedStrainNDMaterial::getCopy(void)
{
  // theMaterial is a typed copy, so its no-argument getCopy() is valid and
  // keeps the type.
  return adoptCopy(theMaterial->getCopy());
}

NDMaterial *StagedStrainNDMaterial::getCopy(const char *type)
{
  // Ask the inner for the view the element requests (PlaneStrain,
  // PlaneStress, AxiSymmetric, ThreeDimensional) and adopt it. This assumes
  // inner.getCopy("ThreeDimensional").getCopy(type) equals
  // inner.getCopy(type), which holds for materials that build typed copies
  // from their parameters (ElasticIsotropic, J2Plasticity). The command checks
  // only the 3D view, so an inner without the requested view is refused here
  // (return 0).
  NDMaterial *innerCopy = theMaterial->getCopy(type);
  if (innerCopy == 0) {
    opserr << "StagedStrainNDMaterial::getCopy - inner material does not support "
              "type '" << type << "'\n";
    return 0;
  }
  return adoptCopy(innerCopy);
}

// =========================================================================== //
//  Parallel / database                                                        //
// =========================================================================== //
int StagedStrainNDMaterial::sendSelf(int cTag, Channel &theChannel)
{
  if (theMaterial == 0) return -1;
  int dbTag = this->getDbTag();
  const int n = eps0.Size();

  static ID dataID(4);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) { matDbTag = theChannel.getDbTag(); theMaterial->setDbTag(matDbTag); }
  dataID(2) = matDbTag;
  dataID(3) = n;                    // eps0 order (0 before the first capture)
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) return -1;

  Vector dataVec(3 + n);            // flags + eps0
  dataVec(0) = active       ? 1.0 : 0.0;
  dataVec(1) = captured     ? 1.0 : 0.0;
  dataVec(2) = eps0Explicit ? 1.0 : 0.0;
  for (int i = 0; i < n; i++) dataVec(3 + i) = eps0(i);
  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) return -2;

  if (theMaterial->sendSelf(cTag, theChannel) < 0) return -3;
  return 0;
}

int StagedStrainNDMaterial::recvSelf(int cTag, Channel &theChannel,
                                     FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID dataID(4);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) return -1;
  this->setTag(dataID(0));
  const int n = dataID(3);

  if (theMaterial == 0) {
    theMaterial = theBroker.getNewNDMaterial(dataID(1));
    if (theMaterial == 0) {
      opserr << "StagedStrainNDMaterial::recvSelf - cannot create inner material classTag "
             << dataID(1) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  Vector dataVec(3 + n);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) return -3;
  active       = (dataVec(0) != 0.0);
  captured     = (dataVec(1) != 0.0);
  eps0Explicit = (dataVec(2) != 0.0);
  if (n > 0) { sizeBuffers(n); for (int i = 0; i < n; i++) eps0(i) = dataVec(3 + i); }

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) return -4;
  // size the buffers to the inner order so "totalStrain" queried before the
  // next setTrialStrain returns a zero vector of the right size
  if (theMaterial->getOrder() > 0) sizeBuffers(theMaterial->getOrder());
  return 0;
}

// =========================================================================== //
//  Misc                                                                       //
// =========================================================================== //
void StagedStrainNDMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"StagedStrainNDMaterial\", ";
    s << "\"inner\": \"" << theMaterial->getTag() << "\", ";
    s << "\"active\": " << (active ? "true" : "false") << ", ";
    s << "\"captured\": " << (captured ? "true" : "false");
    s << "}";
  } else {
    s << "StagedStrainNDMaterial (small-strain staged-activation wrapper), tag: "
      << this->getTag() << endln;
    s << "\tinner material tag: " << theMaterial->getTag() << endln;
    s << "\tactive: " << (active ? "yes" : "no")
      << ", eps0 captured: " << (captured ? "yes" : "no") << endln;
  }
}

int StagedStrainNDMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

int StagedStrainNDMaterial::updateParameter(int /*parameterID*/, Information &/*info*/)
{
  return 0;   // eps0 is not a runtime parameter (-eps0 or capture)
}

// Responses: "eps0" (alias "birthStrain") = the activation strain;
// "totalStrain" = the strain from the original mesh; everything else
// (stress, strain = eps - eps0, tangent, inner-specific output) is delegated
// to the inner.
Response *StagedStrainNDMaterial::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (argc > 0) {
    if (strcmp(argv[0], "eps0") == 0 || strcmp(argv[0], "birthStrain") == 0)
      return new MaterialResponse(this, 201, eps0);
    if (strcmp(argv[0], "totalStrain") == 0)
      return new MaterialResponse(this, 202, totalStrain);
  }
  return theMaterial->setResponse(argv, argc, output);
}

int StagedStrainNDMaterial::getResponse(int responseID, Information &matInfo)
{
  switch (responseID) {
  case 201: return matInfo.setVector(eps0);
  case 202: return matInfo.setVector(totalStrain);
  default:  return theMaterial->getResponse(responseID, matInfo);
  }
}

// =========================================================================== //
//  Sensitivity (DDM), as in InitStrainNDMaterial                              //
// =========================================================================== //
const Vector &StagedStrainNDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  return theMaterial->getStressSensitivity(gradIndex, conditional);
}

int StagedStrainNDMaterial::commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(depsdh, gradIndex, numGrads);
}
