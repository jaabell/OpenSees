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
// InitDefGradNDMaterial: multiplicative staged-activation wrapper for
// finite-strain materials. See InitDefGradNDMaterial.h for the formulation
// and references.

#include <InitDefGradNDMaterial.h>
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
//  3x3 row-major linear-algebra helpers                                       //
// =========================================================================== //
void InitDefGradNDMaterial::setIdentity9(double M[9]) {
  for (int i = 0; i < 9; i++) M[i] = 0.0;
  M[0] = M[4] = M[8] = 1.0;
}

double InitDefGradNDMaterial::det9(const double M[9]) {
  return M[0]*(M[4]*M[8] - M[5]*M[7])
       - M[1]*(M[3]*M[8] - M[5]*M[6])
       + M[2]*(M[3]*M[7] - M[4]*M[6]);
}

int InitDefGradNDMaterial::inv9(const double M[9], double Mi[9]) {
  double d = det9(M);
  if (d == 0.0) return -1;
  double id = 1.0 / d;
  Mi[0] =  (M[4]*M[8] - M[5]*M[7]) * id;
  Mi[1] = -(M[1]*M[8] - M[2]*M[7]) * id;
  Mi[2] =  (M[1]*M[5] - M[2]*M[4]) * id;
  Mi[3] = -(M[3]*M[8] - M[5]*M[6]) * id;
  Mi[4] =  (M[0]*M[8] - M[2]*M[6]) * id;
  Mi[5] = -(M[0]*M[5] - M[2]*M[3]) * id;
  Mi[6] =  (M[3]*M[7] - M[4]*M[6]) * id;
  Mi[7] = -(M[0]*M[7] - M[1]*M[6]) * id;
  Mi[8] =  (M[0]*M[4] - M[1]*M[3]) * id;
  return 0;
}

void InitDefGradNDMaterial::mul9(const double A[9], const double B[9], double C[9]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += A[3*i+k] * B[3*k+j];
      C[3*i+j] = s;
    }
}

// =========================================================================== //
//  Factory:  nDMaterial InitDefGrad $tag $innerTag <-noInitF> <-F0 f11..f33>  //
// =========================================================================== //
void *OPS_InitDefGradNDMaterial(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING invalid args: nDMaterial InitDefGrad $tag $innerTag "
              "<-noInitF> <-F0 f11 f12 f13 f21 f22 f23 f31 f32 f33>\n";
    return 0;
  }
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid ints: nDMaterial InitDefGrad $tag $innerTag\n";
    return 0;
  }

  // The inner must be a finite-strain (F-driven) material: the wrapper
  // re-references the deformation gradient, which only such a material uses.
  NDMaterial *inner = OPS_getNDMaterial(iData[1]);
  if (inner == 0) {
    opserr << "WARNING nDMaterial InitDefGrad " << iData[0]
           << " : inner nDMaterial " << iData[1] << " not found\n";
    return 0;
  }
  if (dynamic_cast<FiniteStrainNDMaterial *>(inner) == 0) {
    opserr << "WARNING nDMaterial InitDefGrad " << iData[0]
           << " : inner nDMaterial " << iData[1]
           << " must be a finite-strain (F-driven) material, e.g. LogStrain. "
              "InitDefGrad is the multiplicative staged-activation wrapper; for a "
              "small-strain element use InitStrain instead.\n";
    return 0;
  }

  // optional flags
  bool   active   = true;
  bool   haveF0   = false;
  double F0in[9];   // populated only when -F0 is given (then fully overwritten)

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-noInitF") == 0) {
      active = false;
    } else if (strcmp(opt, "-F0") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 9) {
        opserr << "WARNING nDMaterial InitDefGrad " << iData[0]
               << " : -F0 needs 9 components (row-major F11..F33)\n";
        return 0;
      }
      int n9 = 9;
      if (OPS_GetDoubleInput(&n9, F0in) != 0) {
        opserr << "WARNING nDMaterial InitDefGrad " << iData[0]
               << " : invalid -F0 components\n";
        return 0;
      }
      haveF0 = true;
    } else {
      opserr << "WARNING nDMaterial InitDefGrad " << iData[0]
             << " : unknown option '" << opt << "'\n";
      return 0;
    }
  }

  // A singular -F0 is refused here; the constructor would otherwise exit.
  if (haveF0) {
    double d = F0in[0]*(F0in[4]*F0in[8] - F0in[5]*F0in[7])
             - F0in[1]*(F0in[3]*F0in[8] - F0in[5]*F0in[6])
             + F0in[2]*(F0in[3]*F0in[7] - F0in[4]*F0in[6]);
    if (d == 0.0) {
      opserr << "WARNING nDMaterial InitDefGrad " << iData[0]
             << " : -F0 is singular (det F0 = 0)\n";
      return 0;
    }
  }

  return new InitDefGradNDMaterial(iData[0], *inner, active, haveF0 ? F0in : 0);
}

// =========================================================================== //
//  Construction                                                               //
// =========================================================================== //
InitDefGradNDMaterial::InitDefGradNDMaterial(int tag, NDMaterial &inner,
                                             bool active_, const double *F0given)
  : FiniteStrainNDMaterial(tag, ND_TAG_InitDefGradNDMaterial),
    theMaterial(0), active(active_), captured(false), f0Explicit(false),
    Frel(3, 3), birthF0(9)
{
  NDMaterial *cp = (strncmp(inner.getType(), "ThreeDimensional", 80) == 0)
                     ? inner.getCopy() : inner.getCopy("ThreeDimensional");
  theMaterial = dynamic_cast<FiniteStrainNDMaterial *>(cp);
  if (theMaterial == 0 || theMaterial->getOrder() != 6) {
    opserr << "InitDefGradNDMaterial: inner must be a 3D (order-6) finite-strain "
              "material\n";
    if (cp != 0) delete cp;
    exit(-1);
  }

  setIdentity9(F0);
  setIdentity9(invF0);
  Frel.Zero(); Frel(0,0) = Frel(1,1) = Frel(2,2) = 1.0;
  birthF0.Zero();

  if (F0given != 0) {                 // explicit F0: captured now
    for (int i = 0; i < 9; i++) F0[i] = F0given[i];
    if (inv9(F0, invF0) < 0) {
      opserr << "InitDefGradNDMaterial: supplied -F0 is singular (det F0 = 0)\n";
      exit(-1);
    }
    captured   = true;
    f0Explicit = true;
  }
}

InitDefGradNDMaterial::InitDefGradNDMaterial()
  : FiniteStrainNDMaterial(0, ND_TAG_InitDefGradNDMaterial),
    theMaterial(0), active(true), captured(false), f0Explicit(false),
    Frel(3, 3), birthF0(9)
{
  setIdentity9(F0);
  setIdentity9(invF0);
  Frel.Zero(); Frel(0,0) = Frel(1,1) = Frel(2,2) = 1.0;
  birthF0.Zero();
}

InitDefGradNDMaterial::~InitDefGradNDMaterial()
{
  if (theMaterial) delete theMaterial;
}

// =========================================================================== //
//  setTrialF: capture F0 at activation, forward F_rel = F F0^-1                //
// =========================================================================== //
int InitDefGradNDMaterial::setTrialF(const Matrix &F)
{
  double Fin[9];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) Fin[3*i+j] = F(i, j);

  if (!active) {                       // pass-through (-noInitF)
    return theMaterial->setTrialF(F);
  }

  if (!captured) {                     // first F seen is the activation F0
    for (int i = 0; i < 9; i++) F0[i] = Fin[i];
    if (inv9(F0, invF0) < 0) {
      opserr << "InitDefGradNDMaterial::setTrialF - birth deformation gradient is "
                "singular (det F0 = " << det9(F0) << "); cannot re-reference\n";
      return -1;
    }
    captured = true;
  }

  double Fr[9];
  mul9(Fin, invF0, Fr);                // F_rel = F F0^-1
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) Frel(i, j) = Fr[3*i+j];

  return theMaterial->setTrialF(Frel);
}

// =========================================================================== //
//  Query: delegated to the inner finite-strain material                       //
// =========================================================================== //
const Vector &InitDefGradNDMaterial::getStress(void)  { return theMaterial->getStress(); }
const Matrix &InitDefGradNDMaterial::getTangent(void) { return theMaterial->getTangent(); }
const Vector &InitDefGradNDMaterial::getStrain(void)  { return theMaterial->getStrain(); }
double        InitDefGradNDMaterial::getRho(void)     { return theMaterial->getRho(); }
double        InitDefGradNDMaterial::getJ(void) const { return theMaterial->getJ(); }

const Matrix &InitDefGradNDMaterial::getInitialTangent(void)
{
  return theMaterial->getInitialTangent();
}

int InitDefGradNDMaterial::getSpatialTangentTensor(double c[3][3][3][3])
{
  return theMaterial->getSpatialTangentTensor(c);
}

// =========================================================================== //
//  State: F0 is fixed after capture; the inner owns its own state             //
// =========================================================================== //
int InitDefGradNDMaterial::commitState(void)
{
  return theMaterial->commitState();
}

int InitDefGradNDMaterial::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int InitDefGradNDMaterial::revertToStart(void)
{
  // Return to the as-constructed state. Auto-capture mode re-arms (a genuine
  // restart re-captures birth at the next setTrialF); an explicitly supplied F0
  // is kept (it is a user-specified property, not trial state).
  if (!f0Explicit) {
    captured = false;
    setIdentity9(F0);
    setIdentity9(invF0);
  }
  return theMaterial->revertToStart();
}

// =========================================================================== //
//  Copies                                                                     //
// =========================================================================== //
NDMaterial *InitDefGradNDMaterial::getCopy(void)
{
  InitDefGradNDMaterial *c =
    new InitDefGradNDMaterial(this->getTag(), *theMaterial, active,
                              f0Explicit ? F0 : 0);
  // carry the live capture state (auto-captured F0 is not passed through the ctor)
  c->captured   = captured;
  c->f0Explicit = f0Explicit;
  for (int i = 0; i < 9; i++) { c->F0[i] = F0[i]; c->invF0[i] = invF0[i]; }
  return c;
}

NDMaterial *InitDefGradNDMaterial::getCopy(const char *type)
{
  if (strncmp(type, "ThreeDimensional", 80) == 0)
    return this->getCopy();
  opserr << "InitDefGradNDMaterial::getCopy - only ThreeDimensional is supported\n";
  return 0;
}

// =========================================================================== //
//  Parallel / database                                                        //
// =========================================================================== //
int InitDefGradNDMaterial::sendSelf(int cTag, Channel &theChannel)
{
  if (theMaterial == 0) return -1;
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) { matDbTag = theChannel.getDbTag(); theMaterial->setDbTag(matDbTag); }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) return -1;

  // flags (3) + F0 (9) = 12; invF0 is rebuilt on recv
  static Vector dataVec(12);
  dataVec(0) = active     ? 1.0 : 0.0;
  dataVec(1) = captured   ? 1.0 : 0.0;
  dataVec(2) = f0Explicit ? 1.0 : 0.0;
  for (int i = 0; i < 9; i++) dataVec(3+i) = F0[i];
  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) return -2;

  if (theMaterial->sendSelf(cTag, theChannel) < 0) return -3;
  return 0;
}

int InitDefGradNDMaterial::recvSelf(int cTag, Channel &theChannel,
                                    FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) return -1;
  this->setTag(dataID(0));

  if (theMaterial == 0) {
    NDMaterial *m = theBroker.getNewNDMaterial(dataID(1));
    theMaterial = dynamic_cast<FiniteStrainNDMaterial *>(m);
    if (theMaterial == 0) {
      opserr << "InitDefGradNDMaterial::recvSelf - cannot create inner finite-strain "
                "material classTag " << dataID(1) << endln;
      if (m != 0) delete m;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(12);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) return -3;
  active     = (dataVec(0) != 0.0);
  captured   = (dataVec(1) != 0.0);
  f0Explicit = (dataVec(2) != 0.0);
  for (int i = 0; i < 9; i++) F0[i] = dataVec(3+i);
  setIdentity9(invF0);
  if (captured && inv9(F0, invF0) < 0) {       // rebuild F0^-1
    opserr << "InitDefGradNDMaterial::recvSelf - received singular F0\n";
    return -4;
  }

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) return -5;
  return 0;
}

// =========================================================================== //
//  Misc                                                                       //
// =========================================================================== //
void InitDefGradNDMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"InitDefGradNDMaterial\", ";
    s << "\"inner\": \"" << theMaterial->getTag() << "\", ";
    s << "\"active\": " << (active ? "true" : "false") << ", ";
    s << "\"captured\": " << (captured ? "true" : "false");
    s << "}";
  } else {
    s << "InitDefGradNDMaterial (multiplicative staged-activation wrapper), tag: "
      << this->getTag() << endln;
    s << "\tinner material tag: " << theMaterial->getTag() << endln;
    s << "\tactive: " << (active ? "yes" : "no")
      << ", F0 captured: " << (captured ? "yes" : "no") << endln;
  }
}

int InitDefGradNDMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

// Responses: "F0" (aliases "birthF", "F0capture") returns the captured
// activation gradient (9, row-major). Every other response (stress = Cauchy,
// strain = Hencky strain relative to activation, tangent, inner-specific
// output) is delegated to the inner material.
Response *InitDefGradNDMaterial::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (argc > 0 &&
      (strcmp(argv[0], "F0") == 0 || strcmp(argv[0], "birthF") == 0 ||
       strcmp(argv[0], "F0capture") == 0)) {
    for (int i = 0; i < 9; i++) birthF0(i) = F0[i];
    return new MaterialResponse(this, 101, birthF0);
  }
  return theMaterial->setResponse(argv, argc, output);
}

int InitDefGradNDMaterial::getResponse(int responseID, Information &matInfo)
{
  if (responseID == 101) {
    for (int i = 0; i < 9; i++) birthF0(i) = F0[i];
    return matInfo.setVector(birthF0);
  }
  return theMaterial->getResponse(responseID, matInfo);
}
