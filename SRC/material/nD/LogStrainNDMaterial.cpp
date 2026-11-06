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
// LogStrainNDMaterial: logarithmic (Hencky) strain-space finite-strain
// material. See LogStrainNDMaterial.h for the formulation, references and
// limitations; the numerical core is LogStrainKernel.h.

#include <LogStrainNDMaterial.h>
#include <LogStrainKernel.h>
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

using namespace logstrain_kernel;

// =========================================================================== //
//  Factory:  nDMaterial LogStrain $tag $innerTag                              //
// =========================================================================== //
void *OPS_LogStrainNDMaterial(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING invalid args: nDMaterial LogStrain $tag $innerTag\n";
    return 0;
  }
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid ints: nDMaterial LogStrain $tag $innerTag\n";
    return 0;
  }
  NDMaterial *inner = OPS_getNDMaterial(iData[1]);
  if (inner == 0) {
    opserr << "WARNING nDMaterial LogStrain " << iData[0]
           << " : inner nDMaterial " << iData[1] << " not found\n";
    return 0;
  }
  // Check that the inner can produce a 3D (order 6) copy before constructing,
  // so invalid input fails the command (return 0) instead of reaching the
  // constructor's exit(-1).
  NDMaterial *probe = (strncmp(inner->getType(), "ThreeDimensional", 80) == 0)
                        ? inner->getCopy() : inner->getCopy("ThreeDimensional");
  if (probe == 0 || probe->getOrder() != 6) {
    opserr << "WARNING nDMaterial LogStrain " << iData[0]
           << " : inner nDMaterial " << iData[1]
           << " must be a 3D (order-6) material\n";
    if (probe != 0) delete probe;
    return 0;
  }
  delete probe;   // the constructor makes its own copy
  return new LogStrainNDMaterial(iData[0], *inner);
}

// =========================================================================== //
//  Construction                                                               //
// =========================================================================== //
void LogStrainNDMaterial::setIdentity(double M[9]) {
  for (int i = 0; i < 9; i++) M[i] = 0.0;
  M[0] = M[4] = M[8] = 1.0;
}

LogStrainNDMaterial::LogStrainNDMaterial(int tag, NDMaterial &inner)
  : FiniteStrainNDMaterial(tag, ND_TAG_LogStrainNDMaterial),
    theMaterial(0), sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), Jdet(1.0)
{
  if (strncmp(inner.getType(), "ThreeDimensional", 80) == 0)
    theMaterial = inner.getCopy();
  else
    theMaterial = inner.getCopy("ThreeDimensional");

  if (theMaterial == 0 || theMaterial->getOrder() != 6) {
    opserr << "LogStrainNDMaterial: inner material must be 3D (order 6)\n";
    exit(-1);
  }
  setIdentity(Fn);
  setIdentity(Be_n);
  setIdentity(Ftrial9);
  setIdentity(BeTrial9);
  setIdentity(Be_trialUpd);
  for (int k = 0; k < 6; k++) { epsFeed_n[k] = 0.0; epsFeedTrial[k] = 0.0; }
  sigmaCauchy.Zero();
  henckyStrain.Zero();
  aTangent.Zero();
}

LogStrainNDMaterial::LogStrainNDMaterial()
  : FiniteStrainNDMaterial(0, ND_TAG_LogStrainNDMaterial),
    theMaterial(0), sigmaCauchy(6), henckyStrain(6), aTangent(6, 6), Jdet(1.0)
{
  setIdentity(Fn);
  setIdentity(Be_n);
  setIdentity(Ftrial9);
  setIdentity(BeTrial9);
  setIdentity(Be_trialUpd);
  for (int k = 0; k < 6; k++) { epsFeed_n[k] = 0.0; epsFeedTrial[k] = 0.0; }
}

LogStrainNDMaterial::~LogStrainNDMaterial()
{
  if (theMaterial) delete theMaterial;
}

// =========================================================================== //
//  setTrialF: Box 14.3 (i), (ii), (iv) and the Section 14.5 tangent           //
// =========================================================================== //
int LogStrainNDMaterial::setTrialF(const Matrix &F)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) Ftrial9[3*i+j] = F(i, j);

  // Reject a non-positive Jacobian: det F <= 0 (inverted or degenerate
  // element) would give a sign-flipped Cauchy stress tau/J and a Be_trial
  // that is not positive definite, whose 1/2 ln has -inf/NaN eigenvalues.
  // A negative return lets the element reject the step.
  Jdet = mat3_det(Ftrial9);
  if (Jdet <= 0.0) {
    opserr << "LogStrainNDMaterial::setTrialF - non-positive det F (" << Jdet
           << "); element inverted/degenerate\n";
    return -1;
  }

  // (i)+(ii) trial elastic left Cauchy-Green tensor and trial Hencky strain
  double Betr[9];
  trial_Be(Ftrial9, Fn, Be_n, Betr);
  for (int i = 0; i < 9; i++) BeTrial9[i] = Betr[i];   // kept for the Section 14.5 tangent
  double epsTr6[6], epsN6[6];
  hencky_voigt(Betr, epsTr6);                          // eps_trial = 1/2 ln Be_trial
  hencky_voigt(Be_n, epsN6);                           // eps_e_n = 1/2 ln Be_n
  for (int k = 0; k < 6; k++) henckyStrain(k) = epsTr6[k];

  // (iii) Plastic inner: cancel the inner's own plastic-strain subtraction.
  // Feed eps_feed = eps_feed_n + (eps_trial - eps_e_n); the inner computes
  // eps_feed - eps_p_inner = eps_trial exactly (the invariant
  // eps_p_inner = eps_feed_n - eps_e_n is maintained on commit). For an
  // elastic inner (eps_p_inner = 0) this reduces to eps_feed = eps_trial.
  static Vector epsFeedV(6);
  for (int k = 0; k < 6; k++) epsFeedV(k) = epsFeed_n[k] + (epsTr6[k] - epsN6[k]);

  if (theMaterial->setTrialStrain(epsFeedV) < 0) {
    opserr << "LogStrainNDMaterial::setTrialF - inner material failed to set "
              "the trial strain\n";
    return -1;
  }
  const Vector &tauV = theMaterial->getStress();   // Kirchhoff stress tau (6)
  const Matrix &D6m  = theMaterial->getTangent();  // d(tau)/d(eps_e) (6x6, elastoplastic)
  double tau6[6], D6[36];
  for (int k = 0; k < 6; k++) tau6[k] = tauV(k);
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) D6[6*I+J] = D6m(I, J);

  // (iv) Cauchy stress tau/J and spatial modulus c = (1/2J)[D:L:B] at
  // Be_trial (Jdet > 0 was checked at entry)
  double sig6[6], c6[36];
  assemble_material(Betr, D6, tau6, Jdet, sig6, c6);
  for (int k = 0; k < 6; k++) sigmaCauchy(k) = sig6[k];
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) aTangent(I, J) = c6[6*I+J];

  // Updated elastic strain eps_e_{n+1} = C_e : tau, with C_e the inverse of
  // the inner initial tangent, then Be_{n+1} = exp(2 eps_e_{n+1}). Exact for
  // an inner with a linear, constant elastic law (tau = D_e : eps_e holds in
  // elastic and plastic steps); see the limitations in the header.
  static Matrix Ce(6, 6);
  Matrix D0(theMaterial->getInitialTangent());
  if (D0.Invert(Ce) < 0) {
    opserr << "LogStrainNDMaterial::setTrialF - inner initial tangent not invertible\n";
    return -1;
  }
  static Vector epsEnp1(6);
  epsEnp1.addMatrixVector(0.0, Ce, tauV, 1.0);     // eps_e_{n+1} = C_e tau (eng. Voigt)
  double epsEnp16[6];
  for (int k = 0; k < 6; k++) epsEnp16[k] = epsEnp1(k);
  be_from_hencky_voigt(epsEnp16, Be_trialUpd);     // Be to commit

  for (int k = 0; k < 6; k++) epsFeedTrial[k] = epsFeedV(k);  // stage feed
  return 0;
}

// =========================================================================== //
//  Query                                                                      //
// =========================================================================== //
const Vector &LogStrainNDMaterial::getStress(void)  { return sigmaCauchy; }
const Matrix &LogStrainNDMaterial::getTangent(void) { return aTangent; }

// Full fourth-order spatial modulus c_ijkl = (1/2J)[D:L:B], the tangent an
// element uses for its consistent stiffness (the 6x6 getTangent() is lossy in
// (k,l)). Evaluated at the trial elastic left Cauchy-Green tensor (Section
// 14.5 uses the trial value, not the returned one; they differ in a plastic
// step), the inner tangent D and J from the most recent setTrialF().
int LogStrainNDMaterial::getSpatialTangentTensor(double c[3][3][3][3])
{
  const Matrix &D6m = theMaterial->getTangent();   // d(tau)/d(eps_e) (6x6)
  double D6[36];
  for (int I = 0; I < 6; I++)
    for (int J = 0; J < 6; J++) D6[6 * I + J] = D6m(I, J);
  spatial_tangent_full(BeTrial9, D6, Jdet, c);
  return 0;
}

const Vector &LogStrainNDMaterial::getStrain(void)  { return henckyStrain; }
double        LogStrainNDMaterial::getRho(void)     { return theMaterial->getRho(); }

const Matrix &LogStrainNDMaterial::getInitialTangent(void)
{
  return theMaterial->getInitialTangent();
}

// =========================================================================== //
//  State (Be_n, F_n and the fed strain are owned here; the inner is wrapped)  //
// =========================================================================== //
int LogStrainNDMaterial::commitState(void)
{
  for (int i = 0; i < 9; i++) { Fn[i] = Ftrial9[i]; Be_n[i] = Be_trialUpd[i]; }
  for (int k = 0; k < 6; k++) epsFeed_n[k] = epsFeedTrial[k];   // protocol state
  return theMaterial->commitState();
}

int LogStrainNDMaterial::revertToLastCommit(void)
{
  // discard the trial state so a commit without a new setTrialF() keeps the
  // committed one
  for (int i = 0; i < 9; i++) {
    Ftrial9[i] = Fn[i]; BeTrial9[i] = Be_n[i]; Be_trialUpd[i] = Be_n[i];
  }
  for (int k = 0; k < 6; k++) epsFeedTrial[k] = epsFeed_n[k];
  return theMaterial->revertToLastCommit();
}

int LogStrainNDMaterial::revertToStart(void)
{
  setIdentity(Fn);
  setIdentity(Be_n);
  setIdentity(Ftrial9);
  setIdentity(BeTrial9);
  setIdentity(Be_trialUpd);
  for (int k = 0; k < 6; k++) { epsFeed_n[k] = 0.0; epsFeedTrial[k] = 0.0; }
  sigmaCauchy.Zero();
  henckyStrain.Zero();
  aTangent.Zero();
  Jdet = 1.0;
  return theMaterial->revertToStart();
}

// =========================================================================== //
//  Copies                                                                     //
// =========================================================================== //
NDMaterial *LogStrainNDMaterial::getCopy(void)
{
  LogStrainNDMaterial *c = new LogStrainNDMaterial(this->getTag(), *theMaterial);
  for (int i = 0; i < 9; i++) { c->Fn[i] = Fn[i]; c->Be_n[i] = Be_n[i]; }
  for (int k = 0; k < 6; k++) c->epsFeed_n[k] = epsFeed_n[k];
  return c;
}

NDMaterial *LogStrainNDMaterial::getCopy(const char *type)
{
  if (strncmp(type, "ThreeDimensional", 80) == 0)
    return this->getCopy();
  opserr << "LogStrainNDMaterial::getCopy - only ThreeDimensional is supported\n";
  return 0;
}

// =========================================================================== //
//  Parallel / database                                                        //
// =========================================================================== //
int LogStrainNDMaterial::sendSelf(int cTag, Channel &theChannel)
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

  static Vector dataVec(24);
  for (int i = 0; i < 9; i++) { dataVec(i) = Fn[i]; dataVec(9+i) = Be_n[i]; }
  for (int k = 0; k < 6; k++) dataVec(18+k) = epsFeed_n[k];
  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) return -2;

  if (theMaterial->sendSelf(cTag, theChannel) < 0) return -3;
  return 0;
}

int LogStrainNDMaterial::recvSelf(int cTag, Channel &theChannel,
                                  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) return -1;
  this->setTag(dataID(0));

  if (theMaterial == 0) {
    theMaterial = theBroker.getNewNDMaterial(dataID(1));
    if (theMaterial == 0) {
      opserr << "LogStrainNDMaterial::recvSelf - cannot create inner material classTag "
             << dataID(1) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(24);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) return -3;
  for (int i = 0; i < 9; i++) { Fn[i] = dataVec(i); Be_n[i] = dataVec(9+i); }
  for (int k = 0; k < 6; k++) epsFeed_n[k] = dataVec(18+k);

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) return -4;
  return 0;
}

// =========================================================================== //
//  Misc                                                                       //
// =========================================================================== //
void LogStrainNDMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LogStrainNDMaterial\", ";
    s << "\"inner\": \"" << theMaterial->getTag() << "\"";
    s << "}";
  } else {
    s << "LogStrainNDMaterial (Hencky log-strain finite-strain adaptor), tag: "
      << this->getTag() << endln;
    s << "\tinner material tag: " << theMaterial->getTag() << endln;
  }
}

int LogStrainNDMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

// Responses: "stress", "strain" and "tangent" return the finite-strain
// measures of this material (Cauchy stress tau/J, Hencky strain, spatial c).
// The inner, queried directly, would report Kirchhoff stress, the fed strain
// and the small-strain D. Other responses (back stress, plastic strain,
// damage) are delegated to the inner.
Response *LogStrainNDMaterial::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (argc > 0) {
    if (strcmp(argv[0], "stress")  == 0 || strcmp(argv[0], "stresses") == 0 ||
        strcmp(argv[0], "Cauchy")  == 0 || strcmp(argv[0], "cauchy")   == 0)
      return new MaterialResponse(this, 1, sigmaCauchy);    // Cauchy stress (6)

    if (strcmp(argv[0], "strain")  == 0 || strcmp(argv[0], "strains")  == 0)
      return new MaterialResponse(this, 2, henckyStrain);   // Hencky strain (6)

    if (strcmp(argv[0], "tangent") == 0)
      return new MaterialResponse(this, 3, aTangent);       // 6x6 projection of c
  }
  return theMaterial->setResponse(argv, argc, output);
}

int LogStrainNDMaterial::getResponse(int responseID, Information &matInfo)
{
  switch (responseID) {
  case 1: return matInfo.setVector(sigmaCauchy);
  case 2: return matInfo.setVector(henckyStrain);
  case 3: return matInfo.setMatrix(aTangent);
  default: return theMaterial->getResponse(responseID, matInfo);
  }
}
