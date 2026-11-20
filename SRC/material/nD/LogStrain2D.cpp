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
// LogStrain2D: plane (2D) logarithmic-strain finite-strain material. See
// LogStrain2D.h for the formulation and references.

#include <LogStrain2D.h>
#include <LogStrainNDMaterial.h>
#include <classTags.h>
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
#include <math.h>

// =========================================================================== //
//  Factory:  nDMaterial LogStrain2D $tag $innerTag <-planeStrain|-planeStress> //
// =========================================================================== //
void *OPS_LogStrain2D(void)
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING invalid args: nDMaterial LogStrain2D $tag $innerTag "
              "<-planeStrain|-planeStress>\n";
    return 0;
  }
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid ints: nDMaterial LogStrain2D $tag $innerTag\n";
    return 0;
  }

  int planeType = LogStrain2D::PlaneStrain;   // default: plane strain
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();
    if (strcmp(opt, "-planeStress") == 0 || strcmp(opt, "-PlaneStress") == 0)
      planeType = LogStrain2D::PlaneStress;
    else if (strcmp(opt, "-planeStrain") == 0 || strcmp(opt, "-PlaneStrain") == 0)
      planeType = LogStrain2D::PlaneStrain;
    else
      opserr << "WARNING nDMaterial LogStrain2D " << iData[0]
             << " : ignoring unknown option '" << opt << "'\n";
  }

  NDMaterial *inner = OPS_getNDMaterial(iData[1]);
  if (inner == 0) {
    opserr << "WARNING nDMaterial LogStrain2D " << iData[0]
           << " : inner nDMaterial " << iData[1] << " not found\n";
    return 0;
  }
  // Validate the inner can produce a 3D (order-6) copy BEFORE constructing, so bad
  // user input fails the command gracefully instead of reaching a hard exit(-1).
  NDMaterial *probe = (strncmp(inner->getType(), "ThreeDimensional", 80) == 0)
                        ? inner->getCopy() : inner->getCopy("ThreeDimensional");
  if (probe == 0 || probe->getOrder() != 6) {
    opserr << "WARNING nDMaterial LogStrain2D " << iData[0]
           << " : inner nDMaterial " << iData[1]
           << " must be a 3D (order-6) material\n";
    if (probe != 0) delete probe;
    return 0;
  }
  delete probe;
  return new LogStrain2D(iData[0], *inner, planeType);
}

// =========================================================================== //
//  Construction                                                               //
// =========================================================================== //
LogStrain2D::LogStrain2D(int tag, NDMaterial &inner3D, int ptype)
  : FiniteStrainND2DMaterial(tag, ND_TAG_LogStrain2D),
    the3D(0), planeType(ptype), lam_n(1.0), lamTrial(1.0),
    sigma3(3), hencky3(3), tangent3(3, 3), Jarea(1.0)
{
  the3D = new LogStrainNDMaterial(tag, inner3D);   // deep-copies the inner
  sigma3.Zero();
  hencky3.Zero();
  tangent3.Zero();
  for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
    for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++) c2[i][j][k][l] = 0.0;
}

LogStrain2D::LogStrain2D()
  : FiniteStrainND2DMaterial(0, ND_TAG_LogStrain2D),
    the3D(0), planeType(PlaneStrain), lam_n(1.0), lamTrial(1.0),
    sigma3(3), hencky3(3), tangent3(3, 3), Jarea(1.0)
{
  sigma3.Zero();
  hencky3.Zero();
  tangent3.Zero();
  for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
    for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++) c2[i][j][k][l] = 0.0;
}

LogStrain2D::~LogStrain2D()
{
  if (the3D) delete the3D;
}

// =========================================================================== //
//  Kinematic lift: build the 3x3 F (F33 given, no out-of-plane shear) and     //
//  drive the internal 3D material.                                            //
// =========================================================================== //
int LogStrain2D::driveLifted(const Matrix &F2, double F33)
{
  static Matrix F3(3, 3);
  F3.Zero();
  F3(0, 0) = F2(0, 0);  F3(0, 1) = F2(0, 1);
  F3(1, 0) = F2(1, 0);  F3(1, 1) = F2(1, 1);
  F3(2, 2) = F33;
  return the3D->setTrialF(F3);
}

// =========================================================================== //
//  setTrialF                                                                  //
// =========================================================================== //
int LogStrain2D::setTrialF(const Matrix &F)
{
  Jarea = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
  if (Jarea <= 0.0) {
    opserr << "LogStrain2D::setTrialF - non-positive det(2x2 F) (" << Jarea
           << "); element inverted/degenerate\n";
    return -1;
  }

  if (planeType == PlaneStrain) {
    if (driveLifted(F, 1.0) < 0) return -1;
    lamTrial = 1.0;
    return buildOutputs(F, /*condense=*/false);
  }

  // ---- Plane stress: solve lambda = F33 from sigma_33(lambda) = 0 (Newton) //
  // The residual is tested on every iterate including the one produced by the
  // final update (loop bound maxit+1, update gated by it < maxit) so a step that
  // converges exactly on the last allowed update is not falsely reported as a
  // non-convergence.
  const int maxit = 30;
  const double tol = 1.0e-12;
  double lam = (lam_n > 0.0) ? lam_n : 1.0;
  bool converged = false;
  for (int it = 0; it <= maxit; it++) {
    if (driveLifted(F, lam) < 0) return -1;
    const Vector &s = the3D->getStress();               // Cauchy stress (6): sigma_33 = s(2)
    double r = s(2);
    double sref = 1.0;
    if (fabs(s(0)) > sref) sref = fabs(s(0));
    if (fabs(s(1)) > sref) sref = fabs(s(1));
    if (fabs(r) <= tol * sref) { converged = true; break; }
    if (it == maxit) break;                             // no updates left to try

    double h = 1.0e-8 * ((fabs(lam) > 1.0) ? fabs(lam) : 1.0);
    if (driveLifted(F, lam + h) < 0) return -1;
    double rp = the3D->getStress()(2);
    double drdl = (rp - r) / h;
    if (drdl == 0.0) {
      opserr << "LogStrain2D::setTrialF - plane-stress Newton: singular d(sigma33)/d(lambda)\n";
      return -1;
    }
    double lamNew = lam - r / drdl;
    if (lamNew <= 0.0) lamNew = 0.5 * lam;               // never cross zero stretch
    lam = lamNew;
  }
  if (!converged) {
    opserr << "LogStrain2D::setTrialF - plane-stress Newton did not converge\n";
    return -1;
  }

  // re-evaluate at the converged lambda so the cached trial state (Be_trial,
  // inner tau and D) is the accepted one before reading stress and tangent
  if (driveLifted(F, lam) < 0) return -1;
  lamTrial = lam;
  return buildOutputs(F, /*condense=*/true);
}

// =========================================================================== //
//  In-plane stress (3), Hencky strain (3), lossy 3x3 c and full c2 tensor    //
//  from the most recent setTrialF of the internal 3D material.               //
// =========================================================================== //
int LogStrain2D::buildOutputs(const Matrix &F2, bool condense)
{
  const Vector &s6 = the3D->getStress();      // Cauchy stress, Voigt {11,22,33,12,23,31}
  const Vector &e6 = the3D->getStrain();      // Hencky strain, engineering Voigt
  sigma3(0) = s6(0); sigma3(1) = s6(1); sigma3(2) = s6(3);
  hencky3(0) = e6(0); hencky3(1) = e6(1); hencky3(2) = e6(3);

  double c3[3][3][3][3];
  the3D->getSpatialTangentTensor(c3);         // full 3D material modulus c = (1/2J)[D:L:B]

  if (!condense) {
    // Plane strain: c2 is the in-plane block of c. The element forms
    // a_ijkl = c2_ijkl - sigma_il delta_jk (i,j,k,l in {0,1}), which equals
    // the in-plane block of the 3D a, since plane strain is 3D with F33 = 1.
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++)
        c2[i][j][k][l] = c3[i][j][k][l];
    voigt3FromTensor(c2, tangent3);
    return 0;
  }

  // Plane stress: condense the out-of-plane direction (index 2) enforcing
  // delta sigma_33 = 0 on the total spatial tangent a = c - sigma_il delta_jk,
  // then add back the in-plane geometric term, so the element's
  // a = c2 - sigma_il delta_jk recovers the condensed total tangent:
  //     a_ijkl  = a3_ijkl - a3_ij22 a3_22kl / a3_2222   (i,j,k,l in {0,1})
  //     c2_ijkl = a_ijkl + sigma_il delta_jk
  double sig[3][3];
  sig[0][0] = s6(0); sig[1][1] = s6(1); sig[2][2] = s6(2);
  sig[0][1] = sig[1][0] = s6(3);
  sig[1][2] = sig[2][1] = s6(4);
  sig[2][0] = sig[0][2] = s6(5);

  // full 3D total spatial tangent a3 = c3 - sigma_il delta_jk
  double a3[3][3][3][3];
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++)
      a3[i][j][k][l] = c3[i][j][k][l] - sig[i][l] * (j == k ? 1.0 : 0.0);

  // Guard the condensation pivot: a_2222 can approach 0 for a nearly
  // perfectly plastic inner under out-of-plane flow, and the division would
  // put NaN into c2. Return -1 so the element rejects the step.
  double denom = a3[2][2][2][2];
  double aScale = 1.0;
  if (fabs(a3[0][0][0][0]) > aScale) aScale = fabs(a3[0][0][0][0]);
  if (fabs(a3[1][1][1][1]) > aScale) aScale = fabs(a3[1][1][1][1]);
  if (fabs(denom) < 1.0e-12 * aScale) {
    opserr << "LogStrain2D::setTrialF - plane-stress condensation: out-of-plane "
              "modulus a_2222 (" << denom << ") ~ 0; cannot condense sigma33 = 0\n";
    return -1;
  }

  for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
    for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++) {
      double a2 = a3[i][j][k][l] - a3[i][j][2][2] * a3[2][2][k][l] / denom;
      c2[i][j][k][l] = a2 + sig[i][l] * (j == k ? 1.0 : 0.0);
    }

  voigt3FromTensor(c2, tangent3);
  return 0;
}

// engineering-Voigt {11,22,12} projection of the (2,2,2,2) modulus (lossy in (k,l))
void LogStrain2D::voigt3FromTensor(const double c[2][2][2][2], Matrix &D)
{
  static const int RI[3] = {0, 1, 0};
  static const int RJ[3] = {0, 1, 1};
  for (int I = 0; I < 3; I++)
    for (int J = 0; J < 3; J++)
      D(I, J) = c[RI[I]][RJ[I]][RI[J]][RJ[J]];
}

// =========================================================================== //
//  Query                                                                      //
// =========================================================================== //
const Vector &LogStrain2D::getStress(void)  { return sigma3; }
const Matrix &LogStrain2D::getTangent(void) { return tangent3; }
const Vector &LogStrain2D::getStrain(void)  { return hencky3; }
double        LogStrain2D::getRho(void)     { return the3D->getRho(); }

int LogStrain2D::getSpatialTangentTensor2D(double c[2][2][2][2])
{
  for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
    for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++)
      c[i][j][k][l] = c2[i][j][k][l];
  return 0;
}

// Initial in-plane tangent: condense the inner's small-strain 6x6 initial modulus
// to the plane face (plane strain = the {11,22,12} sub-block; plane stress =
// static condensation eliminating {33,23,31}).
const Matrix &LogStrain2D::getInitialTangent(void)
{
  static Matrix D2(3, 3);
  const Matrix &D6 = the3D->getInitialTangent();      // inner 6x6 (engineering Voigt)
  const int a[3] = {0, 1, 3};                         // in-plane Voigt {11,22,12}
  if (planeType == PlaneStrain) {
    for (int I = 0; I < 3; I++)
      for (int J = 0; J < 3; J++) D2(I, J) = D6(a[I], a[J]);
    return D2;
  }
  // plane stress: static condensation eliminating b = {33,23,31} (zero stress there)
  const int b[3] = {2, 4, 5};
  Matrix Daa(3, 3), Dab(3, 3), Dba(3, 3), Dbb(3, 3), Dbb_i(3, 3);
  for (int I = 0; I < 3; I++)
    for (int J = 0; J < 3; J++) {
      Daa(I, J) = D6(a[I], a[J]);
      Dab(I, J) = D6(a[I], b[J]);
      Dba(I, J) = D6(b[I], a[J]);
      Dbb(I, J) = D6(b[I], b[J]);
    }
  if (Dbb.Invert(Dbb_i) < 0) { D2 = Daa; return D2; }
  Matrix tmp(3, 3);
  tmp = Dab * Dbb_i * Dba;
  D2 = Daa - tmp;
  return D2;
}

// =========================================================================== //
//  State cycle                                                                //
// =========================================================================== //
int LogStrain2D::commitState(void)
{
  lam_n = lamTrial;
  return the3D->commitState();
}

int LogStrain2D::revertToLastCommit(void)
{
  lamTrial = lam_n;   // drop a diverged trial lambda; setTrialF starts from lam_n
  return the3D->revertToLastCommit();
}

int LogStrain2D::revertToStart(void)
{
  lam_n = 1.0; lamTrial = 1.0;
  sigma3.Zero(); hencky3.Zero(); tangent3.Zero(); Jarea = 1.0;
  for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
    for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++) c2[i][j][k][l] = 0.0;
  return the3D->revertToStart();
}

// =========================================================================== //
//  Copies                                                                     //
// =========================================================================== //
NDMaterial *LogStrain2D::getCopy(void)
{
  LogStrain2D *c = new LogStrain2D();
  c->setTag(this->getTag());
  c->the3D = (LogStrainNDMaterial *)the3D->getCopy();   // deep copy incl. 3D state
  c->planeType = planeType;
  c->lam_n = lam_n;
  c->lamTrial = lamTrial;
  return c;
}

NDMaterial *LogStrain2D::getCopy(const char *type)
{
  if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStress") == 0) {
    LogStrain2D *c = (LogStrain2D *)this->getCopy();
    c->planeType = (strcmp(type, "PlaneStress") == 0) ? PlaneStress : PlaneStrain;
    return c;
  }
  opserr << "LogStrain2D::getCopy - only PlaneStrain / PlaneStress supported "
            "(got '" << type << "')\n";
  return 0;
}

// =========================================================================== //
//  Parallel / database                                                        //
// =========================================================================== //
int LogStrain2D::sendSelf(int cTag, Channel &theChannel)
{
  if (the3D == 0) return -1;
  int dbTag = this->getDbTag();

  static ID dataID(4);
  dataID(0) = this->getTag();
  dataID(1) = the3D->getClassTag();
  int matDbTag = the3D->getDbTag();
  if (matDbTag == 0) { matDbTag = theChannel.getDbTag(); the3D->setDbTag(matDbTag); }
  dataID(2) = matDbTag;
  dataID(3) = planeType;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) return -1;

  static Vector dataVec(1);
  dataVec(0) = lam_n;
  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) return -2;

  if (the3D->sendSelf(cTag, theChannel) < 0) return -3;
  return 0;
}

int LogStrain2D::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID dataID(4);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) return -1;
  this->setTag(dataID(0));
  planeType = dataID(3);

  if (the3D == 0) {
    the3D = (LogStrainNDMaterial *)theBroker.getNewNDMaterial(dataID(1));
    if (the3D == 0) {
      opserr << "LogStrain2D::recvSelf - cannot create inner material classTag "
             << dataID(1) << endln;
      return -2;
    }
  }
  the3D->setDbTag(dataID(2));

  static Vector dataVec(1);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) return -3;
  lam_n = dataVec(0);
  lamTrial = lam_n;

  if (the3D->recvSelf(cTag, theChannel, theBroker) < 0) return -4;
  return 0;
}

// =========================================================================== //
//  Misc                                                                       //
// =========================================================================== //
void LogStrain2D::Print(OPS_Stream &s, int flag)
{
  const char *pt = (planeType == PlaneStress) ? "PlaneStress" : "PlaneStrain";
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"LogStrain2D\", ";
    s << "\"planeType\": \"" << pt << "\"";
    s << "}";
  } else {
    s << "LogStrain2D (2D Hencky log-strain finite-strain material, " << pt
      << "), tag: " << this->getTag() << endln;
  }
}

int LogStrain2D::setParameter(const char **argv, int argc, Parameter &param)
{
  return the3D->setParameter(argv, argc, param);
}

Response *LogStrain2D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (argc > 0) {
    if (strcmp(argv[0], "stress")  == 0 || strcmp(argv[0], "stresses") == 0 ||
        strcmp(argv[0], "Cauchy")  == 0 || strcmp(argv[0], "cauchy")   == 0)
      return new MaterialResponse(this, 1, sigma3);        // in-plane Cauchy stress (3)
    if (strcmp(argv[0], "strain")  == 0 || strcmp(argv[0], "strains")  == 0)
      return new MaterialResponse(this, 2, hencky3);       // in-plane Hencky strain (3)
    if (strcmp(argv[0], "tangent") == 0)
      return new MaterialResponse(this, 3, tangent3);      // 3x3 projection of c
  }
  return the3D->setResponse(argv, argc, output);
}

int LogStrain2D::getResponse(int responseID, Information &matInfo)
{
  switch (responseID) {
  case 1: return matInfo.setVector(sigma3);
  case 2: return matInfo.setVector(hencky3);
  case 3: return matInfo.setMatrix(tangent3);
  default: return the3D->getResponse(responseID, matInfo);
  }
}
