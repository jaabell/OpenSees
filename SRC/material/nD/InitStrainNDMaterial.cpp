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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: Jose Larenas (Universidad de los Andes, Chile)
// Created: Dec 2023

#include <stdlib.h>
#include <math.h>

#include <InitStrainNDMaterial.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>

#include <OPS_Globals.h>

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_InitStrainNDMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;
  NDMaterial *theOtherMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3) {
    opserr << "Want: nDMaterial InitStress tag? otherTag? sig0? <nDim?>" << endln;
  }

  int    iData[2];
  double dData[1];
  int    dim[1];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid nDMaterial InitStrainNDMaterial $tag $otherTag $nDim" << endln;
    return 0;
  }

  theOtherMaterial = OPS_getNDMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "Could not find material with tag: " << iData[1] << "nDMaterial InitStress $tag $otherTag $nDim $sig0" << endln;
    return 0;	
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: nDMaterial InitStress $tag $otherTag $nDim $sig0" << endln;
    return 0;	
  }

  if (numArgs == 4) {
    if (OPS_GetIntInput(&numData, dim) != 0) {
        return 0;
    }
  } else {
    dim[0] = 3;
  }

  Vector sig0(3*dim[0]-3);
  if (dim[0] == 3) {
    sig0(0) = dData[0];
    sig0(1) = dData[0];
    sig0(2) = dData[0];
  } else if (dim[0] == 2) {
    sig0(0) = dData[0];
    sig0(1) = dData[0];
  } else {
    opserr << "nDMaterial InitStress - Invalid number of dimensions: want 2 or 3" << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  if (numArgs == 4) {
    theMaterial = new InitStrainNDMaterial(iData[0], *theOtherMaterial, sig0, dim[0]);
  } else {
    theMaterial = new InitStrainNDMaterial(iData[0], *theOtherMaterial, sig0);
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type InitStrainNDMaterial\n";
    return 0;
  }

  return theMaterial;
}


InitStrainNDMaterial::InitStrainNDMaterial(int tag, NDMaterial &material, const Vector &epsini, int ndim)
  :NDMaterial(tag,ND_TAG_InitStrainNDMaterial), theMaterial(0),
   epsInit(3*ndim-3), sigInit(sigini)
{

  numDim = ndim;
  // get copy of the main material
  if (numDim == 2) {
    theMaterial = material.getCopy("PlaneStrain");
  } else if (numDim == 3) {
    theMaterial = material.getCopy("ThreeDimensional");
  } else {
    opserr << "nDMaterial InitStress - Invalid number of dimensions: want 2 or 3" << endln;
  }

  if (theMaterial == 0) {
    opserr <<  "InitStrainNDMaterial::InitStrainNDMaterial -- failed to get copy of material\n";
    exit(-1);
  }

  epsInit = epsini;
  theMaterial->setTrialStrain(epsInit);
  theMaterial->commitState();
}

InitStrainNDMaterial::InitStrainNDMaterial()
  :NDMaterial(0,ND_TAG_InitStrainNDMaterial), theMaterial(0),
   epsInit(6), sigInit(6)
{

}

InitStrainNDMaterial::~InitStrainNDMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
InitStrainNDMaterial::setTrialStrain(const Vector &strain) 
{
  return theMaterial->setTrialStrain(strain-epsInit);
}

int 
InitStrainNDMaterial::setTrialStrain(const Vector &strain, 
				     const Vector &strainRate)
{
  return theMaterial->setTrialStrain(strain-epsInit, strainRate);
}

int 
InitStrainNDMaterial::setTrialStrainIncr(const Vector &strain) 
{
  return theMaterial->setTrialStrainIncr(strain);
}

int 
InitStrainNDMaterial::setTrialStrainIncr(const Vector &strain, 
					 const Vector &strainRate)
{
  return theMaterial->setTrialStrainIncr(strain, strainRate);
}

const Vector &
InitStrainNDMaterial::getStress(void)
{
  return theMaterial->getStress();
}

const Matrix &
InitStrainNDMaterial::getTangent(void)
{
  return theMaterial->getTangent();  
}

const Matrix &
InitStrainNDMaterial::getInitialTangent(void)
{
  return theMaterial->getInitialTangent();  
}

const Vector & 
InitStrainNDMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

int 
InitStrainNDMaterial::commitState(void)
{	
  return theMaterial->commitState();
}

int 
InitStrainNDMaterial::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int 
InitStrainNDMaterial::revertToStart(void)
{
  int res = 0;
  res = theMaterial->revertToStart();
  res += theMaterial->setTrialStrain(epsInit);
  res += theMaterial->commitState();
  return res;
}

double 
InitStrainNDMaterial::getRho(void) 
{
	return theMaterial->getRho();
}

NDMaterial *
InitStrainNDMaterial::getCopy(void)
{
  InitStrainNDMaterial *theCopy = 
    new InitStrainNDMaterial(this->getTag(), *theMaterial, sigInit, numDim);
        
  return theCopy;
}

NDMaterial *
InitStrainNDMaterial::getCopy(const char *type)
{
  /*if (strcmp(type,"ThreeDimensional") == 0) {
    InitStrainNDMaterial *theCopy = 
      new InitStrainNDMaterial(this->getTag(), *theMaterial, sigInit);

    return theCopy;
  }

  return NDMaterial::getCopy(type);*/
  
  return this->getCopy();
}

const char*
InitStrainNDMaterial::getType(void) const
{
  return theMaterial->getType();
}

int 
InitStrainNDMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if ( matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "InitStrainNDMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(1);
  //dataVec(0) = epsInit;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStrainNDMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "InitStrainNDMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
InitStrainNDMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "InitStrainNDMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "InitStrainNDMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(1);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStrainNDMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  //epsInit = dataVec(0);
  
  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "InitStrainNDMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
InitStrainNDMaterial::Print(OPS_Stream &s, int flag)
{
  s << "InitStrainNDMaterial tag: " << this->getTag() << endln;
  s << "\tMaterial: " << theMaterial->getTag() << endln;
  s << "\tinitital strain: " << epsInit << endln;
}

int 
InitStrainMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"epsInit") == 0) {
    param.setValue(epsInit);
    return param.addObject(1, this);
  }

  // Otherwise, pass it on to the wrapped material
  if (theMaterial)
    return theMaterial->setParameter(argv, argc, param);
  else
    return -1;
}

int
InitStrainMaterial::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    this->epsInit = info.theDouble;
    if (theMaterial) {
      theMaterial->setTrialStrain(localStrain-epsInit);
      theMaterial->commitState();
    } else
      return -1;
  }

  return 0;
}

const Vector &
InitStrainNDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  return theMaterial->getStressSensitivity(gradIndex, conditional);
}

int
InitStrainNDMaterial::commitSensitivity(const Vector &depsdh, 
					int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(depsdh, gradIndex, numGrads);
}
