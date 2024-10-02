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

// Written: Jose A. Abell, UANDES
#ifndef ThermalBoundaryConditionTemperature_CPP
#define ThermalBoundaryConditionTemperature_CPP

                                                                        

#include <ThermalBoundaryConditionTemperature.h>
#include <Vector.h>

Vector ThermalBoundaryConditionTemperature::data(1);

ThermalBoundaryConditionTemperature::ThermalBoundaryConditionTemperature(int tag, int theElementTag, double factor)
  :ElementalLoad(tag, LOAD_TAG_ThermalBoundaryConditionTemperature, theElementTag), m_factor(factor)
{

}

ThermalBoundaryConditionTemperature::ThermalBoundaryConditionTemperature()
  :ElementalLoad(LOAD_TAG_ThermalBoundaryConditionTemperature), m_factor(1.0)
{

}

ThermalBoundaryConditionTemperature::~ThermalBoundaryConditionTemperature()
{

}

const Vector &
ThermalBoundaryConditionTemperature::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_ThermalBoundaryConditionTemperature;
  data(0) = m_factor * loadFactor;

  return data;
}

int 
ThermalBoundaryConditionTemperature::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
ThermalBoundaryConditionTemperature::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
ThermalBoundaryConditionTemperature::Print(OPS_Stream &s, int flag)
{
  s << "ThermalBoundaryConditionTemperature...";
  s << "  element acted on: " << eleTag << endln;;
}

#endif

