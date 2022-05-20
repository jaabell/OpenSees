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

#ifndef PluginUniaxialMaterial_h
#define PluginUniaxialMaterial_h

// Written: Massimo Petracca 
// Created: 02/2020
// Revision: A
//
// Description: This file contains the PluginUniaxialMaterial wrapper

#include <UniaxialMaterial.h>
#include <PluginFrameworkAPI.h>

class PluginMaterialDescriptor;

/**
The PluginUniaxialMaterial class is a wrapper for external materials
using the PluginFrameworkAPI
*/
class PluginUniaxialMaterial : public UniaxialMaterial
{
	friend void* OPS_PluginUniaxialMaterial(void);
	
	// constructor and destructor
public:
	PluginUniaxialMaterial();
	PluginUniaxialMaterial(PluginMaterialDescriptor* descr, PluginMaterialData* d);
	~PluginUniaxialMaterial();

	// from TaggedObject
public:
	void Print(OPS_Stream& s, int flag = 0);

	// from MovableObject
public:
	const char* getClassType() const;
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	int activateParameter(int parameterID);
	int setVariable(const char* variable, Information& info);
	int getVariable(const char* variable, Information& info);

	// from UniaxialMaterial
public:
	int setTrialStrain(double strain, double strainRate = 0);
	int setTrialStrain(double strain, double temperature, double strainRate);
	double getStrain();
	double getStrainRate();
	double getStress();
	double getTangent();
	double getInitialTangent();
	double getDampTangent();
	double getRho();
	int commitState();
	int revertToLastCommit();
	int revertToStart();
	UniaxialMaterial* getCopy();
	Response* setResponse(const char** argv, int argc, OPS_Stream& theOutputStream);
	int getResponse(int responseID, Information& info);
	int getResponseSensitivity(int responseID, int gradIndex, Information& info);
	bool hasFailed();
	double getStressSensitivity(int gradIndex, bool conditional);
	double getStrainSensitivity(int gradIndex);
	double getTangentSensitivity(int gradIndex);
	double getInitialTangentSensitivity(int gradIndex);
	double getDampTangentSensitivity(int gradIndex);
	double getRhoSensitivity(int gradIndex);
	int    commitSensitivity(double strainGradient, int gradIndex, int numGrads);
	double getEnergy();

private:
	PluginMaterialDescriptor* m_descriptor;
	PluginMaterialData* m_data;
	double m_lch;
	bool m_lch_calculated;
};

#endif // PluginUniaxialMaterial_h
