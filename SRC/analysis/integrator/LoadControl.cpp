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
                                                                        
// $Revision: 1.6 $
// $Date: 2007-04-02 23:42:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/LoadControl.cpp,v $
                                                                        
                                                                        
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for LoadControl.
// LoadControl is an algorithmic class for performing a static analysis
// using a load control integration scheme.
//
// What: "@(#) LoadControl.h, revA"



#include <LoadControl.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Node.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Domain.h>
#include<Parameter.h>
#include<ParameterIter.h>
#include<EquiSolnAlgo.h>
#include <elementAPI.h>
#include <iostream>
#include <string.h>

void* OPS_LoadControlIntegrator()
{
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"LoadControl - insufficient arguments\n";
	return 0;
    }

    double lambda;
    int numData = 1;
    if(OPS_GetDoubleInput(&numData,&lambda) < 0) {
	opserr<<"WARNING LoadControl - failed to read double lambda\n";
	return 0;
    }

    // -tangentPredictor is accepted right after lambda or after the
    // optional numIter/minLambda/maxLambda triple
    bool tangentPredictor = false;
    char optBuffer[64];
    if(OPS_GetNumRemainingInputArgs() > 0) {
	const char* opt = OPS_GetStringFromAll(optBuffer, 64);
	if(opt != 0 && strcmp(opt,"-tangentPredictor") == 0)
	    tangentPredictor = true;
	else
	    OPS_ResetCurrentInputArg(-1);
    }

    int numIter = 1;
    double mLambda[2] = {lambda,lambda};
    if(OPS_GetNumRemainingInputArgs() > 2) {
	if(OPS_GetIntInput(&numData,&numIter) < 0) {
	    opserr<<"WARNING LoadControl - failed to read int numIter\n";
	    return 0;
	}
	numData = 2;
	if(OPS_GetDoubleInput(&numData,&mLambda[0]) < 0) {
	    opserr<<"WARNING LoadControl - failed to read double min and max\n";
	    return 0;
	}
    }

    if(!tangentPredictor && OPS_GetNumRemainingInputArgs() > 0) {
	const char* opt = OPS_GetStringFromAll(optBuffer, 64);
	if(opt != 0 && strcmp(opt,"-tangentPredictor") == 0)
	    tangentPredictor = true;
	else
	    OPS_ResetCurrentInputArg(-1);
    }

    LoadControl *theIntegrator = new LoadControl(lambda,numIter,mLambda[0],mLambda[1]);
    if(tangentPredictor)
	theIntegrator->setTangentPredictor(true);

    return theIntegrator;
}

LoadControl::LoadControl(double dLambda, int numIncr, double min, double max, int classtag)
    : StaticIntegrator(classtag),
 deltaLambda(dLambda), 
 specNumIncrStep(numIncr), numIncrLastStep(numIncr),
 dLambdaMin(min), dLambdaMax(max), gradNumber(0), sensitivityFlag(0),
 tangentPredictReq(false), tangentPredict(false), spNotYetEnforced(false),
 fallbackNoted(false), stepLambda(0.0)
{
  // to avoid divide-by-zero error on first update() ensure numIncr != 0
  if (numIncr == 0) {
    opserr << "WARNING LoadControl::LoadControl() - numIncr set to 0, 1 assumed\n";
    specNumIncrStep = 1.0;
    numIncrLastStep = 1.0;
  }
}


LoadControl::~LoadControl()
{
    
}

int 
LoadControl::newStep(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();    
    if (theModel == 0) {
	opserr << "LoadControl::newStep() - no associated AnalysisModel\n";
	return -1;
    }

    // determine delta lambda for this step based on dLambda and #iter of last step
    double factor = specNumIncrStep/numIncrLastStep;
    deltaLambda *=factor;

    if (deltaLambda < dLambdaMin)
      deltaLambda = dLambdaMin;
    else if (deltaLambda > dLambdaMax)
      deltaLambda = dLambdaMax;
    
    double currentLambda = theModel->getCurrentDomainTime();

    currentLambda += deltaLambda;

    if (tangentPredict) {
      // Apply the domain loads only: this assigns every SP_Constraint its
      // value for the new load factor, but the constraint handler does not
      // enforce them, so the elements stay at the committed state. The first
      // iteration then takes the prescribed-motion forcing from -K*du_p
      // (formUnbalance) and the SPs are enforced in update().
      Domain *theModelDomain = theModel->getDomainPtr();
      if (theModelDomain == 0) {
	opserr << "LoadControl::newStep() - no Domain linked to the AnalysisModel\n";
	return -1;
      }
      theModelDomain->applyLoad(currentLambda);
      stepLambda = currentLambda;
      spNotYetEnforced = true;
      numIncrLastStep = 0;
      return 0;
    }

    theModel->applyLoadDomain(currentLambda);

    numIncrLastStep = 0;
   
     return 0;
}
    
int
LoadControl::update(const Vector &deltaU)
{
    AnalysisModel *myModel = this->getAnalysisModel();
    LinearSOE *theSOE = this->getLinearSOE();
    if (myModel == 0 || theSOE == 0) {
	opserr << "WARNING LoadControl::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    myModel->incrDisp(deltaU);    

    // with -tangentPredictor, deltaU of the first iteration is the response
    // of the free dofs to the prescribed increment; enforce the SPs now
    if (tangentPredict && spNotYetEnforced) {
      myModel->applyLoadDomain(stepLambda);
      spNotYetEnforced = false;
    }

    if (myModel->updateDomain() < 0) {
      opserr << "LoadControl::update - model failed to update for new dU\n";
      return -1;
    }

    // Set deltaU for the convergence test
    theSOE->setX(deltaU);

    numIncrLastStep++;

    return 0;
}


int
LoadControl::setDeltaLambda(double newValue)
{
  // we set the #incr at last step = #incr so get newValue incr
  numIncrLastStep = specNumIncrStep;
  deltaLambda = newValue;
  return 0;
}


void
LoadControl::setTangentPredictor(bool onOff)
{
  tangentPredictReq = onOff;
  tangentPredict = onOff;
  spNotYetEnforced = false;
}


int
LoadControl::domainChanged(void)
{
  // a step left half open is abandoned, and a changed domain may now hold
  // non-homogeneous SPs that the fallback in formUnbalance() had not seen
  spNotYetEnforced = false;
  tangentPredict = tangentPredictReq;
  return this->StaticIntegrator::domainChanged();
}


// With -tangentPredictor, between newStep() and the first update() the
// unbalance is P - R(committed) - K*du_p, where du_p are the prescribed
// increments of the non-homogeneous SPs that have been set but not yet
// enforced. If no element supplies such a term (no non-homogeneous SP, or a
// constraint handler other than Transformation) nothing would drive the
// prescribed dofs in the first iteration and a displacement-based test could
// accept an unmoved step, so the SPs are enforced immediately and the
// predictor is switched off until the domain changes.
int
LoadControl::formUnbalance(void)
{
  int res = this->IncrementalIntegrator::formUnbalance();

  if (res < 0 || !tangentPredict || !spNotYetEnforced)
    return res;

  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE *theSOE = this->getLinearSOE();
  if (theModel == 0 || theSOE == 0) {
    opserr << "WARNING LoadControl::formUnbalance() ";
    opserr << "No AnalysisModel or LinearSOE has been set\n";
    return -1;
  }

  FE_EleIter &theEles = theModel->getFEs();
  FE_Element *elePtr;
  int numContributed = 0;
  while ((elePtr = theEles()) != 0) {
    const Vector *spForce = elePtr->getSPTangentForce(this);
    if (spForce != 0) {
      if (theSOE->addB(*spForce, elePtr->getID(), -1.0) < 0) {
	opserr << "WARNING LoadControl::formUnbalance() ";
	opserr << "failed to add the prescribed-motion force to B\n";
	return -1;
      }
      numContributed++;
    }
  }

  if (numContributed > 0)
    return 0;

  if (!fallbackNoted) {
    opserr << "LoadControl -tangentPredictor: no element supplied a K*du_p term ";
    opserr << "(no non-homogeneous sp under constraints Transformation); ";
    opserr << "using the standard predictor\n";
    fallbackNoted = true;
  }
  tangentPredict = false;

  // same state as the standard newStep(): loads applied and SPs enforced
  theModel->applyLoadDomain(stepLambda);
  spNotYetEnforced = false;

  return this->IncrementalIntegrator::formUnbalance();
}


int
LoadControl::sendSelf(int cTag,
		      Channel &theChannel)
{
  Vector data(5);
  data(0) = deltaLambda;
  data(1) = specNumIncrStep;
  data(2) = numIncrLastStep;
  data(3) = dLambdaMin;
  data(4) = dLambdaMax;
  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "LoadControl::sendSelf() - failed to send the Vector\n";
      return -1;
  }
  return 0;
}


int
LoadControl::recvSelf(int cTag,
		      Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(5);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "LoadControl::sendSelf() - failed to send the Vector\n";
      deltaLambda = 0;
      return -1;
  }      
  deltaLambda = data(0);
  specNumIncrStep = data(1);
  numIncrLastStep = data(2);
  dLambdaMin = data(3);
  dLambdaMax = data(4);
  return 0;
}



void
LoadControl::Print(OPS_Stream &s, int flag)
{ 
     AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double currentLambda = theModel->getCurrentDomainTime();
	s << "\t LoadControl - currentLambda: " << currentLambda;
	s << "  deltaLambda: " << deltaLambda;
	if (tangentPredictReq)
	  s << "  -tangentPredictor";
	s << endln;
    } else 
	s << "\t LoadControl - no associated AnalysisModel\n";
    
}

int
LoadControl::formEleResidual(FE_Element* theEle)
{
    if(sensitivityFlag == 0) {  // no sensitivity
	this->StaticIntegrator::formEleResidual(theEle);
    } else {
	theEle->zeroResidual();
	theEle->addResistingForceSensitivity(gradNumber);
    }
    return 0;
}

int
LoadControl::formIndependentSensitivityRHS()
{
    return 0;
}

int
LoadControl::formSensitivityRHS(int passedGradNumber)
{
    sensitivityFlag = 1;

    // Set a couple of data members
    gradNumber = passedGradNumber;

    // get model
    AnalysisModel* theAnalysisModel = this->getAnalysisModel();
    LinearSOE* theSOE = this->getLinearSOE();

    // Loop through elements
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();   
    while((elePtr = theEles()) != 0) {
      theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
    }

    // Loop through the loadPatterns and add the dPext/dh contributions
    static Vector oneDimVectorWithOne(1);
    oneDimVectorWithOne(0) = 1.0;
    static ID oneDimID(1);

    Node *aNode;
    DOF_Group *aDofGroup;
    int nodeNumber, dofNumber, relevantID, i, sizeRandomLoads, numRandomLoads;
    LoadPattern *loadPatternPtr;

    Domain *theDomain = theAnalysisModel->getDomainPtr();
    LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    while((loadPatternPtr = thePatterns()) != 0) {
	const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
	sizeRandomLoads = randomLoads.Size();
	if (sizeRandomLoads == 1) {
	    // No random loads in this load pattern
	}
	else {
	    // Random loads: add contributions to the 'B' vector
	    numRandomLoads = (int)(sizeRandomLoads/2);
	    for (i=0; i<numRandomLoads*2; i=i+2) {
		nodeNumber = (int)randomLoads(i);
		dofNumber = (int)randomLoads(i+1);
		aNode = theDomain->getNode(nodeNumber);
		aDofGroup = aNode->getDOF_GroupPtr();
		const ID &anID = aDofGroup->getID();
		relevantID = anID(dofNumber-1);
		oneDimID(0) = relevantID;
		theSOE->addB(oneDimVectorWithOne, oneDimID);
	    }
	}
    }

    // reset sensitivity flag
    sensitivityFlag = 0;

    return 0;
}

int
LoadControl::saveSensitivity(const Vector &v, int gradNum, int numGrads)
{
    // get model
    AnalysisModel* theAnalysisModel = this->getAnalysisModel();
    
    DOF_GrpIter &theDOFGrps = theAnalysisModel->getDOFs();
    DOF_Group 	*dofPtr;
    
    while ( (dofPtr = theDOFGrps() ) != 0)  {
//	dofPtr->saveSensitivity(v,0,0,gradNum,numGrads);
	dofPtr->saveDispSensitivity(v,gradNum,numGrads);
	
    }
    
    return 0;
}

int 
LoadControl::commitSensitivity(int gradNum, int numGrads)
{
      // get model
    AnalysisModel* theAnalysisModel = this->getAnalysisModel();
    
    // Loop through the FE_Elements and set unconditional sensitivities
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != 0) {
	elePtr->commitSensitivity(gradNum, numGrads);
    }
    
    return 0;
}



// false for LC and true for DC
   bool 
LoadControl::computeSensitivityAtEachIteration()
{

return false;
}





int 
LoadControl::computeSensitivities(void)
{
//  opserr<<" computeSensitivity::start"<<endln; 
    LinearSOE *theSOE = this->getLinearSOE();

    /*
  if (theAlgorithm == 0) {
    opserr << "ERROR the FE algorithm must be defined before ";
    opserr << "the sensitivity algorithm\n";
    return -1;
  }
*/
/*
   // Get pointer to the system of equations (SOE)
	LinearSOE *theSOE = theAlgorithm->getLinearSOEptr();
	if (theSOE == 0) {
	  opserr << "ERROR the FE linearSOE must be defined before ";
	  opserr << "the sensitivity algorithm\n";
	  return -1;
	}

	// Get pointer to incremental integrator
	IncrementalIntegrator *theIncInt = theAlgorithm->getIncrementalIntegratorPtr();
//	IncrementalIntegrator *theIncIntSens=theAlgorithm->getIncrementalIntegratorPtr();//Abbas
	if (theIncInt == 0 ) {
	  opserr << "ERROR the FE integrator must be defined before ";
	  opserr << "the sensitivity algorithm\n";
	  return -1;
	}

	// Form current tangent at converged state
	// (would be nice with an if-statement here in case
	// the current tangent is already formed)
	if (this->formTangent(CURRENT_TANGENT) < 0){
		opserr << "WARNING SensitivityAlgorithm::computeGradients() -";
		opserr << "the Integrator failed in formTangent()\n";
		return -1;
	}
	*/
	// Zero out the old right-hand side of the SOE
	theSOE->zeroB();

	// Form the part of the RHS which are independent of parameter
	this->formIndependentSensitivityRHS();
	AnalysisModel *theModel = this->getAnalysisModel();  
	Domain *theDomain=theModel->getDomainPtr();
	ParameterIter &paramIter = theDomain->getParameters();
	
	Parameter *theParam;
	// De-activate all parameters
	while ((theParam = paramIter()) != 0)
	  theParam->activate(false);
	
	// Now, compute sensitivity wrt each parameter
	int numGrads = theDomain->getNumParameters();
	paramIter = theDomain->getParameters();
	
	while ((theParam = paramIter()) != 0) {

	  // Activate this parameter
	  theParam->activate(true);

	  // Zero the RHS vector
	  theSOE->zeroB();

	  // Get the grad index for this parameter
	  int gradIndex = theParam->getGradIndex();

	  // Form the RHS
	  this->formSensitivityRHS(gradIndex);
         
	  // Solve for displacement sensitivity
	 
	  theSOE->solve();
	  // Save sensitivity to nodes
	  this->saveSensitivity( theSOE->getX(), gradIndex, numGrads );
	 


	  // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
	  this->commitSensitivity(gradIndex, numGrads);
	  
	  // De-activate this parameter for next sensitivity calc
	  theParam->activate(false);
	//  opserr<<"LoadControl::..........ComputeSensitivities. end"<<endln;
	}

	return 0;
}

