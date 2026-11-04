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
                                                                        
// $Revision: 1.5 $
// $Date: 2007-04-02 23:42:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/ArcLength1.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/ArcLength1.C
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for ArcLength1.
// ArcLength1 is an algorithmic class for performing a static analysis
// using the arc length scheme, that is within a load step the following
// constraint is enforced: dU^TdU + alpha^2*dLambda^2 = ArcLength1^2
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and ArcLength1 is a control parameter.
//
// What: "@(#) ArcLength1.C, revA"


#include <ArcLength1.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <stdlib.h>
#include <elementAPI.h>
#include <ContinuationLambda.h>
#include <classTags.h>
#include <string.h>

void* OPS_ArcLength1()
{
    double arcLength;
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING integrator ArcLength arcLength alpha \n";
	return 0;
    }

    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &arcLength) < 0) {
	opserr << "WARNING integrator ArcLength failed to read arc length\n";
	return 0;
    }
    double alpha = 1.0;
    bool haveAlpha = false;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* peek = OPS_GetString();
	bool isOption = (peek != 0 &&
			 (strcmp(peek,"-dt") == 0 || strcmp(peek,"-channel") == 0 ||
			  strcmp(peek,"-duration") == 0 || strcmp(peek,"-target") == 0 ||
			  strcmp(peek,"-signAngle") == 0));
	OPS_ResetCurrentInputArg(-1);
	if (isOption == false) {
	    if (OPS_GetDoubleInput(&numdata, &alpha) < 0) {
		opserr << "WARNING integrator ArcLength1 failed to read alpha\n";
		return 0;
	    }
	    haveAlpha = true;
	}
    }

    // continuous-time (continuation) options; absent => legacy behaviour
    int lambdaChannel = 0;
    double duration = 0.0, target = 0.0, dtFixed = 0.0;
    bool wantContinuation = false;
    bool wantSignAngle = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* arg = OPS_GetString();
	if (arg == 0)
	    break;
	if (strcmp(arg,"-signAngle") == 0) {
	    wantSignAngle = true;
	} else if (strcmp(arg,"-dt") == 0) {
	    numdata = 1;
	    if (OPS_GetNumRemainingInputArgs() < 1 ||
		OPS_GetDoubleInput(&numdata, &dtFixed) < 0) {
		opserr << "WARNING integrator ArcLength1 failed to read -dt\n";
		return 0;
	    }
	    wantContinuation = true;
	} else if (strcmp(arg,"-duration") == 0) {
	    numdata = 1;
	    if (OPS_GetNumRemainingInputArgs() < 1 ||
		OPS_GetDoubleInput(&numdata, &duration) < 0) {
		opserr << "WARNING integrator ArcLength1 failed to read -duration\n";
		return 0;
	    }
	    wantContinuation = true;
	} else if (strcmp(arg,"-target") == 0) {
	    numdata = 1;
	    if (OPS_GetNumRemainingInputArgs() < 1 ||
		OPS_GetDoubleInput(&numdata, &target) < 0) {
		opserr << "WARNING integrator ArcLength1 failed to read -target\n";
		return 0;
	    }
	    wantContinuation = true;
	} else if (strcmp(arg,"-channel") == 0) {
	    numdata = 1;
	    if (OPS_GetNumRemainingInputArgs() < 1 ||
		OPS_GetIntInput(&numdata, &lambdaChannel) < 0) {
		opserr << "WARNING integrator ArcLength1 failed to read -channel\n";
		return 0;
	    }
	    wantContinuation = true;
	}
    }

    ArcLength1* theIntegrator = haveAlpha ? new ArcLength1(arcLength,alpha)
					  : new ArcLength1(arcLength);

    if (wantSignAngle)
	theIntegrator->setSignFromAngle(true);

    if (wantContinuation &&
	theIntegrator->setContinuationTime(lambdaChannel, duration,
					   target, dtFixed) < 0) {
	opserr << "WARNING integrator ArcLength1 - continuation-time options "
	       << "rejected\n";
	delete theIntegrator;
	return 0;
    }

    return theIntegrator;
}

ArcLength1::ArcLength1(double arcLength, double alpha)
:StaticIntegrator(INTEGRATOR_TAGS_ArcLength1),
 arcLength2(arcLength*arcLength), alpha2(alpha*alpha),
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), 
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0),
 signLastDeltaLambdaStep(1),
 useContinuationTime(false), lambdaChannel(0), stageDuration(0.0),
 stageTarget(0.0), dtFixed(0.0), timeStep(0.0), committedLambda(0.0),
 signFromAngle(false), committedUstep(0)
{

}

ArcLength1::~ArcLength1()
{
    // delete any vector object created
    if (deltaUhat != 0)
	delete deltaUhat;
    if (deltaU != 0)
	delete deltaU;
    if (deltaUstep != 0)
	delete deltaUstep;
    if (committedUstep != 0)
	delete committedUstep;
    if (deltaUbar != 0)
	delete deltaUbar;
	if (phat != 0)
	delete phat;
}

double
ArcLength1::continuationDt(void) const
{
    // see ArcLength::continuationDt -- same progress measure
    if (stageTarget != 0.0)
	return stageDuration * sqrt(arcLength2) / fabs(stageTarget);
    return dtFixed;
}

int
ArcLength1::setContinuationTime(int chan, double duration,
				double target, double dt)
{
    if (!OPS_ContinuationLambda::inRange(chan)) {
	opserr << "ArcLength1::setContinuationTime() - lambda channel " << chan
	       << " out of range\n";
	return -1;
    }
    if (target == 0.0 && !(dt > 0.0)) {
	opserr << "ArcLength1::setContinuationTime() - need either -target (the "
	       << "total arc length of the stage) or a positive -dt\n";
	return -1;
    }
    if (target != 0.0 && fabs(target) < sqrt(arcLength2)) {
	opserr << "ArcLength1::setContinuationTime() - -target " << target
	       << " is smaller than one arc-length increment " << sqrt(arcLength2)
	       << "; the stage would overshoot its end time on the first step\n";
	return -1;
    }

    useContinuationTime = true;
    lambdaChannel = chan;
    stageDuration = duration;
    stageTarget = target;
    dtFixed = dt;

    if (OPS_ContinuationLambda::isValid(lambdaChannel) == false)
	OPS_ContinuationLambda::set(lambdaChannel, 0.0);
    committedLambda = OPS_ContinuationLambda::get(lambdaChannel);
    currentLambda = committedLambda;
    OPS_ContinuationLambda::setOwner(lambdaChannel, INTEGRATOR_TAGS_ArcLength1);

    return 0;
}

int
ArcLength1::setSignFromAngle(bool flag)
{
    signFromAngle = flag;
    return 0;
}

int
ArcLength1::setArcLength(double s)
{
    if (!(s > 0.0)) {
	opserr << "ArcLength1::setArcLength() - the arc length must be strictly "
	       << "positive, got " << s << "\n";
	return -1;
    }
    if (useContinuationTime && stageTarget != 0.0 && s > fabs(stageTarget)) {
	opserr << "ArcLength1::setArcLength() - " << s << " exceeds the stage "
	       << "total arc length " << stageTarget << "\n";
	return -1;
    }
    arcLength2 = s*s;
    return 0;
}

int
ArcLength1::commit(void)
{
    int res = this->StaticIntegrator::commit();
    if (res == 0 && useContinuationTime)
	committedLambda = currentLambda;
    if (res == 0 && signFromAngle && committedUstep != 0 && deltaUstep != 0)
	(*committedUstep) = (*deltaUstep);
    return res;
}

int
ArcLength1::revertToLastStep(void)
{
    if (useContinuationTime == false)
	return 0;

    currentLambda = committedLambda;
    deltaLambdaStep = 0.0;
    OPS_ContinuationLambda::set(lambdaChannel, committedLambda);

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	timeStep = theModel->getCurrentDomainTime();
	theModel->applyLoadDomain(timeStep);
	theModel->updateDomain();
    }

    return 0;
}

int
ArcLength1::newStep(void)
{
    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLength1::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    // get the current load factor (see DisplacementControl::newStep)
    if (useContinuationTime == false)
	currentLambda = theModel->getCurrentDomainTime();
    else
	currentLambda = committedLambda;

    if (deltaLambdaStep < 0)
	signLastDeltaLambdaStep = -1;
    else
	signLastDeltaLambdaStep = +1;

    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    // OPT-IN: sign from the angle with the previous converged step; see
    // ArcLength::newStep for the derivation
    if (signFromAngle && committedUstep != 0) {
	double d = (*committedUstep) ^ dUhat;
	if (d < 0.0)
	    signLastDeltaLambdaStep = -1;
	else if (d > 0.0)
	    signLastDeltaLambdaStep = +1;
    }

    // determine delta lambda(1) == dlambda
    double dLambda = sqrt(arcLength2/((dUhat^dUhat)+alpha2));
    dLambda *= signLastDeltaLambdaStep; // base sign of load change
                                        // on what was happening last step
    deltaLambdaStep = dLambda;
    currentLambda += dLambda;

    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);
    if (useContinuationTime == false) {
	theModel->applyLoadDomain(currentLambda);
    } else {
	double dt = this->continuationDt();
	if (!(dt > 0.0)) {
	    opserr << "WARNING ArcLength1::newStep() - continuation dt is " << dt
		   << ", must be strictly positive\n";
	    return -1;
	}
	timeStep = theModel->getCurrentDomainTime() + dt;
	OPS_ContinuationLambda::set(lambdaChannel, currentLambda);
	theModel->applyLoadDomain(timeStep);
    }
    theModel->updateDomain();

    return 0;
}

int
ArcLength1::update(const Vector &dU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLength1::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change

    // determine dUhat    
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();    

    // determine delta lambda(i)
    double a = (*deltaUstep)^(*deltaUbar);
    double b = (*deltaUstep)^(*deltaUhat) + alpha2*deltaLambdaStep;
    if (b == 0) {
      opserr << "ArcLength1::update() - zero denominator,";
      opserr << " alpha was set to 0.0 and zero reference load\n";
      return -1;
    }
    double dLambda = -a/b;

    // determine delta U(i)
    (*deltaU) = (*deltaUbar);    
    deltaU->addVector(1.0, *deltaUhat,dLambda);
    
    // update dU and dlambda
    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda += dLambda;

    // update the model
    theModel->incrDisp(*deltaU);
    if (useContinuationTime == false) {
	theModel->applyLoadDomain(currentLambda);
    } else {
	// lambda changes every iteration; the TIME must not
	OPS_ContinuationLambda::set(lambdaChannel, currentLambda);
	theModel->applyLoadDomain(timeStep);
    }
    theModel->updateDomain();

    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    return 0;
}



int 
ArcLength1::domainChanged(void)
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLength1::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    
    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
	if (deltaUhat != 0)
	    delete deltaUhat;   // delete the old
	deltaUhat = new Vector(size);
	if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
	if (deltaUbar != 0)
	    delete deltaUbar;   // delete the old
	deltaUbar = new Vector(size);
	if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " deltaUbar Vector of size " << size << endln;
	    exit(-1);
	}
    }

    
    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
	if (deltaU != 0)
	    delete deltaU;   // delete the old
	deltaU = new Vector(size);
	if (deltaU == 0 || deltaU->Size() != size) { // check got it
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	if (deltaUstep != 0)
	    delete deltaUstep;  
	deltaUstep = new Vector(size);
	if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (phat == 0 || phat->Size() != size) { 
	if (phat != 0)
	    delete phat;  
	phat = new Vector(size);
	if (phat == 0 || phat->Size() != size) { 
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " phat Vector of size " << size << endln;
	    exit(-1);
	}
    }    

    // only allocated when the angle sign rule is on; see ArcLength
    if (signFromAngle &&
	(committedUstep == 0 || committedUstep->Size() != size)) {
	if (committedUstep != 0)
	    delete committedUstep;
	committedUstep = new Vector(size);
	if (committedUstep == 0 || committedUstep->Size() != size) {
	    opserr << "FATAL ArcLength1::domainChanged() - ran out of memory for";
	    opserr << " committedUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

    // now we have to determine phat
    // do this by incrementing lambda by 1, applying load
    // and getting phat from unbalance.
    if (useContinuationTime) {
	// probe LAMBDA at a frozen time; see ArcLength::domainChanged
	double t = theModel->getCurrentDomainTime();
	double lam = OPS_ContinuationLambda::get(lambdaChannel);

	OPS_ContinuationLambda::set(lambdaChannel, lam + 1.0);
	theModel->applyLoadDomain(t);
	this->formUnbalance();
	(*phat) = theLinSOE->getB();

	OPS_ContinuationLambda::set(lambdaChannel, lam);
	theModel->applyLoadDomain(t);
	currentLambda = lam;
	return 0;
    }

    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);
    this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
    (*phat) = theLinSOE->getB();
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);    
    
    return 0;
}

int
ArcLength1::sendSelf(int cTag,
		    Channel &theChannel)
{
  Vector data(12);
  data(0) = arcLength2;
  data(1) = alpha2;
  data(2) = deltaLambdaStep;
  data(3) = currentLambda;
  data(4)  = signLastDeltaLambdaStep;
  data(5) = useContinuationTime ? 1.0 : 0.0;
  data(6) = lambdaChannel;
  data(7) = stageDuration;
  data(8) = stageTarget;
  data(9) = dtFixed;
  data(10) = committedLambda;
  data(11) = signFromAngle ? 1.0 : 0.0;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLength1::sendSelf() - failed to send the data\n";
      return -1;
  }
  return 0;
}


int
ArcLength1::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(12);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLength1::sendSelf() - failed to send the data\n";
      return -1;
  }

  // set the data
  arcLength2 = data(0);
  alpha2 = data(1);
  deltaLambdaStep = data(2);
  currentLambda = data(3);
  signLastDeltaLambdaStep = data(4);
  useContinuationTime = (data(5) != 0.0);
  lambdaChannel = (int)data(6);
  stageDuration = data(7);
  stageTarget = data(8);
  dtFixed = data(9);
  committedLambda = data(10);
  signFromAngle = (data(11) != 0.0);

  if (useContinuationTime) {
      currentLambda = committedLambda;
      OPS_ContinuationLambda::set(lambdaChannel, committedLambda);
      OPS_ContinuationLambda::setOwner(lambdaChannel, INTEGRATOR_TAGS_ArcLength1);
  }
  return 0;
}

void
ArcLength1::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double cLambda = theModel->getCurrentDomainTime();
	s << "\t ArcLength1 - currentLambda: " << cLambda;
	s << "  ArcLength1: " << sqrt(arcLength2) <<  "  alpha: ";
	s << sqrt(alpha2) << endln;
    } else 
	s << "\t ArcLength1 - no associated AnalysisModel\n";
}








