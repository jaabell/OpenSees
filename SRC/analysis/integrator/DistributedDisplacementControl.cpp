
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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-04-11 23:37:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/DistributedDisplacementControl.cpp,v $
                                                                        
// Written: fmk 
// Created: 07/98
//
// Description: This file contains the class definition for DistributedDisplacementControl.
// DistributedDisplacementControl is an algorithmic class for performing a static analysis
// using the arc length scheme, that is within a load step the following
// constraint is enforced: dU^TdU + alpha^2*dLambda^2 = DistributedDisplacementControl^2
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and DistributedDisplacementControl is a control parameter.
//
// What: "@(#) DistributedDisplacementControl.C, revA"

#include <DistributedDisplacementControl.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <ID.h>
#include <stdlib.h>
#include <ContinuationLambda.h>
#include <classTags.h>

DistributedDisplacementControl::DistributedDisplacementControl(int node, int dof, 
							       double increment, 
							       int numIncr,
							       double min, double max)
:StaticIntegrator(INTEGRATOR_TAGS_DistributedDisplacementControl),
 processID(0), theChannels(0), numChannels(0), 
 theNode(node), theDof(dof), theIncrement(increment),
 theDofID(0), 
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), 
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0),
 specNumIncrStep(numIncr), numIncrLastStep(numIncr),
 minIncrement(min), maxIncrement(max),
 useContinuationTime(false), lambdaChannel(0), stageDuration(0.0),
 stageTarget(0.0), dtFixed(0.0), timeStep(0.0), committedLambda(0.0)
{
  // to avoid divide-by-zero error on first update() ensure numIncr != 0
  if (numIncr == 0) {
    opserr << "WARNING DistributedDisplacementControl::DistributedDisplacementControl() -";
    opserr << " numIncr set to 0, 1 assumed\n";
    specNumIncrStep = 1.0;
    numIncrLastStep = 1.0;
  }
}



DistributedDisplacementControl::DistributedDisplacementControl()
:StaticIntegrator(INTEGRATOR_TAGS_DistributedDisplacementControl),
 processID(0), theChannels(0), numChannels(0), 
 theNode(0), theDof(0), theIncrement(0),
 theDofID(0), deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), 
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0),
 specNumIncrStep(0), numIncrLastStep(0),
 minIncrement(0), maxIncrement(0),
 useContinuationTime(false), lambdaChannel(0), stageDuration(0.0),
 stageTarget(0.0), dtFixed(0.0), timeStep(0.0), committedLambda(0.0)
{

}

DistributedDisplacementControl::~DistributedDisplacementControl()
{
    // delete any vector object created
    if (deltaUhat != 0)
	delete deltaUhat;
    if (deltaU != 0)
	delete deltaU;
    if (deltaUstep != 0)
	delete deltaUstep;
    if (deltaUbar != 0)
	delete deltaUbar;
    if (phat != 0)
	delete phat;
    if (theChannels != 0)
      delete [] theChannels;
}

int
DistributedDisplacementControl::newStep(void)
{
    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING DistributedDisplacementControl::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }


    // determine increment for this iteration
    double factor = specNumIncrStep/numIncrLastStep;
    theIncrement *=factor;

    if (theIncrement < minIncrement)
      theIncrement = minIncrement;
    else if (theIncrement > maxIncrement)
      theIncrement = maxIncrement;

    // get the current load factor.
    // legacy: lambda IS the domain pseudo-time. continuation mode: lambda is an
    // internal accumulator, rolled back by revertToLastStep() -- the legacy path
    // is self-healing on a failed step only because it re-reads the domain time,
    // which Domain::revertToLastCommit() has already restored.
    if (useContinuationTime == false)
      currentLambda = theModel->getCurrentDomainTime();
    else
      currentLambda = committedLambda;

    // determine dUhat
    this->formTangent();



    if (processID == 0)
      theLinSOE->setB(*phat);
    else
      theLinSOE->zeroB();

    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    double dUahat = dUhat(theDofID);
    if (dUahat == 0.0) {
	opserr << "WARNING DistributedDisplacementControl::newStep() ";
	opserr << "dUahat is zero -- zero reference displacement at control node DOF\n";
	return -1;
    }


    // determine delta lambda(1) == dlambda    
    double dLambda = theIncrement/dUahat;

    deltaLambdaStep = dLambda;
    currentLambda += dLambda;
    //    opserr << "DistributedDisplacementControl:newStep " << dUahat  << " " << theDofID << endln;
    //opserr << "DistributedDisplacementControl:newStep " << *phat << endln;
    //   opserr << "DistributedDisplacementControl::newStep() : " << deltaLambdaStep << endln;
    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);
    if (useContinuationTime == false) {
      theModel->applyLoadDomain(currentLambda);
    } else {
      // advance the timeline once per step and freeze it for the iterations
      double dt = this->continuationDt();
      if (!(dt > 0.0)) {
	opserr << "WARNING DistributedDisplacementControl::newStep() - continuation "
	       << "dt is " << dt << ", must be strictly positive; the controlled "
	       << "increment has collapsed to zero or the sign of -target disagrees "
	       << "with the displacement increment\n";
	return -1;
      }
      // dt is a function of stageDuration, stageTarget and theIncrement only --
      // all three are script constants or driven by the iteration count, which
      // is collective. So every rank computes the SAME dt with no communication.
      // At this point in the step the domain time still equals the committed
      // time (nothing has applied a load yet, and a failed step went through
      // Domain::revertToLastCommit which restores it).
      timeStep = theModel->getCurrentDomainTime() + dt;
      OPS_ContinuationLambda::set(lambdaChannel, currentLambda);
      theModel->applyLoadDomain(timeStep);
    }
    if (theModel->updateDomain() < 0) {
      opserr << "DistributedDisplacementControl::newStep - model failed to update for new dU\n";
      return -1;
    }

    numIncrLastStep = 0;

    return 0;
}

int
DistributedDisplacementControl::update(const Vector &dU)
{

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING DistributedDisplacementControl::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change
    double dUabar = (*deltaUbar)(theDofID);
    
    // determine dUhat    
    if (processID == 0)
      theLinSOE->setB(*phat);
    else
      theLinSOE->zeroB();

    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();

    double dUahat = (*deltaUhat)(theDofID);

    if (dUahat == 0.0) {
	opserr << "WARNING DistributedDisplacementControl::update() ";
	opserr << "dUahat is zero -- zero reference displacement at control node DOF\n";
	return -1;
    }

    // determine delta lambda(1) == dlambda    
    double dLambda = -dUabar/dUahat;
    
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
      // lambda changes every iteration; the TIME must not. Freezing it also
      // stops ops_Dt from oscillating under rate-dependent materials.
      OPS_ContinuationLambda::set(lambdaChannel, currentLambda);
      theModel->applyLoadDomain(timeStep);
    }
    if (theModel->updateDomain() < 0) {
      opserr << "DistributedDisplacementControl::update - model failed to update for new dU\n";
      return -1;
    }
	
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    numIncrLastStep++;

    return 0;
}



int 
DistributedDisplacementControl::domainChanged(void)
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING DistributedDisplacementControl::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    
    
    //int size = theModel->getNumEqn(); // ask model in case N+1 space
    const Vector &b = theLinSOE->getB(); // ask model in case N+1 space
    int size = b.Size();

    // first we determine the id of the nodal dof
    Domain *theDomain = theModel->getDomainPtr();
    if (theDomain == 0) {
      opserr << "BUG WARNING DistributedDisplacementControl::domainChanged() - no Domain associated!!";
      return -1;
    }

    theDofID = -1;
    Node *theNodePtr = theDomain->getNode(theNode);
    if (theNodePtr != 0) {
      DOF_Group *theGroup = theNodePtr->getDOF_GroupPtr();
      if (theGroup == 0) {
	opserr << "BUG DistributedDisplacementControl::domainChanged() - no DOF_Group associated with the node!!\n";
	return -1;
      }
      const ID &theID = theGroup->getID();
      if (theDof < 0 || theDof >= theID.Size()) {
	opserr << "DistributedDisplacementControl::domainChanged() - not a valid dof " << theDof << endln;
	return -1;
      }
      theDofID = theID(theDof);
      if (theDofID < 0) {
	opserr << "DistributedDisplacementControl::domainChanged() - constrained dof not a valid a dof\n";;
	return -1;
      }
    }

    static ID data(1);    

    if (processID != 0) {
      Channel *theChannel = theChannels[0];
      data(0) = theDofID;

      theChannel->sendID(0, 0, data);
      theChannel->recvID(0, 0, data);
      theDofID = data(0);
    } 

    // if main domain, collect all theDofID if not -1 then this is the value & send value out
    else {
      
      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->recvID(0, 0, data);
	if (data(0) != -1)
	  theDofID = data(0);
      }
      
      for (int k=0; k<numChannels; k++) {
	Channel *theChannel = theChannels[k];
	data(0) = theDofID;
	theChannel->sendID(0, 0, data);
      }
    }

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
	if (deltaUhat != 0)
	    delete deltaUhat;   // delete the old
	deltaUhat = new Vector(size);
	if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	    opserr << "FATAL DistributedDisplacementControl::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
	if (deltaUbar != 0)
	    delete deltaUbar;   // delete the old
	deltaUbar = new Vector(size);
	if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	    opserr << "FATAL DistributedDisplacementControl::domainChanged() - ran out of memory for";
	    opserr << " deltaUbar Vector of size " << size << endln;
	    exit(-1);
	}
    }
    
    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
	if (deltaU != 0)
	    delete deltaU;   // delete the old
	deltaU = new Vector(size);
	if (deltaU == 0 || deltaU->Size() != size) { // check got it
	    opserr << "FATAL DistributedDisplacementControl::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	if (deltaUstep != 0)
	    delete deltaUstep;  
	deltaUstep = new Vector(size);
	if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	    opserr << "FATAL DistributedDisplacementControl::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (phat == 0 || phat->Size() != size) { 
	if (phat != 0)
	    delete phat;  
	phat = new Vector(size);
	if (phat == 0 || phat->Size() != size) { 
	    opserr << "FATAL DistributedDisplacementControl::domainChanged() - ran out of memory for";
	    opserr << " phat Vector of size " << size << endln;
	    exit(-1);
	}
    }    

    // now we have to determine phat
    // do this by incrementing lambda by 1, applying load
    // and getting phat from unbalance.
    if (useContinuationTime == false) {
      currentLambda = theModel->getCurrentDomainTime();

      currentLambda += 1.0;
      theModel->applyLoadDomain(currentLambda);

      this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
      (*phat) = theLinSOE->getB();

      currentLambda -= 1.0;
      theModel->setCurrentDomainTime(currentLambda);
    } else {
      // Probe LAMBDA, not the time. phat is the derivative of the external load
      // with respect to lambda, so bumping lambda by one at a FROZEN time gives
      // it exactly: every pattern driven by a ContinuationTimeSeries contributes
      // its reference load, and every other pattern contributes zero because its
      // amplitude is evaluated at an unchanged time. Probing the time here (the
      // legacy path) would move no amplitude at all in continuation mode and
      // leave phat == 0.
      double t = theModel->getCurrentDomainTime();
      double lam = OPS_ContinuationLambda::get(lambdaChannel);

      OPS_ContinuationLambda::set(lambdaChannel, lam + 1.0);
      theModel->applyLoadDomain(t);
      this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
      (*phat) = theLinSOE->getB();

      OPS_ContinuationLambda::set(lambdaChannel, lam);
      theModel->applyLoadDomain(t);
      currentLambda = lam;
    }

    // check there is a reference load
    int haveLoad = 0;
    for (int i=0; i<size; i++)
      if ( (*phat)(i) != 0.0 ) {
	haveLoad = 1;
	i = size;
      }

    if (haveLoad == 0) {
      opserr << "WARNING DistributedDisplacementControl::domainChanged() - zero reference load";
      if (useContinuationTime)
	opserr << "\n      in continuation-time mode the reference load pattern must "
	       << "use a ContinuationTimeSeries on lambda channel " << lambdaChannel
	       << "\n";
      return -1;
    }


    if (theDofID == -1) {
      opserr << "DistributedDisplacementControl::setSize() - failed to find valid dof - are the node tag and dof values correct?\n";
      return -1;
    }

    return 0;
}

double
DistributedDisplacementControl::continuationDt(void) const
{
  // progress-normalised: the sum over the stage is exactly stageDuration,
  // because sum|theIncrement| == |stageTarget| by definition of the stage and
  // the controlled displacement cannot reverse sign (that is precisely what
  // lets a continuation method cross a limit point).
  if (stageTarget != 0.0)
    return stageDuration * fabs(theIncrement) / fabs(stageTarget);

  // fixed dt: simple, but drifts off stageDuration when the adaptive stepping
  // changes the number of steps actually taken
  return dtFixed;
}

int
DistributedDisplacementControl::setContinuationTime(int chan, double duration,
						    double target, double dt)
{
  if (!OPS_ContinuationLambda::inRange(chan)) {
    opserr << "DistributedDisplacementControl::setContinuationTime() - lambda "
	   << "channel " << chan << " out of range\n";
    return -1;
  }
  if (target == 0.0 && !(dt > 0.0)) {
    opserr << "DistributedDisplacementControl::setContinuationTime() - need "
	   << "either -target (progress-normalised dt) or a positive -dt\n";
    return -1;
  }

  useContinuationTime = true;
  lambdaChannel = chan;
  stageDuration = duration;
  stageTarget = target;
  dtFixed = dt;

  // inherit a lambda left by a previous continuation stage (the load stays
  // applied); start from zero only if nobody has ever written this channel
  if (OPS_ContinuationLambda::isValid(lambdaChannel) == false)
    OPS_ContinuationLambda::set(lambdaChannel, 0.0);
  committedLambda = OPS_ContinuationLambda::get(lambdaChannel);
  currentLambda = committedLambda;
  OPS_ContinuationLambda::setOwner(lambdaChannel,
				   INTEGRATOR_TAGS_DistributedDisplacementControl);

  return 0;
}

// Make the committed lambda bit-identical on every rank.
//
// dLambda is computed from two components of the solution vector at theDofID,
// so every rank already gets the same value WHEN the distributed SOE returns a
// consistent solution -- which this class requires anyway, since domainChanged()
// hands the same global theDofID to all of them. This exchange therefore does
// not add an assumption; it makes the agreement explicit and turns a silent
// divergence (different load factors per rank => wrong results, no error) into
// a loud one. One scalar per step, negligible against a distributed solve.
//
// Only ever called in continuation mode, so a legacy parallel run is unchanged
// and cannot deadlock on it.
int
DistributedDisplacementControl::syncLambda(void)
{
  if (numChannels == 0 || theChannels == 0)
    return 0; // single process: nothing to agree on

  static Vector lamData(1);

  if (processID != 0) {
    Channel *theChannel = theChannels[0];
    lamData(0) = committedLambda;
    theChannel->sendVector(0, 0, lamData);
    theChannel->recvVector(0, 0, lamData);
    committedLambda = lamData(0); // P0 is authoritative
  }

  else {
    double scale = fabs(committedLambda) > 1.0 ? fabs(committedLambda) : 1.0;
    double tol = 1.0e-10 * scale;

    for (int j = 0; j < numChannels; j++) {
      theChannels[j]->recvVector(0, 0, lamData);
      double diff = fabs(lamData(0) - committedLambda);
      if (diff > tol) {
	opserr << "WARNING DistributedDisplacementControl::commit() - process "
	       << j + 1 << " computed lambda " << lamData(0) << " against "
	       << committedLambda << " on P0 (difference " << diff << "). The "
	       << "distributed SOE is not returning a consistent solution at the "
	       << "control dof; the results are NOT trustworthy. Forcing the P0 "
	       << "value.\n";
      }
    }

    for (int k = 0; k < numChannels; k++) {
      lamData(0) = committedLambda;
      theChannels[k]->sendVector(0, 0, lamData);
    }
  }

  currentLambda = committedLambda;
  OPS_ContinuationLambda::set(lambdaChannel, committedLambda);

  return 0;
}

int
DistributedDisplacementControl::commit(void)
{
  int res = this->StaticIntegrator::commit();
  if (res == 0 && useContinuationTime) {
    committedLambda = currentLambda;
    if (this->syncLambda() < 0)
      return -1;
  }
  return res;
}

int
DistributedDisplacementControl::revertToLastStep(void)
{
  if (useContinuationTime == false)
    return 0;

  // StaticAnalysis calls Domain::revertToLastCommit() BEFORE this, and that
  // re-applies the load while the channel still holds the failed lambda.
  // Restore lambda and re-apply, so the domain is left consistent for whatever
  // the script does next.
  currentLambda = committedLambda;
  deltaLambdaStep = 0.0;
  OPS_ContinuationLambda::set(lambdaChannel, committedLambda);

  AnalysisModel *theModel = this->getAnalysisModel();
  if (theModel != 0) {
    // Domain::revertToLastCommit has already put currentTime back to the
    // committed time
    timeStep = theModel->getCurrentDomainTime();
    theModel->applyLoadDomain(timeStep);
    theModel->updateDomain();
  }

  return 0;
}

int
DistributedDisplacementControl::sendSelf(int cTag, Channel &theChannel)
{
  int sendID =0;

  // if P0 check if already sent. If already sent use old processID; if not allocate a new process 
  // id for remote part of object, enlarge channel * to hold a channel * for this remote object.

  // if not P0, send current processID

  if (processID == 0) {

    // check if already using this object
    bool found = false;
    for (int i=0; i<numChannels; i++)
      if (theChannels[i] == &theChannel) {
	sendID = i+1;
	found = true;
      }

    // if new object, enlarge Channel pointers to hold new channel * & allocate new ID
    if (found == false) {
      int nextNumChannels = numChannels + 1;
      Channel **nextChannels = new Channel *[nextNumChannels];
      if (nextChannels == 0) {
	opserr << "DistributedDisplacementControl::sendSelf() - failed to allocate channel array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	nextChannels[i] = theChannels[i];
      nextChannels[numChannels] = &theChannel;
      
      numChannels = nextNumChannels;
      
      if (theChannels != 0)
	delete [] theChannels;
      
      theChannels = nextChannels;
      
      // allocate new processID for remote object
      sendID = numChannels;
    }

  } else 
    sendID = processID;

  // send remotes processID & info about node, dof and numIncr
  static ID idData(3);
  idData(0) = sendID;
  idData(1) = theNode;
  idData(2) = theDof;
  int res = theChannel.sendID(0, cTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedDisplacementControl::sendSelf() - failed to send data\n";
    return -1;
  }

  static Vector dData(11);
  dData(0) = theIncrement;
  dData(1) = minIncrement;
  dData(2) = maxIncrement;
  dData(3) = specNumIncrStep;
  dData(4) = numIncrLastStep;
  dData(5) = useContinuationTime ? 1.0 : 0.0;
  dData(6) = lambdaChannel;
  dData(7) = stageDuration;
  dData(8) = stageTarget;
  dData(9) = dtFixed;
  dData(10) = committedLambda;
  res = theChannel.sendVector(0, cTag, dData);
  if (res < 0) {
    opserr <<"WARNING DistributedDisplacementControl::recvSelf() - failed to recv vector data\n";
    return -1;
  }	        
  
  return 0;
}


int
DistributedDisplacementControl::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
			      
{
  static ID idData(3);
  int res = theChannel.recvID(0, cTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedDisplacementControl::recvSelf() - failed to recv id data\n";
    return -1;
  }	      
  processID = idData(0);
  theNode = idData(1);
  theDof  = idData(2);

  static Vector dData(11);
  res = theChannel.recvVector(0, cTag, dData);
  if (res < 0) {
    opserr <<"WARNING DistributedDisplacementControl::recvSelf() - failed to recv vector data\n";
    return -1;
  }

  theIncrement = dData(0);
  minIncrement = dData(1);
  maxIncrement = dData(2);
  specNumIncrStep = dData(3);
  numIncrLastStep = dData(4);
  useContinuationTime = (dData(5) != 0.0);
  lambdaChannel = (int)dData(6);
  stageDuration = dData(7);
  stageTarget = dData(8);
  dtFixed = dData(9);
  committedLambda = dData(10);
  currentLambda = committedLambda;

  // the lambda store is a per-process global and is deliberately NOT part of
  // the serialised state of the series; seed it here so the remote copy starts
  // from the same load factor
  if (useContinuationTime) {
    OPS_ContinuationLambda::set(lambdaChannel, committedLambda);
    OPS_ContinuationLambda::setOwner(lambdaChannel,
				     INTEGRATOR_TAGS_DistributedDisplacementControl);
  }

  // set the Channel & processID
  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  return 0;
}

void
DistributedDisplacementControl::Print(OPS_Stream &s, int flag)
{
    // TO FINISH    
}

int
DistributedDisplacementControl::setProcessID(int dTag) 
{
  processID = dTag;
  return 0;
}

int
DistributedDisplacementControl::setChannels(int nChannels, Channel **theC)
{
  numChannels = nChannels;

  if (theChannels != 0)
    delete [] theChannels;

  theChannels = new Channel *[numChannels];
  for (int i=0; i<numChannels; i++)
    theChannels[i] = theC[i];

  return 0;
}
