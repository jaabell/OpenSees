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
// Description: see LadrunoIndirectControl.h. The bookkeeping of dUhat,
// dUbar and deltaUstep follows DisplacementControl, with the control
// component dUhat(theDofID) replaced by the product c . dUhat.

#include <LadrunoIndirectControl.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <math.h>
#include <string.h>
#include <elementAPI.h>
#include <classTags.h>

void *OPS_LadrunoIndirectControl(void)
{
    // Usage: integrator LadrunoIndirectControl $incr
    //          -dof $node $dof $coef <-dof $node $dof $coef ...>
    //          <-iter $numIter $dmin $dmax>
    // crack mouth opening: -dof nA dofA 1.0 -dof nB dofB -1.0
    if (OPS_GetNumRemainingInputArgs() < 4) {
        opserr << "WARNING integrator LadrunoIndirectControl $incr "
                  "-dof $node $dof $coef <-dof ...> <-iter $numIter $dmin $dmax>\n";
        return 0;
    }

    int    numData = 1;
    double incr = 0.0;
    if (OPS_GetDoubleInput(&numData, &incr) < 0) {
        opserr << "WARNING integrator LadrunoIndirectControl failed to read incr\n";
        return 0;
    }

    // collect a variable number of (node, dof, coef) control entries
    static const int MAXC = 64;
    int    nodes[MAXC], dofs[MAXC];
    double coefs[MAXC];
    int    nctrl = 0;

    int    numIter = 1;
    double dmin = incr, dmax = incr;
    bool   haveIter = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-dof") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "WARNING integrator LadrunoIndirectControl -dof expects "
                          "$node $dof $coef\n";
                return 0;
            }
            if (nctrl >= MAXC) {
                opserr << "WARNING integrator LadrunoIndirectControl - too many -dof "
                          "entries (max " << MAXC << ")\n";
                return 0;
            }
            int nd2 = 2, id2[2];
            if (OPS_GetIntInput(&nd2, &id2[0]) < 0) {
                opserr << "WARNING integrator LadrunoIndirectControl -dof failed to "
                          "read node/dof\n";
                return 0;
            }
            int    nd1 = 1; double cf;
            if (OPS_GetDoubleInput(&nd1, &cf) < 0) {
                opserr << "WARNING integrator LadrunoIndirectControl -dof failed to "
                          "read coef\n";
                return 0;
            }
            nodes[nctrl] = id2[0];
            dofs[nctrl]  = id2[1];   // 1-based on input; converted below
            coefs[nctrl] = cf;
            nctrl++;
        } else if (strcmp(opt, "-iter") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING integrator LadrunoIndirectControl -iter expects "
                          "$numIter <$dmin $dmax>\n";
                return 0;
            }
            int nd1 = 1;
            if (OPS_GetIntInput(&nd1, &numIter) < 0) {
                opserr << "WARNING integrator LadrunoIndirectControl -iter failed to "
                          "read numIter\n";
                return 0;
            }
            // optional dmin dmax; tokens that are not numbers are left for
            // the option loop (the Python reader advances past a token it
            // cannot convert, the Tcl reader does not)
            int before = OPS_GetNumRemainingInputArgs();
            if (before >= 2) {
                int nd2 = 2; double mm[2];
                if (OPS_GetDoubleInput(&nd2, mm) == 0) {
                    dmin = mm[0]; dmax = mm[1];
                } else {
                    int consumed = before - OPS_GetNumRemainingInputArgs();
                    if (consumed > 0)
                        OPS_ResetCurrentInputArg(-consumed);
                }
            }
            haveIter = true;
        } else {
            opserr << "WARNING integrator LadrunoIndirectControl - unknown option "
                   << opt << " (ignored)\n";
        }
    }

    if (nctrl < 1) {
        opserr << "WARNING integrator LadrunoIndirectControl - need at least one "
                  "-dof $node $dof $coef control entry\n";
        return 0;
    }
    if (!haveIter) { dmin = incr; dmax = incr; }   // pin increment (non-adaptive)

    // validate control nodes/dofs against the domain
    Domain *theDomain = OPS_GetDomain();
    ID  idNode(nctrl), idDof(nctrl);
    Vector vCoef(nctrl);
    for (int k = 0; k < nctrl; k++) {
        Node *theNode = (theDomain != 0) ? theDomain->getNode(nodes[k]) : 0;
        if (theNode == 0) {
            opserr << "WARNING integrator LadrunoIndirectControl - control node "
                   << nodes[k] << " does not exist\n";
            return 0;
        }
        int numDOF = theNode->getNumberDOF();
        if (dofs[k] <= 0 || dofs[k] > numDOF) {
            opserr << "WARNING integrator LadrunoIndirectControl - invalid dof "
                   << dofs[k] << " at node " << nodes[k] << "\n";
            return 0;
        }
        idNode(k) = nodes[k];
        idDof(k)  = dofs[k] - 1;   // store 0-based
        vCoef(k)  = coefs[k];
    }

    return new LadrunoIndirectControl(idNode, idDof, vCoef, incr, numIter, dmin, dmax);
}

LadrunoIndirectControl::LadrunoIndirectControl(const ID &nodes, const ID &dofs,
                                               const Vector &coefs, double incr,
                                               int numIncr, double min, double max)
 :StaticIntegrator(INTEGRATOR_TAGS_LadrunoIndirectControl),
  ctrlNode(nodes), ctrlDof(dofs), ctrlCoef(coefs),
  theIncrement(incr),
  specNumIncrStep(numIncr > 0 ? numIncr : 1),
  numIncrLastStep(numIncr > 0 ? numIncr : 1),
  minIncrement(min), maxIncrement(max),
  ctrlEqn(nodes.Size()),
  deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), phat(0),
  deltaLambdaStep(0.0), currentLambda(0.0)
{
    ctrlEqn.Zero();
}

LadrunoIndirectControl::LadrunoIndirectControl()
 :StaticIntegrator(INTEGRATOR_TAGS_LadrunoIndirectControl),
  ctrlNode(0), ctrlDof(0), ctrlCoef(0),
  theIncrement(0.0), specNumIncrStep(1), numIncrLastStep(1),
  minIncrement(0.0), maxIncrement(0.0), ctrlEqn(0),
  deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), phat(0),
  deltaLambdaStep(0.0), currentLambda(0.0)
{

}

LadrunoIndirectControl::~LadrunoIndirectControl()
{
    if (deltaUhat  != 0) delete deltaUhat;
    if (deltaUbar  != 0) delete deltaUbar;
    if (deltaU     != 0) delete deltaU;
    if (deltaUstep != 0) delete deltaUstep;
    if (phat       != 0) delete phat;
}

// c . v  over the controlled equations (skips fixed/constrained entries)
double
LadrunoIndirectControl::controlDot(const Vector &v) const
{
    double s = 0.0;
    for (int k = 0; k < ctrlEqn.Size(); k++) {
        int loc = ctrlEqn(k);
        if (loc >= 0) s += ctrlCoef(k) * v(loc);
    }
    return s;
}

int
LadrunoIndirectControl::newStep(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING LadrunoIndirectControl::newStep() - "
                  "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }

    // iteration-count adaptation (Ramm 1981); no-op without -iter (dmin = dmax = incr)
    double factor = double(specNumIncrStep) / (numIncrLastStep > 0 ? numIncrLastStep : 1);
    theIncrement *= factor;
    if (theIncrement < minIncrement)      theIncrement = minIncrement;
    else if (theIncrement > maxIncrement) theIncrement = maxIncrement;

    currentLambda = theModel->getCurrentDomainTime();

    // determine dUhat = K^{-1} phat
    this->formTangent();
    theLinSOE->setB(*phat);
    if (theLinSOE->solve() < 0) {
        opserr << "LadrunoIndirectControl::newStep() - failed in solver\n";
        return -1;
    }
    (*deltaUhat) = theLinSOE->getX();

    double cdotHat = this->controlDot(*deltaUhat);
    if (cdotHat == 0.0) {
        opserr << "WARNING LadrunoIndirectControl::newStep() - control quantity "
                  "c.dUhat is zero (control DOFs not excited by the reference load)\n";
        return -1;
    }

    // dLambda(1) so that c.dU = theIncrement
    double dLambda = theIncrement / cdotHat;
    deltaLambdaStep = dLambda;
    currentLambda  += dLambda;

    (*deltaU) = *deltaUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    theModel->incrDisp(*deltaU);
    theModel->applyLoadDomain(currentLambda);
    if (theModel->updateDomain() < 0) {
        opserr << "LadrunoIndirectControl::newStep() - model failed to update for new dU\n";
        return -1;
    }

    numIncrLastStep = 0;
    return 0;
}

int
LadrunoIndirectControl::update(const Vector &dU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING LadrunoIndirectControl::update() - "
                  "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }

    (*deltaUbar) = dU;   // residual displacement (SOE will be overwritten next)

    // determine dUhat = K^{-1} phat
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();

    double cdotHat = this->controlDot(*deltaUhat);
    if (cdotHat == 0.0) {
        opserr << "WARNING LadrunoIndirectControl::update() - control quantity "
                  "c.dUhat is zero\n";
        return -1;
    }

    // hold the control quantity fixed within the step: c.(dUbar + dLambda dUhat) = 0
    double dLambda = -this->controlDot(*deltaUbar) / cdotHat;

    (*deltaU) = (*deltaUbar);
    deltaU->addVector(1.0, *deltaUhat, dLambda);

    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda   += dLambda;

    theModel->incrDisp(*deltaU);
    theModel->applyLoadDomain(currentLambda);
    if (theModel->updateDomain() < 0) {
        opserr << "LadrunoIndirectControl::update() - model failed to update for new dU\n";
        return -1;
    }

    // set the X soln in the SOE to deltaU for the convergence test
    theLinSOE->setX(*deltaU);
    numIncrLastStep++;
    return 0;
}

int
LadrunoIndirectControl::domainChanged(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING LadrunoIndirectControl::domainChanged() - "
                  "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }
    int size = theModel->getNumEqn();

    if (deltaUhat == 0 || deltaUhat->Size() != size) {
        if (deltaUhat != 0) delete deltaUhat;  deltaUhat  = new Vector(size);
    }
    if (deltaUbar == 0 || deltaUbar->Size() != size) {
        if (deltaUbar != 0) delete deltaUbar;  deltaUbar  = new Vector(size);
    }
    if (deltaU == 0 || deltaU->Size() != size) {
        if (deltaU != 0) delete deltaU;        deltaU     = new Vector(size);
    }
    if (deltaUstep == 0 || deltaUstep->Size() != size) {
        if (deltaUstep != 0) delete deltaUstep; deltaUstep = new Vector(size);
    }
    if (phat == 0 || phat->Size() != size) {
        if (phat != 0) delete phat;            phat       = new Vector(size);
    }
    if (deltaUhat == 0 || deltaUbar == 0 || deltaU == 0 ||
        deltaUstep == 0 || phat == 0 || phat->Size() != size) {
        opserr << "FATAL LadrunoIndirectControl::domainChanged() - ran out of memory "
                  "for Vectors of size " << size << endln;
        return -1;
    }

    // determine phat: increment lambda by 1, apply load, read the unbalance
    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);
    this->formUnbalance();   // assumes unbalance at last was 0
    (*phat) = theLinSOE->getB();
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);

    int haveLoad = 0;
    for (int i = 0; i < size; i++)
        if ((*phat)(i) != 0.0) { haveLoad = 1; break; }
    if (haveLoad == 0) {
        opserr << "WARNING LadrunoIndirectControl::domainChanged() - zero reference load\n";
        return -1;
    }

    // resolve the control equation numbers from the node DOF groups
    Domain *theDomain = theModel->getDomainPtr();
    for (int k = 0; k < ctrlNode.Size(); k++) {
        ctrlEqn(k) = -1;
        Node *theNodePtr = (theDomain != 0) ? theDomain->getNode(ctrlNode(k)) : 0;
        if (theNodePtr == 0) {
            opserr << "WARNING LadrunoIndirectControl::domainChanged() - control node "
                   << ctrlNode(k) << " not found\n";
            return -1;
        }
        DOF_Group *theGroup = theNodePtr->getDOF_GroupPtr();
        if (theGroup == 0) continue;            // no DOF group yet
        const ID &theID = theGroup->getID();
        if (ctrlDof(k) >= 0 && ctrlDof(k) < theID.Size())
            ctrlEqn(k) = theID(ctrlDof(k));     // -1 here means fixed/constrained
    }

    // at least one control DOF must be free
    int haveFree = 0;
    for (int k = 0; k < ctrlEqn.Size(); k++) if (ctrlEqn(k) >= 0) { haveFree = 1; break; }
    if (haveFree == 0) {
        opserr << "WARNING LadrunoIndirectControl::domainChanged() - all control DOFs "
                  "are fixed/constrained; nothing to control\n";
        return -1;
    }
    return 0;
}

int
LadrunoIndirectControl::sendSelf(int cTag, Channel &theChannel)
{
    int n = ctrlNode.Size();
    // header: [nctrl, specNumIncrStep, numIncrLastStep]
    ID header(3);
    header(0) = n;
    header(1) = specNumIncrStep;
    header(2) = numIncrLastStep;
    if (theChannel.sendID(this->getDbTag(), cTag, header) < 0) {
        opserr << "LadrunoIndirectControl::sendSelf() - failed to send header\n";
        return -1;
    }
    // node/dof control table [node_0..node_{n-1}, dof_0..dof_{n-1}]
    ID nd(2 * n);
    for (int k = 0; k < n; k++) { nd(k) = ctrlNode(k); nd(n + k) = ctrlDof(k); }
    if (n > 0 && theChannel.sendID(this->getDbTag(), cTag, nd) < 0) {
        opserr << "LadrunoIndirectControl::sendSelf() - failed to send control table\n";
        return -1;
    }
    // doubles: [coef_0..coef_{n-1}, theIncrement, minIncrement, maxIncrement,
    //           deltaLambdaStep, currentLambda]
    Vector data(n + 5);
    for (int k = 0; k < n; k++) data(k) = ctrlCoef(k);
    data(n)     = theIncrement;
    data(n + 1) = minIncrement;
    data(n + 2) = maxIncrement;
    data(n + 3) = deltaLambdaStep;
    data(n + 4) = currentLambda;
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoIndirectControl::sendSelf() - failed to send data\n";
        return -1;
    }
    return 0;
}

int
LadrunoIndirectControl::recvSelf(int cTag, Channel &theChannel,
                                 FEM_ObjectBroker &theBroker)
{
    ID header(3);
    if (theChannel.recvID(this->getDbTag(), cTag, header) < 0) {
        opserr << "LadrunoIndirectControl::recvSelf() - failed to receive header\n";
        return -1;
    }
    int n = header(0);
    specNumIncrStep = header(1);
    numIncrLastStep = header(2);

    ctrlNode = ID(n);
    ctrlDof  = ID(n);
    ctrlCoef = Vector(n);
    ctrlEqn  = ID(n);
    ctrlEqn.Zero();

    if (n > 0) {
        ID nd(2 * n);
        if (theChannel.recvID(this->getDbTag(), cTag, nd) < 0) {
            opserr << "LadrunoIndirectControl::recvSelf() - failed to receive control table\n";
            return -1;
        }
        for (int k = 0; k < n; k++) { ctrlNode(k) = nd(k); ctrlDof(k) = nd(n + k); }
    }

    Vector data(n + 5);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoIndirectControl::recvSelf() - failed to receive data\n";
        return -1;
    }
    for (int k = 0; k < n; k++) ctrlCoef(k) = data(k);
    theIncrement    = data(n);
    minIncrement    = data(n + 1);
    maxIncrement    = data(n + 2);
    deltaLambdaStep = data(n + 3);
    currentLambda   = data(n + 4);
    return 0;
}

void
LadrunoIndirectControl::Print(OPS_Stream &s, int flag)
{
    s << "\t LadrunoIndirectControl - increment: " << theIncrement
      << "  currentLambda: " << currentLambda << endln;
    s << "\t control quantity c.U over " << ctrlNode.Size() << " dof(s):";
    for (int k = 0; k < ctrlNode.Size(); k++)
        s << " (" << ctrlCoef(k) << " * node " << ctrlNode(k)
          << " dof " << (ctrlDof(k) + 1) << ")";
    s << endln;
}
