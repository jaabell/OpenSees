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
// References:
//   J.O. Hallquist, LS-DYNA Theory Manual, Livermore Software Technology
//   Corporation, 2006, Sections 28.1 and 28.2.
//
// Description: implementation of LadrunoProjectionHandler (see the header).

#include <LadrunoProjectionHandler.h>
#include <LadrunoConstraintProjector.h>
#include <LadrunoProjectionConsumer.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map>
#include <utility>

#include <AnalysisModel.h>
#include <Domain.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <EQ_Constraint.h>
#include <EQ_ConstraintIter.h>
#include <Integrator.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <Subdomain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <classTags.h>

// Rotation of a retained node (rad) above which the small-rotation lever arm
// frozen in C (rigidLink -beam, rigidDiaphragm) is reported as no longer valid.
#define LADRUNO_CCR_STALE_RAD 0.1

// constraints LadrunoProjection <-verbose> <-projectICs> <-icTol $tol>
void *OPS_LadrunoProjectionHandler(void)
{
    bool verbose = false;
    bool projectICs = false;
    double icTol = 1.0e-8;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        if (strcmp(arg, "-verbose") == 0)
            verbose = true;
        else if (strcmp(arg, "-projectICs") == 0)
            projectICs = true;
        else if (strcmp(arg, "-icTol") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int n = 1;
                OPS_GetDoubleInput(&n, &icTol);
            }
        } else {
            opserr << "WARNING LadrunoProjection - unknown option " << arg << " (ignored)\n";
        }
    }
    return new LadrunoProjectionHandler(verbose, projectICs, icTol);
}

LadrunoProjectionHandler::LadrunoProjectionHandler()
  : ConstraintHandler(HANDLER_TAG_LadrunoProjectionHandler),
    vtxNode(), vtxDof(), theGroups(), theProjector(0),
    verbose(false), projectICs(false), icTol(1.0e-8)
{
}

LadrunoProjectionHandler::LadrunoProjectionHandler(bool v, bool p, double t)
  : ConstraintHandler(HANDLER_TAG_LadrunoProjectionHandler),
    vtxNode(), vtxDof(), theGroups(), theProjector(0),
    verbose(v), projectICs(p), icTol(t)
{
}

LadrunoProjectionHandler::~LadrunoProjectionHandler()
{
    if (theProjector != 0) delete theProjector;
}

// PlainHandler assembly without elimination of the constrained DOFs, then
// the constraint groups and their checks.
int
LadrunoProjectionHandler::handle(const ID *nodesLast)
{
    Domain *theDomain = this->getDomainPtr();
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Integrator *theIntegrator = this->getIntegratorPtr();
    if (theDomain == 0 || theModel == 0 || theIntegrator == 0) {
        opserr << "WARNING LadrunoProjectionHandler::handle() - setLinks() not called\n";
        return -1;
    }

    prescribedDOFs.clear();
    prescribedKey.clear();
    derivedSlaves.clear();
    derivedKey.clear();

    // Every SP DOF is removed from the equation set. The non-homogeneous ones
    // (and imposedMotion) are recorded so that applyLoad() imposes their
    // displacement.
    std::multimap<int, SP_Constraint *> allSPs;
    SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP;
    while ((theSP = theSPs()) != 0)
        allSPs.insert(std::make_pair(theSP->getNodeTag(), theSP));

    // DOF_Groups: free = -2, SP = -1. Constrained MP DOFs stay free (-2):
    // they keep their equation and diagonal mass, the projector ties them.
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    int numDOF = 0, countDOF = 0, count3 = 0;
    while ((nodPtr = theNod()) != 0) {
        DOF_Group *dofPtr = new DOF_Group(numDOF++, nodPtr);
        if (dofPtr == 0) {
            opserr << "WARNING LadrunoProjectionHandler::handle() - out of memory (DOF_Group)\n";
            return -4;
        }
        const ID &id = dofPtr->getID();
        for (int j = 0; j < id.Size(); j++) {
            dofPtr->setID(j, -2);
            countDOF++;
        }
        int nodeID = nodPtr->getTag();
        std::multimap<int, SP_Constraint *>::iterator first = allSPs.lower_bound(nodeID);
        std::multimap<int, SP_Constraint *>::iterator last = allSPs.upper_bound(nodeID);
        for (std::multimap<int, SP_Constraint *>::iterator it = first; it != last; it++) {
            int dof = it->second->getDOF_Number();
            const ID &cid = dofPtr->getID();
            if (cid(dof) == -2) {
                dofPtr->setID(dof, -1);
                countDOF--;
                if (it->second->isHomogeneous() == false) {
                    PrescribedDOF p;
                    p.nodeTag = nodeID;
                    p.dof     = dof;
                    p.sp      = it->second;
                    prescribedDOFs.push_back(p);
                    prescribedKey.insert(std::make_pair(nodeID, dof));
                }
            } else {
                opserr << "WARNING LadrunoProjectionHandler::handle() - multiple SPs at DOF "
                       << dof << " node " << nodeID << endln;
            }
        }
        nodPtr->setDOF_GroupPtr(dofPtr);
        theModel->addDOF_Group(dofPtr);
    }
    theModel->setNumEqn(countDOF);

    // subdomain boundary nodes are numbered last (-3), as in PlainHandler
    if (nodesLast != 0)
        for (int i = 0; i < nodesLast->Size(); i++) {
            Node *np = theDomain->getNode((*nodesLast)(i));
            if (np != 0) {
                DOF_Group *dp = np->getDOF_GroupPtr();
                const ID &id = dp->getID();
                for (int j = 0; j < id.Size(); j++)
                    if (id(j) == -2) { dp->setID(j, -3); count3++; }
            }
        }

    // one FE_Element per element, as in PlainHandler
    ElementIter &theEle = theDomain->getElements();
    Element *elePtr;
    int numFe = 0;
    while ((elePtr = theEle()) != 0) {
        if (elePtr->isSubdomain() == true) {
            Subdomain *theSub = (Subdomain *)elePtr;
            if (theSub->doesIndependentAnalysis() == false) {
                FE_Element *fePtr = new FE_Element(numFe++, elePtr);
                if (fePtr == 0) return -5;
                theModel->addFE_Element(fePtr);
                theSub->setFE_ElementPtr(fePtr);
            }
        } else {
            FE_Element *fePtr = new FE_Element(numFe++, elePtr);
            if (fePtr == 0) return -5;
            theModel->addFE_Element(fePtr);
        }
    }

    if (this->buildGroups() < 0)
        return -6;
    // constrained DOFs driven only by prescribed DOFs leave the equation set
    if (this->classifyDerivedSlaves(countDOF) < 0)
        return -8;
    theModel->setNumEqn(countDOF);
    // With no free equation the explicit integrator has nothing to solve and
    // the system would have size zero. Boundary DOFs (count3) are numbered
    // later and count as free.
    if (countDOF == 0 && count3 == 0) {
        opserr << "LadrunoProjectionHandler::handle() - no free equations: every DOF is fixed, "
                  "prescribed (SP or imposedMotion) or tied only to prescribed DOFs. The "
                  "projection handler does not support pure kinematic playback.\n";
        return -9;
    }
    if (this->massGuard() < 0)
        return -7;

    if (verbose)
        opserr << "LadrunoProjectionHandler: " << (int)theGroups.size()
               << " constraint group(s), " << countDOF << " equations.\n";

    return count3;
}

// Union-find over the (node, dof) vertices of all MP and EQ constraints;
// classification of retained and constrained DOFs.
int
LadrunoProjectionHandler::buildGroups(void)
{
    Domain *theDomain = this->getDomainPtr();
    theGroups.clear();
    vtxNode.clear();
    vtxDof.clear();
    rotMonitors.clear();

    std::map<std::pair<int, int>, int> vtxOf;   // (node, dof) -> vertex
    std::vector<int> parent;                     // union-find forest

    // one edge per (constrained i, retained j) with C(i,j) != 0
    struct Edge { int slaveVtx, masterVtx; };
    std::vector<Edge> edges;
    std::vector<int> slaveMPcount;               // number of constraints a vertex is constrained by
    std::vector<char> isMaster, isSlave;

    MP_ConstraintIter &theMPs = theDomain->getMPs();
    MP_Constraint *mpPtr;
    while ((mpPtr = theMPs()) != 0) {
        int cNode = mpPtr->getNodeConstrained();
        int rNode = mpPtr->getNodeRetained();
        const ID &cDofs = mpPtr->getConstrainedDOFs();
        const ID &rDofs = mpPtr->getRetainedDOFs();
        const Matrix &Ccr = mpPtr->getConstraint();

        for (int i = 0; i < cDofs.Size(); i++) {
            std::pair<int, int> cp(cNode, cDofs(i));
            int cv;
            std::map<std::pair<int,int>,int>::iterator itc = vtxOf.find(cp);
            if (itc == vtxOf.end()) {
                cv = (int)vtxNode.size();
                vtxOf[cp] = cv; vtxNode.push_back(cNode); vtxDof.push_back(cDofs(i));
                parent.push_back(cv); slaveMPcount.push_back(0);
                isMaster.push_back(0); isSlave.push_back(0);
            } else cv = itc->second;
            isSlave[cv] = 1;
            slaveMPcount[cv] += 1;

            for (int j = 0; j < rDofs.Size(); j++) {
                double c = Ccr(i, j);
                if (c == 0.0) continue;
                // a retained rotation entering a constrained translation is a
                // frozen lever arm; monitor it
                this->flagRotMonitor(rNode, rDofs(j), cNode, cDofs(i));
                std::pair<int, int> rp(rNode, rDofs(j));
                int rv;
                std::map<std::pair<int,int>,int>::iterator itr = vtxOf.find(rp);
                if (itr == vtxOf.end()) {
                    rv = (int)vtxNode.size();
                    vtxOf[rp] = rv; vtxNode.push_back(rNode); vtxDof.push_back(rDofs(j));
                    parent.push_back(rv); slaveMPcount.push_back(0);
                    isMaster.push_back(0); isSlave.push_back(0);
                } else rv = itr->second;
                isMaster[rv] = 1;
                Edge e; e.slaveVtx = cv; e.masterVtx = rv;
                edges.push_back(e);
            }
        }
    }

    // EQ_Constraint: one constrained DOF and a coefficient vector over
    // retained DOFs, possibly on several nodes (a general row of C)
    EQ_ConstraintIter &theEQs = theDomain->getEQs();
    EQ_Constraint *eqPtr;
    while ((eqPtr = theEQs()) != 0) {
        int cNode = eqPtr->getNodeConstrained();
        int cDof  = eqPtr->getConstrainedDOFs();
        const ID &rNodes = eqPtr->getNodeRetained();
        const ID &rDofs  = eqPtr->getRetainedDOFs();
        const Vector &C  = eqPtr->getConstraint();

        std::pair<int, int> cp(cNode, cDof);
        int cv;
        std::map<std::pair<int,int>,int>::iterator itc = vtxOf.find(cp);
        if (itc == vtxOf.end()) {
            cv = (int)vtxNode.size();
            vtxOf[cp] = cv; vtxNode.push_back(cNode); vtxDof.push_back(cDof);
            parent.push_back(cv); slaveMPcount.push_back(0);
            isMaster.push_back(0); isSlave.push_back(0);
        } else cv = itc->second;
        isSlave[cv] = 1;
        slaveMPcount[cv] += 1;

        for (int k = 0; k < rDofs.Size(); k++) {
            double c = C(k);
            if (c == 0.0) continue;
            std::pair<int, int> rp(rNodes(k), rDofs(k));
            int rv;
            std::map<std::pair<int,int>,int>::iterator itr = vtxOf.find(rp);
            if (itr == vtxOf.end()) {
                rv = (int)vtxNode.size();
                vtxOf[rp] = rv; vtxNode.push_back(rNodes(k)); vtxDof.push_back(rDofs(k));
                parent.push_back(rv); slaveMPcount.push_back(0);
                isMaster.push_back(0); isSlave.push_back(0);
            } else rv = itr->second;
            isMaster[rv] = 1;
            Edge e; e.slaveVtx = cv; e.masterVtx = rv;
            edges.push_back(e);
        }
    }

    int nv = (int)vtxNode.size();
    if (nv == 0)
        return 0;   // no constraints: nothing to project

    // Every constrained or retained node must be in this domain. A missing
    // node means a constraint crossing a partition boundary, which is not
    // supported; doneNumberingDOF() would otherwise drop it as if fixed.
    for (int v = 0; v < nv; v++) {
        if (theDomain->getNode(vtxNode[v]) == 0) {
            opserr << "LadrunoProjectionHandler - constraint references node " << vtxNode[v]
                   << ", which is not in this domain. Constraints crossing partition "
                      "boundaries are not supported; use constraints Transformation for "
                      "the distributed analysis.\n";
            return -1;
        }
    }

    // union of the endpoints of every edge (find with path halving)
    std::vector<int> &par = parent;
    for (int i = 0; i < nv; i++) par[i] = i;
    for (int e = 0; e < (int)edges.size(); e++) {
        int a = edges[e].slaveVtx, b = edges[e].masterVtx;
        while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; }
        while (par[b] != b) { par[b] = par[par[b]]; b = par[b]; }
        if (a != b) par[a] = b;
    }

    // a DOF both retained and constrained (chain), or constrained twice
    for (int v = 0; v < nv; v++) {
        if (isMaster[v] && isSlave[v]) {
            opserr << "LadrunoProjectionHandler - CHAIN: DOF (node " << vtxNode[v]
                   << ", dof " << vtxDof[v] << ") is both a retained and a constrained "
                      "DOF. Constraint chains are not supported; compose the constraints "
                      "or use constraints Transformation.\n";
            return -1;
        }
        if (slaveMPcount[v] > 1) {
            opserr << "LadrunoProjectionHandler - DOUBLE CONSTRAINT: DOF (node "
                   << vtxNode[v] << ", dof " << vtxDof[v] << ") is constrained by more "
                      "than one constraint. Remove the redundant one.\n";
            return -1;
        }
    }

    // one group per root
    std::map<int, int> rootToGroup;
    for (int v = 0; v < nv; v++) {
        int r = v;
        while (par[r] != r) r = par[r];
        if (rootToGroup.find(r) == rootToGroup.end()) {
            rootToGroup[r] = (int)theGroups.size();
            theGroups.push_back(RawGroup());
        }
    }
    for (int v = 0; v < nv; v++) {
        if (!isMaster[v]) continue;
        int r = v; while (par[r] != r) r = par[r];
        theGroups[rootToGroup[r]].masterVtx.push_back(v);
    }

    // constrained rows with their retained DOFs, coefficients and offsets
    MP_ConstraintIter &theMPs2 = theDomain->getMPs();
    while ((mpPtr = theMPs2()) != 0) {
        int cNode = mpPtr->getNodeConstrained();
        int rNode = mpPtr->getNodeRetained();
        const ID &cDofs = mpPtr->getConstrainedDOFs();
        const ID &rDofs = mpPtr->getRetainedDOFs();
        const Matrix &Ccr = mpPtr->getConstraint();
        const Vector &Uc0 = mpPtr->getConstrainedDOFsInitialDisplacement();
        const Vector &Ur0 = mpPtr->getRetainedDOFsInitialDisplacement();
        for (int i = 0; i < cDofs.Size(); i++) {
            int cv = vtxOf[std::make_pair(cNode, cDofs(i))];
            int r = cv; while (par[r] != r) r = par[r];
            RawGroup &g = theGroups[rootToGroup[r]];
            SlaveRec s;
            s.vtx = cv;
            s.delta = (i < Uc0.Size()) ? Uc0(i) : 0.0;
            for (int j = 0; j < rDofs.Size(); j++) {
                double c = Ccr(i, j);
                if (c == 0.0) continue;
                int rv = vtxOf[std::make_pair(rNode, rDofs(j))];
                s.masterVtx.push_back(rv);
                s.coeff.push_back(c);
                if (j < Ur0.Size()) s.delta -= c * Ur0(j);
            }
            g.slaves.push_back(s);
        }
    }

    EQ_ConstraintIter &theEQs2 = theDomain->getEQs();
    while ((eqPtr = theEQs2()) != 0) {
        int cNode = eqPtr->getNodeConstrained();
        int cDof  = eqPtr->getConstrainedDOFs();
        const ID &rNodes = eqPtr->getNodeRetained();
        const ID &rDofs  = eqPtr->getRetainedDOFs();
        const Vector &C  = eqPtr->getConstraint();
        double Uc0 = eqPtr->getConstrainedDOFsInitialDisplacement();
        const Vector &Ur0 = eqPtr->getRetainedDOFsInitialDisplacement();

        int cv = vtxOf[std::make_pair(cNode, cDof)];
        int r = cv; while (par[r] != r) r = par[r];
        RawGroup &g = theGroups[rootToGroup[r]];
        SlaveRec s;
        s.vtx = cv;
        s.delta = Uc0;
        for (int k = 0; k < rDofs.Size(); k++) {
            double c = C(k);
            if (c == 0.0) continue;
            int rv = vtxOf[std::make_pair(rNodes(k), rDofs(k))];
            s.masterVtx.push_back(rv);
            s.coeff.push_back(c);
            if (k < Ur0.Size()) s.delta -= c * Ur0(k);
        }
        g.slaves.push_back(s);
    }
    return 0;
}

// A constrained DOF tied only to prescribed retained DOFs is removed from the
// equation set (-1) and recorded for applyLoad(), where its u, v, a are set
// exactly from the retained DOFs. A constrained DOF tied to both a free and a
// prescribed retained DOF is refused: its displacement tie to the prescribed
// DOF cannot be held by an acceleration projection. A retained DOF fixed by a
// homogeneous SP contributes zero (its -Ur0 term is already in delta).
int
LadrunoProjectionHandler::classifyDerivedSlaves(int &countDOF)
{
    Domain *theDomain = this->getDomainPtr();
    for (int gi = 0; gi < (int)theGroups.size(); gi++) {
        RawGroup &g = theGroups[gi];
        for (int s = 0; s < (int)g.slaves.size(); s++) {
            SlaveRec &sr = g.slaves[s];
            Node *snp = theDomain->getNode(vtxNode[sr.vtx]);
            int ec = (snp != 0) ? snp->getDOF_GroupPtr()->getID()(vtxDof[sr.vtx]) : -1;
            // only a free (-2) constrained DOF can be driven by prescribed
            // DOFs; -1 (fixed or prescribed constrained DOF) is refused in
            // doneNumberingDOF(), -3 is a boundary DOF
            if (ec != -2)
                continue;

            bool hasFree = false, hasPres = false;
            for (int k = 0; k < (int)sr.masterVtx.size(); k++) {
                int mv = sr.masterVtx[k];
                Node *mnp = theDomain->getNode(vtxNode[mv]);
                int me = (mnp != 0) ? mnp->getDOF_GroupPtr()->getID()(vtxDof[mv]) : -1;
                bool mPres = prescribedKey.find(std::make_pair(vtxNode[mv], vtxDof[mv]))
                             != prescribedKey.end();
                if (mPres) hasPres = true;
                else if (me != -1) hasFree = true;
            }

            if (hasFree && hasPres) {
                opserr << "LadrunoProjectionHandler - MIXED RETAINED DOFS: constrained DOF (node "
                       << vtxNode[sr.vtx] << " dof " << vtxDof[sr.vtx] << ") is tied to both "
                          "a free DOF and a prescribed-motion DOF, which the acceleration "
                          "projection cannot enforce. Split the constraint, or use "
                          "constraints Transformation.\n";
                return -1;
            }
            if (!hasPres)
                continue;

            snp->getDOF_GroupPtr()->setID(vtxDof[sr.vtx], -1);
            countDOF--;
            DerivedSlave d;
            d.nodeTag = vtxNode[sr.vtx];
            d.dof     = vtxDof[sr.vtx];
            d.delta   = sr.delta;
            for (int k = 0; k < (int)sr.masterVtx.size(); k++) {
                int mv = sr.masterVtx[k];
                if (prescribedKey.find(std::make_pair(vtxNode[mv], vtxDof[mv])) == prescribedKey.end())
                    continue;
                d.masterNode.push_back(vtxNode[mv]);
                d.masterDof.push_back(vtxDof[mv]);
                d.coeff.push_back(sr.coeff[k]);
            }
            derivedSlaves.push_back(d);
            derivedKey.insert(std::make_pair(d.nodeTag, d.dof));

            if (verbose)
                opserr << "  LadrunoProjection: DOF (node " << d.nodeTag << " dof " << d.dof
                       << ") driven by " << (int)d.masterNode.size()
                       << " prescribed DOF(s), imposed kinematically.\n";
        }
    }
    return 0;
}

// Mass checks on the tied DOFs that keep an equation:
// (1) consistent element mass coupling a tied DOF is refused (the projection
//     assumes the diagonal mass that system Diagonal assembles);
// (2) a tied DOF with zero lumped mass (nodal plus element diagonal) is
//     refused, since it keeps its own equation and the diagonal system would
//     be singular. This is the place where the requirement is enforced for
//     every MP and EQ constraint, independently of the order of the commands.
//     DOFs of nodes whose element mass cannot be mapped are left to the check
//     on the assembled diagonal in LadrunoConstraintProjector::buildMass().
int
LadrunoProjectionHandler::massGuard(void)
{
    Domain *theDomain = this->getDomainPtr();
    if (theGroups.empty()) return 0;

    // tied (node, dof) pairs that keep an equation, with their lumped mass
    std::map<std::pair<int,int>, double> tied;
    for (int g = 0; g < (int)theGroups.size(); g++) {
        std::vector<int> verts(theGroups[g].masterVtx);
        for (int s = 0; s < (int)theGroups[g].slaves.size(); s++)
            verts.push_back(theGroups[g].slaves[s].vtx);
        for (int k = 0; k < (int)verts.size(); k++) {
            int v = verts[k];
            Node *np = theDomain->getNode(vtxNode[v]);
            if (np == 0 || np->getDOF_GroupPtr() == 0) continue;
            if (np->getDOF_GroupPtr()->getID()(vtxDof[v]) != -2) continue;
            tied[std::make_pair(vtxNode[v], vtxDof[v])] = 0.0;
        }
    }
    if (tied.empty()) return 0;

    std::map<int, char> tiedNode;
    for (std::map<std::pair<int,int>,double>::iterator it = tied.begin(); it != tied.end(); ++it)
        tiedNode[it->first.first] = 1;

    // nodal mass
    for (std::map<int,char>::iterator it = tiedNode.begin(); it != tiedNode.end(); ++it) {
        Node *np = theDomain->getNode(it->first);
        if (np == 0) continue;
        const Matrix &Mn = np->getMass();
        for (int d = 0; d < Mn.noRows() && d < Mn.noCols(); d++) {
            std::map<std::pair<int,int>,double>::iterator t = tied.find(std::make_pair(it->first, d));
            if (t != tied.end()) t->second += Mn(d, d);
        }
    }

    std::map<int, char> unmappedNode;
    ElementIter &theEle = theDomain->getElements();
    Element *elePtr;
    while ((elePtr = theEle()) != 0) {
        const ID &enodes = elePtr->getExternalNodes();
        bool incident = false;
        for (int k = 0; k < enodes.Size(); k++)
            if (tiedNode.find(enodes(k)) != tiedNode.end()) { incident = true; break; }
        if (!incident) continue;

        const Matrix &M = elePtr->getMass();
        int n = M.noRows();
        if (n == 0 || M.noCols() != n) continue;

        // local index -> (node, dof), nodes in order with their ndf each
        std::vector<int> idxNode(n), idxDof(n);
        int pos = 0;
        for (int k = 0; k < enodes.Size() && pos < n; k++) {
            Node *np = theDomain->getNode(enodes(k));
            int nd = (np != 0) ? np->getNumberDOF() : 0;
            for (int d = 0; d < nd && pos < n; d++) { idxNode[pos] = enodes(k); idxDof[pos] = d; pos++; }
        }
        if (pos != n) {
            for (int k = 0; k < enodes.Size(); k++) unmappedNode[enodes(k)] = 1;
            continue;
        }

        double maxDiag = 0.0;
        for (int p = 0; p < n; p++) if (fabs(M(p,p)) > maxDiag) maxDiag = fabs(M(p,p));
        if (maxDiag == 0.0) continue;
        double tol = 1.0e-8 * maxDiag;

        for (int p = 0; p < n; p++) {
            std::map<std::pair<int,int>,double>::iterator t =
                tied.find(std::make_pair(idxNode[p], idxDof[p]));
            if (t != tied.end()) t->second += M(p, p);
        }

        for (int p = 0; p < n; p++)
            for (int q = 0; q < n; q++) {
                if (p == q) continue;
                if (fabs(M(p,q)) <= tol) continue;
                bool pTied = tied.find(std::make_pair(idxNode[p], idxDof[p])) != tied.end();
                bool qTied = tied.find(std::make_pair(idxNode[q], idxDof[q])) != tied.end();
                if (pTied || qTied) {
                    opserr << "LadrunoProjectionHandler - element " << elePtr->getTag()
                           << " has a consistent (off-diagonal) mass coupling a tied DOF "
                              "(node " << idxNode[p] << " dof " << idxDof[p] << " and node "
                           << idxNode[q] << " dof " << idxDof[q] << "). The projection handler "
                              "requires lumped mass on tied DOFs; use the lumped mass of the "
                              "element, or constraints Transformation with a banded system.\n";
                    return -1;
                }
            }
    }

    for (std::map<std::pair<int,int>,double>::iterator it = tied.begin(); it != tied.end(); ++it) {
        if (unmappedNode.find(it->first.first) != unmappedNode.end()) continue;
        if (it->second <= 0.0) {
            opserr << "LadrunoProjectionHandler - massless tied DOF (node " << it->first.first
                   << " dof " << it->first.second << "): every DOF of an MP or EQ constraint "
                      "keeps its own equation under constraints LadrunoProjection and needs a "
                      "nonzero lumped mass. Add translational mass, and for rotational ties "
                      "(rigidLink -beam, rigidDiaphragm, equalDOF on rotations) a small "
                      "rotational mass (about 0.01-0.1% of the translational mass of the "
                      "node), or use constraints Transformation or Penalty.\n";
            return -1;
        }
    }
    return 0;
}

// Equation numbers are known: build the projector and pass it to the integrator.
int
LadrunoProjectionHandler::doneNumberingDOF(void)
{
    // this override replaces the base one, which passes the equation numbers
    // to the FE_Elements; call it first
    this->ConstraintHandler::doneNumberingDOF();

    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0)
        return -1;

    if (theProjector != 0) { delete theProjector; theProjector = 0; }
    theProjector = new LadrunoConstraintProjector();
    theProjector->setICPolicy(projectICs, icTol);

    for (int gi = 0; gi < (int)theGroups.size(); gi++) {
        RawGroup &g = theGroups[gi];

        // retained DOFs with an equation become the columns of L; fixed and
        // prescribed retained DOFs are dropped (prescribed ones only drive
        // the constrained DOFs already removed by classifyDerivedSlaves())
        std::map<int, int> colOf;
        std::vector<int> survMasterVtx;
        for (int m = 0; m < (int)g.masterVtx.size(); m++) {
            int v = g.masterVtx[m];
            Node *np = theDomain->getNode(vtxNode[v]);
            int e = (np != 0) ? np->getDOF_GroupPtr()->getID()(vtxDof[v]) : -1;
            if (e >= 0) {
                colOf[v] = (int)survMasterVtx.size();
                survMasterVtx.push_back(v);
            }
        }
        int nRet = (int)survMasterVtx.size();

        // constrained rows with at least one remaining column, and the ones
        // whose retained DOFs are all fixed (acceleration set to zero)
        std::vector<int> rowSlaveEqn;
        std::vector<std::vector<int> > rowCols;
        std::vector<std::vector<double> > rowCoef;
        std::vector<double> rowDelta;
        std::vector<int> fixedEqns;
        for (int s = 0; s < (int)g.slaves.size(); s++) {
            SlaveRec &sr = g.slaves[s];
            Node *np = theDomain->getNode(vtxNode[sr.vtx]);
            int ec = (np != 0) ? np->getDOF_GroupPtr()->getID()(vtxDof[sr.vtx]) : -1;
            if (ec < 0) {
                if (derivedKey.find(std::make_pair(vtxNode[sr.vtx], vtxDof[sr.vtx]))
                    != derivedKey.end())
                    continue;
                if (prescribedKey.find(std::make_pair(vtxNode[sr.vtx], vtxDof[sr.vtx]))
                    != prescribedKey.end())
                    opserr << "LadrunoProjectionHandler - PRESCRIBED CONSTRAINED DOF: DOF "
                              "(node " << vtxNode[sr.vtx] << " dof " << vtxDof[sr.vtx] << ") "
                              "has a prescribed motion (non-homogeneous SP or imposedMotion) "
                              "and is constrained by an MP or EQ constraint. Remove the "
                              "constraint or the prescribed motion.\n";
                else
                    opserr << "LadrunoProjectionHandler - FIXED CONSTRAINED DOF: DOF (node "
                           << vtxNode[sr.vtx] << " dof " << vtxDof[sr.vtx] << ") is constrained "
                              "and also fixed. Remove the fix or the constraint.\n";
                return -2;
            }
            std::vector<int> cols; std::vector<double> coef;
            for (int k = 0; k < (int)sr.masterVtx.size(); k++) {
                std::map<int,int>::iterator it = colOf.find(sr.masterVtx[k]);
                if (it != colOf.end()) { cols.push_back(it->second); coef.push_back(sr.coeff[k]); }
            }
            if (cols.empty()) {
                fixedEqns.push_back(ec);
            } else {
                rowSlaveEqn.push_back(ec);
                rowCols.push_back(cols);
                rowCoef.push_back(coef);
                rowDelta.push_back(sr.delta);
            }
        }

        int nCon = (int)rowSlaveEqn.size();
        if (nRet == 0 && nCon == 0) {
            if (!fixedEqns.empty()) {
                ID fe((int)fixedEqns.size());
                for (int i = 0; i < (int)fixedEqns.size(); i++) fe(i) = fixedEqns[i];
                Matrix L0(0, 0); Vector d0(0); ID a0(0);
                theProjector->addGroup(a0, 0, L0, d0, fe);
            }
            continue;
        }
        if (nRet == 0) {
            // rows without a remaining column were routed to fixedEqns above
            opserr << "LadrunoProjectionHandler - group " << gi
                   << " has constrained rows but no free retained DOF.\n";
            return -3;
        }

        // L = [I; C] with rows retained first, the offsets and the zeroed rows
        int nrows = nRet + nCon;
        ID allEqn(nrows);
        Matrix L(nrows, nRet);
        L.Zero();
        Vector delta(nCon > 0 ? nCon : 1);
        for (int p = 0; p < nRet; p++) {
            int v = survMasterVtx[p];
            allEqn(p) = theDomain->getNode(vtxNode[v])->getDOF_GroupPtr()->getID()(vtxDof[v]);
            L(p, p) = 1.0;
        }
        for (int j = 0; j < nCon; j++) {
            allEqn(nRet + j) = rowSlaveEqn[j];
            for (int k = 0; k < (int)rowCols[j].size(); k++)
                L(nRet + j, rowCols[j][k]) = rowCoef[j][k];
            delta(j) = rowDelta[j];
        }
        ID fe(fixedEqns.empty() ? 0 : (int)fixedEqns.size());
        for (int i = 0; i < (int)fixedEqns.size(); i++) fe(i) = fixedEqns[i];

        theProjector->addGroup(allEqn, nRet, L, delta, fe);

        if (verbose)
            opserr << "  LadrunoProjection group " << gi << ": " << nRet
                   << " retained DOF(s), " << nCon << " constrained, "
                   << fe.Size() << " constrained DOF(s) with fixed retained DOFs.\n";
    }

    // The previous projector was deleted above, so the integrator is always
    // given the new pointer, null when there is no group, never left with a
    // dangling one.
    Integrator *integ = this->getIntegratorPtr();
    LadrunoProjectionConsumer *consumer = dynamic_cast<LadrunoProjectionConsumer *>(integ);
    if (consumer == 0) {
        // without constraints any integrator is acceptable
        if (theProjector->numGroups() > 0) {
            opserr << "LadrunoProjectionHandler::doneNumberingDOF() - this handler requires "
                      "an explicit integrator that applies the projection "
                      "(CentralDifferenceLadruno or ExplicitBathe). Use constraints "
                      "Transformation or Penalty for implicit analyses.\n";
            return -4;
        }
    } else {
        consumer->setConstraintProjector(theProjector->numGroups() > 0 ? theProjector : 0);
    }
    return 0;
}

// A retained rotation entering a constrained translation through C is a
// lever arm frozen at the initial geometry (rigidLink -beam, rigidDiaphragm).
// A rotation tied to a rotation, or a translation tie, has no lever arm.
void
LadrunoProjectionHandler::flagRotMonitor(int masterNode, int masterDof,
                                         int slaveNode, int slaveDof)
{
    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0) return;
    Node *mp = theDomain->getNode(masterNode);
    Node *sp = theDomain->getNode(slaveNode);
    if (mp == 0 || sp == 0) return;
    int ndmM = mp->getCrds().Size();
    int ndfM = mp->getNumberDOF();
    int ndmS = sp->getCrds().Size();
    // 2D: dof 2 is the rotation; 3D with 6 DOFs: dofs 3, 4, 5;
    // translations are 0 .. ndm-1
    bool masterRot  = (ndmM == 2 && masterDof >= 2) || (ndmM == 3 && ndfM >= 6 && masterDof >= 3);
    bool slaveTrans = (slaveDof >= 0 && slaveDof < ndmS);
    if (!masterRot || !slaveTrans) return;
    for (int i = 0; i < (int)rotMonitors.size(); i++)
        if (rotMonitors[i].nodeTag == masterNode) return;
    RotMonitor m;
    m.nodeTag = masterNode;
    m.warned = false;
    if (ndmM == 2) { m.rotDof.push_back(2); }
    else           { m.rotDof.push_back(3); m.rotDof.push_back(4); m.rotDof.push_back(5); }
    const Vector &u = mp->getTrialDisp();
    for (int k = 0; k < (int)m.rotDof.size(); k++) {
        int d = m.rotDof[k];
        m.theta0.push_back((d < u.Size()) ? u(d) : 0.0);
    }
    rotMonitors.push_back(m);
}

// Called by AnalysisModel::updateDomain(time, dt) after Domain::applyLoad(),
// so an ImposedMotionSP has already set the velocity and acceleration of its
// node. The prescribed DOFs have no equation, so the integrator's
// setResponse() leaves the values imposed here unchanged.
int
LadrunoProjectionHandler::applyLoad(void)
{
    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0)
        return 0;
    for (int i = 0; i < (int)prescribedDOFs.size(); i++) {
        PrescribedDOF &p = prescribedDOFs[i];
        if (p.sp == 0)
            continue;
        Node *np = theDomain->getNode(p.nodeTag);
        if (np == 0)
            continue;
        np->setTrialDisp(p.sp->getValue() + p.sp->getInitialValue(), p.dof);
    }

    // constrained DOFs driven by prescribed DOFs, after those DOFs are set:
    //   u_c = sum_k C_k u_k + delta ;  v_c = sum_k C_k v_k ;  a_c = sum_k C_k a_k
    for (int i = 0; i < (int)derivedSlaves.size(); i++) {
        DerivedSlave &d = derivedSlaves[i];
        Node *snp = theDomain->getNode(d.nodeTag);
        if (snp == 0)
            continue;
        double uc = d.delta, vc = 0.0, ac = 0.0;
        for (int k = 0; k < (int)d.masterNode.size(); k++) {
            Node *mnp = theDomain->getNode(d.masterNode[k]);
            if (mnp == 0)
                continue;
            int md = d.masterDof[k];
            const Vector &mu = mnp->getTrialDisp();
            const Vector &mv = mnp->getTrialVel();
            const Vector &ma = mnp->getTrialAccel();
            double c = d.coeff[k];
            if (md >= 0 && md < mu.Size()) uc += c * mu(md);
            if (md >= 0 && md < mv.Size()) vc += c * mv(md);
            if (md >= 0 && md < ma.Size()) ac += c * ma(md);
        }
        int snd = snp->getNumberDOF();
        if (d.dof < 0 || d.dof >= snd)
            continue;
        snp->setTrialDisp(uc, d.dof);
        Vector vv = snp->getTrialVel();  vv(d.dof) = vc;  snp->setTrialVel(vv);
        Vector av = snp->getTrialAccel(); av(d.dof) = ac; snp->setTrialAccel(av);
    }

    // note (once per node) when a frozen lever arm is no longer small-rotation
    for (int i = 0; i < (int)rotMonitors.size(); i++) {
        RotMonitor &m = rotMonitors[i];
        if (m.warned)
            continue;
        Node *np = theDomain->getNode(m.nodeTag);
        if (np == 0)
            continue;
        const Vector &u = np->getTrialDisp();
        double drift = 0.0;
        for (int k = 0; k < (int)m.rotDof.size(); k++) {
            int dd = m.rotDof[k];
            double th = (dd < u.Size()) ? u(dd) : 0.0;
            double a = fabs(th - m.theta0[k]);
            if (a > drift) drift = a;
        }
        if (drift > LADRUNO_CCR_STALE_RAD) {
            opserr << "LadrunoProjectionHandler - NOTE retained node " << m.nodeTag
                   << " has rotated " << drift << " rad since the constraints were built; "
                      "the rigidLink -beam / rigidDiaphragm lever arm is frozen at small "
                      "rotation and the tie is no longer exact. For finite rotation use a "
                      "kinematic coupling element or constraints Transformation.\n";
            m.warned = true;
        }
    }
    return 0;
}

void
LadrunoProjectionHandler::clearAll(void)
{
    Domain *theDomain = this->getDomainPtr();
    if (theDomain != 0) {
        NodeIter &theNod = theDomain->getNodes();
        Node *nodPtr;
        while ((nodPtr = theNod()) != 0)
            nodPtr->setDOF_GroupPtr(0);
    }
    theGroups.clear();
    vtxNode.clear();
    vtxDof.clear();
    prescribedDOFs.clear();
    prescribedKey.clear();
    derivedSlaves.clear();
    derivedKey.clear();
    rotMonitors.clear();
    // the projector stays owned here; doneNumberingDOF() replaces it and
    // rebinds the integrator
}

int
LadrunoProjectionHandler::sendSelf(int cTag, Channel &theChannel)
{
    static ID data(2);
    data(0) = projectICs ? 1 : 0;
    data(1) = verbose ? 1 : 0;
    theChannel.sendID(0, cTag, data);
    static Vector ddata(1);
    ddata(0) = icTol;
    theChannel.sendVector(0, cTag, ddata);
    return 0;
}

int
LadrunoProjectionHandler::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static ID data(2);
    theChannel.recvID(0, cTag, data);
    projectICs = (data(0) == 1);
    verbose = (data(1) == 1);
    static Vector ddata(1);
    theChannel.recvVector(0, cTag, ddata);
    icTol = ddata(0);
    return 0;
}
