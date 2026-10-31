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
// Description: LadrunoContactHandler, see LadrunoContactHandler.h.

#include "LadrunoContactHandler.h"
#include "LadrunoContactFE.h"

#include <LadrunoContactDomain.h>
#include <LadrunoContactSurface.h>
#include <LadrunoContactProjection.h>
#include <LadrunoContactBucketSort.h>
#include <Matrix.h>
#include <Domain.h>
#include <AnalysisModel.h>
#include <Integrator.h>
#include <Node.h>
#include <Vector.h>
#include <NodeIter.h>
#include <DOF_Group.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <EQ_Constraint.h>
#include <EQ_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <FE_Element.h>
#include <Subdomain.h>
#include <ID.h>
#include <classTags.h>
#include <map>
#include <vector>
#include <cmath>

// A declared contact that cannot be built is an error, never a skip: a model
// whose contact silently does not act converges, balances its reactions and
// is wrong. Every refusal in handle() returns through here, and the analysis
// reports the failed domainChanged().
static int
ladrunoContactFatal(void)
{
    return -1;
}

// Pre-flight one contact surface: every node must exist and carry 3D
// coordinates (getCrds()(2) is unchecked in release builds, so a 2D node would
// read out of bounds).
static bool
ladrunoSurfaceNodesOk(Domain *dom, int ctag, const LadrunoContactSurface *s)
{
    const ID &tags = s->getNodeTags();
    for (int i = 0; i < tags.Size(); i++) {
        Node *n = dom->getNode(tags(i));
        if (n == 0) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": node " << tags(i) << " of contactSurface " << s->getTag()
                   << " does not exist\n";
            return false;
        }
        int nd = n->getCrds().Size();
        if (nd != 3) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": node " << tags(i) << " of contactSurface " << s->getTag()
                   << " has " << nd << "D coordinates; node-to-segment contact needs 3D nodes\n";
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// command factory: constraints LadrunoContact
// ----------------------------------------------------------------------------
void *OPS_LadrunoContactHandler(void)
{
    return new LadrunoContactHandler();
}

LadrunoContactHandler::LadrunoContactHandler()
  : ConstraintHandler(HANDLER_TAG_LadrunoContactHandler)
{
}

LadrunoContactHandler::~LadrunoContactHandler()
{
}

// Automatic penalty for ONE (slave, master segment) pair from the initial
// stiffness of the solid element that owns the segment, read through the
// Element interface (getExternalNodes, getNumDOF, getInitialStiff), so any
// solid with 3 DOFs per node works:
//   kn = f_si * mean over the segment nodes of ( n^T K_node n ),
// K_node = the 3x3 diagonal block of the element stiffness at that node and n
// the reference segment normal. This is the element-stiffness form of the
// LS-DYNA penalty k = f_si K A^2 / V (Hallquist 2006, Eq. 26.14a): the diagonal
// stiffness of a solid is of order E A^2 / V, so the geometry factor is
// contained in the assembled matrix. f_si = 0.10, the LS-DYNA SLSFAC default.
// Returns <= 0 when no owning 3-DOF-per-node solid is found or the reference
// normal is ambiguous.
static double
ladrunoResolveAutoKn(Domain *theDomain, Node **segNodes, int nps,
                     const double orientDir[3])
{
    const double F_SI = 0.10;
    // reference (undeformed) segment node coordinates
    double Xref[4][3];
    for (int k = 0; k < nps; k++) {
        const Vector &Xk = segNodes[k]->getCrds();
        for (int d = 0; d < 3; d++) Xref[k][d] = Xk(d);
    }
    // reference outward normal at the segment centre, oriented by orientDir
    double xiC  = (nps == 4) ? 0.0 : (1.0 / 3.0);
    double etaC = (nps == 4) ? 0.0 : (1.0 / 3.0);
    double n[3];
    if (!LadrunoContactProjection::normalOriented(nps, xiC, etaC, Xref, orientDir, n))
        return -1.0;
    int segTags[4];
    for (int k = 0; k < nps; k++) segTags[k] = segNodes[k]->getTag();
    // owning solid = the first element with 3 DOFs per node that contains all
    // the segment nodes
    ElementIter &theEle = theDomain->getElements();
    Element *e;
    while ((e = theEle()) != 0) {
        const ID &en = e->getExternalNodes();
        int nn = en.Size();
        if (nn <= 0 || e->getNumDOF() != 3 * nn) continue;   // not a 3-DOF/node solid
        int loc[4]; bool allIn = true;
        for (int k = 0; k < nps; k++) {
            loc[k] = en.getLocation(segTags[k]);
            if (loc[k] < 0) { allIn = false; break; }
        }
        if (!allIn) continue;
        const Matrix &K = e->getInitialStiff();
        if (K.noRows() != 3 * nn || K.noCols() != 3 * nn) return -1.0;
        double acc = 0.0;
        for (int k = 0; k < nps; k++) {
            int c = 3 * loc[k];
            double Kn[3];
            for (int i = 0; i < 3; i++)
                Kn[i] = K(c + i, c + 0) * n[0] + K(c + i, c + 1) * n[1] + K(c + i, c + 2) * n[2];
            acc += n[0] * Kn[0] + n[1] * Kn[1] + n[2] * Kn[2];
        }
        double kn = F_SI * acc / nps;
        return (kn > 0.0) ? kn : -1.0;
    }
    return -1.0;   // no owning solid element
}

int
LadrunoContactHandler::handle(const ID *nodesLast)
{
    Domain *theDomain = this->getDomainPtr();
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    Integrator *theIntegrator = this->getIntegratorPtr();
    if (theDomain == 0 || theModel == 0 || theIntegrator == 0) {
        opserr << "WARNING LadrunoContactHandler::handle() - setLinks() not called\n";
        return -1;
    }

    // handle() runs again at every domainChanged(); an analysis that removes
    // nodes or elements shrinks the domain under a contact declared once.
    // Prune the removed nodes and retire interactions left without a surface
    // before anything reads a surface.
    if (LadrunoContactDomain *cdPrune = theDomain->getLadrunoContactDomain())
        cdPrune->pruneRemovedNodes(theDomain);

    // MP constraints (equalDOF, rigidLink, rigidDiaphragm) are not enforced
    // by this handler (Plain-style numbering).
    MP_ConstraintIter &theMPcheck = theDomain->getMPs();
    if (theMPcheck() != 0)
        opserr << "WARNING LadrunoContactHandler::handle() - MP constraints (e.g. equalDOF) "
                  "are present but NOT enforced by the contact handler; use a model with "
                  "SP (fix) constraints only.\n";

    // Collect the SPs (domain and load pattern) by node. As in PlainHandler a
    // constrained DOF is removed, so a non-homogeneous (imposed displacement)
    // value cannot be enforced; warn as PlainHandler does. A displacement-driven
    // contact model can use the DisplacementControl integrator instead.
    std::multimap<int, SP_Constraint *> allSPs;
    SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP;
    while ((theSP = theSPs()) != 0) {
        if (theSP->isHomogeneous() == false)
            opserr << "WARNING LadrunoContactHandler::handle() - non-homogeneous SP (imposed "
                      "displacement) at node " << theSP->getNodeTag() << " is NOT enforced by "
                      "the contact handler (homogeneous constraint assumed). Use the "
                      "DisplacementControl integrator for displacement-driven contact.\n";
        allSPs.insert(std::make_pair(theSP->getNodeTag(), theSP));
    }

    // DOF_Groups: free = -2, single-point constrained = -1.
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    int numDOF = 0, countDOF = 0, count3 = 0;
    while ((nodPtr = theNod()) != 0) {
        DOF_Group *dofPtr = new DOF_Group(numDOF++, nodPtr);
        if (dofPtr == 0) {
            opserr << "WARNING LadrunoContactHandler::handle() - out of memory (DOF_Group)\n";
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
            if (dof >= 0 && dof < cid.Size() && cid(dof) == -2) {
                dofPtr->setID(dof, -1);
                countDOF--;
            } else {
                opserr << "WARNING LadrunoContactHandler::handle() - invalid/duplicate SP at DOF "
                       << dof << " node " << nodeID << endln;
            }
        }
        // EQ_Constraints as in PlainHandler: a trivial identity (a single entry
        // equal to 1.0) marks the DOF -4; any other one is reported and not
        // enforced, never dropped without a message.
        {
            EQ_ConstraintIter &theEQs = theDomain->getEQs();
            EQ_Constraint *eqPtr;
            while ((eqPtr = theEQs()) != 0) {
                if (eqPtr->getNodeConstrained() != nodeID)
                    continue;
                if (eqPtr->isTimeVarying() == true)
                    opserr << "WARNING LadrunoContactHandler::handle() - time-varying "
                              "EQ_Constraint for node " << nodeID << "; non-varying assumed\n";
                const Vector &C = eqPtr->getConstraint();
                if (C.Size() > 1 || C(0) != 1.0) {
                    opserr << "WARNING LadrunoContactHandler::handle() - EQ_Constraint at node "
                           << nodeID << " is not a trivial identity: NOT ENFORCED by the "
                              "contact handler (contact together with general equation "
                              "constraints is not supported)\n";
                } else {
                    int dof = eqPtr->getConstrainedDOFs();
                    const ID &cid = dofPtr->getID();
                    if (dof >= 0 && dof < cid.Size() && cid(dof) == -2) {
                        dofPtr->setID(dof, -4);
                        countDOF--;
                    } else {
                        opserr << "WARNING LadrunoContactHandler::handle() - EQ_Constraint DOF "
                               << dof << " already constrained at node " << nodeID << endln;
                    }
                }
            }
        }
        nodPtr->setDOF_GroupPtr(dofPtr);
        theModel->addDOF_Group(dofPtr);
    }
    theModel->setNumEqn(countDOF);

    // nodesLast (subdomain boundary): mark -3, as PlainHandler
    if (nodesLast != 0) {
        for (int i = 0; i < nodesLast->Size(); i++) {
            Node *np = theDomain->getNode((*nodesLast)(i));
            if (np != 0) {
                DOF_Group *dp = np->getDOF_GroupPtr();
                const ID &id = dp->getID();
                for (int j = 0; j < id.Size(); j++)
                    if (id(j) == -2) { dp->setID(j, -3); count3++; }
            }
        }
    }

    // FE_Elements for every Domain element, as PlainHandler. numFe continues
    // into the contact adapter tags.
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

    // --- contact adapters; none without a contact engine (then the handler
    //     is a PlainHandler) ---
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();

    // Node-to-segment penalty contact: for each definition (MASTER_SEGMENTS
    // master, SLAVE_NODES slave) one adapter per (slave node, candidate master
    // segment) pair from the bucket-sort broad phase. A pair that does not
    // project or does not penetrate gives zero force inside the adapter.
    if (cd != 0) {
        cd->clearNtsForce();     // the force readout keys are rebuilt below
        cd->frictionGCBegin();   // friction slots of pairs not rebuilt below are collected

        // A slave node in the slave sets of two contacts gets two adapters and
        // the tractions add up; report it.
        {
            std::multimap<int, int> slaveSeen;   // slaveNodeTag -> contactTag
            for (int c = 0; c < cd->getNumContacts(); c++) {
                const LadrunoContactDomain::Contact &ct = cd->getContact(c);
                if (ct.retired) continue;
                LadrunoContactSurface *ss = cd->getSurface(ct.slaveSurfTag);
                if (ss == 0 || ss->getKind() != LadrunoContactSurface::SLAVE_NODES)
                    continue;
                const ID &st = ss->getNodeTags();
                for (int i = 0; i < st.Size(); i++) {
                    int tag = st(i);
                    if (slaveSeen.find(tag) != slaveSeen.end())
                        opserr << "WARNING LadrunoContactHandler::handle() - slave node "
                               << tag << " appears in more than one contact's slave set; "
                                  "its contact tractions will be ADDED (check the model)\n";
                    slaveSeen.insert(std::make_pair(tag, ct.tag));
                }
            }
        }

        for (int c = 0; c < cd->getNumContacts(); c++) {
            const LadrunoContactDomain::Contact &ct = cd->getContact(c);
            if (ct.retired) continue;
            LadrunoContactSurface *ms = cd->getSurface(ct.masterSurfTag);
            LadrunoContactSurface *ss = cd->getSurface(ct.slaveSurfTag);
            if (ms == 0 || ss == 0) {
                opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": master/slave surface not defined\n";
                return ladrunoContactFatal();
            }
            if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS ||
                ss->getKind() != LadrunoContactSurface::SLAVE_NODES) {
                opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": needs a -master segment surface and a -slave node surface\n";
                return ladrunoContactFatal();
            }
            if (!ladrunoSurfaceNodesOk(theDomain, ct.tag, ms) ||
                !ladrunoSurfaceNodesOk(theDomain, ct.tag, ss))
                return ladrunoContactFatal();
            int nps = ms->getNodesPerSeg();
            if (nps < 3 || nps > 4) {
                opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": nodesPerSeg " << nps << " unsupported (need 3 or 4)\n";
                return ladrunoContactFatal();
            }
            if (!ct.knAuto && ct.kn <= 0.0) {
                // kn == 0 would build adapters that transmit nothing
                opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                       << ": segment contact needs kn > 0 or kn auto (got " << ct.kn << ")\n";
                return ladrunoContactFatal();
            }
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSeg = mTags.Size() / nps;

            // Broad phase: bucket-sort the master segments (reference
            // coordinates) once; each slave queries only its 27 neighbour
            // buckets. The kept set is a superset of every near pair, so the
            // result equals brute force; a huge -cell (one bucket) is brute
            // force. A missing master node (malformed model) is filled with the
            // segment's first valid coordinate, which keeps the superset; the
            // narrow loop refuses such a pair anyway.
            std::vector<double> segCoords((size_t)nSeg * nps * 3, 0.0);
            for (int seg = 0; seg < nSeg; seg++) {
                double *S = &segCoords[(size_t)seg * nps * 3];
                bool haveFirst = false; double first[3] = {0.0, 0.0, 0.0};
                for (int k = 0; k < nps; k++) {
                    Node *mn = theDomain->getNode(mTags(seg * nps + k));
                    if (mn != 0) {
                        const Vector &X = mn->getCrds();
                        S[k*3+0] = X(0); S[k*3+1] = X(1); S[k*3+2] = X(2);
                        if (!haveFirst) { first[0]=X(0); first[1]=X(1); first[2]=X(2); haveFirst = true; }
                    } else {
                        S[k*3+0] = S[k*3+1] = S[k*3+2] = HUGE_VAL;   // mark missing
                    }
                }
                for (int k = 0; k < nps; k++)   // backfill missing with the first valid one
                    if (S[k*3+0] == HUGE_VAL)
                        for (int d = 0; d < 3; d++) S[k*3+d] = first[d];
            }
            LadrunoContactBucketSort::Grid grid(nSeg, nps, segCoords.data(), ct.cellFrac, 1.0);
            std::vector<int> cand(nSeg);

            for (int si = 0; si < sTags.Size(); si++) {
                Node *sn = theDomain->getNode(sTags(si));
                if (sn == 0 || sn->getNumberDOF() < 3) {
                    // the adapter takes the first three DOFs as the translations;
                    // ndf > 3 is allowed (the extra DOFs are not coupled)
                    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                           << " slave node " << sTags(si) << " ndf="
                           << (sn != 0 ? sn->getNumberDOF() : 0)
                           << " < 3 (node-to-segment contact needs the first three DOFs "
                              "to be the translations)\n";
                    return ladrunoContactFatal();
                }
                const Vector &Xs0 = sn->getCrds();
                double slavePt[3] = { Xs0(0), Xs0(1), Xs0(2) };
                int nCand = grid.candidates(slavePt, cand.data());
                for (int ci = 0; ci < nCand; ci++) {
                    int seg = cand[ci];
                    Node *segNodes[4]; int badTag = 0, badNdf = 0;
                    for (int k = 0; k < nps; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * nps + k));
                        if (mn == 0 || mn->getNumberDOF() < 3) {
                            badTag = mTags(seg * nps + k);
                            badNdf = (mn != 0) ? mn->getNumberDOF() : -1;
                            break;
                        }
                        segNodes[k] = mn;
                    }
                    if (badTag != 0) {
                        // a master node shared with a beam or shell passes the
                        // surface pre-flight; refuse it here rather than drop
                        // every pair on its segments
                        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                               << " master segment " << seg << " node " << badTag << ": ndf="
                               << badNdf << " < 3 (or node absent); master segments need "
                                  "three translational DOFs per node\n";
                        return ladrunoContactFatal();
                    }
                    // orientation direction toward the slave's allowed half-space:
                    // the explicit -outward if given, else (slave reference) minus
                    // (segment reference centroid). The derived normal is flipped
                    // to n.orientDir > 0. Use -outward for a slave that starts
                    // penetrated.
                    double orientDir[3];
                    if (ct.hasOutward) {
                        for (int d = 0; d < 3; d++) orientDir[d] = ct.outward[d];
                    } else {
                        double cen[3] = {0.0, 0.0, 0.0};
                        for (int k = 0; k < nps; k++) {
                            const Vector &Xk = segNodes[k]->getCrds();
                            for (int d = 0; d < 3; d++) cen[d] += Xk(d);
                        }
                        const Vector &Xs = sn->getCrds();
                        for (int d = 0; d < 3; d++) orientDir[d] = Xs(d) - cen[d] / nps;
                    }
                    double knUse = ct.kn;
                    if (ct.knAuto) {
                        knUse = ladrunoResolveAutoKn(theDomain, segNodes, nps, orientDir);
                        if (knUse <= 0.0) {
                            // a per-pair skip would leave a hole in the interface
                            opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                   << ct.tag << " slave node " << sTags(si) << " segment "
                                   << seg << ": -kn auto could not size a penalty (no owning "
                                      "solid element, an element without 3 DOFs per node, or "
                                      "an ambiguous normal). Give an explicit kn.\n";
                            return ladrunoContactFatal();
                        }
                    }
                    LadrunoContactFE *fe =
                        new LadrunoContactFE(numFe++, sn, segNodes, nps, knUse, orientDir,
                                             ct.kt, ct.mu, theDomain, ct.tag, seg,
                                             ct.consistentTan, ct.consistentNormal,
                                             ct.cohesion, ct.tauMax, ct.muc);
                    if (fe == 0) return -5;
                    theModel->addFE_Element(fe);
                    if (ct.hasFriction())
                        cd->frictionGCMark(ct.tag, sTags(si), seg);  // live this handle()
                }
            }
        }

        cd->frictionGCEnd();     // erase the friction slots no live pair referenced
    }

    // Rigid analytical planes: one adapter per slave node (connectivity = that
    // node); ndm from the node coordinates.
    if (cd != 0) {
        for (int p = 0; p < cd->getNumRigidPlanes(); p++) {
            const LadrunoContactDomain::RigidPlane &rp = cd->getRigidPlane(p);
            if (rp.retired) continue;
            LadrunoContactSurface *surf = cd->getSurface(rp.slaveSurfTag);
            if (surf == 0) {
                opserr << "FATAL LadrunoContactHandler::handle() - contactPlane " << rp.tag
                       << ": slave surface " << rp.slaveSurfTag << " not defined\n";
                return ladrunoContactFatal();
            }
            const ID &nodeTags = surf->getNodeTags();
            for (int k = 0; k < nodeTags.Size(); k++) {
                Node *sn = theDomain->getNode(nodeTags(k));
                if (sn == 0) {
                    opserr << "FATAL LadrunoContactHandler::handle() - rigid-plane slave node "
                           << nodeTags(k) << " not in domain\n";
                    return ladrunoContactFatal();
                }
                int nd = sn->getCrds().Size();
                if (nd != 3 || sn->getNumberDOF() < nd) {
                    // the adapter couples the first nd translational DOFs; fewer
                    // DOFs than coordinates would leave myID slots unset
                    opserr << "FATAL LadrunoContactHandler::handle() - rigid-plane slave node "
                           << nodeTags(k) << " has ndm=" << nd << ", ndf=" << sn->getNumberDOF()
                           << " (needs 3D coordinates and ndf >= 3)\n";
                    return ladrunoContactFatal();
                }
                LadrunoContactFE *fe = new LadrunoContactFE(numFe++, sn, nd, rp.p0, rp.n, rp.kn,
                                                            theDomain, rp.tag, rp.muc);
                if (fe == 0) return -5;
                theModel->addFE_Element(fe);
            }
        }
    }

    return count3;
}

void
LadrunoContactHandler::clearAll(void)
{
    // reset the node DOF_Group pointers, as PlainHandler; the AnalysisModel owns
    // and deletes the DOF_Groups and FE_Elements, contact adapters included
    Domain *theDomain = this->getDomainPtr();
    if (theDomain == 0)
        return;
    NodeIter &theNod = theDomain->getNodes();
    Node *nodPtr;
    while ((nodPtr = theNod()) != 0)
        nodPtr->setDOF_GroupPtr(0);
}

int
LadrunoContactHandler::sendSelf(int, Channel &)
{
    // The handler holds no state: the contact definitions live on the Domain's
    // contact engine.
    return 0;
}

int
LadrunoContactHandler::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}
