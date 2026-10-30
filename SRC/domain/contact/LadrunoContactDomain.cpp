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
// Description: LadrunoContactDomain, see LadrunoContactDomain.h.

#include "LadrunoContactDomain.h"
#include "LadrunoContactSurface.h"
#include <OPS_Globals.h>
#include <cmath>
#include <algorithm>
#include <Domain.h>
#include <Node.h>
#include <ID.h>

LadrunoContactDomain::LadrunoContactDomain()
  : numCommits(0), numReverts(0)
{
}

LadrunoContactDomain::~LadrunoContactDomain()
{
    for (size_t i = 0; i < theSurfaces.size(); i++)
        delete theSurfaces[i];
    theSurfaces.clear();
    theContacts.clear();
}

int
LadrunoContactDomain::addSurface(LadrunoContactSurface *surf)
{
    if (surf == 0)
        return -1;
    if (getSurface(surf->getTag()) != 0) {
        opserr << "WARNING LadrunoContactDomain::addSurface() - surface tag "
               << surf->getTag() << " already exists\n";
        return -1;
    }
    // The handler derives nSeg = size/nodesPerSeg, so a faceted surface must
    // hold a whole number of segments; a trailing partial segment would be
    // dropped without a message.
    if (surf->getKind() != LadrunoContactSurface::SLAVE_NODES) {
        int nps = surf->getNodesPerSeg();
        int n   = surf->getNodeTags().Size();
        if (nps < 3 || nps > 4 || (n % nps) != 0) {
            opserr << "WARNING LadrunoContactDomain::addSurface() - surface tag "
                   << surf->getTag() << ": " << n << " connectivity entries is not a whole "
                      "number of " << nps << "-node segments (nodesPerSeg must be 3 or 4)\n";
            return -1;
        }
    }
    theSurfaces.push_back(surf);
    return 0;
}

LadrunoContactSurface *
LadrunoContactDomain::getSurface(int tag) const
{
    for (size_t i = 0; i < theSurfaces.size(); i++)
        if (theSurfaces[i]->getTag() == tag)
            return theSurfaces[i];
    return 0;
}

bool
LadrunoContactDomain::warnOnce(int contactTag, int topic)
{
    return theWarnOnce.insert(std::make_pair(contactTag, topic)).second;
}

bool
LadrunoContactDomain::contactTagInUse(int tag) const
{
    for (size_t i = 0; i < theContacts.size(); i++)
        if (theContacts[i].tag == tag) return true;
    for (size_t i = 0; i < theRigidPlanes.size(); i++)
        if (theRigidPlanes[i].tag == tag) return true;
    for (size_t i = 0; i < theMortarContacts.size(); i++)
        if (theMortarContacts[i].tag == tag) return true;
    return false;
}

int
LadrunoContactDomain::addMortarContact(int tag, int masterSurfTag, int slaveSurfTag,
                                       double kn, bool knAuto, double epsN, bool epsNAuto,
                                       double augTol, int maxAug, int ngp,
                                       const double *outward, bool isTie, double mu,
                                       double epsT, bool epsTAuto, double cohesion,
                                       double tauMax, bool consistentTan)
{
    LadrunoContactSurface *ms = getSurface(masterSurfTag);
    LadrunoContactSurface *ss = getSurface(slaveSurfTag);
    if (ms == 0 || ss == 0) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - master/slave surface "
                  "not defined (master " << masterSurfTag << ", slave " << slaveSurfTag << ")\n";
        return -1;
    }
    if (contactTagInUse(tag)) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - contact tag " << tag
               << " already exists (contact tags key the pair state)\n";
        return -1;
    }
    // D is integrated over slave facets and M against master facets, so both
    // surfaces must be faceted
    if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - master surface "
               << masterSurfTag << " must be a -master facet surface\n";
        return -1;
    }
    if (ss->getKind() != LadrunoContactSurface::SLAVE_SEGMENTS) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - slave surface "
               << slaveSurfTag << " must be a -slave-segments facet surface; a -slave node "
                  "set carries no facet, so it cannot give the mortar matrix D\n";
        return -1;
    }
    if (kn < 0.0 || epsN < 0.0) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - kn and epsN must be >= 0 "
                  "(got kn " << kn << ", epsN " << epsN << ")\n";
        return -1;
    }
    if (mu < 0.0 || cohesion < 0.0 || epsT < 0.0) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - mu, cohesion and epsT "
                  "must be >= 0\n";
        return -1;
    }
    if (isTie && (mu > 0.0 || cohesion > 0.0 || tauMax > 0.0)) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - a mesh tie has no "
                  "friction cone; -mu/-cohesion/-tauMax are not allowed with -tie\n";
        return -1;
    }
    if (tauMax > 0.0 && mu <= 0.0 && cohesion <= 0.0) {
        opserr << "WARNING LadrunoContactDomain::addMortarContact() - -tauMax without -mu or "
                  "-cohesion gives the cone min(mu*N + c, tauMax) = 0 (free slip); for a "
                  "shear-capped bond use -cohesion <c>, optionally with -tauMax\n";
        return -1;
    }
    MortarContact m;
    m.tag = tag; m.masterSurfTag = masterSurfTag; m.slaveSurfTag = slaveSurfTag;
    m.kn = kn; m.knAuto = knAuto;
    m.epsN = epsN; m.epsNAuto = epsNAuto;
    m.augTol = (augTol > 0.0) ? augTol : 1e-8;
    m.maxAug = (maxAug > 0) ? maxAug : 20;
    m.ngp = (ngp >= 1) ? ngp : 2;
    m.hasOutward = (outward != 0);
    for (int d = 0; d < 3; d++) m.outward[d] = (outward != 0) ? outward[d] : 0.0;
    m.isTie = isTie;
    m.mu = mu;
    m.epsT = epsT; m.epsTAuto = epsTAuto;
    m.cohesion = cohesion;
    m.tauMax = (tauMax > 0.0) ? tauMax : 0.0;
    m.consistentTan = consistentTan;
    theMortarContacts.push_back(m);
    return 0;
}

LadrunoContactDomain::MortarNormalState &
LadrunoContactDomain::getOrCreateMortarNormalState(int contactTag, int slaveNodeTag)
{
    NodeKey k; k.c = contactTag; k.n = slaveNodeTag;
    return theMortarNormalStates[k];      // created zeroed (lambdaN = 0) when absent
}

LadrunoContactDomain::MortarFrictionState &
LadrunoContactDomain::getOrCreateMortarFrictionState(int contactTag, int slaveNodeTag,
                                                     int slaveFacet, int masterFacet)
{
    MortarPairKey k; k.c = contactTag; k.n = slaveNodeTag; k.sf = slaveFacet; k.mf = masterFacet;
    return theMortarFrictionStates[k];    // created zeroed (engaged = false) when absent
}

void
LadrunoContactDomain::mortarFrictionGCMark(int contactTag, int slaveNodeTag, int slaveFacet,
                                           int masterFacet)
{
    MortarPairKey k; k.c = contactTag; k.n = slaveNodeTag; k.sf = slaveFacet; k.mf = masterFacet;
    liveMortarPairKeys.insert(k);
}

void
LadrunoContactDomain::accumulateMortarGap(int contactTag, int slaveNodeTag, int feTag,
                                          double gtFacet, double aFacet, double epsN)
{
    // replace this adapter's earlier contribution, so the residual may be formed
    // any number of times per iteration
    PairKey fk; fk.c = contactTag; fk.s = slaveNodeTag; fk.g = feTag;
    FacetContrib &fc = theMortarFacetContribs[fk];
    NodeKey nk; nk.c = contactTag; nk.n = slaveNodeTag;
    MortarNormalState &st = theMortarNormalStates[nk];
    st.gtGlobal += gtFacet - fc.gt;
    st.aGlobal += aFacet - fc.a;
    // An automatic epsN is sized per master facet, so a shared slave node sees
    // several values; the largest is kept, which does not depend on the order
    // in which the adapters are evaluated. Reset at every handle().
    if (epsN > st.epsN) st.epsN = epsN;
    fc.gt = gtFacet; fc.a = aFacet;
}

double
LadrunoContactDomain::getMaxMortarPenetration(void) const
{
    double pen = 0.0;
    for (std::map<NodeKey, MortarNormalState>::const_iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it) {
        const MortarNormalState &st = it->second;
        if (st.isTie || st.aGlobal <= 1e-300) continue;  // a tie node has no normal gap
        double gbar = st.gtGlobal / st.aGlobal;          // < 0: penetration
        // only nodes with a compressive augmented pressure count: a node held
        // open by its multiplier carries no pressure whatever its raw gap
        if (st.lambdaN + st.epsN * gbar < 0.0 && -gbar > pen)
            pen = -gbar;
    }
    return pen;
}

void
LadrunoContactDomain::accumulateMortarTie(int contactTag, int slaveNodeTag, int feTag,
                                          const double rFacet[3], double aFacet, double epsTie)
{
    // r_I is linear in the displacements, so the sum over the facet pairs does
    // not depend on the order in which they are evaluated
    PairKey fk; fk.c = contactTag; fk.s = slaveNodeTag; fk.g = feTag;
    FacetContrib &fc = theMortarFacetContribs[fk];
    NodeKey nk; nk.c = contactTag; nk.n = slaveNodeTag;
    MortarNormalState &st = theMortarNormalStates[nk];
    for (int d = 0; d < 3; d++) {
        st.rtGlobal[d] += rFacet[d] - fc.r[d];
        fc.r[d] = rFacet[d];
    }
    st.aGlobal += aFacet - fc.a;
    fc.a = aFacet;
    st.isTie = true;
    if (epsTie > st.epsN) st.epsN = epsTie;   // largest over the facet pairs, as for the gap
}

double
LadrunoContactDomain::getMaxMortarTieResidual(void) const
{
    double res = 0.0;
    for (std::map<NodeKey, MortarNormalState>::const_iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it) {
        const MortarNormalState &st = it->second;
        if (!st.isTie || st.aGlobal <= 1e-300) continue;
        for (int d = 0; d < 3; d++) {
            double rbar = std::fabs(st.rtGlobal[d] / st.aGlobal);
            if (rbar > res) res = rbar;
        }
    }
    return res;
}

void
LadrunoContactDomain::mortarNormalGCBegin(void)
{
    liveNodeKeys.clear();
    liveMortarPairKeys.clear();
}

void
LadrunoContactDomain::mortarNormalGCMark(int contactTag, int slaveNodeTag)
{
    NodeKey k; k.c = contactTag; k.n = slaveNodeTag;
    liveNodeKeys.insert(k);
}

void
LadrunoContactDomain::mortarNormalGCEnd(void)
{
    for (std::map<NodeKey, MortarNormalState>::iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ) {
        if (liveNodeKeys.find(it->first) == liveNodeKeys.end())
            theMortarNormalStates.erase(it++);
        else
            ++it;
    }
    liveNodeKeys.clear();
    for (std::map<MortarPairKey, MortarFrictionState>::iterator it = theMortarFrictionStates.begin();
         it != theMortarFrictionStates.end(); ) {
        if (liveMortarPairKeys.find(it->first) == liveMortarPairKeys.end())
            theMortarFrictionStates.erase(it++);
        else
            ++it;
    }
    liveMortarPairKeys.clear();
    // Adapter tags are reassigned at every handle(), so the contributions are
    // dropped and the sums restart from zero; the adapters fill them again at
    // the next residual evaluation. Only the multipliers carry over.
    theMortarFacetContribs.clear();
    for (std::map<NodeKey, MortarNormalState>::iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it) {
        it->second.gtGlobal = 0.0;
        it->second.aGlobal = 0.0;
        it->second.epsN = 0.0;
        for (int d = 0; d < 3; d++) it->second.rtGlobal[d] = 0.0;
    }
}

int
LadrunoContactDomain::addRigidPlane(int tag, int slaveSurfTag,
                                    const double p0[3], const double n[3], double kn,
                                    double muc)
{
    LadrunoContactSurface *surf = getSurface(slaveSurfTag);
    if (surf == 0) {
        opserr << "WARNING LadrunoContactDomain::addRigidPlane() - slave surface "
               << slaveSurfTag << " not defined\n";
        return -1;
    }
    if (contactTagInUse(tag)) {
        opserr << "WARNING LadrunoContactDomain::addRigidPlane() - contact tag " << tag
               << " already exists (tags are unique across contact/contactPlane)\n";
        return -1;
    }
    if (surf->getKind() != LadrunoContactSurface::SLAVE_NODES) {
        // a segment surface stores repeated connectivity; used as a node set it
        // would count a shared node several times
        opserr << "WARNING LadrunoContactDomain::addRigidPlane() - surface "
               << slaveSurfTag << " is not a slave node-set\n";
        return -1;
    }
    double nrm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (nrm < 1e-14) {
        opserr << "WARNING LadrunoContactDomain::addRigidPlane() - zero normal\n";
        return -1;
    }
    // a normal that is already unit to 1e-12 is stored unchanged (dividing by
    // its recomputed norm would change the last bits)
    if (fabs(nrm - 1.0) < 1e-12)
        nrm = 1.0;
    if (kn <= 0.0) {
        // kn < 0 is attractive with a negative-definite tangent; kn == 0 is inert
        opserr << "WARNING LadrunoContactDomain::addRigidPlane() - penalty kn must be "
               << "> 0 (got " << kn << ")\n";
        return -1;
    }
    RigidPlane rp;
    rp.tag = tag; rp.slaveSurfTag = slaveSurfTag; rp.kn = kn;
    rp.muc = (muc > 0.0) ? muc : 0.0;
    for (int d = 0; d < 3; d++) { rp.p0[d] = p0[d]; rp.n[d] = n[d] / nrm; }
    theRigidPlanes.push_back(rp);
    return 0;
}

int
LadrunoContactDomain::addContact(int tag, int masterSurfTag, int slaveSurfTag,
                                 double kn, double kt, double mu, const double *outward,
                                 bool knAuto, double cellFrac, bool consistentTan,
                                 bool consistentNormal, double cohesion, double tauMax,
                                 double muc)
{
    if (getSurface(masterSurfTag) == 0 || getSurface(slaveSurfTag) == 0) {
        opserr << "WARNING LadrunoContactDomain::addContact() - master/slave surface "
                  "not defined (master " << masterSurfTag << ", slave " << slaveSurfTag << ")\n";
        return -1;
    }
    if (contactTagInUse(tag)) {
        // the tag keys the pair state stores; a duplicate would alias two interactions
        opserr << "WARNING LadrunoContactDomain::addContact() - contact tag " << tag
               << " already exists (contact tags key the pair state)\n";
        return -1;
    }
    if (kt < 0.0) {
        // a negative stick penalty flips the sign of the elastic predictor
        opserr << "WARNING LadrunoContactDomain::addContact() - friction stick penalty kt "
                  "must be >= 0 (got " << kt << ")\n";
        return -1;
    }
    if (mu < 0.0 || cohesion < 0.0) {
        opserr << "WARNING LadrunoContactDomain::addContact() - friction mu and cohesion "
                  "must be >= 0 (got mu " << mu << ", cohesion " << cohesion << ")\n";
        return -1;
    }
    if (tauMax > 0.0 && mu <= 0.0 && cohesion <= 0.0) {
        // cap = min(mu N + c, tauMax) = min(0, tauMax) = 0: free slip, almost
        // certainly not the intent (a capped bond is -cohesion with -tauMax)
        opserr << "WARNING LadrunoContactDomain::addContact() - -tauMax without -mu or "
                  "-cohesion gives the cone min(mu*N + c, tauMax) = 0 (free slip); for a "
                  "shear-capped bond use -cohesion <c>, optionally with -tauMax\n";
        return -1;
    }
    if ((mu > 0.0 || cohesion > 0.0) && kt <= 0.0) {
        opserr << "WARNING LadrunoContactDomain::addContact() - friction needs a "
                  "tangential penalty kt > 0 (contact tag kn kt mu)\n";
        return -1;
    }
    if (kn < 0.0) {
        // kn < 0 is an attractive contact with a negative-definite tangent.
        // kn == 0 is accepted here (a topology-only definition); handle() refuses
        // a segment contact without a positive or automatic penalty.
        opserr << "WARNING LadrunoContactDomain::addContact() - penalty kn must be >= 0 (got "
               << kn << ")\n";
        return -1;
    }
    Contact c;
    c.tag = tag; c.masterSurfTag = masterSurfTag; c.slaveSurfTag = slaveSurfTag;
    c.kn = kn; c.kt = kt; c.mu = mu;
    c.cohesion = (cohesion > 0.0) ? cohesion : 0.0;
    c.tauMax = (tauMax > 0.0) ? tauMax : 0.0;
    c.knAuto = knAuto;
    c.hasOutward = (outward != 0);
    for (int d = 0; d < 3; d++) c.outward[d] = (outward != 0) ? outward[d] : 0.0;
    c.cellFrac = (cellFrac > 0.0) ? cellFrac : 1.0;
    c.consistentNormal = consistentNormal;
    c.consistentTan = consistentTan;
    c.muc = (muc > 0.0) ? muc : 0.0;
    theContacts.push_back(c);
    return 0;
}

void
LadrunoContactDomain::clearNtsForce(void)
{
    theNtsForce.clear();
}

void
LadrunoContactDomain::setNtsForce(int contactTag, int slaveTag, int segIndex, double tn)
{
    PairKey k; k.c = contactTag; k.s = slaveTag; k.g = segIndex;
    theNtsForce[k] = tn;        // overwritten at every residual evaluation
}

double
LadrunoContactDomain::getNtsForce(int slaveTag) const
{
    double f = 0.0;
    for (std::map<PairKey, double>::const_iterator it = theNtsForce.begin();
         it != theNtsForce.end(); ++it)
        if (it->first.s == slaveTag) f += it->second;
    return f;
}

LadrunoContactDomain::FrictionState &
LadrunoContactDomain::getOrCreateFrictionState(int contactTag, int slaveTag, int segIndex)
{
    PairKey k; k.c = contactTag; k.s = slaveTag; k.g = segIndex;
    // operator[] default-constructs a zeroed slot (engaged = false) when absent
    return theFrictionStates[k];
}

void
LadrunoContactDomain::frictionGCBegin(void)
{
    liveKeys.clear();
}

void
LadrunoContactDomain::frictionGCMark(int contactTag, int slaveTag, int segIndex)
{
    PairKey k; k.c = contactTag; k.s = slaveTag; k.g = segIndex;
    liveKeys.insert(k);
}

void
LadrunoContactDomain::frictionGCEnd(void)
{
    // erase every slot whose key the current handle() did not mark live
    for (std::map<PairKey, FrictionState>::iterator it = theFrictionStates.begin();
         it != theFrictionStates.end(); ) {
        if (liveKeys.find(it->first) == liveKeys.end())
            theFrictionStates.erase(it++);     // post-increment keeps the iterator valid
        else
            ++it;
    }
    liveKeys.clear();
}

int
LadrunoContactDomain::commit(bool)
{
    // promote the trial plastic slip and the engagement origin of every slot;
    // the counter lets a test check that the Domain::commit() hook fires
    numCommits++;
    for (std::map<PairKey, FrictionState>::iterator it = theFrictionStates.begin();
         it != theFrictionStates.end(); ++it) {
        FrictionState &st = it->second;
        for (int d = 0; d < 3; d++) {
            st.gpT[d] = st.gpTtrial[d];
            st.gT0committed[d] = st.gT0[d];
        }
        st.engagedCommitted = st.engaged;
    }

    // One Uzawa update of the mortar normal multipliers per commit (Simo and
    // Laursen 1992): lambda_I <- min(0, lambda_I + epsN gbar_I), gbar_I the
    // converged weighted gap summed over the facet pairs of the node. Across
    // load steps this augments with a stock Newton; a held-load loop
    // (ladrunoBeginAugment, analyze at a zero load increment, ladrunoEndAugment)
    // repeats it within a step until ladrunoMortarPenetration is small.
    for (std::map<NodeKey, MortarNormalState>::iterator it = theMortarNormalStates.begin();
         it != theMortarNormalStates.end(); ++it) {
        MortarNormalState &st = it->second;
        if (st.aGlobal <= 1e-300) continue;          // not integrated this step
        if (st.isTie) {
            // mesh tie: an equality, so the update has no clamp and drives the
            // weighted relative displacement to zero at a finite epsTie
            for (int d = 0; d < 3; d++)
                st.lambdaTie[d] += st.epsN * (st.rtGlobal[d] / st.aGlobal);
            continue;
        }
        double gbar = st.gtGlobal / st.aGlobal;
        st.lambdaN = std::min(0.0, st.lambdaN + st.epsN * gbar);
    }

    // mortar friction, per (slave node, facet pair): promote the trial slip,
    // one tangential Uzawa step (lambdaT <- the returned cone-capped traction,
    // which removes the elastic stick creep at a finite epsT), and the
    // engagement copies
    for (std::map<MortarPairKey, MortarFrictionState>::iterator it = theMortarFrictionStates.begin();
         it != theMortarFrictionStates.end(); ++it) {
        MortarFrictionState &st = it->second;
        for (int d = 0; d < 3; d++) {
            st.gpT[d] = st.gpTtrial[d];
            st.lambdaT[d] = st.lambdaTtrial[d];
            st.gT0committed[d] = st.gT0[d];
        }
        st.engagedCommitted = st.engaged;
    }
    return 0;
}

int
LadrunoContactDomain::revertToLastCommit(void)
{
    // drop the trial back to the committed state (a failed implicit step, or a
    // script that calls revertToLastCommit and retries)
    numReverts++;
    for (std::map<PairKey, FrictionState>::iterator it = theFrictionStates.begin();
         it != theFrictionStates.end(); ++it) {
        FrictionState &st = it->second;
        for (int d = 0; d < 3; d++) {
            st.gpTtrial[d] = st.gpT[d];
            st.gT0[d] = st.gT0committed[d];
        }
        st.engaged = st.engagedCommitted;
    }
    // mortar friction: the multipliers lambdaN and lambdaTie are changed only
    // by commit() and need nothing here
    for (std::map<MortarPairKey, MortarFrictionState>::iterator it = theMortarFrictionStates.begin();
         it != theMortarFrictionStates.end(); ++it) {
        MortarFrictionState &st = it->second;
        for (int d = 0; d < 3; d++) {
            st.gpTtrial[d] = st.gpT[d];
            st.lambdaTtrial[d] = st.lambdaT[d];
            st.gT0[d] = st.gT0committed[d];
        }
        st.engaged = st.engagedCommitted;
    }
    return 0;
}

int
LadrunoContactDomain::revertToStart(void)
{
    // Every slot is recreated zeroed at the next handle() or residual
    // evaluation, i.e. in the state of a fresh run.
    theFrictionStates.clear();
    theNtsForce.clear();
    theMortarNormalStates.clear();
    theMortarFacetContribs.clear();
    theMortarFrictionStates.clear();
    return 0;
}

int
LadrunoContactDomain::pruneRemovedNodes(Domain *dom)
{
    if (dom == 0)
        return 0;

    // 1. drop the nodes the analysis removed, keeping facets whole
    for (size_t i = 0; i < theSurfaces.size(); i++) {
        LadrunoContactSurface *s = theSurfaces[i];
        if (s == 0)
            continue;
        int before = s->getNodeTags().Size();
        int dropped = s->pruneMissingNodes(dom);
        if (dropped > 0) {
            opserr << "LadrunoContactDomain: contactSurface " << s->getTag()
                   << " lost " << dropped << " of " << before
                   << " node entries removed by the analysis; "
                   << s->getNodeTags().Size() << " remain\n";
        }
    }

    // 2. retire every interaction left without a usable surface
    int retiredNow = 0;
    for (size_t i = 0; i < theContacts.size(); i++) {
        Contact &c = theContacts[i];
        if (c.retired)
            continue;
        LadrunoContactSurface *ms = getSurface(c.masterSurfTag);
        LadrunoContactSurface *ss = getSurface(c.slaveSurfTag);
        if (ms == 0 || ss == 0 ||
            ms->getNodeTags().Size() == 0 || ss->getNodeTags().Size() == 0) {
            c.retired = true;
            retiredNow++;
            opserr << "LadrunoContactDomain: RETIRING contact " << c.tag
                   << " -- the analysis removed every node of its master or slave "
                      "surface, so the interaction no longer has a body to act on. "
                      "The run continues without it.\n";
        }
    }
    for (size_t i = 0; i < theMortarContacts.size(); i++) {
        MortarContact &m = theMortarContacts[i];
        if (m.retired)
            continue;
        LadrunoContactSurface *ms = getSurface(m.masterSurfTag);
        LadrunoContactSurface *ss = getSurface(m.slaveSurfTag);
        if (ms == 0 || ss == 0 ||
            ms->getNodeTags().Size() == 0 || ss->getNodeTags().Size() == 0) {
            m.retired = true;
            retiredNow++;
            opserr << "LadrunoContactDomain: RETIRING mortar contact " << m.tag
                   << " -- the analysis removed every node of its master or slave "
                      "surface. The run continues without it.\n";
        }
    }
    for (size_t i = 0; i < theRigidPlanes.size(); i++) {
        RigidPlane &r = theRigidPlanes[i];
        if (r.retired)
            continue;
        LadrunoContactSurface *ss = getSurface(r.slaveSurfTag);
        if (ss == 0 || ss->getNodeTags().Size() == 0) {
            r.retired = true;
            retiredNow++;
            opserr << "LadrunoContactDomain: RETIRING contactPlane " << r.tag
                   << " -- the analysis removed every node of its slave surface. "
                      "The run continues without it.\n";
        }
    }
    return retiredNow;
}
