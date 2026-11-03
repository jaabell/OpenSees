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
    return false;
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
    return 0;
}

int
LadrunoContactDomain::revertToStart(void)
{
    // Every slot is recreated zeroed at the next handle() or residual
    // evaluation, i.e. in the state of a fresh run.
    theFrictionStates.clear();
    theNtsForce.clear();
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
