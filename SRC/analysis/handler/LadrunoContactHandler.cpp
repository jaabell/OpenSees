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
#include <LadrunoMortarKernel.h>
#include <LadrunoEdgeKernel.h>
#include <LadrunoContact2DKernel.h>
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
#include <set>
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

// Pre-flight one contact surface: every node must exist and the surface must
// be 2D or 3D as a whole (getCrds()(2) is unchecked in release builds, so a 2D
// node in a 3D lane would read out of bounds). The dimension is returned in
// *dimOut (0 for an empty surface).
static bool
ladrunoSurfaceNodesOk(Domain *dom, int ctag, const LadrunoContactSurface *s, int *dimOut = 0)
{
    const ID &tags = s->getNodeTags();
    int dim = 0;
    for (int i = 0; i < tags.Size(); i++) {
        Node *n = dom->getNode(tags(i));
        if (n == 0) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": node " << tags(i) << " of contactSurface " << s->getTag()
                   << " does not exist\n";
            return false;
        }
        int nd = n->getCrds().Size();
        if (nd != 2 && nd != 3) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": node " << tags(i) << " of contactSurface " << s->getTag()
                   << " has " << nd << "D coordinates; contact needs 2D or 3D nodes\n";
            return false;
        }
        if (dim == 0) {
            dim = nd;
        } else if (nd != dim) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": contactSurface " << s->getTag() << " mixes 2D and 3D nodes (node "
                   << tags(i) << ")\n";
            return false;
        }
    }
    if (dimOut != 0) *dimOut = dim;
    return true;
}

// Both surfaces of one contact must have the same dimension; returns it in
// *dimOut (2 or 3).
static bool
ladrunoContactPairDim(Domain *dom, int ctag, const LadrunoContactSurface *ms,
                      const LadrunoContactSurface *ss, int *dimOut)
{
    int mdim = 0, sdim = 0;
    if (!ladrunoSurfaceNodesOk(dom, ctag, ms, &mdim) ||
        !ladrunoSurfaceNodesOk(dom, ctag, ss, &sdim))
        return false;
    if (mdim != 0 && sdim != 0 && mdim != sdim) {
        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
               << ": the master surface " << ms->getTag() << " is " << mdim
               << "D and the slave surface " << ss->getTag() << " is " << sdim
               << "D; both surfaces of a contact must have one dimension\n";
        return false;
    }
    *dimOut = (mdim == 2 || sdim == 2) ? 2 : 3;
    return true;
}

// Minimum initial length of the 2-node segments of a 2D surface (the
// relative gauge of LadrunoContact2DKernel); <= 0 with a message when a
// segment has zero length or the surface has none.
static double
ladrunoMinSegLength2D(Domain *dom, int ctag, const LadrunoContactSurface *s)
{
    const ID &tags = s->getNodeTags();
    int nSeg = tags.Size() / 2;
    double Lref = 0.0;
    for (int seg = 0; seg < nSeg; seg++) {
        const Vector &A = dom->getNode(tags(2 * seg))->getCrds();
        const Vector &B = dom->getNode(tags(2 * seg + 1))->getCrds();
        double dx = B(0) - A(0), dy = B(1) - A(1);
        double L = std::sqrt(dx * dx + dy * dy);
        if (!(L > 0.0)) {
            opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
                   << ": segment " << seg << " of 2D contactSurface " << s->getTag()
                   << " has zero length\n";
            return -1.0;
        }
        if (Lref <= 0.0 || L < Lref) Lref = L;
    }
    if (nSeg <= 0) {
        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag
               << ": 2D contactSurface " << s->getTag() << " has no segment\n";
        return -1.0;
    }
    return Lref;
}

// Automatic penalty of a 2D master segment: the 2D form of
// ladrunoResolveAutoKn, kn = 0.10 mean over the two nodes of n^T K_node n
// from the initial stiffness of the first 2-DOF-per-node element that holds
// both nodes, n the reference normal sigma perp(t)/|t|.
static double
ladrunoResolveAutoKn2D(Domain *theDomain, Node **segNodes, double sigma)
{
    const double F_SI = 0.10;
    if (segNodes[0] == 0 || segNodes[1] == 0) return -1.0;
    const Vector &XA = segNodes[0]->getCrds();
    const Vector &XB = segNodes[1]->getCrds();
    const double t[2] = { XB(0) - XA(0), XB(1) - XA(1) };
    const double L = std::sqrt(t[0]*t[0] + t[1]*t[1]);
    if (!(L > 0.0)) return -1.0;
    double p[2];
    LadrunoContact2DKernel::perp2(t, p);
    const double n[2] = { sigma * p[0] / L, sigma * p[1] / L };
    int segTags[2] = { segNodes[0]->getTag(), segNodes[1]->getTag() };
    ElementIter &theEle = theDomain->getElements();
    Element *e;
    while ((e = theEle()) != 0) {
        const ID &en = e->getExternalNodes();
        int nn = en.Size();
        if (nn <= 0 || e->getNumDOF() != 2 * nn) continue;   // not a 2-DOF/node solid
        int loc[2]; bool allIn = true;
        for (int k = 0; k < 2; k++) {
            loc[k] = en.getLocation(segTags[k]);
            if (loc[k] < 0) { allIn = false; break; }
        }
        if (!allIn) continue;
        const Matrix &K = e->getInitialStiff();
        if (K.noRows() != 2 * nn || K.noCols() != 2 * nn) return -1.0;
        double acc = 0.0;
        for (int k = 0; k < 2; k++) {
            int c = 2 * loc[k];
            double Kn0 = K(c + 0, c + 0) * n[0] + K(c + 0, c + 1) * n[1];
            double Kn1 = K(c + 1, c + 0) * n[0] + K(c + 1, c + 1) * n[1];
            acc += n[0] * Kn0 + n[1] * Kn1;
        }
        double kn = F_SI * acc / 2.0;
        return (kn > 0.0) ? kn : -1.0;
    }
    return -1.0;
}

// Orientation of a 2D interface, one sign for the whole master surface: the
// direction toward the slave side is -outward when given, else the slave
// centroid minus the master centroid (each node counted once), and every
// master segment votes sign(perp(t).refDir). The interface is refused when
// the reference is degenerate (a flush interface with coincident centroids,
// or a direction tangent to every segment) or when the segments disagree (a
// strongly curved or inconsistently wound master), instead of guessing a
// per-segment sign.
static bool
ladruno2DOrientationVote(Domain *theDomain, int ctag, bool hasOutward, const double outward[3],
                         const ID &mTags, int nSeg, const ID &sTags, double Lref,
                         double &sigmaOut)
{
    double refDir[2] = { 0.0, 0.0 };
    if (hasOutward) {
        refDir[0] = outward[0]; refDir[1] = outward[1];
    } else {
        double mc2[2] = { 0.0, 0.0 }; int mcn = 0;
        std::set<int> mseen;
        for (int i = 0; i < mTags.Size(); i++) {
            if (!mseen.insert(mTags(i)).second) continue;
            const Vector &X = theDomain->getNode(mTags(i))->getCrds();
            mc2[0] += X(0); mc2[1] += X(1); mcn++;
        }
        double sc2[2] = { 0.0, 0.0 }; int scn = 0;
        std::set<int> sseen;
        for (int i = 0; i < sTags.Size(); i++) {
            if (!sseen.insert(sTags(i)).second) continue;
            const Vector &X = theDomain->getNode(sTags(i))->getCrds();
            sc2[0] += X(0); sc2[1] += X(1); scn++;
        }
        if (mcn > 0 && scn > 0) {
            refDir[0] = sc2[0] / scn - mc2[0] / mcn;
            refDir[1] = sc2[1] / scn - mc2[1] / mcn;
        }
    }
    // a centroid offset below 1e-6 Lref is a flush interface
    bool voteMagOk = true;
    if (!hasOutward) {
        const double magFloor = 1.0e-6 * Lref;
        voteMagOk = (refDir[0] * refDir[0] + refDir[1] * refDir[1] > magFloor * magFloor);
    }
    int votePlus = 0, voteMinus = 0;
    for (int seg = 0; seg < nSeg; seg++) {
        const Vector &A = theDomain->getNode(mTags(2 * seg))->getCrds();
        const Vector &B = theDomain->getNode(mTags(2 * seg + 1))->getCrds();
        const double t[2] = { B(0) - A(0), B(1) - A(1) };
        double sg = LadrunoContact2DKernel::sigmaFromRef2D(t, refDir);
        if (sg > 0.0) votePlus++;
        else if (sg < 0.0) voteMinus++;
    }
    if ((votePlus == 0 && voteMinus == 0) || !voteMagOk) {
        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag << ": the "
               << (hasOutward ? "given -outward direction is zero or tangent to the master "
                                "segments"
                              : "slave and master centroids nearly coincide (a flush "
                                "interface)")
               << ", so the side of the 2D interface cannot be decided; give -outward ox oy"
                  " (or -outward winding for node-to-segment contact)\n";
        return false;
    }
    if (votePlus > 0 && voteMinus > 0) {
        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ctag << ": the "
                  "segments of the 2D master surface disagree on the side of the slave ("
               << votePlus << " against " << voteMinus << "; a strongly curved, closed or "
                  "inconsistently wound master); split the contact or give -outward winding\n";
        return false;
    }
    sigmaOut = (votePlus > 0) ? 1.0 : -1.0;
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

// Newell unit normal of a planar tri/quad facet (zero for a degenerate facet)
static void
ladrunoFacetNormalNewell(int npts, const double X[4][3], double n[3])
{
    n[0] = n[1] = n[2] = 0.0;
    for (int i = 0; i < npts; i++) {
        const double *a = X[i]; const double *b = X[(i + 1) % npts];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    double nn = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (nn > 1e-300) { n[0] /= nn; n[1] /= nn; n[2] /= nn; }
}

// --- routing of edge-to-edge contact (-edgeedge), at the reference
//     configuration of the facets ---

// Does x project inside the convex facet X (unit normal nf) with |gap| <= band?
static bool
ladrunoPointInFacet(const double x[3], int npts, const double X[4][3],
                    const double nf[3], double band, double &gap)
{
    double c[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < npts; i++) for (int d = 0; d < 3; d++) c[d] += X[i][d] / npts;
    gap = 0.0;
    for (int d = 0; d < 3; d++) gap += (x[d] - c[d]) * nf[d];
    double proj[3];
    for (int d = 0; d < 3; d++) proj[d] = x[d] - gap * nf[d];
    int sgn = 0;            // same side of every edge (convex facet)
    for (int i = 0; i < npts; i++) {
        const double *a = X[i]; const double *b = X[(i + 1) % npts];
        double e[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
        double w[3] = {proj[0]-a[0], proj[1]-a[1], proj[2]-a[2]};
        double cr[3] = {e[1]*w[2]-e[2]*w[1], e[2]*w[0]-e[0]*w[2], e[0]*w[1]-e[1]*w[0]};
        double dp = cr[0]*nf[0] + cr[1]*nf[1] + cr[2]*nf[2];
        int s = (dp > 1e-300) ? 1 : ((dp < -1e-300) ? -1 : 0);
        if (s != 0) { if (sgn == 0) sgn = s; else if (s != sgn) return false; }
    }
    return std::fabs(gap) <= band;
}

// Distance of two facets: the vertices of each that project inside the other,
// and the closest points of every pair of boundary edges.
static double
ladrunoMinFacetDist(int npS, const double S[4][3], const double nS[3],
                    int npM, const double M[4][3], const double nM[3])
{
    double dmin = HUGE_VAL, gap;
    for (int i = 0; i < npS; i++)
        if (ladrunoPointInFacet(S[i], npM, M, nM, HUGE_VAL, gap)) {
            double g = std::fabs(gap); if (g < dmin) dmin = g;
        }
    for (int i = 0; i < npM; i++)
        if (ladrunoPointInFacet(M[i], npS, S, nS, HUGE_VAL, gap)) {
            double g = std::fabs(gap); if (g < dmin) dmin = g;
        }
    for (int ie = 0; ie < npS; ie++)
        for (int je = 0; je < npM; je++) {
            LadrunoEdgeKernel::ClosestResult cr = LadrunoEdgeKernel::closestPtSegSeg(
                S[ie], S[(ie + 1) % npS], M[je], M[(je + 1) % npM]);
            double w[3] = {cr.c1[0]-cr.c2[0], cr.c1[1]-cr.c2[1], cr.c1[2]-cr.c2[2]};
            double g = std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
            if (g < dmin) dmin = g;
        }
    return dmin;
}

// Routing band of a master facet: -edgeBand when given, else a quarter of the
// mean edge length of the facet.
static double
ladrunoEdgeBand(int npM, const double M[4][3], double given)
{
    if (given > 0.0) return given;
    double el = 0.0;
    for (int k = 0; k < npM; k++) {
        const double *a = M[k]; const double *b = M[(k + 1) % npM];
        double dx = b[0]-a[0], dy = b[1]-a[1], dz = b[2]-a[2];
        el += std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    return 0.25 * (el / npM);
}

// Does the edge-to-edge lane own this (slave facet, master facet) pair? It
// does when the facets are within the band, do not face each other (normals
// more than 50 degrees apart, |nS.nM| < cos 50), no slave vertex lies on the
// master face (that contact is a node-to-face one), and at least one pair of
// boundary edges crosses strictly inside both edges within the band. The
// mortar loop skips a pair that this returns true for and the edge loop
// builds adapters only for such pairs, so the mortar clip and the edge force
// never act on the same pair (in the oblique range the clip is not empty).
static bool
ladrunoEdgeEdgeOwns(int npS, const double S[4][3], int npM, const double M[4][3],
                    double dBand)
{
    const double TAU_FACE = std::cos(50.0 * 3.14159265358979323846 / 180.0);
    double nS[3], nM[3];
    ladrunoFacetNormalNewell(npS, S, nS);
    ladrunoFacetNormalNewell(npM, M, nM);
    if (ladrunoMinFacetDist(npS, S, nS, npM, M, nM) > dBand) return false;   // apart
    double alignFn = std::fabs(nS[0]*nM[0] + nS[1]*nM[1] + nS[2]*nM[2]);
    if (alignFn >= TAU_FACE) return false;                                   // facing: mortar
    double gap;
    for (int k = 0; k < npS; k++)                                           // vertex on face
        if (ladrunoPointInFacet(S[k], npM, M, nM, dBand, gap)) return false;
    for (int ie = 0; ie < npS; ie++)
        for (int je = 0; je < npM; je++) {
            LadrunoEdgeKernel::ClosestResult cr = LadrunoEdgeKernel::closestPtSegSeg(
                S[ie], S[(ie + 1) % npS], M[je], M[(je + 1) % npM]);
            if (cr.status != LadrunoEdgeKernel::EE_OK || !cr.interior) continue;
            double w[3] = {cr.c1[0]-cr.c2[0], cr.c1[1]-cr.c2[1], cr.c1[2]-cr.c2[2]};
            if (std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]) <= dBand) return true;
        }
    return false;
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
        cd->clearVtx2DSide();    // 2D vertex side signs are captured again in each pairing

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
            int pairDim = 3;
            if (!ladrunoContactPairDim(theDomain, ct.tag, ms, ss, &pairDim))
                return ladrunoContactFatal();
            int nps = ms->getNodesPerSeg();
            if (pairDim == 2) {
                // ---- 2D node-to-segment contact (LadrunoContact2DKernel) ----
                if (nps != 2) {
                    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": a 2D master surface needs nodesPerSeg 2 (got " << nps << ")\n";
                    return ladrunoContactFatal();
                }
                if (ct.consistentNormal &&
                    cd->warnOnce(ct.tag, LadrunoContactDomain::WARN_2D_GEOMTAN_NOOP))
                    opserr << "WARNING LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": -geomtan has no effect on a 2D contact (the curvature term "
                              "of the 2D segment tangent is not implemented)\n";
                if (!ct.knAuto && ct.kn <= 0.0) {
                    opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                           << ": segment contact needs kn > 0 or kn auto (got " << ct.kn << ")\n";
                    return ladrunoContactFatal();
                }
                const ID &mTags = ms->getNodeTags();
                const ID &sTags = ss->getNodeTags();
                int nSeg = mTags.Size() / 2;
                double Lref = ladrunoMinSegLength2D(theDomain, ct.tag, ms);
                if (!(Lref > 0.0))
                    return ladrunoContactFatal();
                // Chain integrity: the surface order is the adjacency of the
                // ownership rule and of the vertex pairs, so a node shared by two
                // segments must be the end of segment k and the start of segment
                // k + 1 (a closed loop wraps from the last to the first); a
                // permuted, reversed or forked listing is refused.
                {
                    std::map<int, std::vector<std::pair<int, int> > > occ;
                    for (int seg = 0; seg < nSeg; seg++) {
                        occ[mTags(2 * seg)].push_back(std::make_pair(seg, 0));
                        occ[mTags(2 * seg + 1)].push_back(std::make_pair(seg, 1));
                    }
                    for (std::map<int, std::vector<std::pair<int, int> > >::const_iterator
                             it = occ.begin(); it != occ.end(); ++it) {
                        const std::vector<std::pair<int, int> > &v = it->second;
                        if (v.size() <= 1) continue;
                        bool okLink = false;
                        if (v.size() == 2 && v[0].second != v[1].second) {
                            int tailSeg = (v[0].second == 1) ? v[0].first : v[1].first;
                            int headSeg = (v[0].second == 1) ? v[1].first : v[0].first;
                            if (nSeg > 1 && headSeg == (tailSeg + 1) % nSeg)
                                okLink = true;
                        }
                        if (!okLink) {
                            opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                   << ct.tag << ": node " << it->first << " of 2D master "
                                      "surface " << ct.masterSurfTag << " is shared by "
                                      "segments " << v[0].first << " and " << v[1].first
                                   << " that are not consecutive head to tail; list the "
                                      "segments as n0 n1 n1 n2 ... in chain order\n";
                            return ladrunoContactFatal();
                        }
                    }
                }
                // -outward winding replaces the vote, so a second run wound the
                // other way would not be caught: the master must be one chain
                if (ct.outwardWinding) {
                    for (int seg = 0; seg + 1 < nSeg; seg++) {
                        if (mTags(2 * seg + 1) == mTags(2 * seg + 2)) continue;
                        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                               << ": -outward winding needs 2D master surface "
                               << ct.masterSurfTag << " to be one connected chain, but "
                                  "segment " << seg << " ends at node " << mTags(2 * seg + 1)
                               << " and segment " << seg + 1 << " starts at node "
                               << mTags(2 * seg + 2) << "; declare one contact per chain\n";
                        return ladrunoContactFatal();
                    }
                }
                double sigma = 0.0;
                if (ct.outwardWinding)
                    sigma = 1.0;     // the slave is to the left of the traversal
                else if (!ladruno2DOrientationVote(theDomain, ct.tag, ct.hasOutward, ct.outward,
                                                   mTags, nSeg, sTags, Lref, sigma))
                    return ladrunoContactFatal();
                // concave corners (turn > 0 seen from the contact side), the only
                // vertices whose pair can carry force
                struct LadrunoVtx2D { int segPrev, segNext; Node *prevFar, *vtx, *nextFar; };
                std::vector<LadrunoVtx2D> vtx2d;
                for (int seg = 0; nSeg >= 2 && seg < nSeg; seg++) {
                    int nxt = (seg + 1) % nSeg;
                    if (mTags(2 * seg + 1) != mTags(2 * nxt)) continue;
                    Node *pF = theDomain->getNode(mTags(2 * seg));
                    Node *vN = theDomain->getNode(mTags(2 * seg + 1));
                    Node *nF = theDomain->getNode(mTags(2 * nxt + 1));
                    const Vector &XP = pF->getCrds();
                    const Vector &XV = vN->getCrds();
                    const Vector &XN = nF->getCrds();
                    const double tPrev[2] = { XV(0) - XP(0), XV(1) - XP(1) };
                    const double tNext[2] = { XN(0) - XV(0), XN(1) - XV(1) };
                    double turn = sigma * (tPrev[0] * tNext[1] - tPrev[1] * tNext[0]);
                    if (turn > 0.0) {
                        LadrunoVtx2D v;
                        v.segPrev = seg; v.segNext = nxt;
                        v.prevFar = pF; v.vtx = vN; v.nextFar = nF;
                        vtx2d.push_back(v);
                    }
                }
                // broad phase in the plane z = 0
                std::vector<double> segCoords((size_t)nSeg * 2 * 3, 0.0);
                for (int seg = 0; seg < nSeg; seg++) {
                    double *S = &segCoords[(size_t)seg * 2 * 3];
                    for (int k = 0; k < 2; k++) {
                        const Vector &X = theDomain->getNode(mTags(seg * 2 + k))->getCrds();
                        S[k * 3 + 0] = X(0); S[k * 3 + 1] = X(1); S[k * 3 + 2] = 0.0;
                    }
                }
                LadrunoContactBucketSort::Grid grid(nSeg, 2, segCoords.data(), ct.cellFrac, 1.0);
                std::vector<int> cand(nSeg);
                for (int si = 0; si < sTags.Size(); si++) {
                    Node *sn = theDomain->getNode(sTags(si));
                    // the first two DOFs are the translations; more DOFs (a
                    // rotation, a pressure) are carried along and not coupled
                    if (sn == 0 || sn->getNumberDOF() < 2) {
                        opserr << "FATAL LadrunoContactHandler::handle() - contact " << ct.tag
                               << " slave node " << sTags(si) << " ndf="
                               << (sn != 0 ? sn->getNumberDOF() : 0) << " < 2\n";
                        return ladrunoContactFatal();
                    }
                    const Vector &Xs0 = sn->getCrds();
                    double slavePt[3] = { Xs0(0), Xs0(1), 0.0 };
                    int nCand = grid.candidates(slavePt, cand.data());
                    for (int ci = 0; ci < nCand; ci++) {
                        int seg = cand[ci];
                        Node *segNodes[2];
                        for (int k = 0; k < 2; k++) {
                            Node *mn = theDomain->getNode(mTags(seg * 2 + k));
                            if (mn == 0 || mn->getNumberDOF() < 2) {
                                opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                       << ct.tag << " master segment " << seg << " node "
                                       << mTags(seg * 2 + k) << ": ndf < 2 (or node absent)\n";
                                return ladrunoContactFatal();
                            }
                            segNodes[k] = mn;
                        }
                        double knUse = ct.kn;
                        if (ct.knAuto) {
                            knUse = ladrunoResolveAutoKn2D(theDomain, segNodes, sigma);
                            if (knUse <= 0.0) {
                                opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                       << ct.tag << " slave node " << sTags(si) << " segment "
                                       << seg << ": kn auto could not size a penalty (no "
                                          "owning element with 2 DOFs per node). Give an "
                                          "explicit kn.\n";
                                return ladrunoContactFatal();
                            }
                        }
                        // neighbours in surface order (null at an open end)
                        Node *prevFar = 0, *nextFar = 0;
                        if (nSeg > 1) {
                            int prv = (seg + nSeg - 1) % nSeg;
                            if (prv != seg && mTags(2 * prv + 1) == mTags(2 * seg))
                                prevFar = theDomain->getNode(mTags(2 * prv));
                            int nx2 = (seg + 1) % nSeg;
                            if (nx2 != seg && mTags(2 * seg + 1) == mTags(2 * nx2))
                                nextFar = theDomain->getNode(mTags(2 * nx2 + 1));
                        }
                        LadrunoContactFE *fe = new LadrunoContactFE(
                            numFe++, sn, segNodes, /*nps=*/2, /*ndm2=*/2, knUse, sigma, Lref,
                            prevFar, nextFar, ct.kt, ct.mu, ct.consistentTan, ct.muc,
                            theDomain, ct.tag, seg, ct.cohesion, ct.tauMax);
                        if (fe == 0) return -5;
                        theModel->addFE_Element(fe);
                        if (ct.hasFriction())
                            cd->frictionGCMark(ct.tag, sTags(si), seg);
                        // an end cap at each open end of the chain (pair keys
                        // 2 nSeg + seg and 3 nSeg + seg, apart from the segments
                        // and the concave vertices)
                        for (int end = 0; end < 2; end++) {
                            if ((end == 0 && prevFar != 0) || (end == 1 && nextFar != 0))
                                continue;
                            Node *capNodes[1] = { segNodes[end] };
                            int capIndex = (end == 0 ? 2 : 3) * nSeg + seg;
                            LadrunoContactFE *feCap = new LadrunoContactFE(
                                numFe++, sn, capNodes, /*nps=*/1, /*ndm2=*/2, knUse, sigma,
                                Lref, end == 0 ? 0 : segNodes[0], end == 0 ? segNodes[1] : 0,
                                ct.kt, ct.mu, ct.consistentTan, ct.muc, theDomain, ct.tag,
                                capIndex, ct.cohesion, ct.tauMax);
                            if (feCap == 0) return -5;
                            theModel->addFE_Element(feCap);
                            if (ct.hasFriction())
                                cd->frictionGCMark(ct.tag, sTags(si), capIndex);
                        }
                    }
                    // concave vertex pairs next to the candidate segments (pair
                    // key nSeg + segPrev)
                    for (size_t vi = 0; vi < vtx2d.size(); vi++) {
                        const LadrunoVtx2D &v = vtx2d[vi];
                        bool nearV = false;
                        for (int ci = 0; ci < nCand && !nearV; ci++)
                            if (cand[ci] == v.segPrev || cand[ci] == v.segNext)
                                nearV = true;
                        if (!nearV) continue;
                        double knV = ct.kn;
                        if (ct.knAuto) {
                            Node *prevSeg[2] = { v.prevFar, v.vtx };
                            knV = ladrunoResolveAutoKn2D(theDomain, prevSeg, sigma);
                            if (knV <= 0.0) {
                                opserr << "FATAL LadrunoContactHandler::handle() - contact "
                                       << ct.tag << " vertex " << v.vtx->getTag()
                                       << ": kn auto could not size a penalty. Give an "
                                          "explicit kn.\n";
                                return ladrunoContactFatal();
                            }
                        }
                        Node *vNodes[1] = { v.vtx };
                        LadrunoContactFE *fe = new LadrunoContactFE(
                            numFe++, sn, vNodes, /*nps=*/1, /*ndm2=*/2, knV, sigma, Lref,
                            v.prevFar, v.nextFar, ct.kt, ct.mu, ct.consistentTan, ct.muc,
                            theDomain, ct.tag, nSeg + v.segPrev, ct.cohesion, ct.tauMax);
                        if (fe == 0) return -5;
                        theModel->addFE_Element(fe);
                        if (ct.hasFriction())
                            cd->frictionGCMark(ct.tag, sTags(si), nSeg + v.segPrev);
                    }
                }
                continue;
            }
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
                if ((nd != 2 && nd != 3) || sn->getNumberDOF() < nd) {
                    // the adapter couples the first nd translational DOFs; fewer
                    // DOFs than coordinates would leave myID slots unset
                    opserr << "FATAL LadrunoContactHandler::handle() - rigid-plane slave node "
                           << nodeTags(k) << " has ndm=" << nd << ", ndf=" << sn->getNumberDOF()
                           << " (needs 2D or 3D coordinates and ndf >= ndm)\n";
                    return ladrunoContactFatal();
                }
                LadrunoContactFE *fe = new LadrunoContactFE(numFe++, sn, nd, rp.p0, rp.n, rp.kn,
                                                            theDomain, rp.tag, rp.muc);
                if (fe == 0) return -5;
                theModel->addFE_Element(fe);
            }
        }
    }

    // Mortar segment-to-segment contact: one adapter per (slave facet, master
    // facet) pair. The pairing is exhaustive: the bucket sort of the
    // node-to-segment lane is a superset for point slaves only, and a query at
    // the centroid of an extended slave facet would miss overlaps when the slave
    // facet is coarser than the master buckets. The clip of the adapter rejects
    // a non-overlapping pair at each evaluation (zero force).
    //
    // A mesh tie (-tie) binds the pairs that overlap in the reference
    // configuration, and only those whose mean |gap| over the overlap is at most
    // 0.5 sqrt(A_s), A_s the area of the slave facet. The clip accepts an
    // anti-parallel facet across a closed or folded surface (the far side of a
    // tube projects onto the auxiliary plane of a near facet), and without the
    // proximity condition every slave facet would also be tied to the antipodal
    // master facets; their separation is of the order of the surface diameter,
    // that of the facing pairs of the order of the facet sagitta. This is the
    // exclusion of the mortar generator of LadrunoTie with a tolerance factor 1.
    // An explicit -maxGap replaces it.
    if (cd != 0) {
        cd->mortarNormalGCBegin();   // multiplier slots of nodes not rebuilt below are collected
        for (int c = 0; c < cd->getNumMortarContacts(); c++) {
            const LadrunoContactDomain::MortarContact &mc = cd->getMortarContact(c);
            if (mc.retired) continue;
            LadrunoContactSurface *ms = cd->getSurface(mc.masterSurfTag);
            LadrunoContactSurface *ss = cd->getSurface(mc.slaveSurfTag);
            if (ms == 0 || ss == 0) {
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": master/slave surface not defined\n";
                return ladrunoContactFatal();
            }
            if (ms->getKind() != LadrunoContactSurface::MASTER_SEGMENTS ||
                ss->getKind() != LadrunoContactSurface::SLAVE_SEGMENTS) {
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": needs a -master facet surface and a -slave-segments facet surface\n";
                return ladrunoContactFatal();
            }
            int pairDim = 3;
            if (!ladrunoContactPairDim(theDomain, mc.tag, ms, ss, &pairDim))
                return ladrunoContactFatal();
            int npsM = ms->getNodesPerSeg(), npsS = ss->getNodesPerSeg();
            if (pairDim == 3 && mc.hThickness != 1.0) {
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": -thickness is a 2D option\n";
                return ladrunoContactFatal();
            }
            if (pairDim == 2) {
                // ---- 2D mortar contact (LadrunoContact2DKernel) ----
                if (npsM != 2 || npsS != 2) {
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": 2D surfaces need nodesPerSeg 2\n";
                    return ladrunoContactFatal();
                }
                if (mc.smoothN > 0.0 || mc.smoothT > 0.0 || mc.maxGap > 0.0 ||
                    mc.gapOffset != 0.0 || mc.adjust || mc.edgeEdge) {
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": -smoothN, -smoothT, -maxGap, -gapOffset, -adjust and "
                              "-edgeedge are 3D options (a 2D interface has no edge-to-edge "
                              "primitive)\n";
                    return ladrunoContactFatal();
                }
                if (mc.ngp != 2 && cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_2D_NGP_NOOP))
                    opserr << "WARNING LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": -ngp has no effect on a 2D mortar contact (2-point Gauss on "
                              "the clipped interval)\n";
                const ID &mTags2 = ms->getNodeTags();
                const ID &sTags2 = ss->getNodeTags();
                int nSegM = mTags2.Size() / 2;
                int nSegS = sTags2.Size() / 2;
                double LrefM = ladrunoMinSegLength2D(theDomain, mc.tag, ms);
                double LrefS = ladrunoMinSegLength2D(theDomain, mc.tag, ss);
                if (!(LrefM > 0.0) || !(LrefS > 0.0))
                    return ladrunoContactFatal();
                double Lref = (LrefM < LrefS) ? LrefM : LrefS;   // both surfaces are integrated
                double sigma = 0.0;
                if (!ladruno2DOrientationVote(theDomain, mc.tag, mc.hasOutward, mc.outward,
                                              mTags2, nSegM, sTags2, Lref, sigma))
                    return ladrunoContactFatal();
                // slave nodes covered by a node-to-segment contact on the same
                // master (by node set): a nearly perpendicular pair is then
                // left to it silently
                std::set<int> mcMasterNodes;
                for (int ti = 0; ti < mTags2.Size(); ti++) mcMasterNodes.insert(mTags2(ti));
                std::set<int> ntsSlaveUnion;
                for (int ci = 0; ci < cd->getNumContacts(); ci++) {
                    const LadrunoContactDomain::Contact &ntsCt = cd->getContact(ci);
                    if (ntsCt.retired) continue;
                    LadrunoContactSurface *ntsMaster = cd->getSurface(ntsCt.masterSurfTag);
                    LadrunoContactSurface *ntsSlave = cd->getSurface(ntsCt.slaveSurfTag);
                    if (ntsMaster == 0 || ntsSlave == 0) continue;
                    std::set<int> mset;
                    const ID &nm = ntsMaster->getNodeTags();
                    for (int ti = 0; ti < nm.Size(); ti++) mset.insert(nm(ti));
                    if (mset != mcMasterNodes) continue;
                    const ID &ntsSlaveTags = ntsSlave->getNodeTags();
                    for (int ti = 0; ti < ntsSlaveTags.Size(); ti++)
                        ntsSlaveUnion.insert(ntsSlaveTags(ti));
                }
                bool epsAuto = mc.epsNAuto || mc.knAuto;
                double epsFixed = (mc.epsN > 0.0) ? mc.epsN : mc.kn;
                if (!epsAuto && epsFixed <= 0.0) {
                    opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                           << ": needs a penalty (-epsN value, -epsN auto, or kn > 0)\n";
                    return ladrunoContactFatal();
                }
                std::vector<double> epsAutoPerSeg;
                if (epsAuto) epsAutoPerSeg.assign(nSegM, -1.0);
                double h = mc.hThickness;
                for (int sf = 0; sf < nSegS; sf++) {
                    Node *sNodes[2];
                    for (int k = 0; k < 2; k++) {
                        Node *sn = theDomain->getNode(sTags2(sf * 2 + k));
                        if (sn == 0 || sn->getNumberDOF() < 2) {
                            opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << " slave node " << sTags2(sf * 2 + k)
                                   << ": ndf < 2 (or node absent)\n";
                            return ladrunoContactFatal();
                        }
                        sNodes[k] = sn;
                    }
                    bool ntsCoDeclared = (ntsSlaveUnion.count(sNodes[0]->getTag()) > 0 &&
                                          ntsSlaveUnion.count(sNodes[1]->getTag()) > 0);
                    for (int seg = 0; seg < nSegM; seg++) {
                        Node *mNodes[2];
                        for (int k = 0; k < 2; k++) {
                            Node *mn = theDomain->getNode(mTags2(seg * 2 + k));
                            if (mn == 0 || mn->getNumberDOF() < 2) {
                                opserr << "FATAL LadrunoContactHandler::handle() - mortar "
                                          "contact " << mc.tag << " master segment " << seg
                                       << " node " << mTags2(seg * 2 + k)
                                       << ": ndf < 2 (or node absent)\n";
                                return ladrunoContactFatal();
                            }
                            mNodes[k] = mn;
                        }
                        double epsUse = epsFixed;
                        if (epsAuto) {
                            if (epsAutoPerSeg[seg] < 0.0)
                                epsAutoPerSeg[seg] = ladrunoResolveAutoKn2D(theDomain, mNodes, sigma);
                            epsUse = epsAutoPerSeg[seg];
                            if (epsUse <= 0.0) {
                                opserr << "FATAL LadrunoContactHandler::handle() - mortar "
                                          "contact " << mc.tag << " master segment " << seg
                                       << ": the automatic penalty could not be sized (no "
                                          "owning element with 2 DOFs per node). Give an "
                                          "explicit -epsN.\n";
                                return ladrunoContactFatal();
                            }
                        }
                        bool wantFric = mc.hasFriction();
                        bool epsTFromEpsN = mc.epsTAuto;
                        double epsTuse = mc.epsTAuto ? epsUse : mc.epsT;
                        if (wantFric && epsTuse <= 0.0) {
                            epsTuse = epsUse;
                            epsTFromEpsN = true;
                            if (cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_EPST))
                                opserr << "WARNING LadrunoContactHandler::handle() - mortar "
                                          "contact " << mc.tag << ": friction without -epsT; "
                                          "the tangential penalty is the normal penalty\n";
                        }
                        if (wantFric && mc.consistentTan &&
                            cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_CSL))
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << ": -consistanttan (non-symmetric friction "
                                      "tangent) needs a non-symmetric solver\n";
                        // the one place the thickness enters: penalties and the
                        // cohesion and shear cap are per unit area, the 2D
                        // integrals per unit thickness; an automatic penalty
                        // already carries the element thickness
                        double epsUseH  = epsAuto ? epsUse : epsUse * h;
                        double epsTuseH = epsTFromEpsN ? epsUseH : epsTuse * h;
                        LadrunoContactFE *fe = new LadrunoContactFE(
                            numFe++, sNodes, /*npsS=*/2, mNodes, /*npsM=*/2, /*ndm2=*/2,
                            epsUseH, sigma, Lref, mc.tag, sf, theDomain, mc.isTie, mc.mu,
                            epsTuseH, mc.cohesion * h, mc.tauMax * h, mc.consistentTan, seg,
                            ntsCoDeclared);
                        if (fe == 0) return -5;
                        theModel->addFE_Element(fe);
                        for (int k = 0; k < 2; k++) {
                            cd->mortarNormalGCMark(mc.tag, sNodes[k]->getTag());
                            if (wantFric)
                                cd->mortarFrictionGCMark(mc.tag, sNodes[k]->getTag(), sf, seg);
                        }
                    }
                }
                continue;
            }
            if (npsM < 3 || npsM > 4 || npsS < 3 || npsS > 4) {
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": nodesPerSeg " << npsM << "/" << npsS << " unsupported (need 3 or 4)\n";
                return ladrunoContactFatal();
            }
            // penalty: -epsN if given, else the positional kn; sized per pair when automatic
            bool epsAuto = mc.epsNAuto || mc.knAuto;
            double epsFixed = (mc.epsN > 0.0) ? mc.epsN : mc.kn;
            if (!epsAuto && epsFixed <= 0.0) {
                opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                       << ": needs a penalty (-epsN value, -epsN auto, or kn > 0)\n";
                return ladrunoContactFatal();
            }
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSegM = mTags.Size() / npsM;
            int nSegS = sTags.Size() / npsS;
            int nTieExcluded = 0;

            for (int sf = 0; sf < nSegS; sf++) {
                Node *sNodes[4]; double scen[3] = {0.0, 0.0, 0.0};
                double Sx[4][3];                              // reference coordinates
                for (int k = 0; k < npsS; k++) {
                    Node *sn = theDomain->getNode(sTags(sf * npsS + k));
                    if (sn == 0 || sn->getNumberDOF() < 3) {
                        opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                               << mc.tag << " slave node " << sTags(sf * npsS + k) << " ndf="
                               << (sn != 0 ? sn->getNumberDOF() : -1)
                               << " < 3 (mortar contact needs the first three DOFs to be the "
                                  "translations)\n";
                        return ladrunoContactFatal();
                    }
                    sNodes[k] = sn;
                    const Vector &Xk = sn->getCrds();
                    for (int d = 0; d < 3; d++) { Sx[k][d] = Xk(d); scen[d] += Xk(d); }
                }
                for (int d = 0; d < 3; d++) scen[d] /= npsS;     // slave facet centroid
                // tie proximity tolerance 0.5 sqrt(A_s) (A_s from the facet clipped
                // with itself); < 0: the facet is degenerate, no exclusion
                double tieGapTol = -1.0;
                if (mc.isTie) {
                    double nS[3];
                    ladrunoFacetNormalNewell(npsS, Sx, nS);
                    LadrunoMortarKernel::PairResult prSelf;
                    if (LadrunoMortarKernel::integratePair(npsS, Sx, npsS, Sx, nS, prSelf,
                                                           mc.ngp) == 0 && prSelf.area > 1e-300)
                        tieGapTol = 0.5 * std::sqrt(prSelf.area);
                }

                for (int seg = 0; seg < nSegM; seg++) {
                    Node *mNodes[4]; double mcen[3] = {0.0, 0.0, 0.0};
                    double Mx[4][3];                              // reference coordinates
                    for (int k = 0; k < npsM; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * npsM + k));
                        if (mn == 0 || mn->getNumberDOF() < 3) {
                            opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << " master facet " << seg << " node "
                                   << mTags(seg * npsM + k) << ": ndf="
                                   << (mn != 0 ? mn->getNumberDOF() : -1) << " < 3 (or node "
                                      "absent); master facets need three translational DOFs "
                                      "per node\n";
                            return ladrunoContactFatal();
                        }
                        mNodes[k] = mn;
                        const Vector &Xk = mn->getCrds();
                        for (int d = 0; d < 3; d++) { Mx[k][d] = Xk(d); mcen[d] += Xk(d); }
                    }
                    for (int d = 0; d < 3; d++) mcen[d] /= npsM;
                    // a pair owned by the edge-to-edge lane is not a mortar pair
                    if (mc.edgeEdge &&
                        ladrunoEdgeEdgeOwns(npsS, Sx, npsM, Mx,
                                            ladrunoEdgeBand(npsM, Mx, mc.edgeBand)))
                        continue;
                    // orientation toward the slave's allowed half-space: the explicit
                    // -outward, else (slave centroid - master centroid) at the
                    // reference configuration
                    double orientDir[3];
                    if (mc.hasOutward) {
                        for (int d = 0; d < 3; d++) orientDir[d] = mc.outward[d];
                    } else if (mc.isTie) {
                        // the bond of a tie does not depend on the normal sign
                        for (int d = 0; d < 3; d++) orientDir[d] = scen[d] - mcen[d];
                        double nM[3];
                        ladrunoFacetNormalNewell(npsM, Mx, nM);
                        if (nM[0]*orientDir[0] + nM[1]*orientDir[1] + nM[2]*orientDir[2] == 0.0)
                            for (int d = 0; d < 3; d++) orientDir[d] = nM[d];
                    } else {
                        for (int d = 0; d < 3; d++) orientDir[d] = scen[d] - mcen[d];
                        // The kernel flips the facet normal toward orientDir. On a
                        // coincident pair, and on a zero-gap non-matching pair whose
                        // centroid offset is tangential, the normal component of the
                        // automatic direction vanishes and the sign is decided by
                        // round-off: the contact may then be attractive. Reported once
                        // per contact for a pair closer than two facet sizes.
                        double d2 = orientDir[0]*orientDir[0] + orientDir[1]*orientDir[1]
                                  + orientDir[2]*orientDir[2];
                        double h2 = 0.0;
                        for (int k = 0; k < npsM; k++) {
                            const double *a = Mx[k]; const double *b = Mx[(k + 1) % npsM];
                            double ex = b[0]-a[0], ey = b[1]-a[1], ez = b[2]-a[2];
                            double e2 = ex*ex + ey*ey + ez*ez;
                            if (e2 > h2) h2 = e2;
                        }
                        double nM[3];
                        ladrunoFacetNormalNewell(npsM, Mx, nM);
                        double odn = nM[0]*orientDir[0] + nM[1]*orientDir[1] + nM[2]*orientDir[2];
                        if (odn*odn <= 1e-12 * h2 && d2 <= 4.0 * h2 &&
                            cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_MORTAR_ORIENT)) {
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << ": the automatic orientation is degenerate (the "
                                      "slave-to-master centroid direction has no component along "
                                      "the master facet normal: a coincident or zero-gap "
                                      "non-matching interface); the contact sign then depends on "
                                      "the facet winding and may be attractive. Give -outward.\n";
                        }
                    }
                    // -maxGap: no pair whose slave centroid lies farther than maxGap
                    // from the plane of the master facet, in the configuration of
                    // this handle(); a pairing decision, never switched within a step
                    if (mc.maxGap > 0.0) {
                        double xm[4][3], cs[3] = {0.0, 0.0, 0.0}, cm[3] = {0.0, 0.0, 0.0};
                        for (int k = 0; k < npsS; k++) {
                            const Vector &u = sNodes[k]->getTrialDisp();
                            for (int d = 0; d < 3; d++) cs[d] += (Sx[k][d] + u(d)) / npsS;
                        }
                        for (int k = 0; k < npsM; k++) {
                            const Vector &u = mNodes[k]->getTrialDisp();
                            for (int d = 0; d < 3; d++) {
                                xm[k][d] = Mx[k][d] + u(d);
                                cm[d] += xm[k][d] / npsM;
                            }
                        }
                        double nM[3];
                        ladrunoFacetNormalNewell(npsM, xm, nM);
                        double dist = std::fabs(nM[0] * (cs[0] - cm[0]) + nM[1] * (cs[1] - cm[1])
                                                + nM[2] * (cs[2] - cm[2]));
                        if (dist > mc.maxGap) continue;
                    }
                    if (mc.isTie) {
                        LadrunoMortarKernel::PairResult pr0;
                        if (LadrunoMortarKernel::integratePair(npsS, Sx, npsM, Mx, orientDir,
                                                               pr0, mc.ngp) != 0)
                            continue;                     // no overlap as built: not bonded
                        if (mc.maxGap <= 0.0 && tieGapTol > 0.0 &&
                            pr0.gapL1 > tieGapTol * pr0.area) {
                            nTieExcluded++;               // a distant (antipodal) pair
                            continue;
                        }
                    }
                    double epsUse = epsFixed;
                    if (epsAuto) {
                        epsUse = ladrunoResolveAutoKn(theDomain, mNodes, npsM, orientDir);
                        if (epsUse <= 0.0) {
                            opserr << "FATAL LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << " slave facet " << sf << " master facet " << seg
                                   << ": the automatic penalty could not be sized (no owning "
                                      "solid element, an element without 3 DOFs per node, or an "
                                      "ambiguous normal). Give an explicit -epsN.\n";
                            return ladrunoContactFatal();
                        }
                    }
                    // friction: epsT auto is the normal penalty of the pair; friction
                    // without a tangential penalty would stick with zero force, so
                    // the normal penalty is used then, with a message
                    bool wantFric = mc.hasFriction();
                    double epsTuse = mc.epsTAuto ? epsUse : mc.epsT;
                    if (wantFric && epsTuse <= 0.0) {
                        epsTuse = epsUse;
                        if (cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_EPST))
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << ": friction without -epsT; the tangential "
                                      "penalty is the normal penalty (-epsT auto)\n";
                    }
                    if (wantFric && mc.consistentTan &&
                        cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_CSL))
                        opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                               << mc.tag << ": -consistanttan (non-symmetric friction tangent) "
                                  "needs a non-symmetric solver (system FullGeneral/UmfPack/"
                                  "BandGeneral); symmetric solvers will silently corrupt it.\n";
                    LadrunoContactFE *fe =
                        new LadrunoContactFE(numFe++, sNodes, npsS, mNodes, npsM, epsUse,
                                             orientDir, mc.tag, sf, theDomain, mc.ngp,
                                             mc.isTie, mc.mu, epsTuse, mc.cohesion, mc.tauMax,
                                             mc.consistentTan, seg);
                    if (fe == 0) return -5;
                    if (mc.gapOffset != 0.0 || mc.adjust)
                        fe->setMortarGapShift(mc.gapOffset, mc.adjust, mc.adjustTol);
                    if (mc.smoothN > 0.0 || mc.smoothT > 0.0)
                        fe->setMortarSmoothing(mc.smoothN, mc.smoothT);
                    theModel->addFE_Element(fe);
                    for (int k = 0; k < npsS; k++) {
                        cd->mortarNormalGCMark(mc.tag, sTags(sf * npsS + k));
                        if (wantFric)        // one friction slot per (node, facet pair)
                            cd->mortarFrictionGCMark(mc.tag, sTags(sf * npsS + k), sf, seg);
                    }
                }
            }
            if (nTieExcluded > 0 && cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_TIE_EXCLUDED))
                opserr << "LadrunoContactHandler::handle() - mortar tie " << mc.tag << ": "
                       << nTieExcluded << " overlapping facet pair(s) farther apart than half "
                          "the slave facet size were not tied (the far side of a closed or "
                          "folded surface)\n";
        }
        cd->mortarNormalGCEnd();     // erase the multiplier slots no live pair referenced
    }

    // Edge-to-edge contact (-edgeedge on a mortar contact): over the same
    // exhaustive (slave facet, master facet) enumeration as the mortar lane,
    // every pair owned by the edge lane (ladrunoEdgeEdgeOwns, which the mortar
    // loop above skipped) gets one adapter per pair of boundary edges that
    // cross strictly inside both edges within the band. An edge shared by two
    // facets is built once. The routing uses the reference configuration of
    // this handle(); a pair that slides off its crossing is released by the
    // adapter (large sliding is outside its scope). The mortar loop has already
    // validated the surfaces and their nodes.
    if (cd != 0) {
        cd->edgeGCBegin();
        for (int c = 0; c < cd->getNumMortarContacts(); c++) {
            const LadrunoContactDomain::MortarContact &mc = cd->getMortarContact(c);
            if (mc.retired || !mc.edgeEdge) continue;
            LadrunoContactSurface *ms = cd->getSurface(mc.masterSurfTag);
            LadrunoContactSurface *ss = cd->getSurface(mc.slaveSurfTag);
            if (ms == 0 || ss == 0) continue;
            int npsM = ms->getNodesPerSeg(), npsS = ss->getNodesPerSeg();
            if (npsM < 3 || npsM > 4 || npsS < 3 || npsS > 4) continue;
            bool epsAuto = mc.epsNAuto || mc.knAuto;
            double epsFixed = (mc.epsN > 0.0) ? mc.epsN : mc.kn;
            const ID &mTags = ms->getNodeTags();
            const ID &sTags = ss->getNodeTags();
            int nSegM = mTags.Size() / npsM;
            int nSegS = sTags.Size() / npsS;
            std::set<std::vector<int> > builtEdges;
            for (int sf = 0; sf < nSegS; sf++) {
                Node *sN[4]; bool ok = true; double Scen[3] = {0.0, 0.0, 0.0}, Sx[4][3];
                for (int k = 0; k < npsS; k++) {
                    Node *sn = theDomain->getNode(sTags(sf * npsS + k));
                    if (sn == 0 || sn->getNumberDOF() < 3) { ok = false; break; }
                    sN[k] = sn;
                    const Vector &X = sn->getCrds();
                    for (int d = 0; d < 3; d++) { Sx[k][d] = X(d); Scen[d] += X(d); }
                }
                if (!ok) continue;
                for (int d = 0; d < 3; d++) Scen[d] /= npsS;
                double nS[3]; ladrunoFacetNormalNewell(npsS, Sx, nS);
                for (int seg = 0; seg < nSegM; seg++) {
                    Node *mN[4]; ok = true; double Mcen[3] = {0.0, 0.0, 0.0}, Mx[4][3];
                    for (int k = 0; k < npsM; k++) {
                        Node *mn = theDomain->getNode(mTags(seg * npsM + k));
                        if (mn == 0 || mn->getNumberDOF() < 3) { ok = false; break; }
                        mN[k] = mn;
                        const Vector &X = mn->getCrds();
                        for (int d = 0; d < 3; d++) { Mx[k][d] = X(d); Mcen[d] += X(d); }
                    }
                    if (!ok) continue;
                    for (int d = 0; d < 3; d++) Mcen[d] /= npsM;
                    double dBand = ladrunoEdgeBand(npsM, Mx, mc.edgeBand);
                    if (!ladrunoEdgeEdgeOwns(npsS, Sx, npsM, Mx, dBand)) continue;

                    // reference of the first sign capture: -outward, else the
                    // slave facet normal turned from the master toward the slave
                    // (never zero, unlike a centroid difference of nearly
                    // coincident centroids)
                    double orientDir[3];
                    if (mc.hasOutward) {
                        for (int d = 0; d < 3; d++) orientDir[d] = mc.outward[d];
                    } else {
                        double cdv[3] = {Scen[0]-Mcen[0], Scen[1]-Mcen[1], Scen[2]-Mcen[2]};
                        double dc = nS[0]*cdv[0] + nS[1]*cdv[1] + nS[2]*cdv[2];
                        double sgn = (dc >= 0.0) ? 1.0 : -1.0;
                        for (int d = 0; d < 3; d++) orientDir[d] = sgn * nS[d];
                    }
                    // edge penalty: -edgeKn, else automatic per master facet
                    // (-edgeKn auto, or an automatic mortar penalty), else the
                    // mortar penalty
                    double edgeKnUse = mc.edgeKn;
                    if (edgeKnUse <= 0.0) {
                        if (mc.edgeKnAuto || epsAuto)
                            edgeKnUse = ladrunoResolveAutoKn(theDomain, mN, npsM, orientDir);
                        else
                            edgeKnUse = epsFixed;
                    }
                    if (edgeKnUse <= 0.0) {
                        opserr << "FATAL LadrunoContactHandler::handle() - mortar contact " << mc.tag
                               << " slave facet " << sf << " master facet " << seg
                               << ": the edge-to-edge penalty could not be sized (no owning solid "
                                  "element, an element without 3 DOFs per node, or an ambiguous "
                                  "normal). Give an explicit -edgeKn.\n";
                        return ladrunoContactFatal();
                    }
                    // friction without a tangential penalty would stick with
                    // zero force: the edge penalty is used then, with a message
                    double edgeKtUse = mc.edgeKt;
                    if (mc.hasEdgeFriction() && edgeKtUse <= 0.0) {
                        edgeKtUse = edgeKnUse;
                        if (cd->warnOnce(mc.tag, LadrunoContactDomain::WARN_EDGEKT))
                            opserr << "WARNING LadrunoContactHandler::handle() - mortar contact "
                                   << mc.tag << ": edge-to-edge friction without -edgeKt; the "
                                      "tangential penalty is the edge normal penalty\n";
                    }
                    for (int ie = 0; ie < npsS; ie++) {
                        Node *sa = sN[ie], *sb = sN[(ie + 1) % npsS];
                        for (int je = 0; je < npsM; je++) {
                            Node *ma = mN[je], *mb = mN[(je + 1) % npsM];
                            LadrunoEdgeKernel::ClosestResult cr = LadrunoEdgeKernel::closestPtSegSeg(
                                Sx[ie], Sx[(ie + 1) % npsS], Mx[je], Mx[(je + 1) % npsM]);
                            if (cr.status != LadrunoEdgeKernel::EE_OK || !cr.interior) continue;
                            double w[3] = {cr.c1[0]-cr.c2[0], cr.c1[1]-cr.c2[1], cr.c1[2]-cr.c2[2]};
                            if (std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]) > dBand) continue;
                            int sat = sa->getTag(), sbt = sb->getTag();
                            int mat = ma->getTag(), mbt = mb->getTag();
                            std::vector<int> key(4);
                            key[0] = (sat < sbt) ? sat : sbt; key[1] = (sat < sbt) ? sbt : sat;
                            key[2] = (mat < mbt) ? mat : mbt; key[3] = (mat < mbt) ? mbt : mat;
                            if (!builtEdges.insert(key).second) continue;
                            LadrunoContactFE *fe =
                                new LadrunoContactFE(numFe++, sa, sb, ma, mb, edgeKnUse, orientDir,
                                                     mc.tag, theDomain, mc.edgeMu, edgeKtUse,
                                                     mc.edgeCohesion, mc.edgeTauMax,
                                                     mc.edgeConsistentTan, mc.edgeAlm);
                            if (fe == 0) return -5;
                            theModel->addFE_Element(fe);
                            cd->edgeGCMark(mc.tag, sat, sbt, mat, mbt);
                        }
                    }
                }
            }
        }
        cd->edgeGCEnd();             // erase the edge slots no live pair referenced
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
