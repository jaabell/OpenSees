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
// Description: LadrunoContactFE, see LadrunoContactFE.h.

#include "LadrunoContactFE.h"
#include <cmath>
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <StaticIntegrator.h>
#include <Domain.h>
#include <LadrunoContactDomain.h>
#include <LadrunoContactProjection.h>
#include <LadrunoFrictionKernel.h>
#include <LadrunoMortarKernel.h>
#include <LadrunoEdgeKernel.h>
#include <LadrunoContact2DKernel.h>

// penalty normal traction tn = kn <-gap>+ (Macaulay bracket): zero for gap >= 0,
// so the contact never pulls.
static inline double
ladrunoContactTraction(double kn, double gap)
{
    return (gap < 0.0) ? kn * (-gap) : 0.0;
}

// tributary length below which a 2D mortar slave node carries nothing,
// relative to Lref
static const double MORT2D_TAU_TRIB = 1e-12;

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, int ndm_,
                                   const double p0[3], const double n[3], double kn_,
                                   Domain *dom, int contactTag_, double muc_)
  : FE_Element(tag, /*numDOF_Group=*/1, /*ndof=*/ndm_),
    resid(ndm_), tang(ndm_, ndm_),
    mode(RIGID_PLANE), theSlave(slaveNode), ndm(ndm_), kn(kn_), nps(0),
    theDomain(dom), contactTag(contactTag_), segIndex(0), consistentNormal(false),
    kt(0.0), mu(0.0), cohesion(0.0), tauMax(0.0), consistentTan(false),
    muc(muc_ > 0.0 ? muc_ : 0.0), npsS(0), npsM(0), slaveFacetIndex(0), mortarOrder(2),
    isTie(false), masterFacetIndex(0), hasGapShift(false), gapShiftOffset(0.0),
    gapShiftAdjust(false), gapShiftAdjTol(0.0), gapShiftRefReady(false),
    smoothBand(0.0), smoothSlip(0.0), edgeAlm(false)
{
    for (int i = 0; i < 4; i++) edgeNode[i] = 0;
    // Connectivity = the slave node's DOF_Group; setID() fills myID with its
    // first ndm equation numbers (the translations), as for PenaltySP_FE.
    if (slaveNode != 0) {
        DOF_Group *dg = slaveNode->getDOF_GroupPtr();
        if (dg != 0)
            myDOF_Groups(0) = dg->getTag();
    }
    for (int d = 0; d < 3; d++) { planeP0[d] = p0[d]; planeN[d] = n[d]; orientDir[d] = 0.0; }
    for (int i = 0; i < 4; i++) { segNode[i] = 0; mortarSlave[i] = 0; mortarMaster[i] = 0; }
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes,
                                   int nps_, double kn_, const double odir[3],
                                   double kt_, double mu_, Domain *dom,
                                   int contactTag_, int segIndex_, bool consistentTan_,
                                   bool consistentNormal_, double cohesion_, double tauMax_,
                                   double muc_)
  : FE_Element(tag, /*numDOF_Group=*/1 + nps_, /*ndof=*/3 * (1 + nps_)),
    resid(3 * (1 + nps_)), tang(3 * (1 + nps_), 3 * (1 + nps_)),
    mode(SEGMENT), theSlave(slaveNode), ndm(3), kn(kn_), nps(nps_),
    theDomain(dom), contactTag(contactTag_), segIndex(segIndex_),
    consistentNormal(consistentNormal_),
    kt(kt_), mu(mu_), cohesion(cohesion_), tauMax(tauMax_), consistentTan(consistentTan_),
    muc(muc_ > 0.0 ? muc_ : 0.0), npsS(0), npsM(0), slaveFacetIndex(0), mortarOrder(2),
    isTie(false), masterFacetIndex(0), hasGapShift(false), gapShiftOffset(0.0),
    gapShiftAdjust(false), gapShiftAdjTol(0.0), gapShiftRefReady(false),
    smoothBand(0.0), smoothSlip(0.0), edgeAlm(false)
{
    for (int i = 0; i < 4; i++) edgeNode[i] = 0;
    // Connectivity = slave DOF_Group + each segment-node DOF_Group; setID() fills
    // myID = [slave xyz | seg_1 xyz | ... | seg_nps xyz], the layout B assumes.
    if (slaveNode != 0) {
        DOF_Group *dg = slaveNode->getDOF_GroupPtr();
        if (dg != 0) myDOF_Groups(0) = dg->getTag();
    }
    for (int i = 0; i < nps; i++) {
        segNode[i] = segNodes[i];
        if (segNodes[i] != 0) {
            DOF_Group *dg = segNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(1 + i) = dg->getTag();
        }
    }
    for (int i = nps; i < 4; i++) segNode[i] = 0;
    for (int i = 0; i < 4; i++) { mortarSlave[i] = 0; mortarMaster[i] = 0; }
    for (int d = 0; d < 3; d++) { planeP0[d] = 0.0; planeN[d] = 0.0; }
    // a fixed direction, not the live slave position, so the orientation stays
    // correct after the slave penetrates
    for (int d = 0; d < 3; d++) orientDir[d] = odir[d];
}

LadrunoContactFE::LadrunoContactFE(int tag, Node **slaveNodes, int nps_s,
                                   Node **masterNodes, int nps_m, double epsN,
                                   const double odir[3], int contactTag_, int slaveFacetIndex_,
                                   Domain *dom, int ngp_, bool isTie_, double mu_,
                                   double epsT_, double cohesion_, double tauMax_,
                                   bool consistentTan_, int masterFacetIndex_)
  : FE_Element(tag, /*numDOF_Group=*/nps_s + nps_m, /*ndof=*/3 * (nps_s + nps_m)),
    resid(3 * (nps_s + nps_m)), tang(3 * (nps_s + nps_m), 3 * (nps_s + nps_m)),
    mode(MORTAR), theSlave(0), ndm(3), kn(epsN), nps(0),
    theDomain(dom), contactTag(contactTag_), segIndex(0), consistentNormal(false),
    kt(epsT_), mu(mu_), cohesion(cohesion_), tauMax(tauMax_), consistentTan(consistentTan_),
    muc(0.0), npsS(nps_s), npsM(nps_m), slaveFacetIndex(slaveFacetIndex_),
    mortarOrder(ngp_ >= 1 ? ngp_ : 2), isTie(isTie_), masterFacetIndex(masterFacetIndex_),
    hasGapShift(false), gapShiftOffset(0.0), gapShiftAdjust(false), gapShiftAdjTol(0.0),
    gapShiftRefReady(false), smoothBand(0.0), smoothSlip(0.0), edgeAlm(false)
{
    for (int i = 0; i < 4; i++) edgeNode[i] = 0;
    for (int i = 0; i < 4; i++) gapShiftRef[i] = 0.0;
    // Connectivity = each slave facet DOF_Group, then each master facet
    // DOF_Group; setID() fills myID = [slave_0 xyz | ... | master_0 xyz | ...],
    // the [slave block | master block] layout of the residual and tangent.
    for (int i = 0; i < nps_s; i++) {
        mortarSlave[i] = slaveNodes[i];
        if (slaveNodes[i] != 0) {
            DOF_Group *dg = slaveNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(i) = dg->getTag();
        }
    }
    for (int i = 0; i < nps_m; i++) {
        mortarMaster[i] = masterNodes[i];
        if (masterNodes[i] != 0) {
            DOF_Group *dg = masterNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(nps_s + i) = dg->getTag();
        }
    }
    for (int i = nps_s; i < 4; i++) mortarSlave[i] = 0;
    for (int i = nps_m; i < 4; i++) mortarMaster[i] = 0;
    for (int i = 0; i < 4; i++) segNode[i] = 0;
    for (int d = 0; d < 3; d++) { orientDir[d] = odir[d]; planeP0[d] = 0.0; planeN[d] = 0.0; }
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *sNodeA, Node *sNodeB, Node *mNodeA,
                                   Node *mNodeB, double epsN, const double odir[3],
                                   int contactTag_, Domain *dom, double mu_, double kt_,
                                   double cohesion_, double tauMax_, bool consistentTan_,
                                   bool edgeAlm_)
  : FE_Element(tag, /*numDOF_Group=*/4, /*ndof=*/12),
    resid(12), tang(12, 12),
    mode(EDGE_EDGE), theSlave(0), ndm(3), kn(epsN), nps(0),
    theDomain(dom), contactTag(contactTag_), segIndex(0), consistentNormal(false),
    kt(kt_), mu(mu_), cohesion(cohesion_), tauMax(tauMax_), consistentTan(consistentTan_),
    muc(0.0), npsS(0), npsM(0), slaveFacetIndex(0), mortarOrder(2), isTie(false),
    masterFacetIndex(0), hasGapShift(false), gapShiftOffset(0.0), gapShiftAdjust(false),
    gapShiftAdjTol(0.0), gapShiftRefReady(false), smoothBand(0.0), smoothSlip(0.0),
    edgeAlm(edgeAlm_)
{
    // Connectivity = the four edge-node DOF_Groups; setID() fills
    // myID = [sa xyz | sb xyz | ma xyz | mb xyz], the layout of edgeGeom()
    edgeNode[0] = sNodeA; edgeNode[1] = sNodeB; edgeNode[2] = mNodeA; edgeNode[3] = mNodeB;
    for (int i = 0; i < 4; i++)
        if (edgeNode[i] != 0) {
            DOF_Group *dg = edgeNode[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(i) = dg->getTag();
        }
    for (int i = 0; i < 4; i++) { segNode[i] = 0; mortarSlave[i] = 0; mortarMaster[i] = 0; }
    for (int i = 0; i < 4; i++) gapShiftRef[i] = 0.0;
    for (int d = 0; d < 3; d++) { orientDir[d] = odir[d]; planeP0[d] = 0.0; planeN[d] = 0.0; }
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes, int nps_,
                                   int ndm2, double kn_, double sigma, double LrefIn,
                                   Node *prevFar, Node *nextFar, double kt_, double mu_,
                                   bool consistentTan_, double muc_, Domain *dom,
                                   int contactTag_, int segIndex_, double cohesion_,
                                   double tauMax_)
  : FE_Element(tag, /*numDOF_Group=*/1 + nps_, /*ndof=*/ndm2 * (1 + nps_)),
    resid(ndm2 * (1 + nps_)), tang(ndm2 * (1 + nps_), ndm2 * (1 + nps_)),
    mode(SEGMENT), theSlave(slaveNode), ndm(ndm2), kn(kn_), nps(nps_),
    theDomain(dom), contactTag(contactTag_), segIndex(segIndex_), consistentNormal(false),
    kt(kt_), mu(mu_), cohesion(cohesion_), tauMax(tauMax_), consistentTan(consistentTan_),
    muc(muc_ > 0.0 ? muc_ : 0.0), npsS(0), npsM(0), slaveFacetIndex(0), mortarOrder(2),
    isTie(false), masterFacetIndex(0), hasGapShift(false), gapShiftOffset(0.0),
    gapShiftAdjust(false), gapShiftAdjTol(0.0), gapShiftRefReady(false),
    smoothBand(0.0), smoothSlip(0.0), edgeAlm(false)
{
    // Connectivity = slave DOF_Group + the segment (or vertex) DOF_Groups;
    // setID() fills myID with the first two equations of each node
    if (slaveNode != 0) {
        DOF_Group *dg = slaveNode->getDOF_GroupPtr();
        if (dg != 0) myDOF_Groups(0) = dg->getTag();
    }
    for (int i = 0; i < nps; i++) {
        segNode[i] = segNodes[i];
        if (segNodes[i] != 0) {
            DOF_Group *dg = segNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(1 + i) = dg->getTag();
        }
    }
    for (int i = nps; i < 4; i++) segNode[i] = 0;
    for (int i = 0; i < 4; i++) { mortarSlave[i] = 0; mortarMaster[i] = 0; edgeNode[i] = 0; }
    for (int i = 0; i < 4; i++) gapShiftRef[i] = 0.0;
    for (int d = 0; d < 3; d++) { planeP0[d] = 0.0; planeN[d] = 0.0; orientDir[d] = 0.0; }
    nts2dSigma = sigma;
    nts2dLref  = LrefIn;
    nts2dHasPrevFar = (prevFar != 0);
    nts2dPrevFarTag = (prevFar != 0) ? prevFar->getTag() : 0;
    nts2dHasNextFar = (nextFar != 0);
    nts2dNextFarTag = (nextFar != 0) ? nextFar->getTag() : 0;
}

LadrunoContactFE::LadrunoContactFE(int tag, Node **slaveNodes, int nps_s, Node **masterNodes,
                                   int nps_m, int ndm2, double epsN, double sigma,
                                   double LrefIn, int contactTag_, int slaveFacetIndex_,
                                   Domain *dom, bool isTie_, double mu_, double epsT_,
                                   double cohesion_, double tauMax_, bool consistentTan_,
                                   int masterFacetIndex_, bool ntsCoDeclared_)
  : FE_Element(tag, /*numDOF_Group=*/nps_s + nps_m, /*ndof=*/ndm2 * (nps_s + nps_m)),
    resid(ndm2 * (nps_s + nps_m)), tang(ndm2 * (nps_s + nps_m), ndm2 * (nps_s + nps_m)),
    mode(MORTAR), theSlave(0), ndm(ndm2), kn(epsN), nps(0),
    theDomain(dom), contactTag(contactTag_), segIndex(0), consistentNormal(false),
    kt(epsT_), mu(mu_), cohesion(cohesion_), tauMax(tauMax_), consistentTan(consistentTan_),
    muc(0.0), npsS(nps_s), npsM(nps_m), slaveFacetIndex(slaveFacetIndex_), mortarOrder(2),
    isTie(isTie_), masterFacetIndex(masterFacetIndex_), hasGapShift(false),
    gapShiftOffset(0.0), gapShiftAdjust(false), gapShiftAdjTol(0.0), gapShiftRefReady(false),
    smoothBand(0.0), smoothSlip(0.0), edgeAlm(false)
{
    for (int i = 0; i < nps_s; i++) {
        mortarSlave[i] = slaveNodes[i];
        if (slaveNodes[i] != 0) {
            DOF_Group *dg = slaveNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(i) = dg->getTag();
        }
    }
    for (int i = 0; i < nps_m; i++) {
        mortarMaster[i] = masterNodes[i];
        if (masterNodes[i] != 0) {
            DOF_Group *dg = masterNodes[i]->getDOF_GroupPtr();
            if (dg != 0) myDOF_Groups(nps_s + i) = dg->getTag();
        }
    }
    for (int i = nps_s; i < 4; i++) mortarSlave[i] = 0;
    for (int i = nps_m; i < 4; i++) mortarMaster[i] = 0;
    for (int i = 0; i < 4; i++) { segNode[i] = 0; edgeNode[i] = 0; gapShiftRef[i] = 0.0; }
    for (int d = 0; d < 3; d++) { orientDir[d] = 0.0; planeP0[d] = 0.0; planeN[d] = 0.0; }
    mortar2dSigma = sigma;
    mortar2dLref  = LrefIn;
    mortar2dNtsCoDeclared = ntsCoDeclared_;
}

LadrunoContactFE::~LadrunoContactFE()
{
}

void
LadrunoContactFE::setMortarGapShift(double gapOffset, bool adjust, double adjustTol)
{
    gapShiftOffset = gapOffset;
    gapShiftAdjust = adjust;
    gapShiftAdjTol = (adjustTol > 0.0) ? adjustTol : 0.0;
    hasGapShift = (gapOffset != 0.0) || adjust;
    gapShiftRefReady = false;
    for (int i = 0; i < 4; i++) gapShiftRef[i] = 0.0;
}

void
LadrunoContactFE::setMortarSmoothing(double smoothN, double smoothT)
{
    smoothBand = (smoothN > 0.0) ? smoothN : 0.0;
    smoothSlip = (smoothT > 0.0) ? smoothT : 0.0;
}

void
LadrunoContactFE::mortarSmoothLaw(double pr, double &P, double &dP, double &chi,
                                  double &dchi) const
{
    double x = -pr, S = kn * smoothBand;
    if (x <= -S) { P = 0.0; dP = 0.0; chi = 0.0; dchi = 0.0; return; }
    if (x >= S)  { P = x;   dP = 1.0; chi = 1.0; dchi = 0.0; return; }
    double t = (x + S) / (2.0 * S);
    P = (x + S) * (x + S) / (4.0 * S);
    dP = t;
    chi = t * t * (3.0 - 2.0 * t);
    dchi = 6.0 * t * (1.0 - t) / (2.0 * S);
}

double
LadrunoContactFE::explicitTimeStep(Integrator *)
{
    return 0.0;   // no integrator in this build reports its explicit time step
}

int
LadrunoContactFE::setID(void)
{
    Node *nodes[8]; int nn = 0;
    switch (mode) {
    case RIGID_PLANE:
        nodes[nn++] = theSlave;
        break;
    case SEGMENT:
        nodes[nn++] = theSlave;
        for (int i = 0; i < nps; i++) nodes[nn++] = segNode[i];
        break;
    case MORTAR:
        for (int i = 0; i < npsS; i++) nodes[nn++] = mortarSlave[i];
        for (int i = 0; i < npsM; i++) nodes[nn++] = mortarMaster[i];
        break;
    case EDGE_EDGE:
        for (int i = 0; i < 4; i++) nodes[nn++] = edgeNode[i];
        break;
    default:
        return FE_Element::setID();
    }
    const int per  = ndm;
    const int ndof = myID.Size();
    if (nn != myDOF_Groups.Size() || per * nn != ndof) {
        opserr << "FATAL LadrunoContactFE::setID() - layout mismatch: " << nn << " nodes x "
               << per << " != ndof " << ndof << " (" << myDOF_Groups.Size() << " DOF_Groups)\n";
        return -3;
    }
    int current = 0;
    for (int i = 0; i < nn; i++) {
        DOF_Group *dg = (nodes[i] != 0) ? nodes[i]->getDOF_GroupPtr() : 0;
        if (dg == 0) {
            opserr << "WARNING LadrunoContactFE::setID() - 0 DOF_Group pointer\n";
            return -2;
        }
        const ID &eq = dg->getID();
        if (eq.Size() < per) {
            opserr << "FATAL LadrunoContactFE::setID() - node " << nodes[i]->getTag()
                   << " has " << eq.Size() << " equations < ndm=" << per << "\n";
            return -3;
        }
        for (int j = 0; j < per; j++) myID(current++) = eq(j);
    }
    return 0;
}

bool
LadrunoContactFE::viscousActive(Integrator *theIntegrator) const
{
    if (muc <= 0.0) return false;
    if (dynamic_cast<StaticIntegrator *>(theIntegrator) == 0) return true;   // transient: active
    LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
    if (cd == 0 || cd->warnOnce(contactTag, LadrunoContactDomain::WARN_VISC_STATIC))
        opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": -visc (viscous "
                  "contact stabilization) is disabled under a static integrator; nodal "
                  "velocities are not rates in a static analysis and the damping tangent "
                  "is not assembled there. The dashpot acts again under a transient integrator.\n";
    return false;
}

double
LadrunoContactFE::rigidPlaneGap(void) const
{
    // g = n.(x_s - p0), x_s = X_s + u_s (first ndm components)
    const Vector &X = theSlave->getCrds();
    const Vector &u = theSlave->getTrialDisp();
    double g = 0.0;
    for (int d = 0; d < ndm; d++)
        g += planeN[d] * (X(d) + u(d) - planeP0[d]);
    return g;
}

bool
LadrunoContactFE::segmentActive(double &gap, double n[3], double N[4], double *B,
                                double *gTvec) const
{
    if (mode != SEGMENT || theSlave == 0) return false;
    double Xseg[4][3], xs[3];
    const Vector &Xs = theSlave->getCrds();
    const Vector &us = theSlave->getTrialDisp();
    for (int d = 0; d < 3; d++) xs[d] = Xs(d) + us(d);
    for (int i = 0; i < nps; i++) {
        if (segNode[i] == 0) return false;
        const Vector &Xi = segNode[i]->getCrds();
        const Vector &ui = segNode[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xseg[i][d] = Xi(d) + ui(d);
    }
    if (!LadrunoContactProjection::evalSegment(nps, Xseg, xs, orientDir, gap, n, N))
        return false;
    // gap operator B (1 x ndof) over [u_s | u_1..u_nps]: [ n^T | -N_i n^T ]
    int ndof = 3 * (1 + nps);
    for (int k = 0; k < ndof; k++) B[k] = 0.0;
    for (int d = 0; d < 3; d++) B[d] = n[d];
    for (int i = 0; i < nps; i++)
        for (int d = 0; d < 3; d++) B[3 * (1 + i) + d] = -N[i] * n[d];
    // Tangential relative slip at the contact point. The closest-point
    // projection makes (x_s - xbar) parallel to n, so positions carry no
    // tangential information; the slip is the slave displacement minus the
    // master displacement interpolated at the projection, d = u_s - sum N_i u_i,
    // tangential part (the caller subtracts the engagement origin gT0). Same
    // trial configuration and same shape weights as the gap.
    if (gTvec != 0) {
        double ubar[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < nps; i++) {
            const Vector &ui = segNode[i]->getTrialDisp();
            for (int d = 0; d < 3; d++) ubar[d] += N[i] * ui(d);
        }
        double drel[3] = { us(0)-ubar[0], us(1)-ubar[1], us(2)-ubar[2] };
        LadrunoFrictionKernel::tangentPart(drel, n, gTvec);
    }
    return true;
}

void
LadrunoContactFE::addFrictionTang(double fact, const double n[3], const double N[4],
                                  double tn, const double gTeff[3], const double gpT[3],
                                  bool consistent)
{
    // K_fric = G^T K_ss G, G = [I | -N_i I]: block (a,b) of the ndof x ndof
    // tangent is w_a w_b K_ss, w = [1, -N_0, ..., -N_{nps-1}] over [slave, segment nodes]
    double Kss[3][3];
    LadrunoFrictionKernel::frictionTangentBlock(gTeff, gpT, n, tn, kn, kt, mu, consistent, Kss,
                                                cohesion, tauMax);
    double w[5];
    w[0] = 1.0;
    for (int i = 0; i < nps; i++) w[1 + i] = -N[i];
    int nn = 1 + nps;
    for (int a = 0; a < nn; a++)
        for (int b = 0; b < nn; b++) {
            double wab = fact * w[a] * w[b];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    tang(3 * a + i, 3 * b + j) += wab * Kss[i][j];
        }
}

void
LadrunoContactFE::addNormalGeomTang(double fact)
{
    if (mode != SEGMENT || theSlave == 0) return;
    double Xseg[4][3], xs[3];
    const Vector &Xs = theSlave->getCrds();
    const Vector &us = theSlave->getTrialDisp();
    for (int d = 0; d < 3; d++) xs[d] = Xs(d) + us(d);
    for (int i = 0; i < nps; i++) {
        if (segNode[i] == 0) return;
        const Vector &Xi = segNode[i]->getCrds();
        const Vector &ui = segNode[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xseg[i][d] = Xi(d) + ui(d);
    }
    double Kg[15][15];
    if (!LadrunoContactProjection::normalGeomTangent(nps, Xseg, xs, orientDir, kn, Kg))
        return;                                 // not penetrating or out of bounds
    int ndof = 3 * (1 + nps);
    for (int a = 0; a < ndof; a++)
        for (int b = 0; b < ndof; b++)
            tang(a, b) += fact * Kg[a][b];
}

bool
LadrunoContactFE::mortarActive(double D[4][4], double M[4][4], double g[4], double n[3]) const
{
    if (mode != MORTAR) return false;
    double Xs[4][3], Xm[4][3];
    for (int i = 0; i < npsS; i++) {
        if (mortarSlave[i] == 0) return false;
        const Vector &X = mortarSlave[i]->getCrds();
        const Vector &u = mortarSlave[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xs[i][d] = X(d) + u(d);
    }
    for (int i = 0; i < npsM; i++) {
        if (mortarMaster[i] == 0) return false;
        const Vector &X = mortarMaster[i]->getCrds();
        const Vector &u = mortarMaster[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) Xm[i][d] = X(d) + u(d);
    }
    LadrunoMortarKernel::PairResult pr;
    if (LadrunoMortarKernel::integratePair(npsS, Xs, npsM, Xm, orientDir, pr, mortarOrder) != 0)
        return false;                                // empty or degenerate overlap
    for (int i = 0; i < 4; i++) {
        g[i] = pr.g[i];
        for (int j = 0; j < 4; j++) { D[i][j] = pr.D[i][j]; M[i][j] = pr.M[i][j]; }
    }
    // Gap shift, applied here so that every consumer sees one gap:
    // gt_I <- (gt_I / a_I + s_I) a_I. At the reference configuration the adjust
    // term cancels the gap exactly (the same double minus itself), so an
    // adjusted contact starts with zero pressure, not a round-off residue.
    if (hasGapShift) {
        if (gapShiftAdjust && !gapShiftRefReady) {
            // reference gap of this pair; an empty reference overlap leaves the
            // pair unadjusted (it was not in contact as meshed)
            double Xs0[4][3], Xm0[4][3];
            for (int i = 0; i < npsS; i++) {
                const Vector &X = mortarSlave[i]->getCrds();
                for (int d = 0; d < 3; d++) Xs0[i][d] = X(d);
            }
            for (int i = 0; i < npsM; i++) {
                const Vector &X = mortarMaster[i]->getCrds();
                for (int d = 0; d < 3; d++) Xm0[i][d] = X(d);
            }
            LadrunoMortarKernel::PairResult pr0;
            if (LadrunoMortarKernel::integratePair(npsS, Xs0, npsM, Xm0, orientDir, pr0,
                                                   mortarOrder) == 0) {
                for (int I = 0; I < npsS; I++) {
                    double a0 = 0.0;
                    for (int J = 0; J < npsS; J++) a0 += pr0.D[I][J];
                    if (a0 <= 1e-300) continue;
                    double gbar0 = pr0.g[I] / a0;
                    if (gapShiftAdjTol <= 0.0 || std::fabs(gbar0) <= gapShiftAdjTol)
                        gapShiftRef[I] = -gbar0;
                }
            }
            gapShiftRefReady = true;
        }
        for (int I = 0; I < npsS; I++) {
            double a = 0.0;
            for (int J = 0; J < npsS; J++) a += D[I][J];
            if (a <= 1e-300) continue;
            double gbar = g[I] / a + gapShiftRef[I];
            g[I] = (gbar + gapShiftOffset) * a;
        }
    }
    // the master facet normal oriented by orientDir, the normal the weighted
    // gap was measured along (flat facets)
    return LadrunoMortarKernel::facetNormal(npsM, Xm, orientDir, n);
}

void
LadrunoContactFE::addMortarTang(double fact, bool initialStiff)
{
    double D[4][4], M[4][4], g[4], n[3];
    if (!mortarActive(D, M, g, n)) return;
    int nN = npsS + npsM;
    LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
    if (isTie) {
        // mesh tie: K = epsTie Bt^T diag(1 / a_I) Bt (x) I3, every slave node
        // active and all three components tied; symmetric positive semidefinite,
        // constant (the variation of D and M is not included)
        double Wt[4] = {0.0, 0.0, 0.0, 0.0};
        for (int I = 0; I < npsS; I++) {
            double aFacet = 0.0;
            for (int J = 0; J < npsS; J++) aFacet += D[I][J];
            if (aFacet > 1e-300) Wt[I] = 1.0 / aFacet;
        }
        for (int A = 0; A < nN; A++)
            for (int B = 0; B < nN; B++) {
                double Ks = 0.0;
                for (int I = 0; I < npsS; I++) {
                    double bIA = (A < npsS) ? D[I][A] : -M[I][A - npsS];
                    double bIB = (B < npsS) ? D[I][B] : -M[I][B - npsS];
                    Ks += bIA * Wt[I] * bIB;
                }
                Ks *= kn;
                if (Ks == 0.0) continue;
                for (int d = 0; d < 3; d++)
                    tang(3 * A + d, 3 * B + d) += fact * Ks;
            }
        return;
    }
    // W_I = act_I / a_I with the same per-pair pressure p_I = min(0, lambda_I +
    // epsN gbar_I) as the residual, so the tangent is the derivative of the
    // augmented residual on its active set; lambda_I is constant within the
    // iterations of a step.
    double W[4] = {0.0, 0.0, 0.0, 0.0};
    for (int I = 0; I < npsS; I++) {
        double aFacet = 0.0;
        for (int J = 0; J < npsS; J++) aFacet += D[I][J];
        if (aFacet <= 1e-300) continue;
        double lambdaI = (cd != 0)
            ? cd->getOrCreateMortarNormalState(contactTag, mortarSlave[I]->getTag()).lambdaN
            : 0.0;
        double pr = lambdaI + kn * (g[I] / aFacet);
        if (smoothBand > 0.0) {                       // W_I = P'(x) / a_I, 0..1 over the band
            double P, dP, ch, dch;
            mortarSmoothLaw(pr, P, dP, ch, dch);
            W[I] = dP / aFacet;
            continue;
        }
        if (pr < 0.0) W[I] = 1.0 / aFacet;            // active in compression
        // a shifted node exactly at p = 0 (every adjusted node at the start)
        // takes the closed branch of the kink, where the residual is zero on
        // both sides, so the first iteration sees the interface stiffness
        else if (hasGapShift && pr == 0.0) W[I] = 1.0 / aFacet;
    }
    for (int A = 0; A < nN; A++) {
        for (int B = 0; B < nN; B++) {
            double Ks = 0.0;
            for (int I = 0; I < npsS; I++) {
                double bIA = (A < npsS) ? D[I][A] : -M[I][A - npsS];
                double bIB = (B < npsS) ? D[I][B] : -M[I][B - npsS];
                Ks += bIA * W[I] * bIB;
            }
            Ks *= kn;
            if (Ks == 0.0) continue;
            for (int dA = 0; dA < 3; dA++)
                for (int dB = 0; dB < 3; dB++)
                    tang(3 * A + dA, 3 * B + dB) += fact * Ks * n[dA] * n[dB];
        }
    }

    // Friction tangent K[A][B] += sum_I (b_IA b_IB / a_I) K_ss,I, b the row of
    // Bt = [D, -M], K_ss the 3 x 3 block of LadrunoFrictionKernel at the slip
    // of the residual and the committed state. Symmetric by default; with
    // consistentTan the pressure coupling enters through the normal-gap
    // operator (kn = epsN), non-symmetric. The initial-stiffness path uses the
    // stick block (gTeff = gpT), positive semidefinite.
    if (!hasFriction() || cd == 0)
        return;
    double us[4][3], um[4][3];
    for (int i = 0; i < npsS; i++) {
        const Vector &u = mortarSlave[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) us[i][d] = u(d);
    }
    for (int i = 0; i < npsM; i++) {
        const Vector &u = mortarMaster[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) um[i][d] = u(d);
    }
    double invEpsT = (kt > 0.0) ? 1.0 / kt : 0.0;
    for (int I = 0; I < npsS; I++) {
        double aFacet = 0.0;
        for (int J = 0; J < npsS; J++) aFacet += D[I][J];
        if (aFacet <= 1e-300) continue;
        const LadrunoContactDomain::MortarNormalState &nst =
            cd->getOrCreateMortarNormalState(contactTag, mortarSlave[I]->getTag());
        double pr = nst.lambdaN + kn * (g[I] / aFacet);
        double N_I = -pr, knN = kn, chiI = 1.0, dchiI = 0.0;
        if (smoothBand > 0.0) {
            double dP;
            mortarSmoothLaw(pr, N_I, dP, chiI, dchiI);
            if (N_I <= 0.0) continue;                 // open beyond the band
            knN = kn * dP;                            // dN/dz of the smoothed law
        } else if (pr >= 0.0)
            continue;                                 // friction only in normal contact
        LadrunoContactDomain::MortarFrictionState &st = cd->getOrCreateMortarFrictionState(
            contactTag, mortarSlave[I]->getTag(), slaveFacetIndex, masterFacetIndex);
        double r[3] = {0.0, 0.0, 0.0};
        for (int J = 0; J < npsS; J++)
            for (int d = 0; d < 3; d++) r[d] += D[I][J] * us[J][d];
        for (int K = 0; K < npsM; K++)
            for (int d = 0; d < 3; d++) r[d] -= M[I][K] * um[K][d];
        double rn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
        double gTeff[3];
        for (int d = 0; d < 3; d++)
            gTeff[d] = (r[d] - rn * n[d]) / aFacet - st.gT0[d] + st.lambdaT[d] * invEpsT;
        const double *gtForKss = initialStiff ? st.gpT : gTeff;
        bool useConsistent = consistentTan && !initialStiff;
        double Kss[3][3];
        if (smoothSlip > 0.0)
            LadrunoFrictionKernel::frictionTangentBlockSmooth(gtForKss, st.gpT, n, N_I, knN, kt,
                                                              mu, useConsistent, Kss, cohesion,
                                                              tauMax, smoothSlip);
        else
            LadrunoFrictionKernel::frictionTangentBlock(gtForKss, st.gpT, n, N_I, knN, kt, mu,
                                                        useConsistent, Kss, cohesion, tauMax);
        if (smoothBand > 0.0) {
            // d(chi tF)/du = chi dtF/du + tF (x) dchi/du: the block is scaled by
            // chi, and the consistent path adds the onset coupling
            // kn chi'(x) tF (x) n (non-symmetric; the scatter below supplies
            // dz/du = -b n / a), tF the unscaled traction
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) Kss[i][j] *= chiI;
            if (useConsistent && dchiI > 0.0) {
                double tF[3], gpTr[3];
                if (smoothSlip > 0.0)
                    LadrunoFrictionKernel::frictionReturnMapSmooth(gTeff, st.gpT, N_I, kt, mu, tF,
                                                                   gpTr, cohesion, tauMax,
                                                                   smoothSlip);
                else
                    LadrunoFrictionKernel::frictionReturnMap(gTeff, st.gpT, N_I, kt, mu, tF,
                                                             gpTr, cohesion, tauMax);
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) Kss[i][j] += kn * dchiI * tF[i] * n[j];
            }
        }
        for (int A = 0; A < nN; A++) {
            double bIA = (A < npsS) ? D[I][A] : -M[I][A - npsS];
            if (bIA == 0.0) continue;
            for (int B = 0; B < nN; B++) {
                double bIB = (B < npsS) ? D[I][B] : -M[I][B - npsS];
                double w = fact * bIA * bIB / aFacet;
                if (w == 0.0) continue;
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        tang(3 * A + i, 3 * B + j) += w * Kss[i][j];
            }
        }
    }
}

const Vector &
LadrunoContactFE::getResidual(Integrator *theIntegrator)
{
    resid.Zero();
    if (mode == RIGID_PLANE) {
        double g = rigidPlaneGap();
        if (g < 0.0) {                       // penetration: restoring force along +n
            double tn = ladrunoContactTraction(kn, g);
            for (int d = 0; d < ndm; d++)
                resid(d) = tn * planeN[d];   // drives the slave toward g = 0 (as PenaltySP_FE)
            // viscous normal stabilization f = -muc (n.v_slave) n, opposing the
            // approach rate, while in contact under a transient integrator
            if (theSlave != 0 && viscousActive(theIntegrator)) {
                const Vector &vs = theSlave->getTrialVel();
                double gdot = 0.0;
                for (int d = 0; d < ndm; d++) gdot += planeN[d] * vs(d);
                for (int d = 0; d < ndm; d++) resid(d) += -muc * gdot * planeN[d];
            }
        }
    } else if (mode == SEGMENT && ndm == 2) {
        // 2D node-to-segment: penalty tn = kn <-gap>+ on the pair this adapter
        // owns (segment, concave vertex or end cap), with the vertex side sign
        // committed at the first evaluation that decides it
        double gap, n[3], N2[2], B[6], liveSide = 0.0, gTnow = 0.0;
        const double committedSide = vtx2DCommittedSide();
        bool act2d = segment2DActive(gap, n, N2, B, committedSide, &liveSide, &gTnow);
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (!act2d && cd != 0 && theSlave != 0)
            cd->setNtsForce(contactTag, theSlave->getTag(), segIndex, 0.0);
        // a committed sign that disagrees with the live one means the corner
        // changed type under deformation; the committed sign is kept
        if (nps == 1 && committedSide != 0.0 && liveSide != 0.0 && liveSide != committedSide &&
            cd != 0 && cd->warnOnce(contactTag, LadrunoContactDomain::WARN_VTX2D_SIDE_FLIP))
            opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": the committed "
                      "side sign of a 2D vertex pair differs from the current corner geometry "
                      "(the corner changed type under deformation); the committed sign is kept\n";
        if (act2d) {
            if (nps == 1 && committedSide == 0.0 && liveSide != 0.0 && cd != 0)
                cd->setVtx2DSide(contactTag, theSlave->getTag(), segNode[0]->getTag(), liveSide);
            double tn = ladrunoContactTraction(kn, gap);
            int ndof = 2 * (1 + nps);
            for (int k = 0; k < ndof; k++)
                resid(k) = B[k] * tn;
            if (cd != 0)
                cd->setNtsForce(contactTag, theSlave->getTag(), segIndex, tn);
            if (viscousActive(theIntegrator)) {
                double gdot = 0.0;
                const Vector &vs = theSlave->getTrialVel();
                for (int d = 0; d < 2; d++) gdot += B[d] * vs(d);
                for (int i = 0; i < nps; i++) {
                    const Vector &vm = segNode[i]->getTrialVel();
                    for (int d = 0; d < 2; d++) gdot += B[2 * (1 + i) + d] * vm(d);
                }
                for (int k = 0; k < ndof; k++) resid(k) += -muc * gdot * B[k];
            }
            // scalar friction along t = perp(n) (slot 0 of the friction state)
            if (hasFriction() && cd != 0) {
                LadrunoContactDomain::FrictionState &st =
                    cd->getOrCreateFrictionState(contactTag, theSlave->getTag(), segIndex);
                if (!st.engaged) { st.gT0[0] = gTnow; st.engaged = true; }
                double gTeff = gTnow - st.gT0[0];
                double tFricScalar, gpTtrial0;
                LadrunoFrictionKernel::returnMap1D(gTeff, st.gpT[0], tn, kt, mu, tFricScalar,
                                                   gpTtrial0, cohesion, tauMax);
                st.gpTtrial[0] = gpTtrial0;
                const double th[2] = { -n[1], n[0] };
                resid(0) += tFricScalar * th[0];
                resid(1) += tFricScalar * th[1];
                for (int i = 0; i < nps; i++) {
                    resid(2 * (1 + i) + 0) += -N2[i] * tFricScalar * th[0];
                    resid(2 * (1 + i) + 1) += -N2[i] * tFricScalar * th[1];
                }
            }
        }
    } else if (mode == SEGMENT) {
        double gap, n[3], N[4], B[15], gTvec[3];   // ndof <= 3 (1 + 4) = 15
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (!segmentActive(gap, n, N, B, gTvec)) {
            // a released pair reports zero force; without this the readout kept
            // the last active value until the next handle()
            if (cd != 0 && theSlave != 0)
                cd->setNtsForce(contactTag, theSlave->getTag(), segIndex, 0.0);
        } else {
            double tn = ladrunoContactTraction(kn, gap);   // kn <-gap>+ > 0
            int ndof = 3 * (1 + nps);
            for (int k = 0; k < ndof; k++)
                resid(k) = B[k] * tn;        // r = B^T tn (slave +tn n, master -N_i tn n)

            // normal force readout (ladrunoContactForce); no effect on resid/tang
            if (cd != 0)
                cd->setNtsForce(contactTag, theSlave->getTag(), segIndex, tn);

            // viscous normal stabilization: gdot = B v (normal relative velocity),
            // f = -muc gdot B, opposing the approach, while in contact under a
            // transient integrator (the damping tangent muc B^T B is in addCtoTang)
            if (viscousActive(theIntegrator)) {
                double gdot = 0.0;
                const Vector &vs = theSlave->getTrialVel();
                for (int d = 0; d < 3; d++) gdot += B[d] * vs(d);
                for (int i = 0; i < nps; i++) {
                    const Vector &vm = segNode[i]->getTrialVel();
                    for (int d = 0; d < 3; d++) gdot += B[3 * (1 + i) + d] * vm(d);
                }
                for (int k = 0; k < ndof; k++) resid(k) += -muc * gdot * B[k];
            }

            // --- Coulomb/Tresca friction (the tangent is in addKtToTang) ---
            // A frictionless pair returns before any state is touched.
            if (hasFriction() && cd != 0) {
                LadrunoContactDomain::FrictionState &st =
                    cd->getOrCreateFrictionState(contactTag, theSlave->getTag(), segIndex);
                // capture the tangential origin once, at the first activation, so
                // the tangential drift of a slave before contact does not become a
                // stick traction
                if (!st.engaged) {
                    for (int d = 0; d < 3; d++) st.gT0[d] = gTvec[d];
                    st.engaged = true;
                }
                double gTeff[3];
                for (int d = 0; d < 3; d++) gTeff[d] = gTvec[d] - st.gT0[d];
                double tFric[3], gpTtrial[3];
                // the cone pressure is the current normal force tn
                LadrunoFrictionKernel::frictionReturnMap(gTeff, st.gpT, tn, kt, mu,
                                                         tFric, gpTtrial, cohesion, tauMax);
                // the trial is a function of the committed state only; commit() promotes it
                for (int d = 0; d < 3; d++) st.gpTtrial[d] = gpTtrial[d];
                // as the normal block: slave += tFric, master_i += -N_i tFric
                // (tFric is the force on the slave, opposing its motion)
                for (int d = 0; d < 3; d++) resid(d) += tFric[d];
                for (int i = 0; i < nps; i++)
                    for (int d = 0; d < 3; d++)
                        resid(3 * (1 + i) + d) += -N[i] * tFric[d];
            }
        }
    } else if (mode == MORTAR && ndm == 2) {
        // 2D mortar: the 3D lane on 2-node segments, p_I = min(0, lambda_I +
        // epsN gbar_I) along the master normal, or the tie traction
        double D[2][2], M[2][2], g[2], n[2];
        int refusal = LadrunoContact2DKernel::MORT2D_REF_NONE; double cosT = 0.0;
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (!mortarActive2D(D, M, g, n, refusal, cosT)) {
            // A nearly perpendicular pair is not integrated (the Jacobian
            // 1/|cos| is ill-conditioned). It is reported once, unless a
            // node-to-segment contact on the same surfaces covers it (never
            // for a tie, which a unilateral contact cannot replace).
            if (refusal == LadrunoContact2DKernel::MORT2D_REF_ALIGN && !mortar2dPerpWarned &&
                (isTie || !mortar2dNtsCoDeclared) && cd != 0) {
                mortar2dPerpWarned = true;
                if (cd->warnOnce(contactTag, LadrunoContactDomain::WARN_2D_PERP_NO_NTS)) {
                    opserr << "WARNING LadrunoContactFE - contact " << contactTag
                           << ": 2D mortar slave segment " << slaveFacetIndex
                           << " is nearly perpendicular to its master (|cos| = " << cosT
                           << ") and is not integrated; ";
                    if (isTie)
                        opserr << "a tie there needs a finer mesh or a split surface\n";
                    else
                        opserr << "declare a node-to-segment contact on the same surfaces "
                                  "to cover it\n";
                }
            }
            // withdraw this pair's part of the summed gap / relative displacement
            if (cd != 0) {
                double z[3] = { 0.0, 0.0, 0.0 };
                for (int I = 0; I < 2; I++) {
                    if (isTie)
                        cd->accumulateMortarTie(contactTag, mortarSlave[I]->getTag(),
                                                this->getTag(), z, 0.0, kn);
                    else
                        cd->accumulateMortarGap(contactTag, mortarSlave[I]->getTag(),
                                                this->getTag(), 0.0, 0.0, kn);
                }
            }
            return resid;
        }
        if (isTie) {
            addMortarTie2D(D, M, cd);
            return resid;
        }
        double p[2] = { 0.0, 0.0 };
        for (int I = 0; I < 2; I++) {
            double aFacet = D[I][0] + D[I][1];
            int nodeTag = mortarSlave[I]->getTag();
            if (aFacet <= MORT2D_TAU_TRIB * mortar2dLref) {
                if (cd != 0)
                    cd->accumulateMortarGap(contactTag, nodeTag, this->getTag(), 0.0, 0.0, kn);
                continue;
            }
            double lambdaI = 0.0;
            if (cd != 0) {
                cd->accumulateMortarGap(contactTag, nodeTag, this->getTag(), g[I], aFacet, kn);
                lambdaI = cd->getOrCreateMortarNormalState(contactTag, nodeTag).lambdaN;
            }
            double pr = lambdaI + kn * (g[I] / aFacet);
            p[I] = (pr < 0.0) ? pr : 0.0;
        }
        for (int K = 0; K < 2; K++) {
            double Dp = D[K][0]*p[0] + D[K][1]*p[1];
            resid(2*K + 0) = -Dp * n[0]; resid(2*K + 1) = -Dp * n[1];
        }
        for (int L = 0; L < 2; L++) {
            double Mp = M[0][L]*p[0] + M[1][L]*p[1];
            resid(2*(2 + L) + 0) = Mp * n[0]; resid(2*(2 + L) + 1) = Mp * n[1];
        }
        if (hasFriction() && cd != 0)
            addMortarFriction2D(D, M, n, p, cd);
    } else if (mode == MORTAR) {
        // Augmented Lagrangian normal contact. The nodal pressure of this facet
        // pair is p_I = min(0, lambda_I + epsN gbar_I), gbar_I = gt_I / a_I the
        // weighted gap of this pair. lambda_I is the multiplier of the slave
        // node (committed, on the engine), and gbar_I is the gap of this pair,
        // not the running sum over the pairs of the node: the assembly forms the
        // residual pair by pair, and a running sum would give each pair of a
        // shared node a different pressure, an order-dependent residual. The
        // multiplier term assembles to the global D lambda over the pairs, and
        // the penalty term vanishes as the augmentation converges. The summed
        // gap is reported to the engine for the commit-time Uzawa update and
        // the penetration query only.
        // Force along n, self-equilibrating (the shape functions sum to one):
        // F^s_K = -(D p)_K n, F^m_L = +(M^T p)_L n.
        double D[4][4], M[4][4], g[4], n[3];
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (isTie) {
            if (mortarActive(D, M, g, n)) {
                addMortarTieForce(D, M, cd);
            } else if (cd != 0) {
                // empty overlap: withdraw this pair's part of the summed r
                double z[3] = {0.0, 0.0, 0.0};
                for (int I = 0; I < npsS; I++)
                    if (mortarSlave[I] != 0)
                        cd->accumulateMortarTie(contactTag, mortarSlave[I]->getTag(),
                                                this->getTag(), z, 0.0, kn);
            }
            return resid;
        }
        if (mortarActive(D, M, g, n)) {
            double p[4] = {0.0, 0.0, 0.0, 0.0};
            double chi[4] = {0.0, 0.0, 0.0, 0.0};        // friction onset of the smoothed law
            for (int I = 0; I < npsS; I++) {
                double aFacet = 0.0;                      // a_I = sum_J D_IJ = int N_I over the overlap
                for (int J = 0; J < npsS; J++) aFacet += D[I][J];
                int nodeTag = mortarSlave[I]->getTag();
                if (aFacet <= 1e-300) {                   // node not covered by this overlap
                    if (cd != 0)
                        cd->accumulateMortarGap(contactTag, nodeTag, this->getTag(), 0.0, 0.0, kn);
                    continue;
                }
                double lambdaI = 0.0;
                if (cd != 0) {
                    cd->accumulateMortarGap(contactTag, nodeTag, this->getTag(), g[I], aFacet, kn);
                    lambdaI = cd->getOrCreateMortarNormalState(contactTag, nodeTag).lambdaN;
                }
                double pr = lambdaI + kn * (g[I] / aFacet);
                if (smoothBand > 0.0) {                   // C1 quadratic onset
                    double P, dP, dchi;
                    mortarSmoothLaw(pr, P, dP, chi[I], dchi);
                    p[I] = -P;
                } else
                    p[I] = (pr < 0.0) ? pr : 0.0;         // active in compression
            }
            for (int K = 0; K < npsS; K++) {              // slave block: -(D p)_K n
                double Dp = 0.0;
                for (int I = 0; I < npsS; I++) Dp += D[K][I] * p[I];
                for (int d = 0; d < 3; d++) resid(3 * K + d) = -Dp * n[d];
            }
            for (int L = 0; L < npsM; L++) {              // master block: +(M^T p)_L n
                double Mp = 0.0;
                for (int I = 0; I < npsS; I++) Mp += M[I][L] * p[I];
                for (int d = 0; d < 3; d++) resid(3 * (npsS + L) + d) = Mp * n[d];
            }
            // Coulomb/Tresca friction; a frictionless pair touches no state
            if (hasFriction() && cd != 0)
                addMortarFriction(D, M, n, p, cd, (smoothBand > 0.0) ? chi : 0);
        } else if (cd != 0) {
            // empty overlap at this evaluation: withdraw this pair's part of the
            // summed gap of its slave nodes
            for (int I = 0; I < npsS; I++)
                if (mortarSlave[I] != 0)
                    cd->accumulateMortarGap(contactTag, mortarSlave[I]->getTag(), this->getTag(),
                                            0.0, 0.0, kn);
        }
    } else if (mode == EDGE_EDGE) {
        // Edge-to-edge penalty: f = tN B, tN = -p, p = min(0, lambdaN + epsN gN)
        // (lambdaN = 0 without -edgeAlm, the penalty tN = epsN <-gN>+), with the
        // body-fixed sign of the pair. The master edge takes the opposite
        // reaction (the rows of B sum to zero).
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        LadrunoContactDomain::EdgeEdgeState *st = 0;
        int committedSign = 0;
        if (cd != 0) {
            st = &cd->getOrCreateEdgeEdgeState(contactTag, edgeNode[0]->getTag(),
                                               edgeNode[1]->getTag(), edgeNode[2]->getTag(),
                                               edgeNode[3]->getTag());
            committedSign = st->signN;
        }
        double gN, n[3], s, t, B[4][3]; int usedSign = 0;
        if (edgeGeom(gN, n, s, t, B, committedSign, &usedSign)) {
            // the sign is captured once, at the first engagement, and then kept
            if (st != 0 && st->signN == 0) st->signN = usedSign;
            double lambdaN = (edgeAlm && st != 0) ? st->lambdaN : 0.0;
            double pAug = lambdaN + kn * gN;
            if (edgeAlm && st != 0) { st->gN = gN; st->epsN = kn; }
            if (pAug < 0.0) {
                double tN = -pAug;
                for (int i = 0; i < 4; i++)
                    for (int d = 0; d < 3; d++)
                        resid(3 * i + d) = tN * B[i][d];

                // Coulomb/Tresca friction on the slip of the closest points; a
                // frictionless pair touches no state
                if (st != 0 && hasFriction()) {
                    double gTvec[3]; edgeSlip(s, t, n, gTvec);
                    // the origin is captured at the first engagement, so a
                    // tangential drift before contact does not become a stick traction
                    if (!st->engaged) {
                        for (int d = 0; d < 3; d++) st->gT0[d] = gTvec[d];
                        st->engaged = true;
                    }
                    double gTeff[3];
                    for (int d = 0; d < 3; d++) gTeff[d] = gTvec[d] - st->gT0[d];
                    double tFric[3], gpTtrial[3];
                    LadrunoFrictionKernel::frictionReturnMap(gTeff, st->gpT, tN, kt, mu, tFric,
                                                             gpTtrial, cohesion, tauMax);
                    for (int d = 0; d < 3; d++) st->gpTtrial[d] = gpTtrial[d];
                    // scattered with the weights of the normal force
                    double w[4] = { 1.0 - s, s, -(1.0 - t), -t };
                    for (int i = 0; i < 4; i++)
                        for (int d = 0; d < 3; d++)
                            resid(3 * i + d) += w[i] * tFric[d];
                }
            }
        } else if (edgeAlm && st != 0) {
            // no longer a crossing (the pair slid past an end of an edge or
            // turned parallel): drop the gap so that the next commit does not
            // augment from it and the penetration query ignores the pair;
            // lambdaN is kept until the pair crosses again
            st->gN = 0.0;
            st->epsN = 0.0;
        }
    }
    return resid;
}

void
LadrunoContactFE::addMortarFriction(const double D[4][4], const double M[4][4], const double n[3],
                                    const double p_normal[4], LadrunoContactDomain *cd,
                                    const double *chi)
{
    double us[4][3], um[4][3];
    for (int i = 0; i < npsS; i++) {
        const Vector &u = mortarSlave[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) us[i][d] = u(d);
    }
    for (int i = 0; i < npsM; i++) {
        const Vector &u = mortarMaster[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) um[i][d] = u(d);
    }
    double tFric[4][3] = {{0.0}};
    double invEpsT = (kt > 0.0) ? 1.0 / kt : 0.0;
    for (int I = 0; I < npsS; I++) {
        if (p_normal[I] >= 0.0) continue;             // friction only in normal contact
        double aFacet = 0.0;
        for (int J = 0; J < npsS; J++) aFacet += D[I][J];
        if (aFacet <= 1e-300) continue;
        double N_I = -p_normal[I];
        double r[3] = {0.0, 0.0, 0.0};
        for (int J = 0; J < npsS; J++)
            for (int d = 0; d < 3; d++) r[d] += D[I][J] * us[J][d];
        for (int K = 0; K < npsM; K++)
            for (int d = 0; d < 3; d++) r[d] -= M[I][K] * um[K][d];
        double rn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
        double gbarT[3];
        for (int d = 0; d < 3; d++) gbarT[d] = (r[d] - rn * n[d]) / aFacet;
        LadrunoContactDomain::MortarFrictionState &st = cd->getOrCreateMortarFrictionState(
            contactTag, mortarSlave[I]->getTag(), slaveFacetIndex, masterFacetIndex);
        // the origin is captured at the first engagement, so a tangential drift
        // before contact does not become a stick traction
        if (!st.engaged) {
            for (int d = 0; d < 3; d++) st.gT0[d] = gbarT[d];
            st.engaged = true;
        }
        double gTeff[3];
        for (int d = 0; d < 3; d++) gTeff[d] = (gbarT[d] - st.gT0[d]) + st.lambdaT[d] * invEpsT;
        double tF[3], gpTtrial[3];
        if (smoothSlip > 0.0)                         // rounded stick/slip corner
            LadrunoFrictionKernel::frictionReturnMapSmooth(gTeff, st.gpT, N_I, kt, mu, tF,
                                                           gpTtrial, cohesion, tauMax, smoothSlip);
        else
            LadrunoFrictionKernel::frictionReturnMap(gTeff, st.gpT, N_I, kt, mu, tF, gpTtrial,
                                                     cohesion, tauMax);
        if (chi != 0)                                 // friction onset over the band
            for (int d = 0; d < 3; d++) tF[d] *= chi[I];
        for (int d = 0; d < 3; d++) {
            st.gpTtrial[d] = gpTtrial[d];
            st.lambdaTtrial[d] = -tF[d];
            tFric[I][d] = tF[d];
        }
    }
    for (int K = 0; K < npsS; K++)
        for (int d = 0; d < 3; d++) {
            double s = 0.0;
            for (int I = 0; I < npsS; I++) s += D[K][I] * tFric[I][d];
            resid(3 * K + d) += s;
        }
    for (int L = 0; L < npsM; L++)
        for (int d = 0; d < 3; d++) {
            double s = 0.0;
            for (int I = 0; I < npsS; I++) s += M[I][L] * tFric[I][d];
            resid(3 * (npsS + L) + d) += -s;
        }
}

void
LadrunoContactFE::addMortarTieForce(const double D[4][4], const double M[4][4],
                                    LadrunoContactDomain *cd)
{
    double us[4][3], um[4][3];
    for (int i = 0; i < npsS; i++) {
        const Vector &u = mortarSlave[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) us[i][d] = u(d);
    }
    for (int i = 0; i < npsM; i++) {
        const Vector &u = mortarMaster[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) um[i][d] = u(d);
    }
    double t[4][3] = {{0.0}};
    for (int I = 0; I < npsS; I++) {
        double aFacet = 0.0;
        for (int J = 0; J < npsS; J++) aFacet += D[I][J];
        int nodeTag = mortarSlave[I]->getTag();
        double r[3] = {0.0, 0.0, 0.0};
        if (aFacet <= 1e-300) {                     // node not covered by this overlap
            if (cd != 0)
                cd->accumulateMortarTie(contactTag, nodeTag, this->getTag(), r, 0.0, kn);
            continue;
        }
        for (int J = 0; J < npsS; J++)
            for (int d = 0; d < 3; d++) r[d] += D[I][J] * us[J][d];
        for (int K = 0; K < npsM; K++)
            for (int d = 0; d < 3; d++) r[d] -= M[I][K] * um[K][d];
        double lamTie[3] = {0.0, 0.0, 0.0};
        if (cd != 0) {
            // the summed r feeds the commit-time update and the residual query
            // only; this pair's force uses its own r (order-independent residual)
            cd->accumulateMortarTie(contactTag, nodeTag, this->getTag(), r, aFacet, kn);
            const LadrunoContactDomain::MortarNormalState &st =
                cd->getOrCreateMortarNormalState(contactTag, nodeTag);
            for (int d = 0; d < 3; d++) lamTie[d] = st.lambdaTie[d];
        }
        for (int d = 0; d < 3; d++) t[I][d] = lamTie[d] + kn * (r[d] / aFacet);
    }
    for (int K = 0; K < npsS; K++)
        for (int d = 0; d < 3; d++) {
            double s = 0.0;
            for (int I = 0; I < npsS; I++) s += D[K][I] * t[I][d];
            resid(3 * K + d) = -s;
        }
    for (int L = 0; L < npsM; L++)
        for (int d = 0; d < 3; d++) {
            double s = 0.0;
            for (int I = 0; I < npsS; I++) s += M[I][L] * t[I][d];
            resid(3 * (npsS + L) + d) = s;
        }
}

bool
LadrunoContactFE::edgeGeom(double &gN, double n[3], double &s, double &t, double B[4][3],
                           int committedSign, int *outSign) const
{
    if (mode != EDGE_EDGE) return false;
    double X[4][3];
    for (int i = 0; i < 4; i++) {
        if (edgeNode[i] == 0) return false;
        const Vector &Xi = edgeNode[i]->getCrds();
        const Vector &ui = edgeNode[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) X[i][d] = Xi(d) + ui(d);
    }
    LadrunoEdgeKernel::ClosestResult cr =
        LadrunoEdgeKernel::closestPtSegSeg(X[0], X[1], X[2], X[3]);
    // only a well-conditioned crossing strictly inside both edges is an
    // edge-to-edge contact, where |gN| is the distance of the closest points
    if (cr.status != LadrunoEdgeKernel::EE_OK || !cr.interior) return false;
    s = cr.s; t = cr.t;
    double d1[3], d2[3];
    for (int d = 0; d < 3; d++) { d1[d] = X[1][d] - X[0][d]; d2[d] = X[3][d] - X[2][d]; }
    // A first capture with orientDir perpendicular to n cannot decide the sign
    // (it would be fixed by round-off for the life of the pair): the capture is
    // deferred and the evaluation is treated as no contact. On a symmetric
    // configuration that can last, so it is reported once per contact.
    int chosen = 0;
    gN = LadrunoEdgeKernel::edgeGap(cr.c1, cr.c2, d1, d2, committedSign, orientDir, n, &chosen);
    if (chosen == 0) {
        LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
        if (cd != 0 && cd->warnOnce(contactTag, LadrunoContactDomain::WARN_EDGE_SIGN_DEFER))
            opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": the sign of an "
                      "edge-to-edge pair is not decided yet (the orientation reference is "
                      "perpendicular to the common normal of the edges); the pair carries no "
                      "force until the geometry decides it. Give -outward with a component "
                      "along the expected contact normal.\n";
        return false;
    }
    if (outSign != 0) *outSign = chosen;
    LadrunoEdgeKernel::bOperator(n, s, t, B);
    return true;
}

void
LadrunoContactFE::edgeSlip(double s, double t, const double n[3], double gT[3]) const
{
    const Vector &ua0 = edgeNode[0]->getTrialDisp();
    const Vector &ua1 = edgeNode[1]->getTrialDisp();
    const Vector &ub0 = edgeNode[2]->getTrialDisp();
    const Vector &ub1 = edgeNode[3]->getTrialDisp();
    double drel[3];
    for (int d = 0; d < 3; d++)
        drel[d] = ((1.0 - s) * ua0(d) + s * ua1(d)) - ((1.0 - t) * ub0(d) + t * ub1(d));
    LadrunoFrictionKernel::tangentPart(drel, n, gT);
}

void
LadrunoContactFE::addEdgeTang(double fact, bool initialStiff)
{
    // K = epsN B^T B on the active set of the residual (symmetric, rank one,
    // positive semidefinite) with the normal, s and t frozen. The geometric
    // part (the variation of n, s and t with the displacements) is not
    // included: it vanishes for a pair that does not rotate or slide, and
    // Newton converges linearly otherwise. The active set uses the augmented
    // pressure lambdaN + epsN gN of the residual (lambdaN is constant within
    // the iterations of a step).
    LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
    LadrunoContactDomain::EdgeEdgeState *st = 0;
    int committedSign = 0;
    if (cd != 0) {
        st = &cd->getOrCreateEdgeEdgeState(contactTag, edgeNode[0]->getTag(),
                                           edgeNode[1]->getTag(), edgeNode[2]->getTag(),
                                           edgeNode[3]->getTag());
        committedSign = st->signN;
    }
    double lambdaN = (edgeAlm && st != 0) ? st->lambdaN : 0.0;
    double gN, n[3], s, t, B[4][3]; int usedSign = 0;
    if (!edgeGeom(gN, n, s, t, B, committedSign, &usedSign) || lambdaN + kn * gN >= 0.0)
        return;
    for (int i = 0; i < 4; i++)
        for (int di = 0; di < 3; di++)
            for (int j = 0; j < 4; j++)
                for (int dj = 0; dj < 3; dj++)
                    tang(3 * i + di, 3 * j + dj) += fact * kn * B[i][di] * B[j][dj];
    if (!hasFriction() || (!initialStiff && st == 0))
        return;                       // no engine: the residual carries no friction either
    // friction block w_A w_B K_ss, w = [(1-s), s, -(1-t), -t], at the slip and
    // committed state of the residual; the initial-stiffness path uses the
    // stick block (gTeff = gpT = 0), positive semidefinite and engine-free
    double tN = -(lambdaN + kn * gN);
    double Kss[3][3];
    if (initialStiff) {
        double zero[3] = {0.0, 0.0, 0.0};
        LadrunoFrictionKernel::frictionTangentBlock(zero, zero, n, tN, kn, kt, mu, false, Kss,
                                                    cohesion, tauMax);
    } else {
        double gTvec[3]; edgeSlip(s, t, n, gTvec);
        double gTeff[3];
        for (int d = 0; d < 3; d++) gTeff[d] = st->engaged ? (gTvec[d] - st->gT0[d]) : 0.0;
        LadrunoFrictionKernel::frictionTangentBlock(gTeff, st->gpT, n, tN, kn, kt, mu,
                                                    consistentTan, Kss, cohesion, tauMax);
    }
    double w[4] = { 1.0 - s, s, -(1.0 - t), -t };
    for (int A = 0; A < 4; A++)
        for (int Bn = 0; Bn < 4; Bn++) {
            double wab = fact * w[A] * w[Bn];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    tang(3 * A + i, 3 * Bn + j) += wab * Kss[i][j];
        }
}

Node *
LadrunoContactFE::nts2dPrevFarNode(void) const
{
    if (!nts2dHasPrevFar || theDomain == 0) return 0;
    return theDomain->getNode(nts2dPrevFarTag);
}

Node *
LadrunoContactFE::nts2dNextFarNode(void) const
{
    if (!nts2dHasNextFar || theDomain == 0) return 0;
    return theDomain->getNode(nts2dNextFarTag);
}

bool
LadrunoContactFE::segment2DActive(double &gap, double n[3], double N[2], double B[6],
                                  double committedSide, double *liveSideOut,
                                  double *gTout) const
{
    using namespace LadrunoContact2DKernel;
    gap = 0.0; n[0] = n[1] = n[2] = 0.0; N[0] = N[1] = 0.0;
    for (int k = 0; k < 6; k++) B[k] = 0.0;
    if (liveSideOut != 0) *liveSideOut = 0.0;
    if (gTout != 0) *gTout = 0.0;
    if (mode != SEGMENT || ndm != 2 || theSlave == 0) return false;

    Node *prevFarNode = nts2dPrevFarNode();
    Node *nextFarNode = nts2dNextFarNode();

    const Vector &Xs = theSlave->getCrds();
    const Vector &us = theSlave->getTrialDisp();
    const double xs[2] = { Xs(0) + us(0), Xs(1) + us(1) };

    if (nps == 2) {
        if (segNode[0] == 0 || segNode[1] == 0) return false;
        double X0[2], X1[2];
        {
            const Vector &Xa = segNode[0]->getCrds(); const Vector &ua = segNode[0]->getTrialDisp();
            const Vector &Xb = segNode[1]->getCrds(); const Vector &ub = segNode[1]->getTrialDisp();
            X0[0] = Xa(0) + ua(0); X0[1] = Xa(1) + ua(1);
            X1[0] = Xb(0) + ub(0); X1[1] = Xb(1) + ub(1);
        }
        double xi, g2, n2[2];
        int st = projectSegment2D(X0, X1, xs, nts2dSigma, nts2dLref, xi, g2, n2);
        // past an end: the neighbour, a vertex pair or an end cap owns the slave
        if (st != SEG2D_IN_BOUNDS)
            return false;
        // ordered ownership: the previous segment owns a slave both project on
        if (prevFarNode != 0) {
            const Vector &Xp = prevFarNode->getCrds();
            const Vector &up = prevFarNode->getTrialDisp();
            const double XP[2] = { Xp(0) + up(0), Xp(1) + up(1) };
            double xiP, gP, nP[2];
            if (projectSegment2D(XP, X0, xs, nts2dSigma, nts2dLref, xiP, gP, nP)
                    == SEG2D_IN_BOUNDS)
                return false;
        }
        if (g2 >= 0.0) return false;
        gap = g2; n[0] = n2[0]; n[1] = n2[1];
        shape2D(xi, N);
        bOperatorSegment2D(N[0], N[1], n2, B);
        // scalar slip: (u_s - sum N_i u_i).t, t = perp(n), at the projection
        // of the gap
        if (gTout != 0) {
            const Vector &u0 = segNode[0]->getTrialDisp();
            const Vector &u1 = segNode[1]->getTrialDisp();
            const double ubar[2] = { N[0]*u0(0) + N[1]*u1(0), N[0]*u0(1) + N[1]*u1(1) };
            const double drel[2] = { us(0) - ubar[0], us(1) - ubar[1] };
            const double th[2] = { -n2[1], n2[0] };
            *gTout = drel[0]*th[0] + drel[1]*th[1];
        }
        return true;
    }

    // vertex pair (nps == 1): a concave corner (both neighbours) or an end cap
    // (one neighbour)
    if (nps != 1 || segNode[0] == 0) return false;
    if (prevFarNode == 0 && nextFarNode == 0) return false;

    double XV[2];
    {
        const Vector &Xv = segNode[0]->getCrds();
        const Vector &uv = segNode[0]->getTrialDisp();
        XV[0] = Xv(0) + uv(0); XV[1] = Xv(1) + uv(1);
    }
    double side = 0.0;
    LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;

    if (prevFarNode != 0 && nextFarNode != 0) {
        double XP[2], XN[2];
        {
            const Vector &Xp = prevFarNode->getCrds(); const Vector &up = prevFarNode->getTrialDisp();
            const Vector &Xn = nextFarNode->getCrds(); const Vector &un = nextFarNode->getTrialDisp();
            XP[0] = Xp(0) + up(0); XP[1] = Xp(1) + up(1);
            XN[0] = Xn(0) + un(0); XN[1] = Xn(1) + un(1);
        }
        // the adjacent segments come first
        {
            double xiT, gT, nT[2];
            if (projectSegment2D(XP, XV, xs, nts2dSigma, nts2dLref, xiT, gT, nT) == SEG2D_IN_BOUNDS)
                return false;
            if (projectSegment2D(XV, XN, xs, nts2dSigma, nts2dLref, xiT, gT, nT) == SEG2D_IN_BOUNDS)
                return false;
        }
        const double tPrev[2] = { XV[0] - XP[0], XV[1] - XP[1] };
        const double tNext[2] = { XN[0] - XV[0], XN[1] - XV[1] };
        WedgeResult w = vertexWedge2D(tPrev, tNext, nts2dSigma, xs, XV, nts2dLref);
        // the wedge claim without slack closes the band between the two
        // in-bounds tests
        if (!(w.aPrev > 0.0 && w.aNext < 0.0))
            return false;
        if (liveSideOut != 0) *liveSideOut = w.sideSign;
        side = (committedSide != 0.0) ? committedSide : w.sideSign;
        if (side == 0.0) {
            if (cd != 0 && cd->warnOnce(contactTag, LadrunoContactDomain::WARN_VTX2D_DEFER))
                opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": the side "
                          "sign of a 2D vertex pair is not decided yet (a fold-back corner, or "
                          "a slave on the boundary of the corner wedge); the pair carries no "
                          "force until the geometry decides it\n";
            return false;
        }
    } else {
        // End cap at an open end of the master chain: the radial pair of the
        // terminal vertex, active past the end along the terminal segment and
        // within a reach of 5 percent of Lref, so it does not claim slaves far
        // along the extended line. The side is the side of the extended line,
        // continuous with the segment at the end.
        bool isX0Cap = (prevFarNode == 0);
        Node *far = isX0Cap ? nextFarNode : prevFarNode;
        const Vector &Xf = far->getCrds();
        const Vector &uf = far->getTrialDisp();
        const double XF[2] = { Xf(0) + uf(0), Xf(1) + uf(1) };
        double X0[2], X1[2];
        if (isX0Cap) { X0[0] = XV[0]; X0[1] = XV[1]; X1[0] = XF[0]; X1[1] = XF[1]; }
        else         { X0[0] = XF[0]; X0[1] = XF[1]; X1[0] = XV[0]; X1[1] = XV[1]; }
        double xiChk, gChk, nEnd[2];
        int stChk = projectSegment2D(X0, X1, xs, nts2dSigma, nts2dLref, xiChk, gChk, nEnd);
        if (stChk == SEG2D_DEGENERATE) return false;
        bool onOpenSide = isX0Cap ? (stChk == SEG2D_OUT_LOW) : (stChk == SEG2D_OUT_HIGH);
        if (!onOpenSide) return false;
        static const double ENDCAP_REACH = 5.0e-2;
        const double tSeg[2] = { X1[0] - X0[0], X1[1] - X0[1] };
        const double Lseg = std::sqrt(tSeg[0]*tSeg[0] + tSeg[1]*tSeg[1]);
        double xiPast = isX0Cap ? -xiChk : (xiChk - 1.0);
        if (xiPast * Lseg > ENDCAP_REACH * nts2dLref) return false;
        double liveSide;
        if (std::fabs(gChk) < TAU_PERP_DEFAULT * nts2dLref) liveSide = 0.0;
        else liveSide = (gChk > 0.0) ? 1.0 : -1.0;
        if (liveSideOut != 0) *liveSideOut = liveSide;
        side = (committedSide != 0.0) ? committedSide : liveSide;
        if (side == 0.0) {
            if (cd != 0 && cd->warnOnce(contactTag, LadrunoContactDomain::WARN_VTX2D_DEFER))
                opserr << "WARNING LadrunoContactFE - contact " << contactTag << ": the side "
                          "sign of a 2D end-cap pair is not decided yet (the slave is on the "
                          "extended line of the terminal segment); the pair carries no force "
                          "until the geometry decides it\n";
            return false;
        }
    }

    // the slave on the vertex has no direction; the floor is the conditioning
    // gauge tauPerp Lref
    double g2, n2[2];
    if (vertexEval2D(XV, xs, side, g2, n2, TAU_PERP_DEFAULT * nts2dLref) != VTX2D_OK)
        return false;
    if (g2 >= 0.0) return false;
    gap = g2; n[0] = n2[0]; n[1] = n2[1];
    N[0] = 1.0; N[1] = 0.0;
    bOperatorVertex2D(n2, B);
    if (gTout != 0) {
        const Vector &uv = segNode[0]->getTrialDisp();
        const double drel[2] = { us(0) - uv(0), us(1) - uv(1) };
        const double th[2] = { -n2[1], n2[0] };
        *gTout = drel[0]*th[0] + drel[1]*th[1];
    }
    return true;
}

double
LadrunoContactFE::vtx2DCommittedSide(void) const
{
    if (mode != SEGMENT || ndm != 2 || nps != 1) return 0.0;
    if (theDomain == 0 || theSlave == 0 || segNode[0] == 0) return 0.0;
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
    if (cd == 0) return 0.0;
    return cd->getVtx2DSide(contactTag, theSlave->getTag(), segNode[0]->getTag());
}

void
LadrunoContactFE::addFrictionTang2D(double fact, const double n[3], const double N2[2],
                                    double tn, double gTeff, double gpT, bool consistent)
{
    // the 2 x 2 block from the two scalars of the kernel by the chain rule:
    // d(gTeff)/du_rel = t and d(tn)/du_rel = kn n; dTN already carries kn
    const double th[2] = { -n[1], n[0] };
    double Kss_s = 0.0, dTN_s = 0.0;
    LadrunoFrictionKernel::tangentBlock1D(gTeff, gpT, tn, kn, kt, mu, consistent, Kss_s, dTN_s,
                                          cohesion, tauMax);
    double Kss2[2][2];
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            Kss2[i][j] = Kss_s * th[i] * th[j] + dTN_s * th[i] * n[j];
    double w[3];
    w[0] = 1.0;
    for (int i = 0; i < nps; i++) w[1 + i] = -N2[i];
    int nn = 1 + nps;
    for (int a = 0; a < nn; a++)
        for (int b = 0; b < nn; b++) {
            double wab = fact * w[a] * w[b];
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    tang(2 * a + i, 2 * b + j) += wab * Kss2[i][j];
        }
}

bool
LadrunoContactFE::mortarActive2D(double D[2][2], double M[2][2], double g[2], double n[2],
                                 int &refusalOut, double &cosOut) const
{
    using namespace LadrunoContact2DKernel;
    double Xs[2][2], Xm[2][2];
    for (int i = 0; i < 2; i++) {
        const Vector &X = mortarSlave[i]->getCrds();
        const Vector &u = mortarSlave[i]->getTrialDisp();
        Xs[i][0] = X(0) + u(0); Xs[i][1] = X(1) + u(1);
    }
    for (int i = 0; i < 2; i++) {
        const Vector &X = mortarMaster[i]->getCrds();
        const Vector &u = mortarMaster[i]->getTrialDisp();
        Xm[i][0] = X(0) + u(0); Xm[i][1] = X(1) + u(1);
    }
    PairResult2D pr;
    mortarSegmentPair2D(Xs, Xm, mortar2dSigma, mortar2dLref, pr);
    refusalOut = pr.refusal;
    cosOut = pr.cos_t;
    if (pr.status != MORT2D_OK)
        return false;
    for (int i = 0; i < 2; i++) {
        g[i] = pr.g[i];
        for (int j = 0; j < 2; j++) { D[i][j] = pr.D[i][j]; M[i][j] = pr.M[i][j]; }
    }
    // the current master normal, the one the weighted gap was measured along
    const double tm[2] = { Xm[1][0] - Xm[0][0], Xm[1][1] - Xm[0][1] };
    const double Lm = std::sqrt(tm[0]*tm[0] + tm[1]*tm[1]);
    n[0] = -mortar2dSigma * tm[1] / Lm;
    n[1] =  mortar2dSigma * tm[0] / Lm;
    return true;
}

void
LadrunoContactFE::addMortarFriction2D(const double D[2][2], const double M[2][2],
                                      const double n[2], const double p_normal[2],
                                      LadrunoContactDomain *cd)
{
    double us[2][2], um[2][2];
    for (int i = 0; i < 2; i++) {
        const Vector &u = mortarSlave[i]->getTrialDisp();
        us[i][0] = u(0); us[i][1] = u(1);
    }
    for (int i = 0; i < 2; i++) {
        const Vector &u = mortarMaster[i]->getTrialDisp();
        um[i][0] = u(0); um[i][1] = u(1);
    }
    const double th[2] = { -n[1], n[0] };
    double tFric[2] = { 0.0, 0.0 };
    double invEpsT = (kt > 0.0) ? 1.0 / kt : 0.0;
    for (int I = 0; I < 2; I++) {
        if (p_normal[I] >= 0.0) continue;
        double aFacet = D[I][0] + D[I][1];
        if (aFacet <= MORT2D_TAU_TRIB * mortar2dLref) continue;
        double N_I = -p_normal[I];
        // scalar slip: the tangential part of the weighted relative
        // displacement over a_I, from the engagement origin, offset by the
        // tangential multiplier (slot 0 of the friction state)
        double r[2] = { 0.0, 0.0 };
        for (int J = 0; J < 2; J++) { r[0] += D[I][J]*us[J][0]; r[1] += D[I][J]*us[J][1]; }
        for (int K = 0; K < 2; K++) { r[0] -= M[I][K]*um[K][0]; r[1] -= M[I][K]*um[K][1]; }
        double gT = (r[0]*th[0] + r[1]*th[1]) / aFacet;
        LadrunoContactDomain::MortarFrictionState &st = cd->getOrCreateMortarFrictionState(
            contactTag, mortarSlave[I]->getTag(), slaveFacetIndex, masterFacetIndex);
        if (!st.engaged) { st.gT0[0] = gT; st.engaged = true; }
        double gTeff = (gT - st.gT0[0]) + st.lambdaT[0] * invEpsT;
        double tF_s, spTrial;
        LadrunoFrictionKernel::returnMap1D(gTeff, st.gpT[0], N_I, kt, mu, tF_s, spTrial,
                                           cohesion, tauMax);
        st.gpTtrial[0] = spTrial;
        st.lambdaTtrial[0] = -tF_s;
        tFric[I] = tF_s;
    }
    for (int K = 0; K < 2; K++) {
        double s0 = 0.0, s1 = 0.0;
        for (int I = 0; I < 2; I++) { s0 += D[K][I]*tFric[I]*th[0]; s1 += D[K][I]*tFric[I]*th[1]; }
        resid(2*K + 0) += s0; resid(2*K + 1) += s1;
    }
    for (int L = 0; L < 2; L++) {
        double s0 = 0.0, s1 = 0.0;
        for (int I = 0; I < 2; I++) { s0 += M[I][L]*tFric[I]*th[0]; s1 += M[I][L]*tFric[I]*th[1]; }
        resid(2*(2 + L) + 0) += -s0; resid(2*(2 + L) + 1) += -s1;
    }
}

void
LadrunoContactFE::addMortarTie2D(const double D[2][2], const double M[2][2],
                                 LadrunoContactDomain *cd)
{
    double us[2][2], um[2][2];
    for (int i = 0; i < 2; i++) {
        const Vector &u = mortarSlave[i]->getTrialDisp();
        us[i][0] = u(0); us[i][1] = u(1);
    }
    for (int i = 0; i < 2; i++) {
        const Vector &u = mortarMaster[i]->getTrialDisp();
        um[i][0] = u(0); um[i][1] = u(1);
    }
    double t[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
    for (int I = 0; I < 2; I++) {
        double aFacet = D[I][0] + D[I][1];
        if (aFacet <= MORT2D_TAU_TRIB * mortar2dLref) {
            if (cd != 0) {
                double z[3] = { 0.0, 0.0, 0.0 };
                cd->accumulateMortarTie(contactTag, mortarSlave[I]->getTag(),
                                        this->getTag(), z, 0.0, kn);
            }
            continue;
        }
        double r[2] = { 0.0, 0.0 };
        for (int J = 0; J < 2; J++) { r[0] += D[I][J]*us[J][0]; r[1] += D[I][J]*us[J][1]; }
        for (int K = 0; K < 2; K++) { r[0] -= M[I][K]*um[K][0]; r[1] -= M[I][K]*um[K][1]; }
        double lamTie[2] = { 0.0, 0.0 };
        if (cd != 0) {
            // the three-component accumulators of the 3D tie, third component 0
            int nodeTag = mortarSlave[I]->getTag();
            double rFacet[3] = { r[0], r[1], 0.0 };
            cd->accumulateMortarTie(contactTag, nodeTag, this->getTag(), rFacet, aFacet, kn);
            const LadrunoContactDomain::MortarNormalState &st =
                cd->getOrCreateMortarNormalState(contactTag, nodeTag);
            lamTie[0] = st.lambdaTie[0]; lamTie[1] = st.lambdaTie[1];
        }
        t[I][0] = lamTie[0] + kn * (r[0] / aFacet);
        t[I][1] = lamTie[1] + kn * (r[1] / aFacet);
    }
    for (int K = 0; K < 2; K++) {
        double s0 = 0.0, s1 = 0.0;
        for (int I = 0; I < 2; I++) { s0 += D[K][I]*t[I][0]; s1 += D[K][I]*t[I][1]; }
        resid(2*K + 0) = -s0; resid(2*K + 1) = -s1;
    }
    for (int L = 0; L < 2; L++) {
        double s0 = 0.0, s1 = 0.0;
        for (int I = 0; I < 2; I++) { s0 += M[I][L]*t[I][0]; s1 += M[I][L]*t[I][1]; }
        resid(2*(2 + L) + 0) = s0; resid(2*(2 + L) + 1) = s1;
    }
}

void
LadrunoContactFE::addMortarTang2D(double fact, bool initialStiff)
{
    double D[2][2], M[2][2], g[2], n[2];
    int refusal; double cosT;
    if (!mortarActive2D(D, M, g, n, refusal, cosT)) return;
    const int nN = 4;
    LadrunoContactDomain *cd = (theDomain != 0) ? theDomain->getLadrunoContactDomain() : 0;
    if (isTie) {
        // epsTie Bt^T diag(1/a) Bt (x) I2
        double W[2] = { 0.0, 0.0 };
        for (int I = 0; I < 2; I++) {
            double aFacet = D[I][0] + D[I][1];
            if (aFacet > MORT2D_TAU_TRIB * mortar2dLref) W[I] = 1.0 / aFacet;
        }
        for (int A = 0; A < nN; A++)
            for (int B = 0; B < nN; B++) {
                double Ks = 0.0;
                for (int I = 0; I < 2; I++) {
                    double bIA = (A < 2) ? D[I][A] : -M[I][A - 2];
                    double bIB = (B < 2) ? D[I][B] : -M[I][B - 2];
                    Ks += bIA * W[I] * bIB;
                }
                Ks *= kn;
                if (Ks == 0.0) continue;
                tang(2*A + 0, 2*B + 0) += fact * Ks;
                tang(2*A + 1, 2*B + 1) += fact * Ks;
            }
        return;
    }
    // epsN Bt^T diag(act/a) Bt (x) (n (x) n) on the active set of the residual
    double W[2] = { 0.0, 0.0 };
    for (int I = 0; I < 2; I++) {
        double aFacet = D[I][0] + D[I][1];
        if (aFacet <= MORT2D_TAU_TRIB * mortar2dLref) continue;
        double lambdaI = (cd != 0)
            ? cd->getOrCreateMortarNormalState(contactTag, mortarSlave[I]->getTag()).lambdaN
            : 0.0;
        double pr = lambdaI + kn * (g[I] / aFacet);
        if (pr < 0.0) W[I] = 1.0 / aFacet;
    }
    for (int A = 0; A < nN; A++) {
        for (int B = 0; B < nN; B++) {
            double Ks = 0.0;
            for (int I = 0; I < 2; I++) {
                double bIA = (A < 2) ? D[I][A] : -M[I][A - 2];
                double bIB = (B < 2) ? D[I][B] : -M[I][B - 2];
                Ks += bIA * W[I] * bIB;
            }
            Ks *= kn;
            if (Ks == 0.0) continue;
            for (int dA = 0; dA < 2; dA++)
                for (int dB = 0; dB < 2; dB++)
                    tang(2*A + dA, 2*B + dB) += fact * Ks * n[dA] * n[dB];
        }
    }
    // friction: sum_I (b_IA b_IB / a_I) K2_I, K2 = Kss t (x) t + dTN t (x) n
    // from the scalar kernel at the slip of the residual (stick block on the
    // initial-stiffness path)
    if (!hasFriction() || cd == 0)
        return;
    double us[2][2], um[2][2];
    for (int i = 0; i < 2; i++) {
        const Vector &u = mortarSlave[i]->getTrialDisp();
        us[i][0] = u(0); us[i][1] = u(1);
    }
    for (int i = 0; i < 2; i++) {
        const Vector &u = mortarMaster[i]->getTrialDisp();
        um[i][0] = u(0); um[i][1] = u(1);
    }
    const double th[2] = { -n[1], n[0] };
    double invEpsT = (kt > 0.0) ? 1.0 / kt : 0.0;
    for (int I = 0; I < 2; I++) {
        double aFacet = D[I][0] + D[I][1];
        if (aFacet <= MORT2D_TAU_TRIB * mortar2dLref) continue;
        const LadrunoContactDomain::MortarNormalState &nst =
            cd->getOrCreateMortarNormalState(contactTag, mortarSlave[I]->getTag());
        double pr = nst.lambdaN + kn * (g[I] / aFacet);
        if (pr >= 0.0) continue;
        double N_I = -pr;
        LadrunoContactDomain::MortarFrictionState &st = cd->getOrCreateMortarFrictionState(
            contactTag, mortarSlave[I]->getTag(), slaveFacetIndex, masterFacetIndex);
        double r[2] = { 0.0, 0.0 };
        for (int J = 0; J < 2; J++) { r[0] += D[I][J]*us[J][0]; r[1] += D[I][J]*us[J][1]; }
        for (int K = 0; K < 2; K++) { r[0] -= M[I][K]*um[K][0]; r[1] -= M[I][K]*um[K][1]; }
        double gT = (r[0]*th[0] + r[1]*th[1]) / aFacet;
        double gTeff = (gT - st.gT0[0]) + st.lambdaT[0] * invEpsT;
        double gtForKss = initialStiff ? st.gpT[0] : gTeff;
        bool useConsistent = consistentTan && !initialStiff;
        double Kss_s = 0.0, dTN_s = 0.0;
        LadrunoFrictionKernel::tangentBlock1D(gtForKss, st.gpT[0], N_I, kn, kt, mu,
                                              useConsistent, Kss_s, dTN_s, cohesion, tauMax);
        double K2[2][2];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                K2[i][j] = Kss_s * th[i] * th[j] + dTN_s * th[i] * n[j];
        for (int A = 0; A < nN; A++) {
            double bIA = (A < 2) ? D[I][A] : -M[I][A - 2];
            if (bIA == 0.0) continue;
            for (int B = 0; B < nN; B++) {
                double bIB = (B < 2) ? D[I][B] : -M[I][B - 2];
                double w = fact * bIA * bIB / aFacet;
                if (w == 0.0) continue;
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        tang(2*A + i, 2*B + j) += w * K2[i][j];
            }
        }
    }
}

const Matrix &
LadrunoContactFE::getTangent(Integrator *theIntegrator)
{
    // The integrator decides the K/C/M combination: an explicit integrator that
    // forms only M (and C) leaves the contact stiffness out of its matrix;
    // Newmark calls addKtToTang(c1); a static integrator addKtToTang(1).
    tang.Zero();
    if (theIntegrator != 0)
        theIntegrator->formEleTangent(this);
    return tang;
}

void
LadrunoContactFE::zeroTangent(void)
{
    tang.Zero();
}

void
LadrunoContactFE::addKtToTang(double fact)
{
    // rigid plane: K_c = kn (n (x) n) over the slave DOFs while penetrating
    if (mode == RIGID_PLANE && rigidPlaneGap() < 0.0) {
        for (int i = 0; i < ndm; i++)
            for (int j = 0; j < ndm; j++)
                tang(i, j) += fact * kn * planeN[i] * planeN[j];
    } else if (mode == SEGMENT && ndm == 2) {
        // kn B^T B on the pair this adapter owns (the side sign committed by the
        // residual); the curvature term of -geomtan is not used in 2D
        double gap, n[3], N2[2], B[6], gTnow = 0.0;
        if (segment2DActive(gap, n, N2, B, vtx2DCommittedSide(), 0, &gTnow)) {
            int ndof = 2 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            if (hasFriction() && theDomain != 0) {
                LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
                if (cd != 0) {
                    LadrunoContactDomain::FrictionState &st =
                        cd->getOrCreateFrictionState(contactTag, theSlave->getTag(), segIndex);
                    double gTeff = st.engaged ? (gTnow - st.gT0[0]) : 0.0;
                    double tn = ladrunoContactTraction(kn, gap);
                    addFrictionTang2D(fact, n, N2, tn, gTeff, st.gpT[0], consistentTan);
                }
            }
        }
    } else if (mode == SEGMENT) {
        // K_c = kn B^T B, exact for a flat fixed master where n is constant;
        // with consistentNormal the geometric block kn gN d^2gN/du^2 is added
        // (zero for a flat facet).
        double gap, n[3], N[4], B[15], gTvec[3];
        if (segmentActive(gap, n, N, B, gTvec)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            if (consistentNormal)
                addNormalGeomTang(fact);
            // friction tangent, from the committed gpT so that it is the
            // derivative of the residual evaluated at the same state
            if (hasFriction() && theDomain != 0) {
                LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
                if (cd != 0) {
                    LadrunoContactDomain::FrictionState &st =
                        cd->getOrCreateFrictionState(contactTag, theSlave->getTag(), segIndex);
                    double gTeff[3];
                    for (int d = 0; d < 3; d++)
                        gTeff[d] = st.engaged ? (gTvec[d] - st.gT0[d]) : 0.0;
                    double tn = ladrunoContactTraction(kn, gap);
                    addFrictionTang(fact, n, N, tn, gTeff, st.gpT, consistentTan);
                }
            }
        }
    } else if (mode == MORTAR && ndm == 2) {
        addMortarTang2D(fact);
    } else if (mode == MORTAR) {
        addMortarTang(fact);
    } else if (mode == EDGE_EDGE) {
        addEdgeTang(fact, /*initialStiff=*/false);
    }
}

void
LadrunoContactFE::addKiToTang(double fact)
{
    // Initial-stiffness algorithms (Newton -initial, ModifiedNewton -initial)
    // form the matrix from K_initial; the penalty K_c is used so the contact
    // stiffness is not dropped.
    if (mode == RIGID_PLANE && rigidPlaneGap() < 0.0) {
        for (int i = 0; i < ndm; i++)
            for (int j = 0; j < ndm; j++)
                tang(i, j) += fact * kn * planeN[i] * planeN[j];
    } else if (mode == SEGMENT && ndm == 2) {
        double gap, n[3], N2[2], B[6];
        if (segment2DActive(gap, n, N2, B, vtx2DCommittedSide(), 0)) {
            int ndof = 2 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            if (hasFriction()) {           // the stick block
                double tn = ladrunoContactTraction(kn, gap);
                addFrictionTang2D(fact, n, N2, tn, 0.0, 0.0, false);
            }
        }
    } else if (mode == SEGMENT) {
        double gap, n[3], N[4], B[15];
        if (segmentActive(gap, n, N, B)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * kn * B[i] * B[j];
            // friction: the stick tangent kt G^T P_t G (SPD, so a modified or
            // initial-stiffness Newton iteration contracts); gTeff = gpT = 0
            // forces the stick branch, so no engine state is needed
            if (hasFriction()) {
                double tn = ladrunoContactTraction(kn, gap);
                double zero[3] = {0.0, 0.0, 0.0};
                addFrictionTang(fact, n, N, tn, zero, zero, false);
            }
        }
    } else if (mode == MORTAR && ndm == 2) {
        addMortarTang2D(fact, /*initialStiff=*/true);
    } else if (mode == MORTAR) {
        // the penalty part equals the current tangent (variation of D, M and n
        // not included)
        addMortarTang(fact, /*initialStiff=*/true);
    } else if (mode == EDGE_EDGE) {
        // the penalty part equals the current tangent; friction uses the stick block
        addEdgeTang(fact, /*initialStiff=*/true);
    }
}

void
LadrunoContactFE::addCtoTang(double fact)
{
    // The only contact damping: C = muc B^T B (the normal penalty operator with
    // kn replaced by muc) on the active set of the residual dashpot. Transient
    // integrators call addCtoTang(c2), so implicit Newton gets the consistent
    // damping tangent; a static integrator does not call it.
    if (muc <= 0.0)
        return;
    if (mode == RIGID_PLANE && rigidPlaneGap() < 0.0) {
        for (int i = 0; i < ndm; i++)
            for (int j = 0; j < ndm; j++)
                tang(i, j) += fact * muc * planeN[i] * planeN[j];
    } else if (mode == SEGMENT && ndm == 2) {
        double gap, n[3], N2[2], B[6];
        if (segment2DActive(gap, n, N2, B, vtx2DCommittedSide(), 0)) {
            int ndof = 2 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * muc * B[i] * B[j];
        }
    } else if (mode == SEGMENT) {
        double gap, n[3], N[4], B[15];
        if (segmentActive(gap, n, N, B)) {
            int ndof = 3 * (1 + nps);
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    tang(i, j) += fact * muc * B[i] * B[j];
        }
    }
}

void
LadrunoContactFE::addMtoTang(double)
{
}

const Vector &LadrunoContactFE::getTangForce(const Vector &, double) { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getK_Force(const Vector &, double)   { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getKi_Force(const Vector &, double)  { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getC_Force(const Vector &, double)   { resid.Zero(); return resid; }
const Vector &LadrunoContactFE::getM_Force(const Vector &, double)   { resid.Zero(); return resid; }
