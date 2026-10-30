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

// penalty normal traction tn = kn <-gap>+ (Macaulay bracket): zero for gap >= 0,
// so the contact never pulls.
static inline double
ladrunoContactTraction(double kn, double gap)
{
    return (gap < 0.0) ? kn * (-gap) : 0.0;
}

LadrunoContactFE::LadrunoContactFE(int tag, Node *slaveNode, int ndm_,
                                   const double p0[3], const double n[3], double kn_,
                                   Domain *dom, int contactTag_, double muc_)
  : FE_Element(tag, /*numDOF_Group=*/1, /*ndof=*/ndm_),
    resid(ndm_), tang(ndm_, ndm_),
    mode(RIGID_PLANE), theSlave(slaveNode), ndm(ndm_), kn(kn_), nps(0),
    theDomain(dom), contactTag(contactTag_), segIndex(0), consistentNormal(false),
    kt(0.0), mu(0.0), cohesion(0.0), tauMax(0.0), consistentTan(false),
    muc(muc_ > 0.0 ? muc_ : 0.0), npsS(0), npsM(0), slaveFacetIndex(0), mortarOrder(2),
    isTie(false)
{
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
    isTie(false)
{
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
                                   Domain *dom, int ngp_, bool isTie_)
  : FE_Element(tag, /*numDOF_Group=*/nps_s + nps_m, /*ndof=*/3 * (nps_s + nps_m)),
    resid(3 * (nps_s + nps_m)), tang(3 * (nps_s + nps_m), 3 * (nps_s + nps_m)),
    mode(MORTAR), theSlave(0), ndm(3), kn(epsN), nps(0),
    theDomain(dom), contactTag(contactTag_), segIndex(0), consistentNormal(false),
    kt(0.0), mu(0.0), cohesion(0.0), tauMax(0.0), consistentTan(false), muc(0.0),
    npsS(nps_s), npsM(nps_m), slaveFacetIndex(slaveFacetIndex_),
    mortarOrder(ngp_ >= 1 ? ngp_ : 2), isTie(isTie_)
{
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

LadrunoContactFE::~LadrunoContactFE()
{
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
    // the master facet normal oriented by orientDir, the normal the weighted
    // gap was measured along (flat facets)
    return LadrunoMortarKernel::facetNormal(npsM, Xm, orientDir, n);
}

void
LadrunoContactFE::addMortarTang(double fact, bool)
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
        if (pr < 0.0) W[I] = 1.0 / aFacet;            // active in compression
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
                p[I] = (pr < 0.0) ? pr : 0.0;             // active in compression
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
        } else if (cd != 0) {
            // empty overlap at this evaluation: withdraw this pair's part of the
            // summed gap of its slave nodes
            for (int I = 0; I < npsS; I++)
                if (mortarSlave[I] != 0)
                    cd->accumulateMortarGap(contactTag, mortarSlave[I]->getTag(), this->getTag(),
                                            0.0, 0.0, kn);
        }
    }
    return resid;
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
    } else if (mode == MORTAR) {
        addMortarTang(fact);
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
    } else if (mode == MORTAR) {
        // the penalty part equals the current tangent (variation of D, M and n
        // not included)
        addMortarTang(fact, /*initialStiff=*/true);
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
