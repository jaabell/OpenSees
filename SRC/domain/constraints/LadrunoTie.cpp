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
// References: see LadrunoTie.h.
//
// Description: the surface tie generators (LadrunoTie::generate,
// generateMortar, generateShellSolid) and the command front-end
// OPS_LadrunoTie, registered in the Python and Tcl interpreters.

#include "LadrunoTie.h"

#include <Domain.h>
#include <Node.h>
#include <Element.h>
#include <ElementIter.h>
#include <EQ_Constraint.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <classTags.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

// header-only geometry kernels shared with the contact package
#include <LadrunoContactProjection.h>   // closest-point projection + shape functions
#include <LadrunoContactBucketSort.h>   // bucket-sort broad phase
#include <LadrunoMortarKernel.h>        // clip + Gauss -> mortar D, M, weighted gap

#ifdef HANDLER_TAG_LadrunoProjectionHandler
#include <ConstraintHandler.h>
#endif

#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <string.h>

// --------------------------------------------------------------------------- //
// helpers
// --------------------------------------------------------------------------- //

// Padded 3D reference coords of a node (z=0 for a 2D node). Returns false if the
// node is not in the domain.
static bool ltNodeCoords3(Domain *dom, int tag, double x[3])
{
    Node *nd = dom->getNode(tag);
    if (nd == 0) return false;
    const Vector &X = nd->getCrds();
    x[0] = x[1] = x[2] = 0.0;
    for (int d = 0; d < X.Size() && d < 3; d++) x[d] = X(d);
    return true;
}

// Default tied-DOF set for a slave node: 1..ndm translations, PLUS the 3 rotations
// (4,5,6) when the node is a 3D ndf-6 shell/beam node, so the rotation field is tied
// with the same transfer P (theta_s = sum P_sk theta_{m,k}). A 3D ndf-3 solid node and
// a 2D node tie translations only. -dof overrides the set (e.g. to drop the drilling
// DOF, or to tie translations only).
static void
ltDefaultDofs(Node *s0, ID &dofs)
{
    int ndm = s0->getCrds().Size();
    int ndf = s0->getNumberDOF();
    int nd  = (ndm == 3 && ndf == 6) ? 6 : ndm;      // 3D shell/beam ties its rotations too
    dofs.resize(nd);
    for (int d = 0; d < nd; d++) dofs(d) = d + 1;     // 1-based
}

// Lumped mass on the tied DOFs is a requirement of an explicit projection-type
// constraint handler only: such a handler keeps every tied DOF in a diagonal-mass
// equation set, where a massless or negative-mass tied DOF makes the system singular
// or unstable. Lagrange and Penalty (static or implicit) have no such requirement. The
// check therefore runs only when the constraint handler active at the time of the
// LadrunoTie command is of that type (its class tag defined in classTags.h); with any
// other handler, or none yet, the mass is not inspected.
static bool
ltTiedMassRequired(void)
{
#ifdef HANDLER_TAG_LadrunoProjectionHandler
    ConstraintHandler **theHandler = OPS_GetHandler();
    return (theHandler != 0 && *theHandler != 0 &&
            (*theHandler)->getClassTag() == HANDLER_TAG_LadrunoProjectionHandler);
#else
    return false;
#endif
}

// Per-node set of 1-based DOFs carrying a nonzero mass diagonal, scanned from ELEMENT
// mass matrices (at model-build the nodal mass() is still 0 for solids and shells: their
// mass comes from the element density). The scan is per DOF, so a shell that lumps
// translational mass but neglects rotational inertia (ShellMITC4 and ASDShellQ4 zero
// DOFs 4,5,6 in getMass()) is caught on its tied rotational DOFs. The element mass
// matrix orders each node's rows by getNumberDOF() in external-node order, so global
// DOF d of the node at external position k sits at diagonal (offset + d - 1).
static void
ltScanMassedDOFs(Domain *dom, std::map<int, std::map<int, double> > &massedDOF)
{
    ElementIter &eIter = dom->getElements();
    Element *elePtr;
    while ((elePtr = eIter()) != 0) {
        const Matrix &M = elePtr->getMass();
        int n = M.noRows();
        if (n == 0) continue;
        const ID &en = elePtr->getExternalNodes();
        int off = 0;
        for (int k = 0; k < en.Size(); k++) {
            Node *nd = dom->getNode(en(k));
            int ndf = (nd != 0) ? nd->getNumberDOF() : 0;
            for (int d = 0; d < ndf && (off + d) < n; d++)
                if (std::fabs(M(off + d, off + d)) > 0.0)
                    massedDOF[en(k)][d + 1] += M(off + d, off + d);   // 1-based DOF, SIGNED sum
            off += ndf;
        }
    }
}

// True if slave node `s` carries POSITIVE lumped mass on EVERY tied DOF; else print a
// named refusal and return false. `massedDOF` holds the SIGNED per-DOF element-mass
// diagonal sum (ltScanMassedDOFs); nodal mass is added per DOF. A NEGATIVE assembled
// diagonal (row-sum lumping of a serendipity element gives a hex20 corner -M/8) is
// refused as well.
static bool
ltCheckTiedDofMass(int s, Node *sNode, const ID &dofs,
                   std::map<int, std::map<int, double> > &massedDOF, const char *ctx)
{
    std::map<int, std::map<int, double> >::iterator it = massedDOF.find(s);
    const std::map<int, double> *emass = (it != massedDOF.end()) ? &it->second : 0;
    const Matrix &Mn = sNode->getMass();
    for (int di = 0; di < dofs.Size(); di++) {
        int dof1 = dofs(di);                             // 1-based
        double total = 0.0;
        if (emass != 0) {
            std::map<int, double>::const_iterator jt = emass->find(dof1);
            if (jt != emass->end()) total += jt->second;
        }
        if ((dof1 - 1) < Mn.noRows()) total += Mn(dof1 - 1, dof1 - 1);
        if (total < 0.0) {
            opserr << "WARNING LadrunoTie" << ctx << " - tied slave node " << s
                   << " carries NEGATIVE lumped mass " << total << " on DOF " << dof1
                   << " (row-sum lumping of a serendipity element gives negative corner "
                      "masses, e.g. -M/8 at a hex20 corner). The active constraint handler "
                      "needs a positive diagonal mass: use HRZ lumping (Hinton, Rock and "
                      "Zienkiewicz 1976) or add nodal mass.\n";
            return false;
        }
        if (!(total > 0.0)) {
            opserr << "WARNING LadrunoTie" << ctx << " - tied slave node " << s
                   << " carries no mass on DOF " << dof1;
            if (dof1 >= 4)
                opserr << " (a ROTATIONAL DOF). ShellMITC4 / ASDShellQ4 getMass() neglect "
                          "rotational inertia, so a shell tie needs explicit nodal rotary mass "
                          "(`mass $node 0 0 0 mrx mry mrz`), or tie translations only with "
                          "-dof (e.g. -dof 3 1 2 3).\n";
            else
                opserr << ". Define LadrunoTie AFTER the elements and masses; the active "
                          "constraint handler requires lumped mass on every tied DOF.\n";
            return false;
        }
    }
    return true;
}

// --------------------------------------------------------------------------- //
// Hermite (rotation-consistent) w-theta emission for a shell EDGE tie
// --------------------------------------------------------------------------- //

// Emit ONE EQ_Constraint with mixed retained (node, dof, coef) triples, dropping
// near-zero coefficients. rows: nk retained nodes x 6 DOFs of coefficients.
static bool
ltEmitMixedRow(Domain *dom, int s, int slaveDof0, int nk, const int *mTags,
               const double coef[][6], const char *ctx)
{
    int cnt = 0;
    for (int k = 0; k < nk; k++)
        for (int j = 0; j < 6; j++)
            if (std::fabs(coef[k][j]) > 1e-12) cnt++;
    if (cnt == 0) {
        opserr << "WARNING LadrunoTie" << ctx << " - slave node " << s << " dof "
               << slaveDof0 + 1 << " produced an empty mixed-coefficient row.\n";
        return false;
    }
    ID rNode(cnt), rDOF(cnt);
    Vector Ccr(cnt);
    int p = 0;
    for (int k = 0; k < nk; k++)
        for (int j = 0; j < 6; j++)
            if (std::fabs(coef[k][j]) > 1e-12) {
                rNode(p) = mTags[k];
                rDOF(p)  = j;                       // 0-based, stored
                Ccr(p)   = coef[k][j];
                p++;
            }
    EQ_Constraint *eq = new EQ_Constraint(s, slaveDof0, Ccr, rNode, rDOF);
    if (eq == 0 || dom->addEQ_Constraint(eq) == false) {
        opserr << "WARNING LadrunoTie" << ctx << " - failed to add EQ_Constraint for slave node "
               << s << " dof " << slaveDof0 + 1 << "\n";
        if (eq != 0) delete eq;
        return false;
    }
    return true;
}

// Hermite rows for ONE slave node projecting onto a master facet EDGE (or a
// coincident corner). The shell-NORMAL translation and the slope along the interface
// are reconstructed with the cubic Hermite basis over the edge nodes' (w, theta_t)
// pairs; the in-plane translations and the other rotations keep the linear weights.
// Emits 6 EQ_Constraints (u1..u3, theta4..theta6); returns the count emitted, or -1 on
// a named refusal:
//   u_s     = (I - n n) sum N_k u_k     + n [sum Hw_k (n.u_k)  + Hp_k (e.theta_k)]
//   theta_s = (I - e e) sum N_k theta_k + e [sum dHw_k (n.u_k) + dHp_k (e.theta_k)]
// with t the edge tangent, n the facet normal at the projection, e = t x n, and
// slope = dw/ds = theta.e (Kirchhoff: u = theta x r gives grad w = n x theta, so
// t.grad w = theta.(t x n)). The transfer is invariant under n -> -n (n and e enter in
// even products only).
static int
ltEmitHermiteRows(Domain *dom, int s, Node *sNode,
                  const ID &mfn, int bestF, int nps,
                  const double X[4][3], const double N[4],
                  const double dNx[4], const double dNe[4])
{
    // --- ndf-6 preconditions (shell-to-shell, all 6 DOFs) ---
    if (sNode->getNumberDOF() < 6) {
        opserr << "WARNING LadrunoTie -hermite - slave node " << s << " has ndf "
               << sNode->getNumberDOF() << "; the Hermite w-theta transfer needs ndf-6 "
                  "shell nodes on both sides (it couples translations to rotations).\n";
        return -1;
    }

    // --- which master nodes carry weight: corner (1) / edge (2) / interior (3+) ---
    int nz[4], cnt = 0;
    for (int i = 0; i < nps; i++)
        if (std::fabs(N[i]) > 1e-12) nz[cnt++] = i;

    if (cnt == 1) {
        // coincident with a master corner: the Hermite transfer degenerates to the
        // identity on that node (w_s = w_a, slope_s = slope_a; edge choice irrelevant).
        int mt = mfn(bestF * nps + nz[0]);
        Node *mn = dom->getNode(mt);
        if (mn == 0 || mn->getNumberDOF() < 6) {
            opserr << "WARNING LadrunoTie -hermite - master node " << mt << " has ndf "
                   << (mn ? mn->getNumberDOF() : 0) << " (< 6). The Hermite transfer ties "
                      "shell to shell. For a shell EDGE tied to a solid FACE use -shellSolid "
                      "(the SOLID face nodes are the slaves).\n";
            return -1;
        }
        int emitted = 0;
        for (int d = 0; d < 6; d++) {
            double coef[1][6] = {{0, 0, 0, 0, 0, 0}};
            coef[0][d] = 1.0;
            if (!ltEmitMixedRow(dom, s, d, 1, &mt, coef, " -hermite")) return -1;
            emitted++;
        }
        return emitted;
    }

    if (cnt != 2 || ((nz[1] - nz[0]) != 1 && (nz[1] - nz[0]) != nps - 1)) {
        opserr << "WARNING LadrunoTie -hermite - slave node " << s << " projects into the "
                  "INTERIOR of a master facet (not onto an edge). The Hermite transfer is an "
                  "EDGE (butt-joint) tie: slave nodes must lie on master facet edges. Use "
                  "the standard tie (drop -hermite) for surface-overlap collocation.\n";
        return -1;
    }
    // order the edge (a,b) along the facet winding; for the wrap pair (0, nps-1) the
    // edge runs (nps-1)->0
    int ia = nz[0], ib = nz[1];
    if (ib - ia == nps - 1) { ia = nz[1]; ib = nz[0]; }

    int mA = mfn(bestF * nps + ia), mB = mfn(bestF * nps + ib);
    for (int k = 0; k < 2; k++) {
        int mt = (k == 0) ? mA : mB;
        Node *mn = dom->getNode(mt);
        if (mn == 0 || mn->getNumberDOF() < 6) {
            opserr << "WARNING LadrunoTie -hermite - master node " << mt << " has ndf "
                   << (mn ? mn->getNumberDOF() : 0) << " (< 6). The Hermite transfer ties "
                      "shell to shell. For a shell EDGE tied to a solid FACE use -shellSolid "
                      "(the SOLID face nodes are the slaves).\n";
            return -1;
        }
    }

    // --- the (t, n, e) frame at the projection ---
    double tv[3], hlen = 0.0;
    for (int d = 0; d < 3; d++) { tv[d] = X[ib][d] - X[ia][d]; hlen += tv[d] * tv[d]; }
    hlen = std::sqrt(hlen);
    if (hlen < 1e-300) {
        opserr << "WARNING LadrunoTie -hermite - degenerate (zero-length) master edge "
               << mA << "-" << mB << "\n";
        return -1;
    }
    for (int d = 0; d < 3; d++) tv[d] /= hlen;

    double g1[3] = {0, 0, 0}, g2[3] = {0, 0, 0}, nn[3];
    for (int i = 0; i < nps; i++)
        for (int d = 0; d < 3; d++) { g1[d] += dNx[i] * X[i][d]; g2[d] += dNe[i] * X[i][d]; }
    nn[0] = g1[1]*g2[2] - g1[2]*g2[1];
    nn[1] = g1[2]*g2[0] - g1[0]*g2[2];
    nn[2] = g1[0]*g2[1] - g1[1]*g2[0];
    // strip any tangent component (a warped quad's normal may lean), then normalize
    double nt = nn[0]*tv[0] + nn[1]*tv[1] + nn[2]*tv[2];
    for (int d = 0; d < 3; d++) nn[d] -= nt * tv[d];
    double Ln = std::sqrt(nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2]);
    if (Ln < 1e-300) {
        opserr << "WARNING LadrunoTie -hermite - degenerate master facet normal at slave node "
               << s << "\n";
        return -1;
    }
    for (int d = 0; d < 3; d++) nn[d] /= Ln;
    double ee[3] = { tv[1]*nn[2] - tv[2]*nn[1],
                     tv[2]*nn[0] - tv[0]*nn[2],
                     tv[0]*nn[1] - tv[1]*nn[0] };

    // --- cubic Hermite basis at the edge parameter xi (slope DOFs carry hlen) ---
    double xi = N[ib] / (N[ia] + N[ib]);        // linear edge weights: N_a = 1-xi, N_b = xi
    double x2 = xi * xi, x3 = x2 * xi;
    double Nlin[2] = { 1.0 - xi, xi };
    double Hw[2]  = { 1.0 - 3.0*x2 + 2.0*x3,        3.0*x2 - 2.0*x3 };
    double Hp[2]  = { hlen * (xi - 2.0*x2 + x3),    hlen * (x3 - x2) };
    double dHw[2] = { (-6.0*xi + 6.0*x2) / hlen,    (6.0*xi - 6.0*x2) / hlen };
    double dHp[2] = { 1.0 - 4.0*xi + 3.0*x2,        3.0*x2 - 2.0*xi };

    // --- assemble + emit the 6 mixed rows ---
    int mTags[2] = { mA, mB };
    int emitted = 0;
    for (int d = 0; d < 3; d++) {               // translation rows u_s,d
        double coef[2][6];
        for (int k = 0; k < 2; k++)
            for (int j = 0; j < 3; j++) {
                double dij = (d == j) ? 1.0 : 0.0;
                coef[k][j]     = (dij - nn[d]*nn[j]) * Nlin[k] + nn[d] * Hw[k] * nn[j];
                coef[k][3 + j] = nn[d] * Hp[k] * ee[j];
            }
        if (!ltEmitMixedRow(dom, s, d, 2, mTags, coef, " -hermite")) return -1;
        emitted++;
    }
    for (int r = 0; r < 3; r++) {               // rotation rows theta_s,r
        double coef[2][6];
        for (int k = 0; k < 2; k++)
            for (int j = 0; j < 3; j++) {
                double drj = (r == j) ? 1.0 : 0.0;
                coef[k][j]     = ee[r] * dHw[k] * nn[j];
                coef[k][3 + j] = (drj - ee[r]*ee[j]) * Nlin[k] + ee[r] * dHp[k] * ee[j];
            }
        if (!ltEmitMixedRow(dom, s, 3 + r, 2, mTags, coef, " -hermite")) return -1;
        emitted++;
    }
    return emitted;
}

// --------------------------------------------------------------------------- //
// collocation generator
// --------------------------------------------------------------------------- //

int
LadrunoTie::generate(Domain *dom,
                     const ID &slaves,
                     const ID &mfn, int nps,
                     const ID &dofsIn, double tolFrac,
                     bool hermite)
{
    if (dom == 0) {
        opserr << "WARNING LadrunoTie - domain is not defined\n";
        return -1;
    }
    if (nps != 3 && nps != 4) {
        opserr << "WARNING LadrunoTie - nps must be 3 (tri-3) or 4 (quad-4), got " << nps
               << ". Quadratic master facets (tri-6, quad-8) are supported by -mortar only; "
                  "-hermite is collocation-only and linear-facet-only.\n";
        return -1;
    }
    if (slaves.Size() == 0) {
        opserr << "WARNING LadrunoTie - no slave nodes given\n";
        return -1;
    }
    if (mfn.Size() == 0 || (mfn.Size() % nps) != 0) {
        opserr << "WARNING LadrunoTie - master facet node count " << mfn.Size()
               << " is not a positive multiple of nps " << nps << "\n";
        return -1;
    }
    const int nf = mfn.Size() / nps;

    // tied DOFs: default = 1..ndm translations (+ rotations 4,5,6 for a 3D ndf-6 node).
    ID dofs = dofsIn;
    if (dofs.Size() == 0) {
        Node *s0 = dom->getNode(slaves(0));
        if (s0 == 0) {
            opserr << "WARNING LadrunoTie - slave node " << slaves(0) << " not in domain\n";
            return -1;
        }
        ltDefaultDofs(s0, dofs);
    }

    // The Hermite transfer couples translations to rotations, so a partial DOF set
    // would emit a one-sided coupling: require the full 1..6 tie (the default for a
    // 3D ndf-6 node).
    if (hermite) {
        bool full6 = (dofs.Size() == 6);
        for (int di = 0; full6 && di < 6; di++)
            if (dofs(di) != di + 1) full6 = false;
        if (!full6) {
            opserr << "WARNING LadrunoTie -hermite - the Hermite w-theta transfer couples "
                      "translations to rotations and needs the FULL 6-DOF tie (1..6, the "
                      "default for a 3D ndf-6 node). Drop -dof, or drop -hermite for a "
                      "partial tie.\n";
            return -1;
        }
    }

    // --- master facets: flat reference coords (nf*nps*3), node-disjoint set, sizes ---
    std::vector<double> seg((size_t)nf * nps * 3, 0.0);
    std::set<int> masterSet;
    std::vector<double> facetSize(nf, 0.0);
    for (int f = 0; f < nf; f++) {
        double Xf[4][3];
        for (int i = 0; i < nps; i++) {
            int tag = mfn(f * nps + i);
            double x[3];
            if (!ltNodeCoords3(dom, tag, x)) {
                opserr << "WARNING LadrunoTie - master facet node " << tag << " not in domain\n";
                return -1;
            }
            masterSet.insert(tag);
            for (int d = 0; d < 3; d++) { seg[((size_t)f * nps + i) * 3 + d] = x[d]; Xf[i][d] = x[d]; }
        }
        // characteristic facet size = max pairwise node distance (the diagonal)
        double dmax = 0.0;
        for (int i = 0; i < nps; i++)
            for (int j = i + 1; j < nps; j++) {
                double dx = Xf[i][0]-Xf[j][0], dy = Xf[i][1]-Xf[j][1], dz = Xf[i][2]-Xf[j][2];
                double d = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (d > dmax) dmax = d;
            }
        facetSize[f] = dmax;
    }

    // --- broad phase over the master facets (reference configuration) ---
    LadrunoContactBucketSort::Grid grid(nf, nps, seg.data());
    std::vector<int> cand(nf);

    // --- per-DOF element mass, only when the active handler requires it ---
    const bool massRequired = ltTiedMassRequired();
    std::map<int, std::map<int, double> > massedDOF;
    if (massRequired)
        ltScanMassedDOFs(dom, massedDOF);

    std::set<int> doneSlaves;
    int emitted = 0;

    for (int si = 0; si < slaves.Size(); si++) {
        int s = slaves(si);

        // a slave listed twice would be doubly constrained
        if (doneSlaves.count(s)) {
            opserr << "WARNING LadrunoTie - slave node " << s << " is listed more than once. Each "
                      "slave node is tied to exactly one facet (no double constraints). Remove the "
                      "duplicate.\n";
            return -1;
        }
        doneSlaves.insert(s);

        // slave and master surfaces must be node-disjoint (no chained constraints)
        if (masterSet.count(s)) {
            opserr << "WARNING LadrunoTie - node " << s << " is BOTH a slave and a master facet node. "
                      "The slave and master surfaces must be node-disjoint (a slave retained by "
                      "another tie row would chain the constraints).\n";
            return -1;
        }

        Node *sNode = dom->getNode(s);
        double xs[3];
        if (sNode == 0 || !ltNodeCoords3(dom, s, xs)) {
            opserr << "WARNING LadrunoTie - slave node " << s << " not in domain\n";
            return -1;
        }

        if (massRequired && !ltCheckTiedDofMass(s, sNode, dofs, massedDOF, ""))
            return -1;

        // --- nearest in-bounds facet (broad phase first, brute-force fallback) ---
        int    bestF = -1;
        double bestDist = 0.0, bestXi = 0.0, bestEta = 0.0;
        for (int pass = 0; pass < 2 && bestF < 0; pass++) {
            int nc = (pass == 0) ? grid.candidates(xs, cand.data()) : nf;
            for (int ci = 0; ci < nc; ci++) {
                int f = (pass == 0) ? cand[ci] : ci;
                double X[4][3];
                for (int i = 0; i < nps; i++)
                    for (int d = 0; d < 3; d++) X[i][d] = seg[((size_t)f * nps + i) * 3 + d];
                double xi, eta;
                if (LadrunoContactProjection::project(nps, X, xs, xi, eta, 1e-12, 30) != 0)
                    continue;   // out-of-bounds or degenerate facet
                double N[4], dNx[4], dNe[4], xb[3];
                LadrunoContactProjection::shape(nps, xi, eta, N, dNx, dNe);
                LadrunoContactProjection::interp(nps, N, X, xb);
                double dd = std::sqrt((xs[0]-xb[0])*(xs[0]-xb[0]) +
                                      (xs[1]-xb[1])*(xs[1]-xb[1]) +
                                      (xs[2]-xb[2])*(xs[2]-xb[2]));
                if (bestF < 0 || dd < bestDist) { bestF = f; bestDist = dd; bestXi = xi; bestEta = eta; }
            }
        }
        if (bestF < 0) {
            opserr << "WARNING LadrunoTie - slave node " << s << " does not project in-bounds onto any "
                      "master facet. Check the master surface covers the slave surface (nps/facets correct).\n";
            return -1;
        }

        // conforming geometry: the slave must lie on the master surface
        double tol = tolFrac * (facetSize[bestF] > 0.0 ? facetSize[bestF] : 1.0);
        if (bestDist > tol) {
            opserr << "WARNING LadrunoTie - slave node " << s << " projects " << bestDist
                   << " off the master surface (tol " << tol << " = " << tolFrac
                   << " x facet size). The tie requires conforming interface geometry: the slave "
                      "must lie on the master surface (it is not moved onto it). Move the node "
                      "onto the interface, or relax -tol if the offset is intentional.\n";
            return -1;
        }

        // --- shape weights at the projection; emit one EQ_Constraint per tied DOF ---
        double X[4][3];
        for (int i = 0; i < nps; i++)
            for (int d = 0; d < 3; d++) X[i][d] = seg[((size_t)bestF * nps + i) * 3 + d];
        double N[4], dNx[4], dNe[4];
        LadrunoContactProjection::shape(nps, bestXi, bestEta, N, dNx, dNe);

        // the Hermite w-theta transfer replaces the per-DOF linear rows entirely (mixed
        // u <- (u, theta) rows over the projected master EDGE)
        if (hermite) {
            int e6 = ltEmitHermiteRows(dom, s, sNode, mfn, bestF, nps, X, N, dNx, dNe);
            if (e6 < 0) return -1;
            emitted += e6;
            continue;
        }

        // drop near-zero weights so a corner/edge collocation does not couple master
        // nodes that carry no share of this slave
        int cnt = 0;
        for (int i = 0; i < nps; i++) if (std::fabs(N[i]) > 1e-12) cnt++;

        for (int di = 0; di < dofs.Size(); di++) {
            int dof1 = dofs(di);                 // 1-based
            if (dof1 < 1 || dof1 > sNode->getNumberDOF()) {
                opserr << "WARNING LadrunoTie - slave node " << s << " has no DOF " << dof1
                       << " (numberDOF " << sNode->getNumberDOF() << ")\n";
                return -1;
            }
            // like-to-like only: every retained master node must carry this DOF (an ndf-6
            // shell slave tied to ndf-3 solid master nodes would reference nonexistent
            // master rotations)
            for (int i = 0; i < nps; i++)
                if (std::fabs(N[i]) > 1e-12) {
                    Node *mn = dom->getNode(mfn(bestF * nps + i));
                    if (mn == 0 || dof1 > mn->getNumberDOF()) {
                        opserr << "WARNING LadrunoTie - slave node " << s << " tied DOF " << dof1
                               << " references master node " << mfn(bestF * nps + i)
                               << " which has no DOF " << dof1 << " (ndf mismatch). This mode "
                                  "ties like to like (same ndf). For a shell EDGE tied to a "
                                  "solid FACE use -shellSolid (the SOLID face nodes are the "
                                  "slaves, -masterEdge is the shell edge), or restrict -dof to "
                                  "the DOFs both sides share (e.g. -dof 3 1 2 3).\n";
                        return -1;
                    }
                }
            ID rNode(cnt), rDOF(cnt);
            Vector Ccr(cnt);
            int p = 0;
            for (int i = 0; i < nps; i++)
                if (std::fabs(N[i]) > 1e-12) {
                    rNode(p) = mfn(bestF * nps + i);
                    rDOF(p)  = dof1 - 1;          // 0-based, stored
                    Ccr(p)   = N[i];              // u_s = sum N_i u_{m,i}  (Ccr = N, no negation)
                    p++;
                }
            EQ_Constraint *eq = new EQ_Constraint(s, dof1 - 1, Ccr, rNode, rDOF);
            if (eq == 0 || dom->addEQ_Constraint(eq) == false) {
                opserr << "WARNING LadrunoTie - failed to add EQ_Constraint for slave node " << s
                       << " dof " << dof1 << "\n";
                if (eq != 0) delete eq;
                return -1;
            }
            emitted++;
        }
    }

    opserr << "LadrunoTie: emitted " << emitted << " EQ_Constraint(s) tying " << (int)doneSlaves.size()
           << " slave node(s) to " << nf << " master facet(s). Enforce with constraints Lagrange "
              "or Penalty.\n";
    return emitted;
}

// --------------------------------------------------------------------------- //
// integral mortar generator
// --------------------------------------------------------------------------- //

int
LadrunoTie::generateMortar(Domain *dom,
                           const ID &sfn, int npsS,
                           const ID &mfn, int npsM,
                           const ID &dofsIn, double tolFrac, const double *outward,
                           bool dual)
{
    if (dom == 0) { opserr << "WARNING LadrunoTie - domain is not defined\n"; return -1; }
    // quadratic serendipity facets: quad8 on both sides, tri6 on the MASTER side only
    if (npsS == 6) {
        opserr << "WARNING LadrunoTie -mortar - tri6 SLAVE facets are refused: the integrals of "
                  "the tri6 corner shape functions over the facet are exactly zero, so the dual "
                  "scaling and the coverage measures divide by zero. Swap master and slave (tri6 "
                  "MASTER facets are accepted), or use quad8 faces on the slave side.\n";
        return -1;
    }
    if (npsS != 3 && npsS != 4 && npsS != 8) {
        opserr << "WARNING LadrunoTie -mortar - slave nps must be 3, 4, or 8, got " << npsS << "\n"; return -1; }
    if (npsM != 3 && npsM != 4 && npsM != 6 && npsM != 8) {
        opserr << "WARNING LadrunoTie -mortar - master nps must be 3, 4, 6, or 8, got " << npsM << "\n"; return -1; }
    if (sfn.Size() == 0 || (sfn.Size() % npsS) != 0) {
        opserr << "WARNING LadrunoTie -mortar - slave facet node count " << sfn.Size()
               << " not a positive multiple of npsS " << npsS << "\n"; return -1; }
    if (mfn.Size() == 0 || (mfn.Size() % npsM) != 0) {
        opserr << "WARNING LadrunoTie -mortar - master facet node count " << mfn.Size()
               << " not a positive multiple of npsM " << npsM << "\n"; return -1; }
    const int nfS = sfn.Size() / npsS;
    const int nfM = mfn.Size() / npsM;

    // --- unique slave / master node lists (tag <-> local index) ---
    std::map<int,int> sIdx, mIdx;
    std::vector<int> sTag, mTag;
    for (int i = 0; i < sfn.Size(); i++)
        if (sIdx.find(sfn(i)) == sIdx.end()) { sIdx[sfn(i)] = (int)sTag.size(); sTag.push_back(sfn(i)); }
    for (int i = 0; i < mfn.Size(); i++)
        if (mIdx.find(mfn(i)) == mIdx.end()) { mIdx[mfn(i)] = (int)mTag.size(); mTag.push_back(mfn(i)); }
    const int Ns = (int)sTag.size(), Nm = (int)mTag.size();

    // node-disjoint slave / master surfaces (no chained constraints)
    for (int i = 0; i < Ns; i++)
        if (mIdx.find(sTag[i]) != mIdx.end()) {
            opserr << "WARNING LadrunoTie -mortar - node " << sTag[i] << " is on BOTH the slave and "
                      "master surface. They must be node-disjoint (no chained constraints).\n";
            return -1;
        }

    // tied DOFs default = 1..ndm (+ rotations 4,5,6 for a 3D ndf-6 node)
    ID dofs = dofsIn;
    if (dofs.Size() == 0) {
        Node *s0 = dom->getNode(sTag[0]);
        if (s0 == 0) { opserr << "WARNING LadrunoTie -mortar - slave node " << sTag[0] << " not in domain\n"; return -1; }
        ltDefaultDofs(s0, dofs);
    }

    // --- master facet reference coords (flat nfM*npsM*3) + average outward normal ---
    const int ncM = LadrunoMortarKernel::ncorner(npsM);   // corners-first ordering
    std::vector<double> segM((size_t)nfM * npsM * 3, 0.0);
    double refDir[3] = {0.0, 0.0, 0.0};
    for (int f = 0; f < nfM; f++) {
        double X[8][3];
        for (int i = 0; i < npsM; i++) {
            double x[3];
            if (!ltNodeCoords3(dom, mfn(f * npsM + i), x)) {
                opserr << "WARNING LadrunoTie -mortar - master facet node " << mfn(f*npsM+i) << " not in domain\n";
                return -1;
            }
            for (int d = 0; d < 3; d++) { segM[((size_t)f*npsM+i)*3+d] = x[d]; X[i][d] = x[d]; }
        }
        // accumulate the (unoriented) CORNER-polygon normal for a default refDir
        double a1[3], a2[3], nrm[3];
        for (int d = 0; d < 3; d++) { a1[d] = X[1][d]-X[0][d]; a2[d] = X[ncM-1][d]-X[0][d]; }
        nrm[0]=a1[1]*a2[2]-a1[2]*a2[1]; nrm[1]=a1[2]*a2[0]-a1[0]*a2[2]; nrm[2]=a1[0]*a2[1]-a1[1]*a2[0];
        double L = std::sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]+nrm[2]*nrm[2]);
        if (L > 1e-300) for (int d = 0; d < 3; d++) refDir[d] += nrm[d]/L;
    }
    if (outward != 0) { for (int d = 0; d < 3; d++) refDir[d] = outward[d]; }
    double Lr = std::sqrt(refDir[0]*refDir[0]+refDir[1]*refDir[1]+refDir[2]*refDir[2]);
    if (Lr < 1e-12) {
        opserr << "WARNING LadrunoTie -mortar - could not infer an outward direction (master normals "
                  "cancel; the surface may be closed or curved). Pass -outward ox oy oz, a direction "
                  "not perpendicular to any facet.\n";
        return -1;
    }
    for (int d = 0; d < 3; d++) refDir[d] /= Lr;

    // --- refuse a SELF-OVERLAPPING master surface ---
    // The per-facet coverage check below sums clip areas over the master facets, so two
    // masters covering the same region double-count coverage and can MASK a genuine
    // hole (a doubly-listed master half exactly cancels a same-size slave protrusion).
    // Adjacent facets sharing an edge clip to a sliver (skipped); a finite COINCIDENT
    // overlap (near-zero mean gap: duplicate or doubly-listed facets, overlapping
    // patches) is a modeling error. The gap condition keeps a CURVED or closed master
    // surface (distant facets can overlap on each other's auxiliary planes, but with a
    // large mean gap) out of the refusal.
    for (int fa = 0; fa < nfM; fa++) {
        double Xa[8][3];
        for (int i = 0; i < npsM; i++)
            for (int d = 0; d < 3; d++) Xa[i][d] = segM[((size_t)fa*npsM+i)*3+d];
        LadrunoMortarKernel::PairResult prA;
        int rcA = LadrunoMortarKernel::integratePair(npsM, Xa, npsM, Xa, refDir, prA, 2);
        if (rcA == -2) {
            opserr << "WARNING LadrunoTie -mortar - master facet " << fa << " (first node "
                   << mfn(fa*npsM) << ") is an INVALID quadratic facet: a midside node is off "
                      "its edge midpoint (curved or distorted face). Quadratic facets must be "
                      "flat and straight-edged (midside nodes at the edge midpoints).\n";
            return -1;
        }
        double areaA = (rcA == 0) ? prA.area : 0.0;
        if (areaA <= 1e-300) continue;                  // degenerate facet: caught elsewhere
        for (int fb = fa + 1; fb < nfM; fb++) {
            double Xb[8][3];
            for (int i = 0; i < npsM; i++)
                for (int d = 0; d < 3; d++) Xb[i][d] = segM[((size_t)fb*npsM+i)*3+d];
            LadrunoMortarKernel::PairResult pr;
            if (LadrunoMortarKernel::integratePair(npsM, Xa, npsM, Xb, refDir, pr, 2) != 0)
                continue;                               // disjoint / edge-sliver / degenerate
            if (pr.area <= 1e-6 * areaA)
                continue;                               // numerical dust
            double gbar = pr.gapL1 / pr.area;           // mean separation over the overlap
            if (gbar <= tolFrac * 0.5 * std::sqrt(areaA)) {
                opserr << "WARNING LadrunoTie -mortar - master facets " << fa << " and " << fb
                       << " (first nodes " << mfn(fa*npsM) << ", " << mfn(fb*npsM)
                       << ") OVERLAP coincidently (shared area " << pr.area << ", mean gap "
                       << gbar << "). The master surface must not overlap itself: a "
                          "duplicated or doubly-listed facet double-counts the coverage check "
                          "and can mask an uncovered slave region. Remove the duplicate/overlap.\n";
                return -1;
            }
        }
    }

    // --- assemble global D (Ns x Ns), M (Ns x Nm) + per-facet checks ---
    // All slave x master facet pairs are integrated: the clip returns at once for a
    // non-overlapping pair, and the exhaustive sweep captures EVERY overlap.
    //
    // The checks are per slave facet and sign-free, for every facet order (the
    // quad8 corner row sums of D are NEGATIVE, -A/12, so nodal row-sum measures are
    // meaningless there):
    //   coverage  = sum over master facets of the overlap area, against the full facet
    //               area from the facet clipped with itself;
    //   gap       = sum of the L1 gap integrals / covered area (no cancellation).
    // Pair proximity: a pair whose own mean |gap| exceeds the conforming tolerance
    // tolFrac * 0.5 * sqrt(A) is EXCLUDED from D and M. On a closed or folded surface
    // (a tube, a box) the far side projects onto the auxiliary plane of a near facet
    // with a separation of the order of the surface diameter; integrating that pair
    // would tie each slave facet to the antipodal master facets as well. A normal
    // orientation test cannot detect those pairs (the facets are parallel), the
    // separation can. Tolerance 0.1 on a tube of radius 0.5 with 16 master and 24
    // slave facets: near pairs separate by at most 0.0096, antipodal pairs by 0.98.
    // Serendipity products need the degree-6 rule; linear pairs keep order 4.
    const int gpOrder = (npsS > 4 || npsM > 4) ? 6 : 4;
    Matrix D(Ns, Ns), M(Ns, Nm);
    D.Zero(); M.Zero();
    std::vector<double> pairGapTol(nfS, -1.0);   // per slave facet; < 0 = no exclusion
    double totalArea = 0.0;
    int nExcluded = 0;
    for (int fs = 0; fs < nfS; fs++) {
        double Xs[8][3];
        for (int i = 0; i < npsS; i++) {
            double x[3];
            if (!ltNodeCoords3(dom, sfn(fs*npsS+i), x)) {
                opserr << "WARNING LadrunoTie -mortar - slave facet node " << sfn(fs*npsS+i) << " not in domain\n";
                return -1;
            }
            for (int d = 0; d < 3; d++) Xs[i][d] = x[d];
        }
        // full facet area via self-clip (the facet against itself = the whole facet)
        LadrunoMortarKernel::PairResult prSelf;
        int rcSelf = LadrunoMortarKernel::integratePair(npsS, Xs, npsS, Xs, refDir, prSelf, gpOrder);
        if (rcSelf == -2) {
            opserr << "WARNING LadrunoTie -mortar - slave facet " << fs << " (first node "
                   << sfn(fs*npsS) << ") is an INVALID quadratic facet: a midside node is off "
                      "its edge midpoint (curved or distorted face). Quadratic facets must be "
                      "flat and straight-edged (midside nodes at the edge midpoints).\n";
            return -1;
        }
        double areaFull = (rcSelf == 0) ? prSelf.area : 0.0;
        if (areaFull > 1e-300)
            pairGapTol[fs] = tolFrac * 0.5 * std::sqrt(areaFull);
        double areaCov = 0.0, gapL1 = 0.0;
        double areaExcl = 0.0, gapExclMin = 0.0;
        for (int fm = 0; fm < nfM; fm++) {
            double Xm[8][3];
            for (int i = 0; i < npsM; i++)
                for (int d = 0; d < 3; d++) Xm[i][d] = segM[((size_t)fm*npsM+i)*3+d];
            LadrunoMortarKernel::PairResult pr;
            int rc = LadrunoMortarKernel::integratePair(npsS, Xs, npsM, Xm, refDir, pr, gpOrder);
            if (rc == -2) {
                opserr << "WARNING LadrunoTie -mortar - master facet " << fm << " (first node "
                       << mfn(fm*npsM) << ") is an INVALID quadratic facet: a midside node is off "
                          "its edge midpoint (curved or distorted face). Quadratic facets must be "
                          "flat and straight-edged (midside nodes at the edge midpoints).\n";
                return -1;
            }
            if (rc != 0)
                continue;                              // empty/degenerate overlap
            if (pairGapTol[fs] >= 0.0 && pr.gapL1 > pairGapTol[fs] * pr.area) {
                // pair farther apart than the conforming tolerance: not part of the bond
                double g = pr.gapL1 / pr.area;
                if (areaExcl <= 0.0 || g < gapExclMin) gapExclMin = g;
                areaExcl += pr.area;
                nExcluded++;
                continue;
            }
            areaCov += pr.area;
            gapL1   += pr.gapL1;
            for (int a = 0; a < npsS; a++) {
                int I = sIdx[sfn(fs*npsS+a)];
                for (int b = 0; b < npsS; b++)
                    D(I, sIdx[sfn(fs*npsS+b)]) += pr.D[a][b];
                for (int k = 0; k < npsM; k++)
                    M(I, mIdx[mfn(fm*npsM+k)]) += pr.M[a][k];
            }
        }
        totalArea += areaCov;
        // A slave facet must be covered over essentially its WHOLE area, else the slave
        // surface protrudes beyond the master there (a partial, extrapolated bond).
        // 0.1% slack for clip/quadrature noise.
        if (areaFull > 1e-300 && areaCov < (1.0 - 1.0e-3) * areaFull) {
            if (areaExcl > 0.0)
                opserr << "WARNING LadrunoTie -mortar - slave facet " << fs << " (first node "
                       << sfn(fs*npsS) << ") sits a mean |gap| of at least " << gapExclMin
                       << " off the master surface (tol " << pairGapTol[fs] << ") over "
                       << (1.0 - areaCov / areaFull) * 100.0 << "% of its area. The tie requires "
                          "conforming interface geometry (coincident surfaces). Move the slave "
                          "onto the interface, or relax -tol.\n";
            else
                opserr << "WARNING LadrunoTie -mortar - slave facet " << fs << " (first node "
                       << sfn(fs*npsS) << ") is only " << (areaCov / areaFull) * 100.0
                       << "% covered by the master surface (the slave surface protrudes beyond the "
                          "master). A tie needs the master surface to cover the whole slave "
                          "surface. Extend the master surface or trim the slave.\n";
            return -1;
        }
        if (areaCov > 1e-300) {
            // 0.5*sqrt(A) ~ the tributary length of a quad node (sqrt(A/4))
            double facetScale = 0.5 * std::sqrt((areaFull > 1e-300) ? areaFull : areaCov);
            double gbar = gapL1 / areaCov;             // mean |gap| over the bonded area
            if (gbar > tolFrac * facetScale) {
                opserr << "WARNING LadrunoTie -mortar - slave facet " << fs << " (first node "
                       << sfn(fs*npsS) << ") sits a mean |gap| " << gbar << " off the master "
                          "surface (tol " << tolFrac * facetScale << "). The tie requires "
                          "conforming interface geometry (coincident surfaces). Move the slave "
                          "onto the interface, or relax -tol.\n";
                return -1;
            }
        }
    }
    if (totalArea <= 1e-300) {
        opserr << "WARNING LadrunoTie -mortar - the slave and master surfaces do not overlap "
                  "(zero integrated area). Check the facets are coincident and the nps/winding correct.\n";
        return -1;
    }
    if (nExcluded > 0)
        opserr << "LadrunoTie -mortar - NOTE: " << nExcluded << " overlapping facet pair(s) "
                  "separated by more than the conforming tolerance were excluded (the far side "
                  "of a closed or folded surface).\n";

    // --- lumped mass on every tied slave DOF, only when the active handler requires it ---
    if (ltTiedMassRequired()) {
        std::map<int, std::map<int, double> > massedDOF;
        ltScanMassedDOFs(dom, massedDOF);
        for (int I = 0; I < Ns; I++) {
            Node *nd = dom->getNode(sTag[I]);
            if (nd == 0) { opserr << "WARNING LadrunoTie -mortar - slave node " << sTag[I] << " not in domain\n"; return -1; }
            if (!ltCheckTiedDofMass(sTag[I], nd, dofs, massedDOF, " -mortar")) {
                if (npsS > 4)
                    opserr << "HINT LadrunoTie -mortar - serendipity (hex20/quad8) slave surfaces "
                              "need HRZ lumped mass (Hinton, Rock and Zienkiewicz 1976): row-sum "
                              "lumping gives zero or NEGATIVE corner masses.\n";
                return -1;
            }
        }
    }

    // --- condense P (= D^{-1} M) ------------------------------------------------
    Matrix P(Ns, Nm);
    P.Zero();
    if (!dual) {
        // STANDARD basis: one global dense solve. P is dense (every slave couples to
        // every master through D^{-1}).
        if (D.Solve(M, P) < 0) {
            opserr << "WARNING LadrunoTie -mortar - failed to solve D P = M (the slave-interface mass D is "
                      "singular; an uncovered slave node likely passed the coverage check).\n";
            return -1;
        }
    } else {
        // DUAL (biorthogonal) basis (Wohlmuth 2000): replace the slave TEST functions
        // N_I^s with a dual basis psi_I such that INT psi_I N_J^s = delta_IJ INT N_I, so
        // D_dual is DIAGONAL and P = D_dual^{-1} M_dual is LOCAL (slave I ties only to the
        // masters under its own facet support). psi = A^e N is a per-slave-FACET linear
        // transform with
        //        A^e = diag(c^e) (D^e)^{-1},  c^e_a = sum_b D^e_ab = INT_e N_a  (facet row sum),
        // so D^e_dual = A^e D^e = diag(c^e) EXACTLY, for any (even distorted) facet. Built
        // from the same integratePair D^e/M^e in a second, setup-time integration pass.
        std::vector<double> Ddual(Ns, 0.0);          // diagonal slave mass
        Matrix Mdual(Ns, Nm);
        Mdual.Zero();
        for (int fs = 0; fs < nfS; fs++) {
            double Xs[8][3];
            for (int i = 0; i < npsS; i++) {
                double x[3];
                if (!ltNodeCoords3(dom, sfn(fs*npsS+i), x)) return -1;   // (already validated above)
                for (int d = 0; d < 3; d++) Xs[i][d] = x[d];
            }
            // this facet's D^e (npsS x npsS) and M^e (npsS x Nm) over the master facets it overlaps
            Matrix De(npsS, npsS), Me(npsS, Nm);
            De.Zero(); Me.Zero();
            for (int fm = 0; fm < nfM; fm++) {
                double Xm[8][3];
                for (int i = 0; i < npsM; i++)
                    for (int d = 0; d < 3; d++) Xm[i][d] = segM[((size_t)fm*npsM+i)*3+d];
                LadrunoMortarKernel::PairResult pr;
                if (LadrunoMortarKernel::integratePair(npsS, Xs, npsM, Xm, refDir, pr, gpOrder) != 0)
                    continue;                          // (-2 already refused in the first pass)
                if (pairGapTol[fs] >= 0.0 && pr.gapL1 > pairGapTol[fs] * pr.area)
                    continue;                          // excluded far pair (same test as above)
                for (int a = 0; a < npsS; a++) {
                    for (int b = 0; b < npsS; b++) De(a,b) += pr.D[a][b];
                    for (int k = 0; k < npsM; k++) Me(a, mIdx[mfn(fm*npsM+k)]) += pr.M[a][k];
                }
            }
            // c^e = row sums of D^e; skip a facet with no overlap at all. For a quad8 slave
            // the CORNER row sums are NEGATIVE (-A/12); the sign cancels in P = Mdual/Ddual.
            double cRow[8] = {0,0,0,0,0,0,0,0};
            double deNorm = 0.0;
            for (int a = 0; a < npsS; a++)
                for (int b = 0; b < npsS; b++) { cRow[a] += De(a,b); deNorm += std::fabs(De(a,b)); }
            if (deNorm < 1e-300) continue;
            // Y = (D^e)^{-1} M^e ; then M^e_dual = diag(c^e) Y ; scatter into Mdual, Ddual.
            Matrix Y(npsS, Nm);
            if (De.Solve(Me, Y) < 0) {
                opserr << "WARNING LadrunoTie -mortar -dual - a slave facet mass D^e is singular "
                          "(near-degenerate facet). Check the slave facet geometry.\n";
                return -1;
            }
            for (int a = 0; a < npsS; a++) {
                int I = sIdx[sfn(fs*npsS+a)];
                Ddual[I] += cRow[a];
                for (int k = 0; k < Nm; k++) Mdual(I,k) += cRow[a] * Y(a,k);
            }
        }
        // Quad8 corner dual masses are negative by construction; only a (near-)ZERO dual
        // mass means an uncovered or cancellation-degenerate node. The threshold is
        // RELATIVE to the largest dual mass: the partition-of-unity check below cannot
        // catch a near-cancelled mass (the dual row sums are identically 1).
        double maxAbsDdual = 0.0;
        for (int I = 0; I < Ns; I++)
            if (std::fabs(Ddual[I]) > maxAbsDdual) maxAbsDdual = std::fabs(Ddual[I]);
        for (int I = 0; I < Ns; I++) {
            if (std::fabs(Ddual[I]) <= 1e-12 * maxAbsDdual) {
                opserr << "WARNING LadrunoTie -mortar -dual - slave node " << sTag[I]
                       << " has (near-)zero dual mass " << Ddual[I] << " (uncovered, or a "
                          "corner/midside cancellation). Extend the master surface or check "
                          "the slave facet geometry.\n";
                return -1;
            }
            for (int k = 0; k < Nm; k++) P(I,k) = Mdual(I,k) / Ddual[I];
        }
    }
    // Every transfer row must satisfy partition of unity (sum_k P_Ik = 1), which is exact
    // for P = D^{-1}M (D 1 = M 1). A row that drifts from 1 means an ill-conditioned solve
    // (a near-singular D accepted by the factorization): refuse rather than emit a tie
    // that does not transfer a constant field.
    for (int I = 0; I < Ns; I++) {
        double rs = 0.0;
        for (int k = 0; k < Nm; k++) rs += P(I, k);
        if (std::fabs(rs - 1.0) > 1e-6) {
            opserr << "WARNING LadrunoTie -mortar - the mortar transfer for slave node " << sTag[I]
                   << " violates partition of unity (Sum P = " << rs << " != 1): the interface matrix D "
                      "is ill-conditioned. Check for sliver overlaps / near-degenerate slave facets.\n";
            return -1;
        }
    }

    // --- emit one EQ_Constraint per slave node per tied DOF: u_s = sum_k P_sk u_{m,k} ---
    int emitted = 0;
    for (int I = 0; I < Ns; I++) {
        int s = sTag[I];
        Node *sNode = dom->getNode(s);
        if (sNode == 0) { opserr << "WARNING LadrunoTie -mortar - slave node " << s << " not in domain\n"; return -1; }

        // master columns with a non-negligible weight (drop true zeros only; P is dense)
        int cnt = 0;
        for (int k = 0; k < Nm; k++) if (std::fabs(P(I,k)) > 1e-12) cnt++;
        if (cnt == 0) {
            opserr << "WARNING LadrunoTie -mortar - slave node " << s << " produced an empty transfer row.\n";
            return -1;
        }

        for (int di = 0; di < dofs.Size(); di++) {
            int dof1 = dofs(di);
            if (dof1 < 1 || dof1 > sNode->getNumberDOF()) {
                opserr << "WARNING LadrunoTie -mortar - slave node " << s << " has no DOF " << dof1 << "\n";
                return -1;
            }
            // like-to-like only: every retained master node must carry this DOF
            for (int k = 0; k < Nm; k++)
                if (std::fabs(P(I,k)) > 1e-12) {
                    Node *mn = dom->getNode(mTag[k]);
                    if (mn == 0 || dof1 > mn->getNumberDOF()) {
                        opserr << "WARNING LadrunoTie -mortar - slave node " << s << " tied DOF " << dof1
                               << " references master node " << mTag[k] << " which has no DOF " << dof1
                               << " (ndf mismatch). The mortar mode ties like to like (same ndf). "
                                  "For a shell EDGE tied to a solid FACE use -shellSolid "
                                  "(collocation only; the SOLID face nodes are the slaves), or "
                                  "restrict -dof to the shared DOFs (e.g. -dof 3 1 2 3).\n";
                        return -1;
                    }
                }
            ID rNode(cnt), rDOF(cnt);
            Vector Ccr(cnt);
            int p = 0;
            for (int k = 0; k < Nm; k++)
                if (std::fabs(P(I,k)) > 1e-12) {
                    rNode(p) = mTag[k];
                    rDOF(p)  = dof1 - 1;
                    Ccr(p)   = P(I,k);                 // u_s = sum_k P_sk u_{m,k}  (Ccr = P row)
                    p++;
                }
            EQ_Constraint *eq = new EQ_Constraint(s, dof1 - 1, Ccr, rNode, rDOF);
            if (eq == 0 || dom->addEQ_Constraint(eq) == false) {
                opserr << "WARNING LadrunoTie -mortar - failed to add EQ_Constraint for slave node " << s
                       << " dof " << dof1 << "\n";
                if (eq != 0) delete eq;
                return -1;
            }
            emitted++;
        }
    }

    opserr << "LadrunoTie (mortar" << (dual ? ", dual" : "") << "): emitted " << emitted
           << " EQ_Constraint(s) tying " << Ns << " slave node(s) to " << Nm << " master node(s) over "
           << nfS << " slave facet(s). Enforce with constraints Lagrange or Penalty.\n";
    return emitted;
}

// --------------------------------------------------------------------------- //
// shell-to-solid plane-section generator (-shellSolid)
// --------------------------------------------------------------------------- //

// Closest-point projection of a slave onto the master edge POLYLINE (closed-form
// point-to-segment, brute force over the segments; edge lists are short).
// seg = flat reference coords, 6 doubles per segment (a.xyz b.xyz). In-bounds
// means the UNCLAMPED parameter lies in [0,1] within a small parametric tol on
// at least one segment; among those the closest wins. Returns false if the slave
// projects outside the polyline everywhere (a named refusal at the call site).
static bool
ltProjectEdge(int nseg, const double *seg, const double xs[3],
              int &bestSeg, double &bestXi, double foot[3])
{
    const double xitol = 1e-9;
    bestSeg = -1;
    double bestDist = 0.0;
    for (int k = 0; k < nseg; k++) {
        const double *a = seg + (size_t)k * 6;
        const double *b = a + 3;
        double ab[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
        double L2 = ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2];
        if (L2 <= 0.0) continue;                 // degenerate (refused at setup)
        double xi = ((xs[0]-a[0])*ab[0] + (xs[1]-a[1])*ab[1] + (xs[2]-a[2])*ab[2]) / L2;
        if (xi < -xitol || xi > 1.0 + xitol) continue;
        double xic = (xi < 0.0) ? 0.0 : ((xi > 1.0) ? 1.0 : xi);
        double f[3] = { a[0]+xic*ab[0], a[1]+xic*ab[1], a[2]+xic*ab[2] };
        double dd = std::sqrt((xs[0]-f[0])*(xs[0]-f[0]) + (xs[1]-f[1])*(xs[1]-f[1]) +
                              (xs[2]-f[2])*(xs[2]-f[2]));
        if (bestSeg < 0 || dd < bestDist) {
            bestSeg = k; bestDist = dd; bestXi = xic;
            foot[0] = f[0]; foot[1] = f[1]; foot[2] = f[2];
        }
    }
    return bestSeg >= 0;
}

int
LadrunoTie::generateShellSolid(Domain *dom,
                               const ID &slaves,
                               const ID &edgeNodes,
                               double thickness, double tolFrac)
{
    if (dom == 0) {
        opserr << "WARNING LadrunoTie -shellSolid - domain is not defined\n";
        return -1;
    }
    if (slaves.Size() == 0) {
        opserr << "WARNING LadrunoTie -shellSolid - no slave (solid face) nodes given\n";
        return -1;
    }
    if (edgeNodes.Size() == 0 || (edgeNodes.Size() % 2) != 0) {
        opserr << "WARNING LadrunoTie -shellSolid - master edge node count " << edgeNodes.Size()
               << " is not a positive multiple of 2 (segment node pairs a1 b1 a2 b2 ..)\n";
        return -1;
    }
    const int nseg = edgeNodes.Size() / 2;

    // --- master edge segments: reference coords, ndf-6 check, node-disjoint set ---
    std::vector<double> seg((size_t)nseg * 6, 0.0);
    std::set<int> masterSet;
    for (int k = 0; k < nseg; k++) {
        for (int i = 0; i < 2; i++) {
            int tag = edgeNodes(k * 2 + i);
            Node *mn = dom->getNode(tag);
            double x[3];
            if (mn == 0 || !ltNodeCoords3(dom, tag, x)) {
                opserr << "WARNING LadrunoTie -shellSolid - master edge node " << tag
                       << " not in domain\n";
                return -1;
            }
            if (mn->getNumberDOF() < 6) {
                opserr << "WARNING LadrunoTie -shellSolid - master edge node " << tag
                       << " has ndf " << mn->getNumberDOF() << " (< 6). The plane-section "
                          "rows put moments on the shell edge; the master edge must be "
                          "ndf-6 shell nodes (the SOLID face nodes are the slaves).\n";
                return -1;
            }
            masterSet.insert(tag);
            for (int d = 0; d < 3; d++) seg[((size_t)k * 2 + i) * 3 + d] = x[d];
        }
        const double *a = seg.data() + (size_t)k * 6;
        const double *b = a + 3;
        double L2 = (b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) + (b[2]-a[2])*(b[2]-a[2]);
        if (L2 <= 1e-300) {
            opserr << "WARNING LadrunoTie -shellSolid - degenerate (zero-length) master edge "
                      "segment " << edgeNodes(k*2) << "-" << edgeNodes(k*2+1) << "\n";
            return -1;
        }
    }

    // without -thickness the arm check is skipped (with a note)
    if (thickness <= 0.0)
        opserr << "LadrunoTie -shellSolid - NOTE: -thickness not given; the arm check "
                  "|d| <= (0.5+tol)*t is SKIPPED. A slave far off the shell mid-surface "
                  "would get a spurious long lever arm; pass -thickness t to catch that.\n";

    // --- tied DOFs are the solid translations; element mass only when the active
    //     handler requires it (a solid with a density satisfies it). ---
    const bool massRequired = ltTiedMassRequired();
    std::map<int, std::map<int, double> > massedDOF;
    if (massRequired)
        ltScanMassedDOFs(dom, massedDOF);
    ID dofs(3);
    for (int d = 0; d < 3; d++) dofs(d) = d + 1;

    std::set<int> doneSlaves;
    int emitted = 0;

    for (int si = 0; si < slaves.Size(); si++) {
        int s = slaves(si);

        // a slave listed twice would be doubly constrained
        if (doneSlaves.count(s)) {
            opserr << "WARNING LadrunoTie -shellSolid - slave node " << s << " is listed more "
                      "than once. Each solid face node gets exactly one plane-section row set. "
                      "Remove the duplicate.\n";
            return -1;
        }
        doneSlaves.insert(s);

        // slave and master edge must be node-disjoint (no chained constraints)
        if (masterSet.count(s)) {
            opserr << "WARNING LadrunoTie -shellSolid - node " << s << " is BOTH a slave and a "
                      "master edge node. The solid face and the shell edge must be "
                      "node-disjoint (no chained constraints).\n";
            return -1;
        }

        Node *sNode = dom->getNode(s);
        double xs[3];
        if (sNode == 0 || !ltNodeCoords3(dom, s, xs)) {
            opserr << "WARNING LadrunoTie -shellSolid - slave node " << s << " not in domain\n";
            return -1;
        }
        if (sNode->getNumberDOF() < 3) {
            opserr << "WARNING LadrunoTie -shellSolid - slave node " << s << " has ndf "
                   << sNode->getNumberDOF() << " (< 3); the tie constrains the 3 solid "
                      "translations.\n";
            return -1;
        }

        if (massRequired && !ltCheckTiedDofMass(s, sNode, dofs, massedDOF, " -shellSolid"))
            return -1;

        // --- closest-point edge projection; FROZEN arm d = x_s - xbar(xi) ---
        int bestSeg_;
        double xi, foot[3];
        if (!ltProjectEdge(nseg, seg.data(), xs, bestSeg_, xi, foot)) {
            opserr << "WARNING LadrunoTie -shellSolid - slave node " << s << " does not project "
                      "in-bounds onto any master edge segment (parametric xi outside [0,1] "
                      "everywhere). The shell edge polyline must span the solid face; extend "
                      "-masterEdge or check the node lists.\n";
            return -1;
        }
        double dv[3] = { xs[0]-foot[0], xs[1]-foot[1], xs[2]-foot[2] };
        double dlen = std::sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);

        // arm check: a slave farther than half the thickness (+slack) off the shell
        // mid-surface is not this shell's through-thickness material. The slack floors
        // at 2% of the thickness: the projection tolerance (default 1e-6) would refuse a
        // legitimate outer-surface node at |d| = 0.5t nudged a few ppm past it by
        // rounding or curvature. A larger -tol widens it.
        double armSlack = (tolFrac > 0.02) ? tolFrac : 0.02;
        if (thickness > 0.0 && dlen > (0.5 + armSlack) * thickness) {
            opserr << "WARNING LadrunoTie -shellSolid - slave node " << s << " sits " << dlen
                   << " off the master edge, beyond (0.5+slack)*t = "
                   << (0.5 + armSlack) * thickness << " (t = " << thickness << "). It is not "
                      "this shell's through-thickness material and would get a spurious long "
                      "lever arm. Trim the slave list, fix -thickness, or relax -tol.\n";
            return -1;
        }

        // --- the 3 plane-section rows: N_j on translations, N_j*(-[d]x) on rotations.
        //     -[d]x annihilates its own axis, so the drilling (theta.n ~ theta.d) column
        //     is ~zero and the 1e-12 filter in ltEmitMixedRow drops it: the shell
        //     drilling DOF stays free (a |d| = 0 slave gets identity translation rows). ---
        int mTags[2] = { edgeNodes(bestSeg_ * 2), edgeNodes(bestSeg_ * 2 + 1) };
        double Nw[2] = { 1.0 - xi, xi };
        double B[3][3] = { { 0.0,    dv[2], -dv[1] },     // -skew(d) = d(theta x d)/d(theta)
                           { -dv[2], 0.0,    dv[0] },
                           { dv[1], -dv[0],  0.0   } };
        for (int i = 0; i < 3; i++) {                     // slave translation row u_s,i
            double coef[2][6];
            for (int k = 0; k < 2; k++)
                for (int j = 0; j < 3; j++) {
                    coef[k][j]     = (i == j) ? Nw[k] : 0.0;
                    coef[k][3 + j] = Nw[k] * B[i][j];
                }
            if (!ltEmitMixedRow(dom, s, i, 2, mTags, coef, " -shellSolid")) return -1;
            emitted++;
        }
    }

    opserr << "LadrunoTie (shellSolid): emitted " << emitted << " EQ_Constraint(s) tying "
           << (int)doneSlaves.size() << " solid face node(s) to the shell edge ("
           << nseg << " segment(s), " << (int)masterSet.size() << " node(s)). Enforce with "
              "constraints Lagrange or Penalty.\n";
    return emitted;
}

// --------------------------------------------------------------------------- //
// command front-end
//
//  collocation (default):
//   LadrunoTie -slaveNodes <ns> s1.. -masterFacets <npsM> <nfM> m..
//              [-dof <nd> d1..] [-tol <frac>] [-hermite]
//   (-hermite: rotation-consistent cubic Hermite w-theta transfer for ndf-6 shell
//    EDGE ties; requires the full 6-DOF tie.)
//
//  integral mortar:
//   LadrunoTie -mortar -slaveFacets <npsS> <nfS> s.. -masterFacets <npsM> <nfM> m..
//              [-dof <nd> d1..] [-tol <frac>] [-outward ox oy oz] [-dual]
//   (-dual: biorthogonal basis, sparse P. npsS in {3,4,8}, npsM in {3,4,6,8};
//    quadratic facets corners first, flat, midside nodes at the edge midpoints.)
//
//  shell-to-solid:
//   LadrunoTie -shellSolid -slaveNodes <ns> s1.. -masterEdge <nseg> a1 b1 a2 b2 ..
//              [-thickness <t>] [-tol <frac>]
//   (solid FACE nodes = slaves, shell EDGE polyline = master; 3 plane-section
//    rows per slave, u_s = sum N_j (u_j + theta_j x d); tied DOFs {1,2,3};
//    -mortar/-hermite/-dual/-dof/-outward are refused.)
//
// Flat token streams (as equationConstraint): Tcl and Python both pass node tags
// as separate arguments.
// --------------------------------------------------------------------------- //

// read "<n> v1..vn" into id; returns false on a parse/count error.
static bool ltReadCountedID(ID &id, const char *what)
{
    int one = 1, n = 0;
    if (OPS_GetIntInput(&one, &n) < 0 || n < 1) {
        opserr << "WARNING LadrunoTie " << what << " - need a positive count\n";
        return false;
    }
    id.resize(n);
    for (int i = 0; i < n; i++) {
        int v;
        if (OPS_GetIntInput(&one, &v) < 0) {
            opserr << "WARNING LadrunoTie " << what << " - could not read entry " << i + 1 << "\n";
            return false;
        }
        id(i) = v;
    }
    return true;
}

// read "<nseg> (2*nseg tags)" into id (edge segment node pairs); false on error.
static bool ltReadEdge(ID &id, const char *what)
{
    int one = 1, nseg = 0;
    if (OPS_GetIntInput(&one, &nseg) < 0 || nseg < 1) {
        opserr << "WARNING LadrunoTie " << what << " - need a positive segment count\n";
        return false;
    }
    id.resize(nseg * 2);
    for (int i = 0; i < nseg * 2; i++) {
        int v;
        if (OPS_GetIntInput(&one, &v) < 0) {
            opserr << "WARNING LadrunoTie " << what << " - could not read edge node " << i + 1
                   << " of " << nseg * 2 << "\n";
            return false;
        }
        id(i) = v;
    }
    return true;
}

// read "<nps> <nf> (nps*nf tags)" into id; sets nps; returns false on error.
static bool ltReadFacets(ID &id, int &nps, const char *what)
{
    int one = 1, nf = 0;
    if (OPS_GetIntInput(&one, &nps) < 0 || (nps != 3 && nps != 4 && nps != 6 && nps != 8)) {
        opserr << "WARNING LadrunoTie " << what << " - need nps (3, 4, 6, or 8; quadratic "
                  "facets are -mortar only)\n";
        return false;
    }
    if (OPS_GetIntInput(&one, &nf) < 0 || nf < 1) {
        opserr << "WARNING LadrunoTie " << what << " - need a positive facet count\n";
        return false;
    }
    id.resize(nf * nps);
    for (int i = 0; i < nf * nps; i++) {
        int v;
        if (OPS_GetIntInput(&one, &v) < 0) {
            opserr << "WARNING LadrunoTie " << what << " - could not read facet node " << i + 1
                   << " of " << nf * nps << "\n";
            return false;
        }
        id(i) = v;
    }
    return true;
}

int OPS_LadrunoTie()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) {
        opserr << "WARNING LadrunoTie - domain is not defined\n";
        return -1;
    }
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "WARNING want - LadrunoTie -slaveNodes ns s1.. -masterFacets npsM nf m.. "
                  "<-dof nd d1..> <-tol frac> <-hermite>   |   LadrunoTie -mortar -slaveFacets npsS nf s.. "
                  "-masterFacets npsM nf m.. <-dof..> <-tol..> <-outward ox oy oz> <-dual>   |   "
                  "LadrunoTie -shellSolid -slaveNodes ns s1.. -masterEdge nseg a1 b1 .. "
                  "<-thickness t> <-tol frac>\n";
        return -1;
    }

    ID slaveNodes(0), slaveFacetNodes(0), masterFacetNodes(0), masterEdgeNodes(0), dofs(0);
    int npsM = 0, npsS = 0;
    double tolFrac = 1.0e-6;
    double outward[3] = {0.0, 0.0, 0.0};
    double thickness = -1.0;
    bool mortar = false, haveSlaveNodes = false, haveSlaveFacets = false;
    bool haveMaster = false, haveOutward = false, dual = false, hermite = false;
    bool shellSolid = false, haveEdge = false, haveThickness = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-mortar") == 0) {
            mortar = true;
        }
        else if (strcmp(opt, "-dual") == 0) {   // biorthogonal basis, sparse P (mortar only)
            dual = true;
        }
        else if (strcmp(opt, "-hermite") == 0) { // Hermite w-theta shell-edge transfer (collocation only)
            hermite = true;
        }
        else if (strcmp(opt, "-shellSolid") == 0) { // shell EDGE <- solid FACE plane-section tie
            shellSolid = true;
        }
        else if (strcmp(opt, "-masterEdge") == 0) {
            if (!ltReadEdge(masterEdgeNodes, "-masterEdge")) return -1;
            haveEdge = true;
        }
        else if (strcmp(opt, "-thickness") == 0) {
            int one = 1;
            if (OPS_GetDoubleInput(&one, &thickness) < 0 || thickness <= 0.0) {
                opserr << "WARNING LadrunoTie -thickness - need a positive shell thickness\n";
                return -1;
            }
            haveThickness = true;
        }
        else if (strcmp(opt, "-slaveNodes") == 0 || strcmp(opt, "-slave") == 0) {
            if (!ltReadCountedID(slaveNodes, "-slaveNodes")) return -1;
            haveSlaveNodes = true;
        }
        else if (strcmp(opt, "-slaveFacets") == 0) {
            if (!ltReadFacets(slaveFacetNodes, npsS, "-slaveFacets")) return -1;
            haveSlaveFacets = true;
        }
        else if (strcmp(opt, "-masterFacets") == 0 || strcmp(opt, "-master") == 0) {
            if (!ltReadFacets(masterFacetNodes, npsM, "-masterFacets")) return -1;
            haveMaster = true;
        }
        else if (strcmp(opt, "-dof") == 0) {
            if (!ltReadCountedID(dofs, "-dof")) return -1;
        }
        else if (strcmp(opt, "-tol") == 0) {
            int one = 1;
            if (OPS_GetDoubleInput(&one, &tolFrac) < 0 || tolFrac < 0.0) {
                opserr << "WARNING LadrunoTie -tol - need a non-negative number\n";
                return -1;
            }
        }
        else if (strcmp(opt, "-outward") == 0) {
            int three = 3;
            if (OPS_GetDoubleInput(&three, outward) < 0) {
                opserr << "WARNING LadrunoTie -outward - need ox oy oz\n";
                return -1;
            }
            haveOutward = true;
        }
        else {
            opserr << "WARNING LadrunoTie - unknown option " << opt << "\n";
            return -1;
        }
    }

    // -masterEdge / -thickness belong to -shellSolid only
    if (!shellSolid && (haveEdge || haveThickness)) {
        opserr << "WARNING LadrunoTie - -masterEdge/-thickness apply only to the -shellSolid "
                  "mode (the shell EDGE <- solid FACE plane-section tie). Add -shellSolid, or "
                  "use -masterFacets for the facet-based modes.\n";
        return -1;
    }

    int rc;
    if (shellSolid) {
        // collocation only, and the tied-DOF set is fixed by the operator (solid
        // translations {1,2,3}): named refusals for the other options
        if (mortar || dual) {
            opserr << "WARNING LadrunoTie -shellSolid - the plane-section shell-to-solid tie is "
                      "COLLOCATION only (a weak form would need the plane-section basis inside "
                      "the mortar integral). Drop -mortar/-dual.\n";
            return -1;
        }
        if (hermite) {
            opserr << "WARNING LadrunoTie -shellSolid - -hermite is the shell-to-SHELL edge "
                      "transfer; -shellSolid has its own plane-section transfer. Drop one.\n";
            return -1;
        }
        if (dofs.Size() > 0) {
            opserr << "WARNING LadrunoTie -shellSolid - the tied DOFs are FIXED to the solid "
                      "slave translations {1,2,3} by the plane-section operator (shell drilling "
                      "stays free); -dof is not selectable here.\n";
            return -1;
        }
        if (haveMaster || haveSlaveFacets) {
            opserr << "WARNING LadrunoTie -shellSolid - the master is the shell EDGE polyline "
                      "(-masterEdge nseg a1 b1 ..) and the slaves are the solid FACE nodes "
                      "(-slaveNodes); -masterFacets/-slaveFacets do not apply.\n";
            return -1;
        }
        if (haveOutward) {
            opserr << "WARNING LadrunoTie -shellSolid - -outward orients the mortar surface "
                      "normal and does not apply here; the plane-section arm direction comes "
                      "from each slave's closest-point projection onto the shell edge. Drop "
                      "-outward.\n";
            return -1;
        }
        if (!haveSlaveNodes || !haveEdge) {
            opserr << "WARNING LadrunoTie -shellSolid - both -slaveNodes (solid face nodes) and "
                      "-masterEdge (shell edge segment pairs) are required\n";
            return -1;
        }
        rc = LadrunoTie::generateShellSolid(theDomain, slaveNodes, masterEdgeNodes,
                                            haveThickness ? thickness : -1.0, tolFrac);
    } else if (mortar) {
        if (hermite) {
            opserr << "WARNING LadrunoTie -hermite - the Hermite w-theta transfer applies only to "
                      "the collocation mode (a weak-form Hermite tie would need the Hermite basis "
                      "inside the mortar integral). Drop -mortar.\n";
            return -1;
        }
        if (!haveSlaveFacets || !haveMaster) {
            opserr << "WARNING LadrunoTie -mortar - both -slaveFacets and -masterFacets are required\n";
            return -1;
        }
        rc = LadrunoTie::generateMortar(theDomain, slaveFacetNodes, npsS, masterFacetNodes, npsM,
                                        dofs, tolFrac, haveOutward ? outward : 0, dual);
    } else {
        if (dual) {
            opserr << "WARNING LadrunoTie -dual - the dual/biorthogonal basis applies only to the "
                      "integral mortar mode; add -mortar (collocation has no D to diagonalize).\n";
            return -1;
        }
        if (haveSlaveFacets) {   // refuse rather than silently drop the facets
            opserr << "WARNING LadrunoTie - -slaveFacets belongs to -mortar (collocation ties "
                      "NODES: -slaveNodes). Add -mortar, or list the nodes instead.\n";
            return -1;
        }
        if (!haveSlaveNodes || !haveMaster) {
            opserr << "WARNING LadrunoTie - both -slaveNodes and -masterFacets are required "
                      "(or -mortar with -slaveFacets)\n";
            return -1;
        }
        rc = LadrunoTie::generate(theDomain, slaveNodes, masterFacetNodes, npsM, dofs, tolFrac,
                                  hermite);
    }
    return (rc < 0) ? -1 : 0;
}
