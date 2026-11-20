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
//   Wriggers, P. (2006). Computational Contact Mechanics, 2nd ed. Springer,
//     Berlin. doi:10.1007/978-3-540-32609-0 (closest-point projection,
//     normal gap and its linearization for node-to-segment contact).
//   Laursen, T.A. (2002). Computational Contact and Impact Mechanics.
//     Springer, Berlin. doi:10.1007/978-3-662-04864-1
//   Wriggers, P. and Simo, J.C. (1985). A note on tangent stiffness for fully
//     nonlinear contact problems. Communications in Applied Numerical
//     Methods, 1(5), 199-203. doi:10.1002/cnm.1630010503
//
// Description: LadrunoContactProjection, header-only closest-point projection
// geometry for node-to-segment (NTS) and segment-to-segment (mortar) contact
// on 3-node triangle and 4-node bilinear quadrilateral master segments. It uses
// no OpenSees types (raw double arrays and <cmath> only), so the NTS adapter
// and the mortar Gauss-point loop share one projection.
//
//   evalSegment()     NTS one-shot: projection, oriented normal and gap; true
//                     only for an in-bounds, penetrating slave.
//   evalSegmentSmooth() as evalSegment(), with the normal interpolated from an
//                     averaged nodal-normal field (LadrunoContactNormalField).
//   projectFull()     mortar: the geometry of any converged projection (gap as
//                     data, no penetration or in-bounds test), with the
//                     covariant tangents and the surface metric g_ab = g_a.g_b.
//   normalGeomTangent() the geometric part of the consistent NTS normal tangent.
//
// The projection is a bounded closest-point Newton iteration (maxIt iterations;
// a segment with |det K| < 1e-14 |g1|^2 |g2|^2 is degenerate). The outward
// normal is derived from the segment geometry and oriented by a reference
// direction, never taken from the node winding, so reversing the winding gives
// the same normal.

#ifndef LadrunoContactProjection_h
#define LadrunoContactProjection_h

#include <cmath>
#include "LadrunoContactNormalField.h"   // smoothNormal() blend for evalSegmentSmooth()

namespace LadrunoContactProjection {

// ----------------------------------------------------------------- small vec3
inline double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline void cross3(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}
inline double norm3(const double a[3]) { return std::sqrt(dot3(a, a)); }

// ------------------------------------------------------- shape fns + derivs
// nps = 3 (linear tri, (xi,eta) area coords) or 4 (bilinear quad on [-1,1]^2).
// Fills N[nps], dNxi[nps], dNeta[nps].
inline void shape(int nps, double xi, double eta,
                  double N[4], double dNxi[4], double dNeta[4]) {
    if (nps == 4) {
        N[0] = 0.25*(1-xi)*(1-eta);  N[1] = 0.25*(1+xi)*(1-eta);
        N[2] = 0.25*(1+xi)*(1+eta);  N[3] = 0.25*(1-xi)*(1+eta);
        dNxi[0] = -0.25*(1-eta);  dNxi[1] =  0.25*(1-eta);
        dNxi[2] =  0.25*(1+eta);  dNxi[3] = -0.25*(1+eta);
        dNeta[0] = -0.25*(1-xi);  dNeta[1] = -0.25*(1+xi);
        dNeta[2] =  0.25*(1+xi);  dNeta[3] =  0.25*(1-xi);
    } else { // tri-3
        N[0] = 1 - xi - eta;  N[1] = xi;  N[2] = eta;  N[3] = 0.0;
        dNxi[0] = -1.0;  dNxi[1] = 1.0;  dNxi[2] = 0.0;  dNxi[3] = 0.0;
        dNeta[0] = -1.0; dNeta[1] = 0.0; dNeta[2] = 1.0; dNeta[3] = 0.0;
    }
}

// xbar = sum N_i X_i ; tangents g1 = d xbar/d xi, g2 = d xbar/d eta.  X is nps x 3.
inline void interp(int nps, const double N[4], const double X[4][3], double xbar[3]) {
    for (int d = 0; d < 3; d++) {
        xbar[d] = 0.0;
        for (int i = 0; i < nps; i++) xbar[d] += N[i] * X[i][d];
    }
}
inline void tangents(int nps, const double dNxi[4], const double dNeta[4],
                     const double X[4][3], double g1[3], double g2[3]) {
    for (int d = 0; d < 3; d++) {
        g1[d] = 0.0; g2[d] = 0.0;
        for (int i = 0; i < nps; i++) { g1[d] += dNxi[i]*X[i][d]; g2[d] += dNeta[i]*X[i][d]; }
    }
}

// in-bounds test on the parent domain (small tolerance).
inline bool inBounds(int nps, double xi, double eta) {
    const double t = 1e-9;
    if (nps == 4) return (xi >= -1 - t && xi <= 1 + t && eta >= -1 - t && eta <= 1 + t);
    return (xi >= -t && eta >= -t && xi + eta <= 1 + t);
}

// ----------------------------------------------------- bounded projection
// Closest-point Newton iteration on (x_s - xbar).g_alpha = 0, capped at maxIt
// iterations. The Gauss-Newton metric K = [[g1.g1, g1.g2], [g2.g1, g2.g2]]
// drops the O(gap) curvature correction and is SPD for a non-degenerate
// segment. The segment is rejected as degenerate when |det K| < 1e-14 K00 K11
// (sin^2 of the angle between g1 and g2 below 1e-14; dimensionless, so
// independent of the coordinate scale).
//
// Two convergence tests, either suffices:
//   (a) |R| < tolR, an absolute residual test. R = d.g_alpha has units of
//       length^2 and a round-off floor of about eps |X| |g|: for coordinates
//       far from the origin (x = 5e4 with 500-unit facets the floor is about
//       3e-10) an absolute 1e-12 cannot be reached.
//   (b) |dxi| + |deta| < tolP after the update, a scale-free test on the
//       parent coordinates (which are O(1)); 1e-8 bounds the footpoint error
//       by about 1e-8 h. On a flat facet the residual collapses in one step and
//       (b) exits with the same (xi, eta) one iteration earlier; on a warped
//       facet Gauss-Newton converges linearly and the footpoint differs from
//       the residual-converged one by about tolP in parent coordinates.
// returns: 0 converged in-bounds, 1 converged out-of-bounds, -1 no valid projection.
inline int project(int nps, const double X[4][3], const double xs[3],
                   double &xi, double &eta, double tolR, int maxIt,
                   double tolP = 1e-8) {
    xi = (nps == 4) ? 0.0 : (1.0/3.0);
    eta = (nps == 4) ? 0.0 : (1.0/3.0);
    double N[4], dNxi[4], dNeta[4], xbar[3], g1[3], g2[3];
    bool converged = false;
    for (int it = 0; it < maxIt; it++) {
        shape(nps, xi, eta, N, dNxi, dNeta);
        interp(nps, N, X, xbar);
        tangents(nps, dNxi, dNeta, X, g1, g2);
        double d[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
        double R0 = dot3(d, g1), R1 = dot3(d, g2);
        if (std::sqrt(R0*R0 + R1*R1) < tolR) { converged = true; break; }
        double K00 = dot3(g1,g1), K01 = dot3(g1,g2), K11 = dot3(g2,g2);
        double detK = K00*K11 - K01*K01;
        if (std::fabs(detK) < 1e-14 * (K00*K11 + 1e-300)) return -1;  // degenerate (angle test)
        double dxi  = ( K11*R0 - K01*R1) / detK;
        double deta = (-K01*R0 + K00*R1) / detK;
        xi += dxi; eta += deta;
        if (std::fabs(dxi) + std::fabs(deta) < tolP) { converged = true; break; }
    }
    if (!converged) return -1;
    return inBounds(nps, xi, eta) ? 0 : 1;
}

// Clamp (xi,eta) to the parent facet (quad [-1,1]^2; tri {xi>=0, eta>=0,
// xi+eta<=1}), so an out-of-bounds Newton result maps to a point on the facet.
// Used by voteSignRobust() to get a footpoint for a slave that projects past a
// facet edge; a sign vote needs a footpoint on the right facet, not the exact
// boundary-closest point.
inline void clampParam(int nps, double &xi, double &eta) {
    if (nps == 4) {
        if (xi < -1.0) xi = -1.0; else if (xi > 1.0) xi = 1.0;
        if (eta < -1.0) eta = -1.0; else if (eta > 1.0) eta = 1.0;
    } else {
        if (xi < 0.0) xi = 0.0;
        if (eta < 0.0) eta = 0.0;
        double s = xi + eta;
        if (s > 1.0) {                        // project onto the hypotenuse xi+eta=1, then re-clamp >= 0
            double h = 0.5 * (s - 1.0);
            xi -= h; eta -= h;
            if (xi < 0.0) { eta += xi; xi = 0.0; }
            if (eta < 0.0) { xi += eta; eta = 0.0; }
        }
    }
}

// Global outward sign of a master surface by a distance-weighted majority vote
// of the slaves. For each slave: find its nearest master facet (closest-point
// projection clamped to the facet), take that facet's coherent (sigma applied)
// unit normal n, and vote
//   w = n . (slave - surfaceCentroid)
// where surfaceCentroid is the average of the facet node slots, an interior
// reference for an open convex patch. The sign is sign(sum w), each vote
// weighted by |w|. Returns the sign (+1/-1); *conf = |sum w| / sum |w| in [0,1]
// (1: all slaves agree; 0: an ambiguous two-sided cloud); *nVoted = number of
// slaves that voted (0: nothing projected, the caller falls back to
// LadrunoContactNormalField::voteSign).
//
// The local facet normal is used because an aggregate normal sum nearly
// cancels on a curved (domed) master, and when the slave cloud grazes the
// master edge-on the chord between centroids is nearly tangent to it, so the
// aggregate vote is ill-conditioned. The local normal at an edge-grazing slave
// carries the lateral component, so n.chord stays well conditioned.
//
// The surface centroid, rather than the local footpoint, is the reference
// because a slave placed slightly penetrating the surface in the reference
// configuration has a footpoint separation pointing inward; the interior
// centroid keeps slave - centroid pointing outward. A non-convex open patch
// whose centroid is not interior needs an explicit outward direction.
//
// segCoords and slaveCoords must be the same configuration; the reference
// configuration makes the vote independent of the deformation.
inline double voteSignRobust(int nps, int nSeg, const double *segCoords, const int *sigma,
                             const double *slaveCoords, int nSlaves,
                             double *conf, int *nVoted) {
    if (nVoted != 0) *nVoted = 0;
    long np = (long)nSeg * nps;
    if (np <= 0) { if (conf != 0) *conf = 0.0; return +1.0; }
    // interior reference: the surface centroid (average of the facet node slots)
    double cen[3] = {0.0, 0.0, 0.0};
    for (long i = 0; i < np; i++) for (int d = 0; d < 3; d++) cen[d] += segCoords[i*3 + d];
    for (int d = 0; d < 3; d++) cen[d] /= (double)np;

    double acc = 0.0, wsum = 0.0;
    int voted = 0;
    for (int p = 0; p < nSlaves; p++) {
        const double *xs = slaveCoords + (size_t)p * 3;
        // nearest facet by clamped closest-point distance
        double bestDist2 = 0.0; int bestSeg = -1;
        for (int s = 0; s < nSeg; s++) {
            double X[4][3];
            for (int k = 0; k < nps; k++)
                for (int d = 0; d < 3; d++) X[k][d] = segCoords[((size_t)s * nps + k) * 3 + d];
            double xi, eta;
            if (project(nps, X, xs, xi, eta, 1e-12, 10) < 0) continue;   // degenerate facet: skip
            clampParam(nps, xi, eta);
            double N[4], dNx[4], dNe[4], foot[3];
            shape(nps, xi, eta, N, dNx, dNe);
            interp(nps, N, X, foot);
            double dd[3] = { xs[0]-foot[0], xs[1]-foot[1], xs[2]-foot[2] };
            double dist2 = dot3(dd, dd);
            if (bestSeg < 0 || dist2 < bestDist2) { bestDist2 = dist2; bestSeg = s; }
        }
        if (bestSeg < 0) continue;                                       // no valid projection
        double fn[3];
        LadrunoContactNormalField::newellAreaNormal(nps, segCoords + (size_t)bestSeg * nps * 3, fn);
        double sgn = (double)sigma[bestSeg];
        double nf[3] = { sgn*fn[0], sgn*fn[1], sgn*fn[2] };              // coherent normal (G not applied)
        double nfn = norm3(nf);
        if (nfn < 1e-300) continue;                                     // sliver facet: no vote
        double d2c[3] = { xs[0]-cen[0], xs[1]-cen[1], xs[2]-cen[2] };   // chord from the centroid
        double w = dot3(nf, d2c) / nfn;   // signed distance of the slave from the centroid along n
        acc += w; wsum += (w >= 0.0 ? w : -w);
        voted++;
    }
    if (nVoted != 0) *nVoted = voted;
    if (voted == 0 || wsum < 1e-300) { if (conf != 0) *conf = 0.0; return +1.0; }
    if (conf != 0) *conf = (acc >= 0.0 ? acc : -acc) / wsum;
    return (acc < 0.0) ? -1.0 : +1.0;
}

// Derived outward normal, oriented so n.refDir > 0 (refDir points toward the
// slave's allowed half-space). Orienting by a fixed direction rather than the
// live slave position keeps the sign after the slave penetrates; reversing the
// node order flips the cross product and the refDir test flips it back.
// refDir need not be unit.
inline bool normalOriented(int nps, double xi, double eta, const double X[4][3],
                           const double refDir[3], double n[3]) {
    double N[4], dNxi[4], dNeta[4], xbar[3], g1[3], g2[3], raw[3];
    shape(nps, xi, eta, N, dNxi, dNeta);
    interp(nps, N, X, xbar);
    tangents(nps, dNxi, dNeta, X, g1, g2);
    cross3(g1, g2, raw);
    double j = norm3(raw);
    if (j < 1e-300) return false;
    for (int d = 0; d < 3; d++) n[d] = raw[d] / j;
    // If refDir is numerically perpendicular to n the outward sense is
    // ambiguous (e.g. the slave reference point lies in the segment plane):
    // refuse the pair rather than guess a sign; the caller should pass an
    // explicit -outward direction.
    double proj = dot3(n, refDir);
    if (std::fabs(proj) < 1e-12 * (norm3(refDir) + 1e-300)) return false;
    if (proj < 0.0) for (int d = 0; d < 3; d++) n[d] = -n[d];
    return true;
}

// One-shot segment evaluation: project xs onto the segment, derive the oriented
// normal (toward refDir) and the gap. Returns the active flag (penetrating and
// in-bounds). When active, fills gap (< 0), n[3] and N[nps] (shape functions at
// the projection).
inline bool evalSegment(int nps, const double X[4][3], const double xs[3],
                        const double refDir[3], double &gap, double n[3], double N[4],
                        double tolR = 1e-12, int maxIt = 10) {
    double xi, eta;
    int st = project(nps, X, xs, xi, eta, tolR, maxIt);
    if (st != 0) return false;                       // out of bounds or no valid projection
    if (!normalOriented(nps, xi, eta, X, refDir, n)) return false;
    double dN1[4], dN2[4], xbar[3];
    shape(nps, xi, eta, N, dN1, dN2);
    interp(nps, N, X, xbar);
    double d[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
    gap = dot3(n, d);
    return (gap < 0.0);                              // penetrating only
}

// Smoothed-normal one-shot segment evaluation. Identical to evalSegment()
// except that the normal is the averaged nodal-normal field interpolated at the
// projection, n = normalize(sum N_i n_i), instead of the faceted normal. The
// closest-point projection on the facet is unchanged; only the force direction
// changes. nodalNorm = the segment's nps nodal normals
// (LadrunoContactNormalField::nodalNormals). On a degenerate blend the faceted
// normalOriented(refDir) is used.
//
// Facet-ownership guard (sharedEdge != 0): at a sharp convex ridge a slave
// clearly over one facet can have its closest point on the adjacent facet land
// on their shared edge, reading a penetration proportional to the slave's
// lateral distance from the ridge (a large spurious ejecting force). A
// projection that lands on a shared edge and reads a penetration larger than
// edgeGapFrac h (h about the facet size) is rejected. A near-edge projection
// with a small gap (a slave at the ridge, where the true owner also projects
// near the edge) is kept, so the apex still holds. A free (boundary) edge and
// an interior penetration are never rejected. sharedEdge == 0: no guard.
// Returns the active flag (penetrating, in-bounds and owned).
inline bool evalSegmentSmooth(int nps, const double X[4][3], const double xs[3],
                              const double nodalNorm[4][3], const double refDir[3],
                              double &gap, double n[3], double N[4],
                              const int *sharedEdge = 0, double edgeGapFrac = 0.05,
                              double edgeTol = 0.1, double tolR = 1e-12, int maxIt = 10) {
    double xi, eta;
    int st = project(nps, X, xs, xi, eta, tolR, maxIt);
    if (st != 0) return false;                       // same projection test as evalSegment
    double dNxi[4], dNeta[4], xbar[3];
    shape(nps, xi, eta, N, dNxi, dNeta);
    // smoothed normal from the nodal field; faceted oriented normal on a degenerate blend
    if (!LadrunoContactNormalField::smoothNormal(nps, N, nodalNorm, n)) {
        if (!normalOriented(nps, xi, eta, X, refDir, n)) return false;
    }
    interp(nps, N, X, xbar);
    double d[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
    gap = dot3(n, d);
    if (gap >= 0.0) return false;                    // not penetrating
    // ownership: drop a large-gap projection on a shared edge (the neighbour owns it)
    if (sharedEdge != 0 &&
        LadrunoContactNormalField::onSharedInteriorEdge(nps, xi, eta, sharedEdge, edgeTol)) {
        double g1[3], g2[3];
        tangents(nps, dNxi, dNeta, X, g1, g2);
        double h = norm3(g1) + norm3(g2);            // about the facet extent
        if (-gap > edgeGapFrac * h) return false;
    }
    return true;                                     // penetrating, in-bounds, owned
}

// ------------------------------------------------ rich projection (mortar GP loop)
// Everything a segment-to-segment Gauss-point loop needs at one slave point.
// Unlike evalSegment(), projectFull() does not test penetration or bounds: the
// mortar overlap clip has already restricted the Gauss point to the overlap,
// so gap and status are returned as data.
struct Projection {
    int    status;       // 0 conv in-bounds, 1 conv out-of-bounds, -1 no valid projection
    double xi, eta;      // master parametric coords at the closest point
    double gap;          // gN = n.(xs - xbar)  (signed; < 0 penetration)
    double n[3];         // oriented unit normal (toward refDir)
    double t1[3], t2[3]; // covariant tangents g1 = d xbar/d xi, g2 = d xbar/d eta (not unit)
    double g[2][2];      // covariant surface metric g_ab = g_a . g_b
    double phi[4];       // master shape fns phi_i at the projection (phi[3]=0 for tri-3)
};

// returns status (also in p.status): 0 conv in-bounds, 1 conv out-of-bounds,
// -1 no valid projection (degenerate segment) or ambiguous normal (refDir
// perpendicular to n).
inline int projectFull(int nps, const double X[4][3], const double xs[3],
                       const double refDir[3], Projection &p,
                       double tolR = 1e-12, int maxIt = 10) {
    // zero the struct so a failed projection never leaves stale geometry
    p.gap = 0.0;
    for (int d = 0; d < 3; d++) { p.n[d] = 0.0; p.t1[d] = 0.0; p.t2[d] = 0.0; }
    p.g[0][0] = p.g[0][1] = p.g[1][0] = p.g[1][1] = 0.0;
    for (int i = 0; i < 4; i++) p.phi[i] = 0.0;

    int st = project(nps, X, xs, p.xi, p.eta, tolR, maxIt);
    p.status = st;
    if (st < 0) return st;                            // no valid projection

    double dNxi[4], dNeta[4], g1[3], g2[3], xbar[3];
    shape(nps, p.xi, p.eta, p.phi, dNxi, dNeta);
    interp(nps, p.phi, X, xbar);
    tangents(nps, dNxi, dNeta, X, g1, g2);
    if (!normalOriented(nps, p.xi, p.eta, X, refDir, p.n)) { p.status = -1; return -1; }

    for (int d = 0; d < 3; d++) { p.t1[d] = g1[d]; p.t2[d] = g2[d]; }
    p.g[0][0] = dot3(g1, g1);  p.g[0][1] = dot3(g1, g2);
    p.g[1][0] = p.g[0][1];     p.g[1][1] = dot3(g2, g2);
    double d[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
    p.gap = dot3(p.n, d);
    return st;                                        // 0 in-bounds, 1 out-of-bounds
}

// ----------------------------------------- geometric part of the NTS normal tangent
// Second derivatives of the shape functions (constant in xi, eta for tri-3 and
// quad-4): tri-3 all zero; quad-4 only the cross term N_{i,xi eta}.
inline void shapeD2(int nps, double Nxx[4], double Nee[4], double Nxe[4]) {
    for (int i = 0; i < 4; i++) { Nxx[i] = 0.0; Nee[i] = 0.0; Nxe[i] = 0.0; }
    if (nps == 4) {
        // d^2/(d xi d eta) of 0.25(1 +- xi)(1 +- eta): +-0.25 with the corner sign pattern
        Nxe[0] =  0.25; Nxe[1] = -0.25; Nxe[2] =  0.25; Nxe[3] = -0.25;
    }
    // tri-3: zero (a planar facet, zero curvature)
}

// Geometric block of the consistent frictionless NTS tangent:
//   K_geom = kn gN H,   H = dB/du = d^2 gN/du^2,
// the variation of the normal and of the shape functions at the moving
// projection point (curvature and large-sliding terms) that the main term
// kn B^T B omits (Wriggers and Simo 1985; Wriggers 2006). The block is
// symmetric, and its slave rows vanish for a flat facet, so for a flat fixed
// master the total tangent equals kn B^T B.
//
// Fills Kgeom[ndof][ndof] (ndof = 3 (1 + nps) <= 15) and returns true when the
// pair is active (in-bounds and penetrating); the caller adds fact Kgeom to its
// tangent. X = master segment current coordinates (nps x 3), xs = slave current
// position, refDir orients the normal (n.refDir > 0). The function re-projects
// the slave; the projection is deterministic, so it reproduces the xbar, n and
// gap of the residual.
inline bool normalGeomTangent(int nps, const double X[4][3], const double xs[3],
                              const double refDir[3], double kn, double Kgeom[15][15]) {
    const int ndof = 3 * (1 + nps);
    for (int a = 0; a < ndof; a++)
        for (int b = 0; b < ndof; b++) Kgeom[a][b] = 0.0;

    double xi, eta;
    if (project(nps, X, xs, xi, eta, 1e-12, 10) != 0) return false;   // out of bounds / no projection

    double N[4], dNxi[4], dNeta[4], Nxx[4], Nee[4], Nxe[4];
    shape(nps, xi, eta, N, dNxi, dNeta);
    shapeD2(nps, Nxx, Nee, Nxe);
    double xbar[3], a1[3], a2[3];
    interp(nps, N, X, xbar);
    tangents(nps, dNxi, dNeta, X, a1, a2);
    double raw[3]; cross3(a1, a2, raw);
    double jn = norm3(raw);
    if (jn < 1e-300) return false;
    double n[3]; for (int d = 0; d < 3; d++) n[d] = raw[d] / jn;
    double pr = dot3(n, refDir);
    if (std::fabs(pr) < 1e-12 * (norm3(refDir) + 1e-300)) return false;  // ambiguous
    if (pr < 0.0) for (int d = 0; d < 3; d++) n[d] = -n[d];
    double dvec[3] = { xs[0]-xbar[0], xs[1]-xbar[1], xs[2]-xbar[2] };
    double gN = dot3(n, dvec);
    if (gN >= 0.0) return false;                       // not penetrating

    // covariant metric m_ab, its inverse, contravariant tangents a^a = m^{ab} a_b
    double m00 = dot3(a1,a1), m01 = dot3(a1,a2), m11 = dot3(a2,a2);
    double detm = m00*m11 - m01*m01;
    if (std::fabs(detm) < 1e-300) return false;
    double mi00 = m11/detm, mi01 = -m01/detm, mi11 = m00/detm;
    double au1[3], au2[3];
    for (int d = 0; d < 3; d++) {
        au1[d] = mi00*a1[d] + mi01*a2[d];
        au2[d] = mi01*a1[d] + mi11*a2[d];
    }
    // curvature tensor kappa_ab = n.a_{a,b}, a_{a,b} = sum N_{i,ab} X_i
    double axx[3] = {0,0,0}, aee[3] = {0,0,0}, axe[3] = {0,0,0};
    for (int i = 0; i < nps; i++)
        for (int d = 0; d < 3; d++) {
            axx[d] += Nxx[i]*X[i][d]; aee[d] += Nee[i]*X[i][d]; axe[d] += Nxe[i]*X[i][d];
        }
    double k00 = dot3(n,axx), k01 = dot3(n,axe), k11 = dot3(n,aee);   // kappa (k10 = k01)
    // projection-sensitivity metric H2_ab = m_ab - gN kappa_ab and its inverse
    double H200 = m00 - gN*k00, H201 = m01 - gN*k01, H211 = m11 - gN*k11;
    double detH2 = H200*H211 - H201*H201;
    if (std::fabs(detH2) < 1e-300) return false;
    double Hi00 = H211/detH2, Hi01 = -H201/detH2, Hi11 = H200/detH2;

    // per-DOF projected components Q_b[l], P_b[l] (l = slave xyz, then master nodes)
    double Q0[15], Q1[15], P0[15], P1[15];
    for (int l = 0; l < ndof; l++) { Q0[l]=Q1[l]=P0[l]=P1[l]=0.0; }
    for (int d = 0; d < 3; d++) { Q0[d] = a1[d]; Q1[d] = a2[d]; }      // slave: P = 0
    for (int j = 0; j < nps; j++)
        for (int d = 0; d < 3; d++) {
            int l = 3*(1+j) + d;
            Q0[l] = -N[j]*a1[d]; Q1[l] = -N[j]*a2[d];
            P0[l] = dNxi[j]*n[d]; P1[l] = dNeta[j]*n[d];
        }
    // d xi^a/d u_l = H2^{ab}(Q_b + gN P_b); d n/d u_l = -(c0 a^1 + c1 a^2),
    // c_g = P_g[l] + kappa_g0 zeta0 + kappa_g1 zeta1. Assemble H[:,l] scaled by kn gN.
    for (int l = 0; l < ndof; l++) {
        double r0 = Q0[l] + gN*P0[l], r1 = Q1[l] + gN*P1[l];
        double Z0 = Hi00*r0 + Hi01*r1, Z1 = Hi01*r0 + Hi11*r1;
        double c0 = P0[l] + k00*Z0 + k01*Z1;
        double c1 = P1[l] + k01*Z0 + k11*Z1;
        double dn[3];
        for (int d = 0; d < 3; d++) dn[d] = -(c0*au1[d] + c1*au2[d]);
        double w = kn * gN;
        for (int d = 0; d < 3; d++) Kgeom[d][l] = w * dn[d];           // slave rows
        for (int i = 0; i < nps; i++) {
            double dNil = dNxi[i]*Z0 + dNeta[i]*Z1;
            for (int d = 0; d < 3; d++)
                Kgeom[3*(1+i)+d][l] = w * (-dNil*n[d] - N[i]*dn[d]);   // master-i rows
        }
    }
    return true;
}

} // namespace LadrunoContactProjection

#endif
