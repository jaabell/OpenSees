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
//   Ericson, C. (2005). Real-Time Collision Detection. Morgan Kaufmann,
//     San Francisco, Section 5.1.9. doi:10.1201/b14581
//   Wriggers, P. (2006). Computational Contact Mechanics, 2nd ed. Springer,
//     Berlin, Section 4.3. doi:10.1007/978-3-540-32609-0
//   Zavarise, G. and Wriggers, P. (2000). Contact with friction between
//     beams in 3-D space. International Journal for Numerical Methods in
//     Engineering, 49(8), 977-1006.
//     doi:10.1002/1097-0207(20001120)49:8<977::AID-NME986>3.0.CO;2-C
//
// Description: LadrunoEdgeKernel, header-only segment-to-segment closest-point
// geometry for edge-to-edge contact. Two facets that meet edge-on (their
// normals nearly perpendicular) have no in-plane overlap for the mortar clip
// to integrate, and no slave node of either projects onto the other face, so
// neither the mortar nor the node-to-segment formulation sees the contact
// between the interiors of the two crossing edges. The kernel uses no OpenSees
// types (raw double arrays and <cmath>).
//
//   closestPtSegSeg()  closest points of two segments by the clamped branch
//                      sequence of Ericson (2005), Section 5.1.9, over all
//                      Voronoi regions of the parameter square; returns the
//                      parameters (s, t), the closest points, the
//                      conditioning measure denom and a strict-interior flag
//                      (both parameters inside [delta, 1 - delta]).
//   edgeNormal()       common-perpendicular unit normal n = (d1 x d2)/|d1 x d2|,
//                      defined for non-parallel edges only.
//   edgeGap()          signed gap gN = w.n with a body-fixed sign (see below).
//   bOperator()        first variation dgN = B [dp1, dq1, dp2, dq2], exact at
//                      the closest point (the derivatives of s and t drop out
//                      because w is orthogonal to both edges there).
//
// Conditioning: denom = a e - b^2 = |d1|^2 |d2|^2 sin^2(theta), theta the
// angle between the edges, governs both the uniqueness of the closest pair
// and the conditioning of n (which grows as 1/sin theta). A pair with
// sin^2(theta) < TAU_PAR_DEFAULT = sin^2(15 deg) is refused as parallel; a
// pair with an edge shorter than 1e-6 times the other is refused as
// degenerate (a relative test, so the result does not depend on the length
// unit).
//
// Sign: the cross-product normal has no intrinsic orientation, and orienting
// it by w (n.w >= 0) makes the gap |w| >= 0 always, so penetration through
// the edge would never be detected. The sign is therefore fixed once, from a
// reference direction at the first engagement, and then held for the life of
// the pair.

#ifndef LadrunoEdgeKernel_h
#define LadrunoEdgeKernel_h

#include <cmath>

namespace LadrunoEdgeKernel {

// status codes
enum Status { EE_OK = 0, EE_PARALLEL = -1, EE_DEGENERATE = -2 };

static const double TAU_PAR_DEFAULT = 0.0669872981077807;  // sin^2(15 deg): parallel refusal
static const double TAU_LEN_REL     = 1e-12;                // squared length ratio of a collapsed edge
static const double EE_MARGIN       = 1e-6;                 // strict-interior parametric margin

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
inline double clamp01(double x) { return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }

// ------------------------------------------------------ closest point result
struct ClosestResult {
    int    status;       // EE_OK / EE_PARALLEL / EE_DEGENERATE
    double s, t;         // closest parameters on seg1 (slave) and seg2 (master), each in [0,1]
    double c1[3], c2[3]; // closest points c1 = p1 + s d1 (slave), c2 = p2 + t d2 (master)
    double denom;        // a e - b^2 = |d1|^2 |d2|^2 sin^2(theta)
    double a, e;         // squared edge lengths
    bool   interior;     // s and t both in [EE_MARGIN, 1 - EE_MARGIN]: a true edge crossing
};

// Closest points of seg1 = p1->q1 (slave) and seg2 = p2->q2 (master). On
// EE_PARALLEL the s = 0 fallback of the branch sequence is still returned.
inline ClosestResult closestPtSegSeg(const double p1[3], const double q1[3],
                                     const double p2[3], const double q2[3],
                                     double tau_par = TAU_PAR_DEFAULT) {
    double d1[3], d2[3], r[3];
    for (int k = 0; k < 3; k++) { d1[k] = q1[k]-p1[k]; d2[k] = q2[k]-p2[k]; r[k] = p1[k]-p2[k]; }
    double a = dot3(d1, d1);   // |seg1|^2
    double e = dot3(d2, d2);   // |seg2|^2

    ClosestResult o;
    o.status = EE_OK; o.s = 0.0; o.t = 0.0; o.denom = 0.0; o.a = a; o.e = e; o.interior = false;
    for (int k = 0; k < 3; k++) { o.c1[k] = p1[k]; o.c2[k] = p2[k]; }

    // collapsed edge, relative to the longer edge of the pair
    double lenRef = (a > e) ? a : e;
    if (lenRef <= 0.0 || a <= TAU_LEN_REL * lenRef || e <= TAU_LEN_REL * lenRef) {
        o.status = EE_DEGENERATE; return o;
    }

    double c = dot3(d1, r);
    double b = dot3(d1, d2);
    double f = dot3(d2, r);
    double denom = a*e - b*b;          // >= 0
    o.denom = denom;

    // the clamped branch sequence of Ericson (2005), Section 5.1.9
    double s = (denom > 1e-300) ? clamp01((b*f - c*e) / denom) : 0.0;
    double t = (b*s + f) / e;
    if (t < 0.0)      { t = 0.0; s = clamp01(-c / a); }
    else if (t > 1.0) { t = 1.0; s = clamp01((b - c) / a); }

    o.s = s; o.t = t;
    for (int k = 0; k < 3; k++) { o.c1[k] = p1[k] + s*d1[k]; o.c2[k] = p2[k] + t*d2[k]; }
    o.interior = (s >= EE_MARGIN && s <= 1.0 - EE_MARGIN &&
                  t >= EE_MARGIN && t <= 1.0 - EE_MARGIN);

    // parallel or ill-conditioned normal: sin^2(theta) < tau_par
    if (denom < tau_par * a * e) o.status = EE_PARALLEL;
    return o;
}

// Common-perpendicular unit normal n = (d1 x d2)/|d1 x d2|; false if the
// edges are parallel or of zero length.
inline bool edgeNormal(const double d1[3], const double d2[3], double n[3]) {
    double m[3];
    cross3(d1, d2, m);
    double nrm = norm3(m);
    if (nrm < 1e-300) return false;
    for (int k = 0; k < 3; k++) n[k] = m[k] / nrm;
    return true;
}

// Signed edge gap gN = w.n, w = c1 - c2 (slave minus master), with a
// body-fixed sign:
//   committedSign = +1 or -1: applied as given (a pair already engaged);
//   committedSign = 0 and signRef != 0: n is oriented so that n.signRef > 0
//     and the chosen sign is returned in *outSign for the caller to keep.
//     When signRef is numerically perpendicular to n
//     (|n.signRef| < 1e-12 |signRef|) the sign cannot be decided: *outSign = 0
//     and gN = 0 are returned, and the caller treats the evaluation as no
//     contact and decides at a later configuration;
//   committedSign = 0 and signRef = 0: the raw cross-product orientation.
// The oriented normal is written to n.
inline double edgeGap(const double c1[3], const double c2[3], const double d1[3],
                      const double d2[3], int committedSign,
                      const double *signRef, double n[3], int *outSign = 0) {
    if (!edgeNormal(d1, d2, n)) { if (outSign) *outSign = 0; return 0.0; }
    double w[3];
    for (int k = 0; k < 3; k++) w[k] = c1[k] - c2[k];
    int sign;
    if (committedSign == 1 || committedSign == -1) {
        sign = committedSign;
    } else if (signRef != 0) {
        double proj = dot3(n, signRef);
        if (std::fabs(proj) < 1e-12 * (norm3(signRef) + 1e-300)) {
            if (outSign) *outSign = 0;
            return 0.0;
        }
        sign = (proj >= 0.0) ? 1 : -1;
    } else {
        sign = 1;
    }
    for (int k = 0; k < 3; k++) n[k] *= sign;
    if (outSign) *outSign = sign;
    return dot3(w, n);
}

// First variation dgN = B . [dp1, dq1, dp2, dq2]; fills B[4][3] (rows p1, q1,
// p2, q2). The penalty residual tN B and the main tangent epsN B^T B are built
// on it.
inline void bOperator(const double n[3], double s, double t, double B[4][3]) {
    for (int k = 0; k < 3; k++) {
        B[0][k] = (1.0 - s) * n[k];   // p1 (slave node 0)
        B[1][k] =        s  * n[k];   // q1 (slave node 1)
        B[2][k] = -(1.0 - t) * n[k];  // p2 (master node 0)
        B[3][k] =       -t  * n[k];   // q2 (master node 1)
    }
}

} // namespace LadrunoEdgeKernel

#endif
