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
//   Puso, M.A. and Laursen, T.A. (2004). A mortar segment-to-segment contact
//     method for large deformation solid mechanics. Computer Methods in
//     Applied Mechanics and Engineering, 193(6-8), 601-629.
//     doi:10.1016/j.cma.2003.10.010
//   Popp, A., Gee, M.W. and Wall, W.A. (2009). A finite deformation mortar
//     contact formulation using a primal-dual active set strategy.
//     International Journal for Numerical Methods in Engineering, 79(11),
//     1354-1391. doi:10.1002/nme.2614
//   Wriggers, P. (2006). Computational Contact Mechanics, 2nd ed. Springer,
//     Berlin. doi:10.1007/978-3-540-32609-0
//   Sutherland, I.E. and Hodgman, G.W. (1974). Reentrant polygon clipping.
//     Communications of the ACM, 17(1), 32-42. doi:10.1145/360767.360802
//   Dunavant, D.A. (1985). High degree efficient symmetrical Gaussian
//     quadrature rules for the triangle. International Journal for Numerical
//     Methods in Engineering, 21(6), 1129-1148. doi:10.1002/nme.1620210612
//
// Description: LadrunoMortarKernel, header-only segment-to-segment mortar
// integration for one slave/master facet pair: overlap clip, sub-triangle
// Gauss quadrature, the mortar matrices D_IJ and M_IK and the weighted gap
// g_I. It uses no OpenSees types (raw double arrays, <cmath> and
// LadrunoContactProjection).
//
// The overlap of slave and master is clipped on an auxiliary plane of the
// master facet (Sutherland-Hodgman) and fan-triangulated. Each Gauss point is
// projected onto the master facet with LadrunoContactProjection::projectFull.
//
// Standard (not dual) basis:
//   D_IJ += w J N_I^s N_J^s                  (slave-slave mortar matrix)
//   M_IK += w J N_I^s phi_K^m                (slave-master mortar matrix)
//   g_I  += w J N_I^s ((x_s - x_m).n)        (weighted gap, normal inside the integral)
// with J = A_aux / |n_s.n0| (aux-plane cell area to slave-surface area; exact
// and constant for a flat slave facet). The kernel returns D, M, g and the
// per-pair measures only; tractions, multipliers and friction belong to the
// caller.
//
// Scope and limitations:
//   - Flat, convex tri/quad facets. A non-convex (concave or bow-tie) facet is
//     refused (integratePair returns -1); it would otherwise be integrated over
//     its convex hull. A warped slave facet is not refused, but the constant
//     J = A_aux/|n_s.n0| then biases the area measure (about 0.7% for an
//     out-of-plane lift of 0.3 times the edge); partition of unity and the
//     patch test still hold because they reduce to sum phi = 1.
//   - Serendipity facets, quad8 (nps 8) and tri6 (nps 6) with corner nodes
//     first (the 20-node brick face convention), are accepted. All geometry
//     (aux plane, clip, inverse map, projection) runs on the corner polygon;
//     only the basis evaluation uses the full serendipity functions. This is
//     exact for a flat facet whose midside nodes sit at the edge midpoints,
//     because the serendipity map then reduces to the corner map; that
//     condition is checked, and a curved or displaced-midside facet returns
//     status -2 (invalid geometry: the caller must stop with an error; -1
//     means an empty overlap and may be skipped).
//
// Verified properties: clip areas of known cases (including the 2(sqrt(2)-1)
// octagon); partition of unity sum_K M_IK = sum_J D_IJ to 1e-12 on matching
// and non-matching meshes; the constant-pressure patch test on a non-matching
// two-block mesh (P = D^-1 M transfers a uniform field to 1e-6); the weighted
// gap is exact for uniform and tilted-plane gaps.

#ifndef LadrunoMortarKernel_h
#define LadrunoMortarKernel_h

#include <cmath>
#include "LadrunoContactProjection.h"

namespace LadrunoMortarKernel {

namespace LP = LadrunoContactProjection;   // projection geometry (projectFull, shape, ...)

// limits: a facet has <= 8 nodes (serendipity) but its corner polygon has <= 4; the
// clip runs on corner polygons, so a convex tri/quad ∩ tri/quad clip still has ≤8
// vertices and ≤8 fan sub-triangles; the GP rule has ≤12 points per sub-triangle
// (degree-6 Dunavant rule for quadratic pairs).
static const int MAXN  = 8;    // max facet nodes (serendipity quad8)
static const int MAXV  = 8;    // max clip-polygon vertices
static const int MAXGP = 12;   // max Gauss points per sub-triangle

// corner count of a facet: tri3/quad4 are all corners; tri6/quad8 are corners-first.
inline int ncorner(int nps) { return (nps <= 4) ? nps : nps / 2; }

// ------------------------------------------------------------------ result of one pair
struct PairResult {
    int    status;     // 0 accumulated, -1 empty/degenerate overlap (skip),
                       // -2 INVALID quadratic facet geometry (caller must hard-error)
    int    ngp;        // number of Gauss points integrated (0 ⇒ no contribution)
    double area;       // ∫ dΓ over the overlap (slave-surface measure)
    double gapL1;      // ∫ |g_N| dΓ over the overlap (a sign-free gap measure)
    double D[MAXN][MAXN]; // slave–slave mortar  D_IJ   (rows/cols = slave nodes)
    double M[MAXN][MAXN]; // slave–master mortar M_IK   (rows = slave, cols = master nodes)
    double g[MAXN];    // weighted gap g̃_I           (one per slave node)
};

// =============================================================== small 2D helpers
inline double signedArea2(const double P[][2], int n) {
    double a = 0.0;
    for (int i = 0; i < n; i++) { int j = (i + 1) % n; a += P[i][0]*P[j][1] - P[j][0]*P[i][1]; }
    return 0.5 * a;
}

// reverse the winding of an n-vertex polygon in place (to enforce CCW).
inline void reverse2(double P[][2], int n) {
    for (int i = 0, j = n - 1; i < j; i++, j--) {
        double tx = P[i][0], ty = P[i][1];
        P[i][0] = P[j][0]; P[i][1] = P[j][1];
        P[j][0] = tx;      P[j][1] = ty;
    }
}
inline void ensureCCW(double P[][2], int n) { if (signedArea2(P, n) < 0.0) reverse2(P, n); }

// convexity test for a CCW polygon: every consecutive turn must be a (weak) left turn.
// The mortar path requires convex projected facets: Sutherland-Hodgman needs a convex
// CLIP polygon, and a non-convex SUBJECT would be silently integrated over its convex
// hull (wrong area) or truncate the clip buffer. A flat,
// convex tri/quad facet always projects to a convex polygon; a concave/bow-tie or
// strongly-warped facet does not, and the caller refuses the pair.
inline bool isConvex2(const double P[][2], int n, double tol = 1e-12) {
    if (n < 3) return false;
    double scale = 0.0;                          // edge-length scale for a relative tol
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        double ex = P[j][0]-P[i][0], ey = P[j][1]-P[i][1];
        double e = std::sqrt(ex*ex + ey*ey);
        if (e > scale) scale = e;
    }
    if (scale < 1e-300) return false;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n, k = (i + 2) % n;
        double ax = P[j][0]-P[i][0], ay = P[j][1]-P[i][1];
        double bx = P[k][0]-P[j][0], by = P[k][1]-P[j][1];
        if (ax*by - ay*bx < -tol * scale * scale) return false;   // a right turn ⇒ concave
    }
    return true;
}

// intersection of segment S->E with the infinite line A->B (Sutherland-Hodgman edge).
inline void iline(const double S[2], const double E[2], const double A[2],
                  const double B[2], double out[2]) {
    double d1x = E[0]-S[0], d1y = E[1]-S[1];
    double d2x = B[0]-A[0], d2y = B[1]-A[1];
    double den = d1x*d2y - d1y*d2x;
    if (std::fabs(den) < 1e-300) { out[0] = S[0]; out[1] = S[1]; return; }
    double t = ((A[0]-S[0])*d2y - (A[1]-S[1])*d2x) / den;
    out[0] = S[0] + t*d1x;  out[1] = S[1] + t*d1y;
}

// Sutherland-Hodgman: subject ∩ clip, both CCW, CLIP CONVEX (tri/quad ⇒ convex). Writes
// the (convex) intersection into out[][2], returns its vertex count (≤ ns+nc, ≤ MAXV).
inline int clipPolygon(const double subj[][2], int ns, const double clip[][2], int nc,
                       double out[][2]) {
    double buf[2][MAXV + 4][2];
    int cur = 0, n = ns;
    for (int i = 0; i < ns; i++) { buf[cur][i][0] = subj[i][0]; buf[cur][i][1] = subj[i][1]; }
    for (int e = 0; e < nc; e++) {
        const double *A = clip[e];
        const double *B = clip[(e + 1) % nc];
        double ex = B[0]-A[0], ey = B[1]-A[1];
        int nxt = 1 - cur, m = 0;
        if (n == 0) { cur = nxt; break; }
        const double *S = buf[cur][n - 1];
        double sS = ex*(S[1]-A[1]) - ey*(S[0]-A[0]);   // >0 ⇒ inside (left of A->B)
        for (int i = 0; i < n && m < MAXV + 3; i++) {
            const double *E = buf[cur][i];
            double sE = ex*(E[1]-A[1]) - ey*(E[0]-A[0]);
            if (sE >= 0.0) {
                if (sS < 0.0) { iline(S, E, A, B, buf[nxt][m]); m++; }
                buf[nxt][m][0] = E[0]; buf[nxt][m][1] = E[1]; m++;
            } else if (sS >= 0.0) {
                iline(S, E, A, B, buf[nxt][m]); m++;
            }
            S = E; sS = sE;
        }
        cur = nxt; n = m;
    }
    if (n > MAXV) n = MAXV;
    for (int i = 0; i < n; i++) { out[i][0] = buf[cur][i][0]; out[i][1] = buf[cur][i][1]; }
    return n;
}

// merge near-coincident consecutive vertices; returns the surviving count.
inline int dedupe(double P[][2], int n, double tol = 1e-12) {
    double tmp[MAXV][2]; int m = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        double dx = P[i][0]-P[j][0], dy = P[i][1]-P[j][1];
        if (std::sqrt(dx*dx + dy*dy) > tol) { tmp[m][0] = P[i][0]; tmp[m][1] = P[i][1]; m++; }
    }
    for (int i = 0; i < m; i++) { P[i][0] = tmp[i][0]; P[i][1] = tmp[i][1]; }
    return m;
}

// =============================================================== facet geometry
// unit normal of a (flat) facet, oriented so n·refDir > 0. quad: cross of the diagonals.
inline bool facetNormal(int nps, const double X[4][3], const double refDir[3], double n[3]) {
    double a[3], b[3], raw[3];
    if (nps == 3) {
        for (int d = 0; d < 3; d++) { a[d] = X[1][d]-X[0][d]; b[d] = X[2][d]-X[0][d]; }
    } else {
        for (int d = 0; d < 3; d++) { a[d] = X[2][d]-X[0][d]; b[d] = X[3][d]-X[1][d]; }
    }
    LP::cross3(a, b, raw);
    double j = LP::norm3(raw);
    if (j < 1e-300) return false;
    for (int d = 0; d < 3; d++) n[d] = raw[d] / j;
    if (LP::dot3(n, refDir) < 0.0) for (int d = 0; d < 3; d++) n[d] = -n[d];
    return true;
}

// auxiliary plane from the MASTER facet: origin x0=centroid, normal n0, orthonormal axes.
inline bool auxPlane(int nps, const double X[4][3], const double refDir[3],
                     double x0[3], double n0[3], double e1[3], double e2[3]) {
    for (int d = 0; d < 3; d++) {
        x0[d] = 0.0;
        for (int i = 0; i < nps; i++) x0[d] += X[i][d];
        x0[d] /= nps;
    }
    if (!facetNormal(nps, X, refDir, n0)) return false;
    double t[3];
    for (int d = 0; d < 3; d++) t[d] = X[1][d] - X[0][d];
    double tn = LP::dot3(t, n0);
    for (int d = 0; d < 3; d++) e1[d] = t[d] - tn * n0[d];
    double e1n = LP::norm3(e1);
    if (e1n < 1e-300) return false;
    for (int d = 0; d < 3; d++) e1[d] /= e1n;
    LP::cross3(n0, e1, e2);
    return true;
}

inline void to2d(const double X[3], const double x0[3], const double e1[3],
                 const double e2[3], double uv[2]) {
    double d[3] = { X[0]-x0[0], X[1]-x0[1], X[2]-x0[2] };
    uv[0] = LP::dot3(d, e1);  uv[1] = LP::dot3(d, e2);
}

// recover ξ_s s.t. Σ N_I(ξ_s)·UV_I = q  (UV = slave nodes in aux 2D). tri-3: affine;
// quad: bounded inverse-bilinear Newton (same pattern as LadrunoContactProjection).
// Returns TRUE on convergence (residual < tolR / non-degenerate affine map); FALSE lets
// the caller skip that GP instead of integrating a wrong root. With the
// convex-facet test in integratePair this does not fail, but, as in
// LadrunoContactProjection::project, a last iterate is never returned silently.
inline bool inverseIsomap2D(int nps, const double UV[4][2], const double q[2],
                            double &xi, double &eta, double tolR = 1e-13, int maxIt = 20) {
    if (nps == 3) {
        double a00 = UV[1][0]-UV[0][0], a01 = UV[2][0]-UV[0][0];
        double a10 = UV[1][1]-UV[0][1], a11 = UV[2][1]-UV[0][1];
        double det = a00*a11 - a01*a10;
        double scale = std::sqrt((a00*a00+a10*a10)*(a01*a01+a11*a11));
        if (std::fabs(det) < 1e-14 * (scale + 1e-300)) { xi = eta = 0.0; return false; }
        double rx = q[0]-UV[0][0], ry = q[1]-UV[0][1];
        xi  = ( a11*rx - a01*ry) / det;
        eta = (-a10*rx + a00*ry) / det;
        return true;
    }
    xi = 0.0; eta = 0.0;
    double N[4], dNxi[4], dNeta[4];
    bool conv = false;
    for (int it = 0; it < maxIt; it++) {
        LP::shape(4, xi, eta, N, dNxi, dNeta);
        double rx = -q[0], ry = -q[1], j00 = 0, j01 = 0, j10 = 0, j11 = 0;
        for (int i = 0; i < 4; i++) {
            rx += N[i]*UV[i][0];  ry += N[i]*UV[i][1];
            j00 += dNxi[i]*UV[i][0];  j01 += dNeta[i]*UV[i][0];
            j10 += dNxi[i]*UV[i][1];  j11 += dNeta[i]*UV[i][1];
        }
        if (std::sqrt(rx*rx + ry*ry) < tolR) { conv = true; break; }
        double det = j00*j11 - j01*j10;
        if (std::fabs(det) < 1e-300) break;
        double dxi  = ( j11*rx - j01*ry) / det;
        double deta = (-j10*rx + j00*ry) / det;
        xi  -= dxi;
        eta -= deta;
        // As in LadrunoContactProjection::project: the absolute tolR
        // (1e-13) operates on a length-unit residual (aux-plane UV ~ facet size h) whose
        // noise floor is ~eps·h; for facet edges ≳2000 length units it is unreachable
        // and the Gauss point would be skipped (378 of 400 points at h = 5e4). The parametric step
        // is scale-free (parent coords are O(1)); step-converged returns the identical
        // (ξ,η) as residual-converged one iteration later for the affine-dominated maps
        // the convex-facet test admits.
        if (std::fabs(dxi) + std::fabs(deta) < 1e-10) { conv = true; break; }
    }
    return conv;
}

// =============================================================== full (serendipity) basis
// Evaluate the full facet basis at parent (ξ,η): tri3/quad4 delegate to LP::shape
// (identical values on the linear path); tri6/quad8 are the standard serendipity
// functions, corners first. N must hold MAXN entries. Derivatives are not provided: all
// geometry (tangents, normals, inverse maps) runs on the corner polygon by design.
inline void shapeFull(int nps, double xi, double eta, double N[MAXN]) {
    if (nps <= 4) {
        double dNxi[4], dNeta[4];
        LP::shape(nps, xi, eta, N, dNxi, dNeta);
        for (int i = nps; i < MAXN; i++) N[i] = 0.0;
        return;
    }
    if (nps == 8) {
        const double xs[4] = {-1.0, 1.0, 1.0, -1.0};
        const double es[4] = {-1.0, -1.0, 1.0, 1.0};
        for (int i = 0; i < 4; i++)
            N[i] = 0.25 * (1 + xi*xs[i]) * (1 + eta*es[i]) * (xi*xs[i] + eta*es[i] - 1);
        N[4] = 0.5 * (1 - xi*xi) * (1 - eta);    // mid 0-1
        N[5] = 0.5 * (1 + xi) * (1 - eta*eta);   // mid 1-2
        N[6] = 0.5 * (1 - xi*xi) * (1 + eta);    // mid 2-3
        N[7] = 0.5 * (1 - xi) * (1 - eta*eta);   // mid 3-0
        return;
    }
    // tri6: area coords λ0=1-ξ-η, λ1=ξ, λ2=η; corners λ(2λ-1), mids 4λaλb.
    double l0 = 1 - xi - eta, l1 = xi, l2 = eta;
    N[0] = l0 * (2*l0 - 1);  N[1] = l1 * (2*l1 - 1);  N[2] = l2 * (2*l2 - 1);
    N[3] = 4*l0*l1;          N[4] = 4*l1*l2;          N[5] = 4*l2*l0;
    N[6] = 0.0;              N[7] = 0.0;
}

// A quadratic facet is admissible only when every midside node sits at
// its edge midpoint (3D distance ≤ tol · max corner-edge length; this catches a midside
// lifted off its straight edge, i.e. a genuinely curved EDGE). It deliberately does NOT
// check corner coplanarity: the serendipity map reduces identically to the corner
// (bilinear) map for midpoint midsides even with warped corners (checked
// 2026-08-04, verified numerically to ~1e-15 over random non-planar corner sets), so a
// warped quad8 degrades exactly like a warped quad4: the ~0.7%
// constant-J area bias noted in the file header. Linear facets are
// trivially valid. A violation must surface as status -2 (hard error at the caller),
// never a silent skip.
inline bool validQuadraticFacet(int nps, const double X[MAXN][3], double tol = 1e-6) {
    int npc = ncorner(nps);
    if (npc == nps) return true;
    double scale = 0.0;
    for (int a = 0; a < npc; a++) {
        int b = (a + 1) % npc;
        double e[3] = { X[b][0]-X[a][0], X[b][1]-X[a][1], X[b][2]-X[a][2] };
        double L = LP::norm3(e);
        if (L > scale) scale = L;
    }
    if (scale < 1e-300) return false;
    for (int k = 0; k < npc; k++) {                     // midside k lies on edge (k, k+1)
        int a = k, b = (k + 1) % npc;
        double d[3] = { X[npc+k][0] - 0.5*(X[a][0]+X[b][0]),
                        X[npc+k][1] - 0.5*(X[a][1]+X[b][1]),
                        X[npc+k][2] - 0.5*(X[a][2]+X[b][2]) };
        if (LP::norm3(d) > tol * scale) return false;
    }
    return true;
}

// barycentric triangle Gauss rule, weights SUM TO 1 (area applied separately). order<=2
// ⇒ 3-pt (degree 2); order<=4 ⇒ the 6-pt degree-4 rule; else the 12-pt degree-6
// Dunavant rule (serendipity N·φ products reach degree 6). Returns the count.
inline int triQuadrature(int order, double bary[MAXGP][3], double w[MAXGP]) {
    if (order <= 2) {
        const double b[3][3] = {{2./3,1./6,1./6},{1./6,2./3,1./6},{1./6,1./6,2./3}};
        for (int i = 0; i < 3; i++) { for (int k = 0; k < 3; k++) bary[i][k] = b[i][k]; w[i] = 1./3; }
        return 3;
    }
    if (order <= 4) {
        const double a1 = 0.445948490915965, b1 = 0.108103018168070, w1 = 0.223381589678011;
        const double a2 = 0.091576213509771, b2 = 0.816847572980459, w2 = 0.109951743655322;
        const double b[6][3] = {{b1,a1,a1},{a1,b1,a1},{a1,a1,b1},{b2,a2,a2},{a2,b2,a2},{a2,a2,b2}};
        const double ww[6]   = {w1,w1,w1,w2,w2,w2};
        for (int i = 0; i < 6; i++) { for (int k = 0; k < 3; k++) bary[i][k] = b[i][k]; w[i] = ww[i]; }
        return 6;
    }
    // Dunavant degree-6, 12 points, weights sum to 1 (exact for monomials with p+q≤6).
    const double p1 = 0.501426509658179, q1 = 0.249286745170910, w1 = 0.116786275726379;
    const double p2 = 0.873821971016996, q2 = 0.063089014491502, w2 = 0.050844906370207;
    const double p3 = 0.053145049844817, q3 = 0.310352451033784;
    const double r3 = 1.0 - p3 - q3,     w3 = 0.082851075618374;
    const double b[12][3] = {
        {p1,q1,q1},{q1,p1,q1},{q1,q1,p1},
        {p2,q2,q2},{q2,p2,q2},{q2,q2,p2},
        {p3,q3,r3},{p3,r3,q3},{q3,p3,r3},{q3,r3,p3},{r3,p3,q3},{r3,q3,p3}};
    const double ww[12] = {w1,w1,w1,w2,w2,w2,w3,w3,w3,w3,w3,w3};
    for (int i = 0; i < 12; i++) { for (int k = 0; k < 3; k++) bary[i][k] = b[i][k]; w[i] = ww[i]; }
    return 12;
}

// =============================================================== the mortar pair
// Integrate ONE slave/master facet pair (the caller's broad phase supplies the pairing).
// Accumulates D, M, g̃, area, gapL1 into `out` (zeroed here). Returns 0 if a contribution
// was added, -1 if the overlap is empty / a sliver / degenerate (caller simply skips it),
// -2 if a QUADRATIC facet violates the midside-at-midpoint condition (the caller
// MUST hard-error with a named message; silently skipping -2 would drop the bond).
//
// Quadratic (tri6/quad8, corners-first) facets are admitted; all geometry below
// runs on the CORNER polygon (counts ncs/ncm) and only the basis evaluation (shapeFull)
// sees the full node set, exact because the check pins midsides to midpoints.
inline int integratePair(int nps_s, const double Xs[MAXN][3],
                         int nps_m, const double Xm[MAXN][3],
                         const double refDir[3], PairResult &out, int order = 2) {
    out.status = -1; out.ngp = 0; out.area = 0.0; out.gapL1 = 0.0;
    for (int i = 0; i < MAXN; i++) { out.g[i] = 0.0;
        for (int j = 0; j < MAXN; j++) { out.D[i][j] = 0.0; out.M[i][j] = 0.0; } }

    if (!validQuadraticFacet(nps_s, Xs) || !validQuadraticFacet(nps_m, Xm)) {
        out.status = -2;
        return -2;                                      // invalid quadratic geometry
    }
    const int ncs = ncorner(nps_s), ncm = ncorner(nps_m);

    double x0[3], n0[3], e1[3], e2[3], n_s[3];
    if (!auxPlane(ncm, Xm, refDir, x0, n0, e1, e2)) return -1;
    if (!facetNormal(ncs, Xs, refDir, n_s)) return -1;
    double cos_t = std::fabs(LP::dot3(n_s, n0));
    if (cos_t < 1e-12) return -1;                       // slave ⟂ aux plane: ill-posed

    // slave corners in aux 2D: node order (for the GP back-map) and a CCW copy (clip subject)
    double UVs[4][2], subj[4][2], clip[4][2];
    for (int i = 0; i < ncs; i++) { to2d(Xs[i], x0, e1, e2, UVs[i]);
        subj[i][0] = UVs[i][0]; subj[i][1] = UVs[i][1]; }
    for (int i = 0; i < ncm; i++) to2d(Xm[i], x0, e1, e2, clip[i]);
    ensureCCW(subj, ncs);
    ensureCCW(clip, ncm);
    // Scope: flat, convex facets (corner polygons). A concave/bow-tie or strongly-
    // warped facet projects to a non-convex polygon ⇒ Sutherland-Hodgman would silently
    // integrate the slave over its convex hull and truncate the clip buffer.
    // Refuse the pair rather than return a wrong area.
    if (!isConvex2(subj, ncs) || !isConvex2(clip, ncm)) return -1;

    double poly[MAXV][2];
    int np = clipPolygon(subj, ncs, clip, ncm, poly);
    np = dedupe(poly, np);
    if (np < 3) return -1;

    double A_s = std::fabs(signedArea2(subj, ncs));
    double A_m = std::fabs(signedArea2(clip, ncm));
    double A_min = (A_s < A_m) ? A_s : A_m;
    double A_poly = std::fabs(signedArea2(poly, np));
    if (A_poly < 1e-12 * A_min) return -1;              // sliver guard

    double cen[2] = {0.0, 0.0};
    for (int i = 0; i < np; i++) { cen[0] += poly[i][0]; cen[1] += poly[i][1]; }
    cen[0] /= np; cen[1] /= np;

    double bary[MAXGP][3], w[MAXGP];
    int ngpr = triQuadrature(order, bary, w);

    // fan from the centroid → sub-triangles; Gauss each
    for (int s = 0; s < np; s++) {
        const double *V1 = poly[s];
        const double *V2 = poly[(s + 1) % np];
        double triP[3][2] = { {cen[0],cen[1]}, {V1[0],V1[1]}, {V2[0],V2[1]} };
        double A_aux = std::fabs(signedArea2(triP, 3));
        if (A_aux < 1e-12 * A_min) continue;
        double A_phys = A_aux / cos_t;                  // aux-plane area → slave-surface area
        for (int gp = 0; gp < ngpr; gp++) {
            double q[2] = {
                bary[gp][0]*triP[0][0] + bary[gp][1]*triP[1][0] + bary[gp][2]*triP[2][0],
                bary[gp][0]*triP[0][1] + bary[gp][1]*triP[1][1] + bary[gp][2]*triP[2][1] };
            double xi_s, eta_s;
            if (!inverseIsomap2D(ncs, UVs, q, xi_s, eta_s)) continue;    // back-map failed
            double Ns[MAXN];
            shapeFull(nps_s, xi_s, eta_s, Ns);          // full slave basis at the mapped GP
            double x_s[3] = {0,0,0};
            for (int i = 0; i < nps_s; i++) for (int d = 0; d < 3; d++) x_s[d] += Ns[i]*Xs[i][d];

            LP::Projection p;
            int st = LP::projectFull(ncm, Xm, x_s, refDir, p);
            if (st < 0) continue;                        // no valid projection: skip GP
            // master basis at the projection: linear masters use p.phi verbatim (the
            // linear path); quadratic masters evaluate the serendipity
            // basis at the same corner-map (ξ,η), exact under the midside condition.
            double phi[MAXN];
            if (nps_m <= 4) { for (int k = 0; k < 4; k++) phi[k] = p.phi[k]; }
            else            shapeFull(nps_m, p.xi, p.eta, phi);
            double x_m[3] = {0,0,0};
            for (int k = 0; k < nps_m; k++) for (int d = 0; d < 3; d++) x_m[d] += phi[k]*Xm[k][d];
            double dx[3] = { x_s[0]-x_m[0], x_s[1]-x_m[1], x_s[2]-x_m[2] };
            double gN = LP::dot3(p.n, dx);

            double wJ = w[gp] * A_phys;
            out.area += wJ;
            out.gapL1 += wJ * std::fabs(gN);
            for (int i = 0; i < nps_s; i++) {
                out.g[i] += wJ * Ns[i] * gN;
                for (int j = 0; j < nps_s; j++) out.D[i][j] += wJ * Ns[i] * Ns[j];
                for (int k = 0; k < nps_m; k++) out.M[i][k] += wJ * Ns[i] * phi[k];
            }
            out.ngp++;
        }
    }
    out.status = (out.ngp > 0) ? 0 : -1;
    return out.status;
}

} // namespace LadrunoMortarKernel

#endif
