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
//     Berlin. doi:10.1007/978-3-540-32609-0
//   Laursen, T.A. (2002). Computational Contact and Impact Mechanics.
//     Springer, Berlin. doi:10.1007/978-3-662-04864-1
//   Puso, M.A. and Laursen, T.A. (2004). A mortar segment-to-segment contact
//     method for large deformation solid mechanics. Computer Methods in Applied
//     Mechanics and Engineering, 193(6-8), 601-629. doi:10.1016/j.cma.2003.10.010
//
// Description: LadrunoContact2DKernel, header-only geometry of 2D (-ndm 2)
// contact: a slave point against a 2-node master line segment, the corner
// (vertex) pair that covers what both adjacent segments refuse, and the
// segment-to-segment mortar integral of two 2-node segments. It uses no
// OpenSees types (raw double arrays and <cmath>) and holds no state: the
// orientation of a master surface, the committed side sign of a vertex pair
// and the ownership decisions are the caller's.
//
//   projectSegment2D()   closed-form projection of a point on a segment, the
//                        oriented normal and the gap, with a refusal status
//                        (no clamping: a point past an end belongs to the
//                        neighbouring segment or to a vertex pair).
//   vertexWedge2D()      classification at a shared vertex: is the slave in
//                        the vertex's Voronoi wedge, is the corner convex or
//                        concave seen from the contact side, and on which side
//                        of the corner is the slave.
//   vertexEval2D()       the vertex pair: radial normal and signed gap from a
//                        committed side sign.
//   sigmaFromRef2D()     the orientation sign that turns a segment normal
//                        toward a reference direction (0 when ambiguous).
//   mortarSegmentPair2D() mortar matrices D, M and the weighted gap of a slave
//                        and a master segment.
//
// Corners. A convex corner's Voronoi wedge lies outside the master body, so a
// slave in it has a positive gap and the vertex pair never carries force; the
// region inside the corner is covered by both adjacent segments. A concave
// (re-entrant) corner is the opposite: the region behind it, inside the
// material, is refused by both adjacent segments and is entirely penetrating,
// a coverage hole of finite measure. Vertex pairs are therefore built at
// concave corners only, and at the open ends of a master chain (end caps).
//
// Tolerances are relative or dimensionless: tolIn is a parametric slack on xi;
// tauSeg refuses a segment shorter than tauSeg Lref, Lref the minimum initial
// segment length of the surface (a caller contract; Lref <= 0 is refused);
// tolLen (vertexEval2D) is a length and is passed relative to Lref; tauPerp is
// the conditioning gate of the two ambiguity tests (a reference direction
// tangent to the segment, a fold-back corner).

#ifndef LadrunoContact2DKernel_h
#define LadrunoContact2DKernel_h

#include <cmath>

namespace LadrunoContact2DKernel {

// projectSegment2D status; OUT_LOW and OUT_HIGH are refusals, split because
// the routing at a vertex needs to know which end the slave ran past
enum SegStatus {
    SEG2D_IN_BOUNDS  =  0,   // xi in [-tolIn, 1 + tolIn]
    SEG2D_OUT_LOW    =  1,   // xi < -tolIn (past X0)
    SEG2D_OUT_HIGH   =  2,   // xi > 1 + tolIn (past X1)
    SEG2D_DEGENERATE = -1    // zero-length segment or Lref <= 0
};

enum VtxStatus {
    VTX2D_OK         =  0,
    VTX2D_DEGENERATE = -1    // the slave is on the vertex (no direction)
};

// corner type seen from the contact side; corner == sign(turnSin)
enum CornerType {
    C2D_CONVEX   = -1,
    C2D_STRAIGHT =  0,
    C2D_CONCAVE  = +1
};

static const double TOL_IN_BOUNDS_DEFAULT = 1e-9;    // parametric slack on xi
static const double TAU_SEG_DEFAULT       = 1e-8;    // zero-length refusal, relative to Lref
static const double TAU_PERP_DEFAULT      = 1e-12;   // ambiguity gate

inline double dot2(const double a[2], const double b[2]) { return a[0]*b[0] + a[1]*b[1]; }
inline double cross2(const double a[2], const double b[2]) { return a[0]*b[1] - a[1]*b[0]; }
inline double norm2(const double a[2]) { return std::sqrt(dot2(a, a)); }
// perp(t) = (-t_y, t_x), the +90 degree rotation; n = sigma perp(t) / L
inline void perp2(const double a[2], double p[2]) { p[0] = -a[1]; p[1] = a[0]; }

// linear 2-node segment x(xi) = (1 - xi) X0 + xi X1, xi in [0, 1]
inline void shape2D(double xi, double N[2]) { N[0] = 1.0 - xi; N[1] = xi; }

// Projection of xs on the segment X0 -> X1:
//   t = X1 - X0, xi = (xs - X0).t / (t.t), n = sigma perp(t) / |t|,
//   gap = (xs - x(xi)).n  (< 0: penetration).
// On OUT_LOW / OUT_HIGH xi, n and gap are filled (the gap is the distance to
// the infinite line, not a contact gap); on DEGENERATE all outputs are zero.
inline int projectSegment2D(const double X0[2], const double X1[2], const double xs[2],
                            double sigma, double Lref,
                            double &xi, double &gap, double n[2],
                            double tolIn  = TOL_IN_BOUNDS_DEFAULT,
                            double tauSeg = TAU_SEG_DEFAULT)
{
    xi = 0.0; gap = 0.0; n[0] = 0.0; n[1] = 0.0;

    const double t[2] = { X1[0] - X0[0], X1[1] - X0[1] };
    const double L2 = dot2(t, t);
    const double L  = std::sqrt(L2);

    if (!(Lref > 0.0))        return SEG2D_DEGENERATE;   // also traps NaN
    if (L < tauSeg * Lref)    return SEG2D_DEGENERATE;

    const double d0[2] = { xs[0] - X0[0], xs[1] - X0[1] };
    xi = dot2(d0, t) / L2;

    n[0] = sigma * (-t[1]) / L;
    n[1] = sigma * ( t[0]) / L;

    // through the shape functions, so that the residual and the B operator see
    // the same x(xi)
    double N[2];
    shape2D(xi, N);
    const double xbar[2] = { N[0]*X0[0] + N[1]*X1[0],
                             N[0]*X0[1] + N[1]*X1[1] };
    const double d[2] = { xs[0] - xbar[0], xs[1] - xbar[1] };
    gap = dot2(d, n);

    if (xi < -tolIn)       return SEG2D_OUT_LOW;
    if (xi > 1.0 + tolIn)  return SEG2D_OUT_HIGH;
    return SEG2D_IN_BOUNDS;
}

// Result of vertexWedge2D (classification only).
struct WedgeResult {
    int    inWedge;    // 1: xs is past both adjacent segments toward the vertex
    int    corner;     // CornerType, sign(turnSin)
    double sideSign;   // +1 exterior (separated) side, -1 interior (penetrating)
                       // side, 0 undetermined (fold-back corner, xs at the vertex
                       // or within the conditioning gate): the caller defers
    double turnSin;    // sigma cross(t_prev, t_next) / (|t_prev| |t_next|)
    double aPrev;      // (xs - Xv).t_prev / |t_prev|^2 = xi_prev - 1
    double aNext;      // (xs - Xv).t_next / |t_next|^2 = xi_next
};

// Corner classifier at the vertex Xv shared by the segment ending there
// (tPrev = Xv - X_prev) and the one starting there (tNext = X_next - Xv).
// The vertex owns xs when xs ran off both segments toward the vertex,
// aPrev > tolIn and aNext < -tolIn. Because xi_prev is formed from X_prev it
// resolves the overshoot only to ulp(1), while aPrev resolves it exactly, so
// the caller decides ownership by: (1) the previous segment owns if in
// bounds, (2) else the next one, (3) else the vertex if aPrev > 0 and
// aNext < 0 (no slack), (4) else nobody. The side sign is taken from the
// bisector of the two outward normals, which is correct for every turn angle
// below pi (a single normal is not past a right angle, and sign(turnSin) is a
// coin flip at a nearly straight or a fold-back corner).
inline WedgeResult vertexWedge2D(const double tPrev[2], const double tNext[2],
                                 double sigma, const double xs[2], const double Xv[2],
                                 double Lref,
                                 double tolIn   = TOL_IN_BOUNDS_DEFAULT,
                                 double tauSeg  = TAU_SEG_DEFAULT,
                                 double tauPerp = TAU_PERP_DEFAULT)
{
    WedgeResult w;
    w.inWedge = 0; w.corner = C2D_STRAIGHT; w.sideSign = 0.0;
    w.turnSin = 0.0; w.aPrev = 0.0; w.aNext = 0.0;

    const double LP2 = dot2(tPrev, tPrev), LN2 = dot2(tNext, tNext);
    const double LP  = std::sqrt(LP2),     LN  = std::sqrt(LN2);

    // the refusal of projectSegment2D on the same gauge, so that the
    // ownership map stays consistent
    if (!(Lref > 0.0))                                   return w;
    if (LP < tauSeg * Lref || LN < tauSeg * Lref)         return w;

    const double r[2] = { xs[0] - Xv[0], xs[1] - Xv[1] };
    w.aPrev = dot2(r, tPrev) / LP2;
    w.aNext = dot2(r, tNext) / LN2;
    w.inWedge = (w.aPrev > tolIn && w.aNext < -tolIn) ? 1 : 0;

    w.turnSin = sigma * (tPrev[0]*tNext[1] - tPrev[1]*tNext[0]) / (LP * LN);
    if (w.turnSin < 0.0)      w.corner = C2D_CONVEX;
    else if (w.turnSin > 0.0) w.corner = C2D_CONCAVE;
    else                      w.corner = C2D_STRAIGHT;

    const double nP[2] = { sigma*(-tPrev[1])/LP, sigma*(tPrev[0])/LP };
    const double nN[2] = { sigma*(-tNext[1])/LN, sigma*(tNext[0])/LN };
    const double bis[2] = { nP[0] + nN[0], nP[1] + nN[1] };
    const double bl = norm2(bis);
    const double rl = norm2(r);
    // a fold-back corner (antiparallel normals, vanishing bisector) or a slave
    // within the gate of the wedge boundary leaves sideSign = 0
    if (bl > 2.0*tauPerp && rl > 0.0) {
        const double p = dot2(r, bis);
        if (std::fabs(p) > tauPerp * rl) w.sideSign = (p > 0.0) ? 1.0 : -1.0;
    }
    return w;
}

// Vertex pair: n = sigmaSide (xs - Xv)/|xs - Xv|, gap = sigmaSide |xs - Xv|,
// B = [n | -n]. |xs - Xv| is a distance, so the sign must come from the side
// classification; it is committed at the first capture and passed here. On
// the wedge boundaries the pair equals the adjacent segment at its end.
// tolLen (> 0, relative to the surface scale) refuses a slave on the vertex.
inline int vertexEval2D(const double Xv[2], const double xs[2], double sigmaSide,
                        double &gap, double n[2], double tolLen)
{
    gap = 0.0; n[0] = 0.0; n[1] = 0.0;

    const double r[2] = { xs[0] - Xv[0], xs[1] - Xv[1] };
    const double d = norm2(r);
    if (!(tolLen > 0.0)) return VTX2D_DEGENERATE;
    if (d < tolLen)      return VTX2D_DEGENERATE;

    n[0] = sigmaSide * r[0] / d;
    n[1] = sigmaSide * r[1] / d;
    gap  = sigmaSide * d;
    return VTX2D_OK;
}

// The sign that turns the segment normal perp(t)/L toward refDir:
// sign(perp(t).refDir / L), and 0 when refDir is numerically tangent to the
// segment (|proj| < tauPerp |refDir|), where the caller must refuse rather
// than guess.
inline double sigmaFromRef2D(const double t[2], const double refDir[2],
                             double tauPerp = TAU_PERP_DEFAULT)
{
    const double L = norm2(t);
    const double R = norm2(refDir);
    if (!(L > 0.0) || !(R > 0.0)) return 0.0;
    double p[2];
    perp2(t, p);
    const double proj = dot2(p, refDir) / L;
    if (std::fabs(proj) < tauPerp * R) return 0.0;
    return (proj > 0.0) ? 1.0 : -1.0;
}

// gap operator of a segment pair over [slave xy | X0 xy | X1 xy]
inline void bOperatorSegment2D(double N0, double N1, const double n[2], double B[6])
{
    B[0] =  n[0];      B[1] =  n[1];
    B[2] = -N0*n[0];   B[3] = -N0*n[1];
    B[4] = -N1*n[0];   B[5] = -N1*n[1];
}

// gap operator of a vertex pair over [slave xy | vertex xy]
inline void bOperatorVertex2D(const double n[2], double B[4])
{
    B[0] =  n[0];   B[1] =  n[1];
    B[2] = -n[0];   B[3] = -n[1];
}

// ---------------------------------------------------------------- 2D mortar
// Segment-to-segment mortar of two 2-node segments: the slave end points are
// projected on the affine parametrization of the master, the interval
// [xi_m(s0), xi_m(s1)] is clipped against [0, 1] and integrated with a 2-point
// Gauss rule on the slave-trace measure dGamma_s = (L_m / |t_s.t_m|) dxi_m
// (the 2D form of the auxiliary-plane area over the cosine of the 3D kernel).
// The gap is evaluated inside the integral.
enum Mort2DStatus { MORT2D_OK = 0, MORT2D_SKIP = -1 };
enum Mort2DRefusal {
    MORT2D_REF_NONE  = 0,
    MORT2D_REF_EMPTY = 1,   // clipped overlap below tauOv min(Ls, Lm)/Lm (sliver or disjoint)
    MORT2D_REF_DEGEN = 2,   // a zero-length segment (tauSeg Lref)
    MORT2D_REF_ALIGN = 3    // |t_s.t_m| < tauAlign: nearly perpendicular segments,
                            // where the Jacobian 1/|cos| is ill-conditioned
};

// Both relative. tauOv must stay far below tauAlign: for a fully contained
// slave the overlap is Ls |cos| / Lm against tauOv Ls / Lm, so the sliver
// test fires exactly when |cos| < tauOv, and the alignment test must classify
// such a pair first (ALIGN is tested before EMPTY).
static const double MORT2D_TAU_OV_DEFAULT    = 1e-9;
static const double MORT2D_TAU_ALIGN_DEFAULT = 1e-6;

struct PairResult2D {
    int    status;    // MORT2D_OK or MORT2D_SKIP
    int    refusal;   // Mort2DRefusal (MORT2D_REF_NONE when status is MORT2D_OK)
    double cos_t;     // t_s.t_m (unit tangents), set before the ALIGN test
    int    ngp;       // Gauss points integrated (0 or 2)
    double len;       // integral of dGamma_s over the overlap
    double gapL1;     // integral of |gN| dGamma_s
    double D[2][2];   // D_ij = int N_i^s N_j^s dGamma_s
    double M[2][2];   // M_ik = int N_i^s N_k^m dGamma_s
    double g[2];      // weighted gap g_i = int N_i^s gN dGamma_s
};

// Xs, Xm: segment end points (the caller's configuration); sigma: the
// orientation of the interface; Lref: the minimum initial segment length.
// Refusal order DEGEN, ALIGN, EMPTY.
inline int mortarSegmentPair2D(const double Xs[2][2], const double Xm[2][2],
                               double sigma, double Lref, PairResult2D &out,
                               double tauOv    = MORT2D_TAU_OV_DEFAULT,
                               double tauAlign = MORT2D_TAU_ALIGN_DEFAULT,
                               double tauSeg   = TAU_SEG_DEFAULT)
{
    out.status = MORT2D_SKIP; out.refusal = MORT2D_REF_EMPTY; out.cos_t = 0.0;
    out.ngp = 0; out.len = 0.0; out.gapL1 = 0.0;
    for (int i = 0; i < 2; i++) {
        out.g[i] = 0.0;
        for (int j = 0; j < 2; j++) { out.D[i][j] = 0.0; out.M[i][j] = 0.0; }
    }

    const double ts[2] = { Xs[1][0] - Xs[0][0], Xs[1][1] - Xs[0][1] };
    const double tm[2] = { Xm[1][0] - Xm[0][0], Xm[1][1] - Xm[0][1] };
    const double Ls = norm2(ts), Lm = norm2(tm);
    if (!(Lref > 0.0))               { out.refusal = MORT2D_REF_DEGEN; return out.status; }
    if (Ls <= tauSeg * Lref || Lm <= tauSeg * Lref) {
        out.refusal = MORT2D_REF_DEGEN; return out.status;
    }

    const double tsHat[2] = { ts[0] / Ls, ts[1] / Ls };
    const double tmHat[2] = { tm[0] / Lm, tm[1] / Lm };
    const double cosT = dot2(tsHat, tmHat);
    out.cos_t = cosT;
    if (std::fabs(cosT) < tauAlign) {
        out.refusal = MORT2D_REF_ALIGN;
        return out.status;
    }

    // slave end points on the master parametrization, clipped to [0, 1]
    const double Lm2 = dot2(tm, tm);
    const double d0[2] = { Xs[0][0] - Xm[0][0], Xs[0][1] - Xm[0][1] };
    const double d1[2] = { Xs[1][0] - Xm[0][0], Xs[1][1] - Xm[0][1] };
    const double al0 = dot2(d0, tm) / Lm2;
    const double al1 = dot2(d1, tm) / Lm2;
    const double lo = (al0 < al1) ? al0 : al1;
    const double hi = (al0 < al1) ? al1 : al0;
    const double a  = (lo > 0.0) ? lo : 0.0;
    const double b  = (hi < 1.0) ? hi : 1.0;
    const double thr = tauOv * ((Ls < Lm) ? Ls : Lm) / Lm;
    if (b - a <= thr) {
        out.refusal = MORT2D_REF_EMPTY;
        return out.status;
    }

    const double nx = -sigma * tm[1] / Lm;
    const double ny =  sigma * tm[0] / Lm;

    const double J = Lm / std::fabs(cosT);
    const double dal = al1 - al0;
    // 2-point Gauss-Legendre on [0, 1]
    static const double GP2_U0 = 0.21132486540518711775;
    static const double GP2_U1 = 0.78867513459481288225;
    static const double GP2_W  = 0.5;
    const double u[2] = { GP2_U0, GP2_U1 };
    const double span = b - a;

    for (int gp = 0; gp < 2; gp++) {
        const double xim = a + span * u[gp];
        const double wgt = GP2_W * span;
        const double xis = (xim - al0) / dal;
        const double Ns[2] = { 1.0 - xis, xis };
        const double Nm[2] = { 1.0 - xim, xim };
        const double xsx = Xs[0][0] + xis * ts[0];
        const double xsy = Xs[0][1] + xis * ts[1];
        const double xmx = Xm[0][0] + xim * tm[0];
        const double xmy = Xm[0][1] + xim * tm[1];
        const double gN = (xsx - xmx) * nx + (xsy - xmy) * ny;
        const double wJ = wgt * J;
        out.len   += wJ;
        out.gapL1 += wJ * std::fabs(gN);
        for (int i = 0; i < 2; i++) {
            out.g[i] += wJ * Ns[i] * gN;
            for (int j = 0; j < 2; j++) out.D[i][j] += wJ * Ns[i] * Ns[j];
            for (int k = 0; k < 2; k++) out.M[i][k] += wJ * Ns[i] * Nm[k];
        }
        out.ngp++;
    }
    out.status = MORT2D_OK;
    out.refusal = MORT2D_REF_NONE;
    return out.status;
}

} // namespace LadrunoContact2DKernel

#endif
