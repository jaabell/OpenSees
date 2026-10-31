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
//     Berlin. doi:10.1007/978-3-540-32609-0 (penalty regularization of
//     Coulomb friction, return mapping and its consistent linearization).
//   Laursen, T.A. (2002). Computational Contact and Impact Mechanics.
//     Springer, Berlin. doi:10.1007/978-3-662-04864-1
//   Giannakopoulos, A.E. (1989). The return mapping method for the
//     integration of friction constitutive relations. Computers & Structures,
//     32(1), 157-167.
//   Dassault Systemes Simulia. Abaqus Theory Guide, Section 5.2.3 (frictional
//     behavior: Coulomb friction with a shear stress limit).
//
// Description: LadrunoFrictionKernel, header-only Coulomb/Tresca friction for
// penalty contact: the return map of the tangential traction and its
// consistent tangent. It uses no OpenSees types (raw double arrays and
// <cmath> only).
//
// Friction cone with an optional cohesion and shear cap:
//   cap = min( mu N + c , tauMax )      N = contact pressure >= 0, c = cohesion
//   - Coulomb:             c = 0, tauMax <= 0 (no cap)   -> cap = mu N
//   - Coulomb + cohesion:  c > 0                         -> cap = mu N + c
//   - Coulomb + shear cap: tauMax > 0                    -> cap = min(mu N + c, tauMax)
//   - Tresca:              mu = 0, c = tau, tauMax <= 0  -> cap = c
// The cone is zero when the pair is separated (N <= 0). The pressure coupling
// of the tangent, d t_T / d g_N, is non-zero only where the cap depends on N:
// mu on the Coulomb branch, 0 when the tauMax cap binds or when separated.

#ifndef LadrunoFrictionKernel_h
#define LadrunoFrictionKernel_h

#include <cmath>

namespace LadrunoFrictionKernel {

inline double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline double norm3(const double a[3]) { return std::sqrt(dot3(a, a)); }

// tangential component of v with respect to the unit normal n: v - (v.n) n
inline void tangentPart(const double v[3], const double n[3], double out[3]) {
    double vn = dot3(v, n);
    for (int d = 0; d < 3; d++) out[d] = v[d] - vn * n[d];
}

// critical shear cap of the cone: capC = mu N + c in contact (N > 0), 0 when
// separated; cap = min(capC, tauMax) when tauMax > 0
inline double frictionCap(double N, double mu, double cohesion, double tauMax) {
    double capC = (N > 0.0) ? (mu * N + cohesion) : 0.0;
    return (tauMax > 0.0 && tauMax < capC) ? tauMax : capC;
}

// Return map of the tangential traction (a one-dimensional plasticity problem
// with yield radius cap):
//
//   tT*  = kt (gTeff - gpT)        trial elastic tangential traction
//   stick (|tT*| <= cap): tT = tT*,  gpT unchanged
//   slip  (|tT*| >  cap): nh = tT*/|tT*|, dlambda = (|tT*| - cap)/kt,
//                         tT = cap nh,  gpT += dlambda nh
//
// tT* follows the slave's tangential motion, so the force on the slave is
// returned already negated, tFric = -tT, and opposes the motion. gpTtrial is a
// function of the committed state only (assigned, not accumulated), so a
// repeated residual evaluation gives the same result. Returns true when
// slipping.
inline bool frictionReturnMap(const double gTeff[3], const double gpT[3],
                              double N, double kt, double mu,
                              double tFric[3], double gpTtrial[3],
                              double cohesion = 0.0, double tauMax = 0.0) {
    double tTtr[3];
    for (int d = 0; d < 3; d++) tTtr[d] = kt * (gTeff[d] - gpT[d]);
    double cap  = frictionCap(N, mu, cohesion, tauMax);
    double norm = norm3(tTtr);
    if (cap <= 0.0) {
        // zero cone radius: free slip, zero tangential traction, the slip takes
        // the whole tangential motion (never the elastic trial traction)
        for (int d = 0; d < 3; d++) { tFric[d] = 0.0; gpTtrial[d] = gTeff[d]; }
        return true;
    }
    if (norm <= cap) {                               // stick
        for (int d = 0; d < 3; d++) { tFric[d] = -tTtr[d]; gpTtrial[d] = gpT[d]; }
        return false;
    }
    double inv  = 1.0 / norm;                        // norm > cap >= 0
    double dlam = (norm - cap) / kt;
    for (int d = 0; d < 3; d++) {
        double nh   = tTtr[d] * inv;
        tFric[d]    = -cap * nh;
        gpTtrial[d] = gpT[d] + dlam * nh;
    }
    return true;
}

// Slave block of the friction tangent, K_ss = D_TT P_t + d_TN (x) n (3 x 3),
// for the assembled implicit tangent K_fric = G^T K_ss G, G = [I | -N_i I].
//   stick: K_ss = kt P_t, P_t = I - n (x) n
//   slip:  K_ss = (cap kt / |tT*|) (P_t - nh (x) nh)
// consistent = false drops the coupling d_TN (x) n and keeps the tangent
// symmetric, which is safe with any linear solver (a symmetric solver reads
// one triangle only and would silently lose a non-symmetric term). consistent
// = true adds d_TN = -(d cap/d N) kn nh, the full consistent tangent
// (quadratic convergence; needs a non-symmetric solver such as FullGeneral,
// UmfPack or BandGeneral). d cap/d N = mu on the Coulomb branch, 0 when the
// tauMax cap binds or when separated.
inline void frictionTangentBlock(const double gTeff[3], const double gpT[3],
                                 const double n[3], double N, double kn, double kt,
                                 double mu, bool consistent, double Kss[3][3],
                                 double cohesion = 0.0, double tauMax = 0.0) {
    double tTtr[3];
    for (int d = 0; d < 3; d++) tTtr[d] = kt * (gTeff[d] - gpT[d]);
    double capC   = (N > 0.0) ? (mu * N + cohesion) : 0.0;
    bool   capped = (tauMax > 0.0 && tauMax < capC);   // tauMax binds: d cap/d N = 0
    double cap    = capped ? tauMax : capC;
    double nrm = norm3(tTtr);
    double Pt[3][3];                                 // tangent-plane projector I - n (x) n
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) Pt[i][j] = (i == j ? 1.0 : 0.0) - n[i]*n[j];
    if (cap <= 0.0) {                                // free slip: K_ss = 0, as in the return map
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] = 0.0;
        return;
    }
    if (nrm <= cap) {                                // stick: K_ss = kt P_t
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] = kt * Pt[i][j];
        return;
    }
    double inv = 1.0 / nrm, nh[3];                   // slip
    for (int d = 0; d < 3; d++) nh[d] = tTtr[d] * inv;
    double s = cap * kt * inv;                        // cap kt / |tT*|
    for (int i = 0; i < 3; i++)                       // (cap kt / |tT*|)(P_t - nh (x) nh)
        for (int j = 0; j < 3; j++) Kss[i][j] = s * (Pt[i][j] - nh[i]*nh[j]);
    if (consistent) {                                // + d_TN (x) n, d_TN = -(d cap/d N) kn nh
        double dCap_dN = (capped || N <= 0.0) ? 0.0 : mu;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] += (-dCap_dN * kn * nh[i]) * n[j];
    }
}

// Rounded stick/slip corner (mortar -smoothT r). The traction magnitude of the
// return map above is min(rho, cap), rho = |tT*|, with a kink at rho = cap that
// a moving slip front puts into the residual at every node it crosses. Here the
// corner is rounded over |rho - cap| < delta = r cap (0 < r < 1) by a C1
// quadratic blend:
//   rho <= cap - delta: phi = rho                            (stick)
//   rho >= cap + delta: phi = cap                            (slip)
//   otherwise:          phi = rho - (rho - cap + delta)^2 / (4 delta)
// T = phi nh, and the slip that keeps T = kt (gTeff - gpTtrial) is
// gpTtrial = gpT + (rho - phi)/kt nh, a function of the committed state only.
// r -> 0 recovers the return map above. dphi_dcap is the total derivative
// (delta moves with cap), needed by the pressure coupling.
inline void smoothCorner(double rho, double cap, double r,
                         double &phi, double &dphi_drho, double &dphi_dcap) {
    double delta = r * cap;
    if (rho <= cap - delta) { phi = rho; dphi_drho = 1.0; dphi_dcap = 0.0; return; }
    if (rho >= cap + delta) { phi = cap; dphi_drho = 0.0; dphi_dcap = 1.0; return; }
    double e = rho - cap + delta;                    // in (0, 2 delta)
    phi = rho - e * e / (4.0 * delta);
    dphi_drho = 1.0 - e / (2.0 * delta);
    dphi_dcap = e / (2.0 * delta) + r * (e * e / (4.0 * delta * delta) - e / (2.0 * delta));
}

inline bool frictionReturnMapSmooth(const double gTeff[3], const double gpT[3],
                                    double N, double kt, double mu,
                                    double tFric[3], double gpTtrial[3],
                                    double cohesion, double tauMax, double r) {
    double tTtr[3];
    for (int d = 0; d < 3; d++) tTtr[d] = kt * (gTeff[d] - gpT[d]);
    double cap = frictionCap(N, mu, cohesion, tauMax);
    double rho = norm3(tTtr);
    if (cap <= 0.0) {                                // free slip, as the return map above
        for (int d = 0; d < 3; d++) { tFric[d] = 0.0; gpTtrial[d] = gTeff[d]; }
        return true;
    }
    double phi, dr, dc;
    smoothCorner(rho, cap, r, phi, dr, dc);
    if (rho <= cap - r * cap) {                      // stick (also rho == 0)
        for (int d = 0; d < 3; d++) { tFric[d] = -tTtr[d]; gpTtrial[d] = gpT[d]; }
        return false;
    }
    double inv = 1.0 / rho, dlam = (rho - phi) / kt;
    for (int d = 0; d < 3; d++) {
        double nh   = tTtr[d] * inv;
        tFric[d]    = -phi * nh;
        gpTtrial[d] = gpT[d] + dlam * nh;
    }
    return true;
}

// K_ss = dT/dgTeff = kt [phi'(rho) nh (x) nh + (phi/rho)(P_t - nh (x) nh)];
// consistent adds -(dphi/dcap)(dcap/dN) kn nh (x) n (kn = dN/dz, the normal
// slope of the caller).
inline void frictionTangentBlockSmooth(const double gTeff[3], const double gpT[3],
                                       const double n[3], double N, double kn, double kt,
                                       double mu, bool consistent, double Kss[3][3],
                                       double cohesion, double tauMax, double r) {
    double tTtr[3];
    for (int d = 0; d < 3; d++) tTtr[d] = kt * (gTeff[d] - gpT[d]);
    double capC   = (N > 0.0) ? (mu * N + cohesion) : 0.0;
    bool   capped = (tauMax > 0.0 && tauMax < capC);
    double cap    = capped ? tauMax : capC;
    double rho = norm3(tTtr);
    double Pt[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) Pt[i][j] = (i == j ? 1.0 : 0.0) - n[i]*n[j];
    if (cap <= 0.0) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] = 0.0;
        return;
    }
    double phi, dr, dc;
    smoothCorner(rho, cap, r, phi, dr, dc);
    if (rho <= cap - r * cap) {                      // stick
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] = kt * Pt[i][j];
        return;
    }
    double inv = 1.0 / rho, nh[3];
    for (int d = 0; d < 3; d++) nh[d] = tTtr[d] * inv;
    double s = phi * kt * inv;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Kss[i][j] = s * (Pt[i][j] - nh[i]*nh[j]) + kt * dr * nh[i]*nh[j];
    if (consistent) {
        double dCap_dN = (capped || N <= 0.0) ? 0.0 : mu;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Kss[i][j] += (-dc * dCap_dN * kn * nh[i]) * n[j];
    }
}

// --- scalar forms for 2D contact. The tangent space of a 2D contact is one
//     dimensional, t = perp(n), so the slip, the plastic slip and the traction
//     are scalars along t. They are kept as scalars (not as 2-vectors
//     differenced componentwise), so that a rotation of the contact normal
//     does not turn a normal displacement into a spurious slip.

inline double sgn1D(double x) { return (x < 0.0) ? -1.0 : 1.0; }

// Scalar return map: tT* = kt (s - sp); stick |tT*| <= cap: tFric = -tT*;
// slip: tFric = -cap sign(tT*), sp += (|tT*| - cap)/kt sign(tT*). A zero cap
// (separated, no cohesion) gives free slip. Returns true when slipping.
inline bool returnMap1D(double s, double sp, double N, double kt, double mu,
                        double &tFric, double &spTrial,
                        double cohesion = 0.0, double tauMax = 0.0)
{
    const double tTtr = kt * (s - sp);
    const double cap  = frictionCap(N, mu, cohesion, tauMax);
    const double atr  = std::fabs(tTtr);
    if (cap <= 0.0) {
        tFric = 0.0; spTrial = s;
        return true;
    }
    if (atr <= cap) {
        tFric = -tTtr; spTrial = sp;
        return false;
    }
    const double sg   = sgn1D(tTtr);
    const double dlam = (atr - cap) / kt;
    tFric   = -cap * sg;
    spTrial = sp + dlam * sg;
    return true;
}

// Scalar tangent of returnMap1D: Kss = d(-tFric)/ds (kt in stick, 0 in slip)
// and, with consistent, dTN = d(-tFric)/d(normal gap) = -(dcap/dN) kn sign(tT*)
// in slip (the pressure coupling); kn is already included in dTN.
inline void tangentBlock1D(double s, double sp, double N, double kn, double kt, double mu,
                           bool consistent, double &Kss, double &dTN,
                           double cohesion = 0.0, double tauMax = 0.0)
{
    const double tTtr   = kt * (s - sp);
    const double capC   = (N > 0.0) ? (mu * N + cohesion) : 0.0;
    const bool   capped = (tauMax > 0.0 && tauMax < capC);
    const double cap    = capped ? tauMax : capC;
    const double atr    = std::fabs(tTtr);
    dTN = 0.0;
    if (cap <= 0.0) {
        Kss = 0.0;
        return;
    }
    if (atr <= cap) {
        Kss = kt;
        return;
    }
    Kss = 0.0;
    if (consistent) {
        const double dCap_dN = (capped || N <= 0.0) ? 0.0 : mu;
        const double sg = sgn1D(tTtr);
        dTN = -dCap_dN * kn * sg;
    }
}

} // namespace LadrunoFrictionKernel

#endif
