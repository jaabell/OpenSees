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

} // namespace LadrunoFrictionKernel

#endif
