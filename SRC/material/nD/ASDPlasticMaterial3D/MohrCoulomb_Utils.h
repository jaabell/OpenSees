/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
** ****************************************************************** */

// Mohr-Coulomb function shared by MohrCoulomb_YF (with phi) and
// MohrCoulomb_PF (with psi):
//     F = K(theta) sqrt(J2) + I1 sin(phi) / 3 - c cos(phi)
//     K(theta) = cos(theta) - sin(theta) sin(phi) / sqrt(3)       |theta| <= theta_T
//     K(theta) = A - B sin(3 theta)                                |theta| >  theta_T
// with the Lode angle theta in [-30, 30] deg (sin 3theta = -3 sqrt(3) J3 / (2 J2^1.5)).
// The exact surface has corners at theta = +-30 deg (the triaxial compression and
// extension meridians), where its gradient does not exist. Above theta_T the
// corners are rounded as in Abbo & Sloan (1995), "A smooth hyperbolic
// approximation to the Mohr-Coulomb yield criterion", Computers & Structures
// 54(3):427-441: A and B make K and dK/dtheta continuous at theta_T, and dK/dtheta
// vanishes at +-30 deg, so the gradient is smooth everywhere. The rounded surface
// lies inside the exact one; at phi = 30 deg, K grows on the compression
// meridian by about 0.1% for theta_T = 29.9 deg (0.45% for 29.5, 1.8% for 28).
//
// theta_T_deg: <= 0 selects the default 29.9 deg, >= 30 the exact (sharp) surface.

#ifndef MohrCoulomb_Utils_H
#define MohrCoulomb_Utils_H

#include <cmath>
#include "EigenAPI.h"

namespace MohrCoulombUtils
{

constexpr double DEFAULT_THETA_T_DEG = 29.9;

// A and B of the rounded part, K = A - B sin(3 theta), for the sign sg of
// theta. They depend only on theta_T and phi, so the last pair is kept (the
// same material calls this many times per integration).
inline void rounding_AB(double theta_T_deg, double sin_phi, double sg, double &A, double &B)
{
    static double last_tT = -1.0, last_sp = 2.0, Ap = 0, Bp = 0, Am = 0, Bm = 0;
    if (theta_T_deg != last_tT || sin_phi != last_sp) {
        const double tT = theta_T_deg * M_PI / 180.0;
        for (int k = 0; k < 2; k++) {
            const double s = k == 0 ? 1.0 : -1.0;
            const double a = std::cos(tT) / 3.0 * (3.0 + std::tan(tT) * std::tan(3.0 * tT)
                             + s / std::sqrt(3.0) * (std::tan(3.0 * tT) - 3.0 * std::tan(tT)) * sin_phi);
            const double b = 1.0 / (3.0 * std::cos(3.0 * tT)) * (s * std::sin(tT) + sin_phi * std::cos(tT) / std::sqrt(3.0));
            if (k == 0) { Ap = a; Bp = b; } else { Am = a; Bm = b; }
        }
        last_tT = theta_T_deg;
        last_sp = sin_phi;
    }
    if (sg > 0) { A = Ap; B = Bp; } else { A = Am; B = Bm; }
}

inline double lode_K(double theta, double sin_phi, double theta_T_deg)
{
    const double deg = M_PI / 180.0;
    if (!(theta_T_deg > 0.0))
        theta_T_deg = DEFAULT_THETA_T_DEG;
    if (theta_T_deg >= 30.0 || std::abs(theta) <= theta_T_deg * deg)
        return std::cos(theta) - std::sin(theta) * sin_phi / std::sqrt(3.0);

    double A, B;
    rounding_AB(theta_T_deg, sin_phi, theta > 0 ? 1.0 : -1.0, A, B);
    return A - B * std::sin(3.0 * theta);
}

// F(sigma) for friction (or dilatancy) angle phi_rad and cohesion c
inline double function(const VoigtVector& sigma, double phi_rad, double c, double theta_T_deg)
{
    const double sin_phi = std::sin(phi_rad);
    const double J2 = sigma.getJ2();
    const double K = lode_K(sigma.lodeAngle(), sin_phi, theta_T_deg);
    return K * std::sqrt(J2) + sigma.getI1() * sin_phi / 3.0 - c * std::cos(phi_rad);
}

// dF/dsigma in closed form, in tensor components (as gradient() below), by the
// chain rule through I1, J2 and the Lode angle (Sloan & Booker 1986):
//     dF/dsigma = sin(phi)/3 delta + K/(2 sqrt(J2)) s + sqrt(J2) K'(theta) dtheta/dsigma
//     dtheta/dsigma = -sqrt(3) / (2 cos(3 theta) J2^1.5) (dJ3/dsigma - 1.5 J3/J2 s)
//     dJ3/dsigma = s.s - 2/3 J2 delta
// On the rounded part K' = -3 B cos(3 theta), so cos(3 theta) cancels.
// Returns false where the gradient is not defined: at the cone axis (J2 = 0)
// and on the sharp surface (theta_T >= 30) at the corners. The caller then
// falls back to central differences.
inline bool gradient_closed_form(const VoigtVector& sigma, double phi_rad, double theta_T_deg,
                                 VoigtVector& out)
{
    const double deg = M_PI / 180.0;
    const double sqrt3 = std::sqrt(3.0);
    const double sin_phi = std::sin(phi_rad);

    const double p = sigma.getI1() / 3.0;
    const double s11 = sigma(0) - p, s22 = sigma(1) - p, s33 = sigma(2) - p;
    const double s12 = sigma(3), s23 = sigma(4), s13 = sigma(5);

    const double J2 = 0.5 * (s11 * s11 + s22 * s22 + s33 * s33) + s12 * s12 + s23 * s23 + s13 * s13;
    const double scale = sigma(0) * sigma(0) + sigma(1) * sigma(1) + sigma(2) * sigma(2)
                         + 2.0 * (s12 * s12 + s23 * s23 + s13 * s13);
    if (!(J2 > 1e-20 * scale) || !(J2 > 0.0))
        return false;
    const double J3 = s11 * s22 * s33 - s11 * s23 * s23 - s22 * s13 * s13 - s33 * s12 * s12
                      + 2.0 * s12 * s23 * s13;
    const double sqJ2 = std::sqrt(J2);

    double x = -3.0 * sqrt3 * J3 / (2.0 * J2 * sqJ2);
    if (x > 1.0) x = 1.0;
    else if (x < -1.0) x = -1.0;
    const double theta = std::asin(x) / 3.0;

    if (!(theta_T_deg > 0.0))
        theta_T_deg = DEFAULT_THETA_T_DEG;

    // K and K'(theta) / cos(3 theta)
    double K, Kp_over_c3;
    if (theta_T_deg >= 30.0 || std::abs(theta) <= theta_T_deg * deg) {
        const double ct = std::cos(theta), st = std::sin(theta);
        K = ct - st * sin_phi / sqrt3;
        const double c3 = std::sqrt(std::max(0.0, 1.0 - x * x));   // cos(3 theta) >= 0
        if (!(c3 > 1e-8))
            return false;
        Kp_over_c3 = (-st - ct * sin_phi / sqrt3) / c3;
    } else {
        double A, B;
        rounding_AB(theta_T_deg, sin_phi, theta > 0 ? 1.0 : -1.0, A, B);
        K = A - B * x;                       // sin(3 theta) = x
        Kp_over_c3 = -3.0 * B;
    }

    // s.s
    const double t11 = s11 * s11 + s12 * s12 + s13 * s13;
    const double t22 = s12 * s12 + s22 * s22 + s23 * s23;
    const double t33 = s13 * s13 + s23 * s23 + s33 * s33;
    const double t12 = s11 * s12 + s12 * s22 + s13 * s23;
    const double t23 = s12 * s13 + s22 * s23 + s23 * s33;
    const double t13 = s11 * s13 + s12 * s23 + s13 * s33;

    const double a_d = sin_phi / 3.0;                       // on delta
    const double a_s = K / (2.0 * sqJ2);                    // on s
    const double a_t = -0.5 * sqrt3 * Kp_over_c3 / J2;      // on (dJ3 - 1.5 J3/J2 s)
    const double r = 1.5 * J3 / J2;
    const double d2 = 2.0 / 3.0 * J2;

    out(0) = a_d + a_s * s11 + a_t * (t11 - d2 - r * s11);
    out(1) = a_d + a_s * s22 + a_t * (t22 - d2 - r * s22);
    out(2) = a_d + a_s * s33 + a_t * (t33 - d2 - r * s33);
    out(3) = a_s * s12 + a_t * (t12 - r * s12);
    out(4) = a_s * s23 + a_t * (t23 - r * s23);
    out(5) = a_s * s13 + a_t * (t13 - r * s13);
    return true;
}

// dF/dsigma by central differences, in tensor components (see strain_like):
// perturbing storage entry i >= 3 moves both sigma_ij and sigma_ji, so it is halved.
// ds: perturbation relative to |sigma| (default 1e-6 when <= 0)
template <class F>
inline VoigtVector gradient(const VoigtVector& sigma, double ds, F && fun)
{
    if (!(ds > 0.0))
        ds = 1e-6;
    const double h = std::max(ds, ds * sigma.norm());
    VoigtVector result;
    for (int i = 0; i < 6; ++i) {
        VoigtVector s1 = sigma, s2 = sigma;
        s1(i) += h;
        s2(i) -= h;
        result(i) = (fun(s1) - fun(s2)) / (2.0 * h);
        if (i >= 3) result(i) *= 0.5;
    }
    return result;
}

}

#endif
