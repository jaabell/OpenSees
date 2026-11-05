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

inline double lode_K(double theta, double sin_phi, double theta_T_deg)
{
    const double deg = M_PI / 180.0;
    if (!(theta_T_deg > 0.0))
        theta_T_deg = DEFAULT_THETA_T_DEG;
    if (theta_T_deg >= 30.0 || std::abs(theta) <= theta_T_deg * deg)
        return std::cos(theta) - std::sin(theta) * sin_phi / std::sqrt(3.0);

    const double tT = theta_T_deg * deg;
    const double sg = theta > 0 ? 1.0 : -1.0;
    const double A = std::cos(tT) / 3.0 * (3.0 + std::tan(tT) * std::tan(3.0 * tT)
                     + sg / std::sqrt(3.0) * (std::tan(3.0 * tT) - 3.0 * std::tan(tT)) * sin_phi);
    const double B = 1.0 / (3.0 * std::cos(3.0 * tT)) * (sg * std::sin(tT) + sin_phi * std::cos(tT) / std::sqrt(3.0));
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
