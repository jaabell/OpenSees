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

// Hoek-Brown shear function shared by HoekBrown_YF (with mb) and HoekBrown_PF
// (with mb_psi), written in invariants so that its edges can be rounded:
//     F = (sigma1 - sigma3) - sigma_ci * max(mb sigma3 / sigma_ci + s, 0)^a
// with sigma1 >= sigma3 the major/minor principal stresses, compression positive.
// For a tension-positive stress with mean sigma_m, J2 and Lode angle theta in
// [-30, 30] deg (sin 3theta = -3 sqrt(3) J3 / (2 J2^1.5)):
//     sigma1 - sigma3 = 2 sqrt(J2) cos(theta)
//     sigma3          = -sigma_m - 2/sqrt(3) sqrt(J2) sin(theta + 2 pi/3)
// Both trigonometric factors have a nonzero slope at theta = +-30 deg (the
// triaxial compression and extension meridians), which are the edges of the
// surface. Above theta_T each factor is replaced by A - B sin(3 theta), continuous
// with its slope at theta_T and flat at +-30 deg (the construction of Abbo &
// Sloan 1995 for Mohr-Coulomb, applied to each factor). The rounded surface lies
// slightly inside the exact one. Near zero confinement the strength is very
// sensitive to sigma3 (s is small), so theta_T must be close to 30 deg: for
// GSI 50, mi 10 and theta_T = 29.9 deg the uniaxial strength drops 0.7%, the
// triaxial compression strength at 5 MPa confinement 0.15% and the extension
// strength 0.05% (6.9%, 1.5% and 0.5% for theta_T = 29 deg).
//
// theta_T_deg: <= 0 selects the default 29.9 deg, >= 30 the exact (sharp) surface.

#ifndef HoekBrown_Invariants_H
#define HoekBrown_Invariants_H

#include <cmath>
#include <algorithm>
#include "EigenAPI.h"

namespace HoekBrownInvariants
{

constexpr double DEFAULT_THETA_T_DEG = 29.9;

// value of f(theta), replaced by A - B sin(3 theta) above theta_T (f, df: the function and its slope)
template <class Fn, class DFn>
inline double rounded(double theta, double theta_T_deg, Fn && f, DFn && df)
{
    const double deg = M_PI / 180.0;
    if (!(theta_T_deg > 0.0))
        theta_T_deg = DEFAULT_THETA_T_DEG;
    if (theta_T_deg >= 30.0 || std::abs(theta) <= theta_T_deg * deg)
        return f(theta);
    const double t = (theta > 0 ? 1.0 : -1.0) * theta_T_deg * deg;
    const double B = -df(t) / (3.0 * std::cos(3.0 * t));
    const double A = f(t) + B * std::sin(3.0 * t);
    return A - B * std::sin(3.0 * theta);
}

// F for a tension-positive stress
inline double shear_function(const VoigtVector& sigma, double sigma_ci, double mb, double s, double a, double theta_T_deg)
{
    const double sm = sigma.getI1() / 3.0;
    const double r = std::sqrt(std::max(sigma.getJ2(), 0.0));
    const double theta = sigma.lodeAngle();
    const double c = rounded(theta, theta_T_deg,
        [](double x) { return std::cos(x); }, [](double x) { return -std::sin(x); });
    const double t = rounded(theta, theta_T_deg,
        [](double x) { return std::sin(x + 2.0 * M_PI / 3.0); }, [](double x) { return std::cos(x + 2.0 * M_PI / 3.0); });
    const double sigma1_minus_sigma3 = 2.0 * r * c;
    const double sigma3 = -sm - 2.0 / std::sqrt(3.0) * r * t;
    const double arg = std::max(mb * sigma3 / sigma_ci + s, 0.0);
    return sigma1_minus_sigma3 - sigma_ci * std::pow(arg, a);
}

// dF/dsigma by central differences, in tensor components (see strain_like):
// perturbing storage entry i >= 3 moves both sigma_ij and sigma_ji, so it is halved.
template <class F>
inline VoigtVector gradient(const VoigtVector& sigma, double ds, F && fun)
{
    if (!(ds > 0.0))
        ds = 1e-8;
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
