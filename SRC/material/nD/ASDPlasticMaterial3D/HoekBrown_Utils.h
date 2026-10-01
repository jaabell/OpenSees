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

// Implementation: José Abell (UANDES), Massimo Petracca (ASDEA)
//
// HoekBrown_Utils.h
//
// Utility functions for computing Hoek-Brown rock mass parameters
// Based on Hoek & Brown (2018) - "The Hoek-Brown failure criterion and GSI - 2018 edition"

#ifndef HoekBrown_Utils_H
#define HoekBrown_Utils_H

#include <cmath>
#include <tuple>
#include <stdexcept>

namespace HoekBrownUtils {

/**
 * Compute the reduced material constant mb for rock mass
 * 
 * mb = mi * exp((GSI - 100) / (28 - 14*D))
 * 
 * @param mi   Material constant for intact rock (dimensionless, typically 4-35)
 * @param GSI  Geological Strength Index (0-100)
 * @param D    Disturbance factor (0-1)
 * @return     Reduced material constant mb
 */
inline double compute_mb(double mi, double GSI, double D) {
    if (GSI < 0 || GSI > 100) {
        throw std::invalid_argument("GSI must be between 0 and 100");
    }
    if (D < 0 || D > 1) {
        throw std::invalid_argument("D must be between 0 and 1");
    }
    
    double denominator = 28.0 - 14.0 * D;
    if (std::abs(denominator) < 1e-10) {
        denominator = 1e-10;  // Prevent division by zero
    }
    
    return mi * std::exp((GSI - 100.0) / denominator);
}

/**
 * Compute the rock mass parameter s
 * 
 * s = exp((GSI - 100) / (9 - 3*D))
 * 
 * For intact rock (GSI=100): s = 1
 * For very poor rock (low GSI): s → 0
 * 
 * @param GSI  Geological Strength Index (0-100)
 * @param D    Disturbance factor (0-1)
 * @return     Rock mass parameter s (0 to 1)
 */
inline double compute_s(double GSI, double D) {
    if (GSI < 0 || GSI > 100) {
        throw std::invalid_argument("GSI must be between 0 and 100");
    }
    if (D < 0 || D > 1) {
        throw std::invalid_argument("D must be between 0 and 1");
    }
    
    double denominator = 9.0 - 3.0 * D;
    if (std::abs(denominator) < 1e-10) {
        denominator = 1e-10;  // Prevent division by zero
    }
    
    return std::exp((GSI - 100.0) / denominator);
}

/**
 * Compute the rock mass parameter a
 * 
 * a = 0.5 + (1/6) * (exp(-GSI/15) - exp(-20/3))
 * 
 * For intact rock or high GSI: a ≈ 0.5
 * For very poor rock mass: a → 0.65
 * 
 * @param GSI  Geological Strength Index (0-100)
 * @return     Rock mass parameter a (approximately 0.5 to 0.65)
 */
inline double compute_a(double GSI) {
    if (GSI < 0 || GSI > 100) {
        throw std::invalid_argument("GSI must be between 0 and 100");
    }
    
    return 0.5 + (1.0/6.0) * (std::exp(-GSI / 15.0) - std::exp(-20.0 / 3.0));
}

/**
 * Compute all three derived Hoek-Brown parameters at once
 * 
 * @param mi   Material constant for intact rock
 * @param GSI  Geological Strength Index (0-100)
 * @param D    Disturbance factor (0-1)
 * @return     Tuple of (mb, s, a)
 */
inline std::tuple<double, double, double> compute_rock_mass_parameters(double mi, double GSI, double D) {
    double mb = compute_mb(mi, GSI, D);
    double s = compute_s(GSI, D);
    double a = compute_a(GSI);
    
    return std::make_tuple(mb, s, a);
}

/**
 * Compute the tensile strength of the rock mass
 * 
 * σt = -s * σci / mb (approximate, exact for a = 0.5)
 * 
 * @param sigma_ci  Unconfined compressive strength
 * @param mb        Reduced material constant
 * @param s         Rock mass parameter s
 * @return          Tensile strength (negative value, tension is negative)
 */
inline double compute_tensile_strength(double sigma_ci, double mb, double s) {
    if (std::abs(mb) < 1e-15) {
        return 0.0;  // No tensile strength if mb is zero
    }
    return -s * sigma_ci / mb;
}

/**
 * Compute the uniaxial compressive strength of the rock mass
 * 
 * σcm = σci * s^a
 * 
 * @param sigma_ci  Unconfined compressive strength of intact rock
 * @param s         Rock mass parameter s
 * @param a         Rock mass parameter a
 * @return          Uniaxial compressive strength of rock mass
 */
inline double compute_rock_mass_ucs(double sigma_ci, double s, double a) {
    if (s <= 0) {
        return 0.0;  // Zero strength if s is zero
    }
    return sigma_ci * std::pow(s, a);
}

/**
 * Compute rock mass deformation modulus using the Hoek-Diederichs equation (2006)
 * 
 * E_rm = E_i * (0.02 + (1 - D/2) / (1 + exp((60 + 15*D - GSI) / 11)))
 * 
 * @param E_i   Intact rock Young's modulus
 * @param GSI   Geological Strength Index (0-100)
 * @param D     Disturbance factor (0-1)
 * @return      Rock mass deformation modulus
 */
inline double compute_rock_mass_modulus(double E_i, double GSI, double D) {
    double numerator = 1.0 - D / 2.0;
    double denominator = 1.0 + std::exp((60.0 + 15.0 * D - GSI) / 11.0);
    
    return E_i * (0.02 + numerator / denominator);
}

/**
 * Estimate mi from rock type (approximate values from Hoek 2018)
 * 
 * @param rock_type  String identifying rock type
 * @return           Estimated mi value
 */
inline double estimate_mi_from_rock_type(const std::string& rock_type) {
    // These are approximate values from Table 1 in Hoek & Brown (2018)
    // Actual values should be determined from triaxial testing
    
    // Very weak rocks (mi = 4-8)
    if (rock_type == "chalk" || rock_type == "rocksalt") return 6.0;
    if (rock_type == "claystone") return 4.0;
    
    // Weak sedimentary rocks (mi = 8-12)
    if (rock_type == "shale" || rock_type == "siltstone") return 7.0;
    if (rock_type == "mudstone") return 5.0;
    
    // Medium strength sedimentary rocks (mi = 12-18)
    if (rock_type == "sandstone") return 17.0;
    if (rock_type == "limestone") return 12.0;
    if (rock_type == "dolomite") return 9.0;
    
    // Strong metamorphic rocks (mi = 18-25)
    if (rock_type == "marble") return 9.0;
    if (rock_type == "slate") return 11.0;
    if (rock_type == "gneiss") return 28.0;
    if (rock_type == "schist") return 12.0;
    if (rock_type == "quartzite") return 20.0;
    
    // Very strong igneous rocks (mi = 25-35)
    if (rock_type == "granite") return 32.0;
    if (rock_type == "basalt") return 25.0;
    if (rock_type == "diorite") return 25.0;
    if (rock_type == "andesite") return 25.0;
    if (rock_type == "rhyolite") return 25.0;
    if (rock_type == "gabbro") return 27.0;
    
    // Default for unknown rock type
    return 15.0;
}

/**
 * Compute equivalent Mohr-Coulomb parameters for a given stress range
 * 
 * This function fits a linear Mohr-Coulomb envelope to the nonlinear
 * Hoek-Brown criterion over a specified stress range.
 * 
 * @param sigma_ci  Unconfined compressive strength
 * @param mb        Reduced material constant
 * @param s         Rock mass parameter s
 * @param a         Rock mass parameter a
 * @param sigma3_max Maximum confining stress for the fit
 * @return          Tuple of (cohesion c, friction angle phi in degrees)
 */
inline std::tuple<double, double> compute_equivalent_mohr_coulomb(
    double sigma_ci, double mb, double s, double a, double sigma3_max) {
    
    // Using the relationships from Hoek et al. (2002)
    // These are simplified approximations
    
    double sigma_cm = sigma_ci * std::pow(s, a);  // Rock mass UCS
    
    // Compute sigma'_3n = sigma3_max / sigma_ci
    double sigma3n = sigma3_max / sigma_ci;
    
    // Compute friction angle (simplified)
    // sin(phi) ≈ (6 * a * mb * (s + mb * sigma3n)^(a-1)) / 
    //            (2 * (1 + a) * (2 + a) + 6 * a * mb * (s + mb * sigma3n)^(a-1))
    double term = std::pow(s + mb * sigma3n, a - 1.0);
    double numerator = 6.0 * a * mb * term;
    double denominator = 2.0 * (1.0 + a) * (2.0 + a) + numerator;
    double sin_phi = numerator / denominator;
    
    // Ensure sin_phi is in valid range
    sin_phi = std::max(-1.0, std::min(1.0, sin_phi));
    double phi = std::asin(sin_phi) * 180.0 / M_PI;  // Convert to degrees
    
    // Compute cohesion
    // c = sigma_ci * ((1 + 2*a)*s + (1-a)*mb*sigma3n) * (s + mb*sigma3n)^(a-1) /
    //     ((1+a)*(2+a)*sqrt(1 + (6*a*mb*(s+mb*sigma3n)^(a-1))/((1+a)*(2+a))))
    double term1 = (1.0 + 2.0*a)*s + (1.0 - a)*mb*sigma3n;
    double sqrt_term = std::sqrt(1.0 + numerator / ((1.0+a)*(2.0+a)));
    double c = sigma_ci * term1 * term / ((1.0+a)*(2.0+a)*sqrt_term);
    
    return std::make_tuple(c, phi);
}

} // namespace HoekBrownUtils

#endif // HoekBrown_Utils_H
