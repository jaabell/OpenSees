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

// Original implementation: José Abell (UANDES), Massimo Petracca (ASDEA)
//
// ASDPlasticMaterial3D
//
// Fully general templated material class for plasticity modeling

#ifndef StiffSoil_EL_H
#define StiffSoil_EL_H

#include "../ElasticityBase.h"
#include "../AllASDModelParameterTypes.h"

#include <iostream>
#include <cmath>
#include <algorithm>

/**
 * @brief StiffSoil_EL - Stress-dependent elasticity for the Hardening Soil model
 * 
 * The elastic modulus for unloading/reloading is stress-dependent (PLAXIS Eq. 15.3):
 * 
 *   E_ur = E_ur^ref * ((c*cos(φ) + σ3*sin(φ)) / (c*cos(φ) + p_ref*sin(φ)))^m
 * 
 * where:
 *   - E_ur^ref : Reference unloading/reloading modulus at reference pressure
 *   - c        : Cohesion
 *   - φ        : Friction angle
 *   - σ3       : Minor principal stress (most compressive, positive in compression)
 *   - p_ref    : Reference pressure
 *   - m        : Power exponent (typically 0.5 < m < 1.0)
 *   - ν        : Poisson's ratio (typically 0.2 for unloading/reloading)
 * 
 * A minimum stress limit (p_lim) is applied to prevent zero or negative stiffness
 * when σ3 approaches zero or becomes tensile. This is set to ~10% of p_ref.
 */
class StiffSoil_EL : public ElasticityBase<StiffSoil_EL>
{
public:

    static constexpr const char* NAME = "StiffSoil_EL";

    StiffSoil_EL(): ElasticityBase<StiffSoil_EL>::ElasticityBase()
    {
    }

    ELASTICITY_MATRIX
    {
        using namespace std;

        // Get parameters
        double Eur_ref = GET_PARAMETER_VALUE(SS_Eur_ref);  // Reference E for unload/reload
        double nu = GET_PARAMETER_VALUE(PoissonsRatio);    // Poisson's ratio
        double pref = GET_PARAMETER_VALUE(SS_pref);        // Reference pressure
        double m = GET_PARAMETER_VALUE(SS_m);              // Power exponent
        double phi = GET_PARAMETER_VALUE(MC_phi) * M_PI / 180.0;  // Friction angle in radians
        double c = GET_PARAMETER_VALUE(MC_c);              // Cohesion

        // Convert to geotechnical convention (compression positive)
        VoigtVector sig_geo = -stress;

        // Get minor principal stress σ3 (most compressive in geotechnical convention)
        // principalStresses() returns (σ3, σ2, σ1) in ascending order
        auto [sigma3, sigma2, sigma1] = sig_geo.principalStresses();

        // Apply minimum stress limit to prevent zero/negative stiffness
        // Set p_lim to ~10% of reference pressure (as recommended in PLAXIS documentation)
        double p_lim = 0.1 * pref;
        sigma3 = std::max(sigma3, p_lim);

        // Compute stress-dependent elastic modulus (Eq. 15.3)
        double sin_phi = sin(phi);
        double cos_phi = cos(phi);

        double numerator = c * cos_phi + sigma3 * sin_phi;
        double denominator = c * cos_phi + pref * sin_phi;

        // Ensure denominator is positive
        if (denominator < 1e-12) {
            denominator = 1e-12;
        }

        // Ensure numerator is positive (for stability)
        numerator = std::max(numerator, 0.1 * denominator);

        double E_ur = Eur_ref * pow(numerator / denominator, m);

        // Ensure minimum stiffness
        E_ur = std::max(E_ur, 0.01 * Eur_ref);

        // Build isotropic elastic stiffness matrix
        double lambda = (nu * E_ur) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        double mu = E_ur / (2.0 * (1.0 + nu));

        EE_MATRIX.setZero();

        EE_MATRIX(0, 0) = EE_MATRIX(1, 1) = EE_MATRIX(2, 2) = 2.0 * mu + lambda;
        EE_MATRIX(0, 1) = EE_MATRIX(1, 0) = lambda;
        EE_MATRIX(0, 2) = EE_MATRIX(2, 0) = lambda;
        EE_MATRIX(1, 2) = EE_MATRIX(2, 1) = lambda;
        EE_MATRIX(3, 3) = mu;
        EE_MATRIX(4, 4) = mu;
        EE_MATRIX(5, 5) = mu;

        return EE_MATRIX;
    }

    using parameters_t = std::tuple<SS_Eur_ref, PoissonsRatio, SS_pref, SS_m, MC_phi, MC_c>;

};

#endif
