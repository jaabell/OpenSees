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

#ifndef StiffSoilCap_YF_H
#define StiffSoilCap_YF_H

#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>
#include "../AllASDModelParameterTypes.h"

/**
 * @brief StiffSoilCap_YF - Elliptical cap yield surface for the Hardening Soil model
 * 
 * The cap yield surface is defined as (PLAXIS Eq. 15.11):
 * 
 *   F_c = (q* / alpha)^2 + p^2 - p_c^2 = 0
 * 
 * where:
 *   - q* = q / f(theta) is the scaled deviatoric stress
 *   - f(theta) = (3 - sin(phi)) / (2 * (sqrt(3)*cos(theta) - sin(theta)*sin(phi)))
 *   - alpha is the shape factor for the elliptical cap
 *   - p_c is the cap hardening variable (internal variable)
 *   - p is the mean stress (positive in compression for geotechnical convention)
 * 
 * The hardening law is (PLAXIS Eq. 15.13):
 *   eps_v^{p-cap} = (beta / (1-m)) * (p_c / p_ref)^{1-m}
 * 
 * @tparam CapPressureType The internal variable type for cap pressure p_c
 */
template<class CapPressureType>
class StiffSoilCap_YF : public YieldFunctionBase<StiffSoilCap_YF<CapPressureType>> // CRTP
{
public:

    static constexpr const char* NAME = "StiffSoilCap_YF";

    StiffSoilCap_YF():
        YieldFunctionBase<StiffSoilCap_YF<CapPressureType>>::YieldFunctionBase() 
        {}

    YIELD_FUNCTION 
    {
        using namespace std;

        // Get parameters
        double phi = GET_PARAMETER_VALUE(MC_phi) * M_PI / 180.0;  // friction angle in radians
        double alpha = GET_PARAMETER_VALUE(SS_alpha);              // cap shape factor

        // Get internal variable (cap pressure)
        auto pc = GET_TRIAL_INTERNAL_VARIABLE(CapPressureType);

        // Convert to geotechnical convention (compression positive)
        VoigtVector sig_geo = -sigma;
        
        // Mean stress p (positive in compression)
        double p = sig_geo.meanStress();
        
        // Deviatoric stress q
        double q = sqrt(3.0 * sig_geo.getJ2());  // q = sqrt(3*J2)
        
        // Lode angle theta
        double theta = sig_geo.lodeAngle();
        
        // Scaling function f(theta) for Mohr-Coulomb shape in deviatoric plane
        // f(theta) = (3 - sin(phi)) / (2 * (sqrt(3)*cos(theta) - sin(theta)*sin(phi)))
        double sin_phi = sin(phi);
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);
        
        double f_theta_denom = sqrt(3.0) * cos_theta - sin_theta * sin_phi;
        
        // Avoid division by zero
        double f_theta;
        if (abs(f_theta_denom) > 1e-12) {
            f_theta = (3.0 - sin_phi) / (2.0 * f_theta_denom);
        } else {
            f_theta = 1.0;  // fallback for numerical stability
        }
        
        // Scaled deviatoric stress q*
        double q_star = q / f_theta;
        
        // Cap yield function: F_c = (q*/alpha)^2 + p^2 - p_c^2
        double Fc = (q_star / alpha) * (q_star / alpha) + p * p - pc.value() * pc.value();

        return Fc;
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {  
        using namespace std;

        double ds = GET_PARAMETER_VALUE(MC_ds);
        double sigma_norm = sigma.norm();

        // Perturbation to smooth the YF
        ds = ds * std::max(1.0, sigma_norm);

        // Helper lambda for numerical differentiation
        auto computeNumericalDerivative = [this, &internal_variables_storage, &parameters_storage](const VoigtVector& sig, double perturbation) -> VoigtVector {
            VoigtVector result;
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sig;
                VoigtVector SIG2 = sig;

                SIG1(i) += perturbation;
                SIG2(i) -= perturbation;

                double yf1 = YF(SIG1);
                double yf2 = YF(SIG2);

                result(i) = (yf1 - yf2) / (2 * perturbation);
            }
            return result;
        };

        vv_out = computeNumericalDerivative(sigma, ds);

        return vv_out;
    }

    YIELD_FUNCTION_HARDENING
    {
        // df/dp_c = -2 * p_c
        // The hardening contribution is df/dp_c * dp_c/d_lambda = df/dp_c * h_pc
        auto pc = GET_TRIAL_INTERNAL_VARIABLE(CapPressureType);
        
        double df_dpc = -2.0 * pc.value();
        
        double dbl_result = (df_dpc * GET_INTERNAL_VARIABLE_HARDENING(CapPressureType)).value();
        
        return dbl_result;
    }

    using internal_variables_t = std::tuple<CapPressureType>;

    using parameters_t = std::tuple<MC_phi, MC_ds, SS_alpha, SS_pref, SS_m, SS_beta>;

private:

    static VoigtVector vv_out; //For returning VoigtVector's
};

template <class CapPressureType>
VoigtVector StiffSoilCap_YF<CapPressureType>::vv_out;


#endif
