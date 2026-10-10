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

#ifndef StiffSoilCap_PF_H
#define StiffSoilCap_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
#include "../AllASDModelParameterTypes.h"
#include <cmath>

/**
 * @brief StiffSoilCap_PF - Plastic flow direction for the cap (volumetric) mechanism
 * 
 * The cap mechanism uses ASSOCIATED flow, meaning the plastic potential equals
 * the yield function:
 * 
 *   G_c = F_c = (q* / α)² + p² - p_c²
 * 
 * where:
 *   - q* = q / f(θ) is the scaled deviatoric stress
 *   - f(θ) = (3 - sin(φ)) / (2 * (√3*cos(θ) - sin(θ)*sin(φ)))
 *   - α is the cap shape factor
 *   - p_c is the cap hardening variable
 *   - p is the mean stress (positive in compression)
 * 
 * The flow direction is computed as dG_c/dσ using numerical differentiation.
 * 
 * @tparam CapPressureType Internal variable type for cap pressure p_c
 */
template<class CapPressureType>
class StiffSoilCap_PF : public PlasticFlowBase<StiffSoilCap_PF<CapPressureType>>
{
public:

    static constexpr const char* NAME = "StiffSoilCap_PF";

    StiffSoilCap_PF():
        PlasticFlowBase<StiffSoilCap_PF<CapPressureType>>::PlasticFlowBase()
        {}

    /**
     * @brief Plastic potential function for cap mechanism (associated flow)
     * 
     * G_c = (q* / α)² + p² - p_c²
     * 
     * where q* = q / f(θ) with f(θ) for Mohr-Coulomb shape in deviatoric plane
     */
    template<typename StorageType, typename ParameterStorageType>
    double plasticPotential(const VoigtVector& sigma, 
                            const StorageType& internal_variables_storage,
                            const ParameterStorageType& parameters_storage) const
    {
        using namespace std;

        // Get parameters
        double phi = GET_PARAMETER_VALUE(MC_phi) * M_PI / 180.0;
        double alpha = GET_PARAMETER_VALUE(SS_alpha);

        // Get internal variable (cap pressure)
        auto pc = GET_TRIAL_INTERNAL_VARIABLE(CapPressureType);

        // Convert to geotechnical convention (compression positive)
        VoigtVector sig_geo = -sigma;

        // Mean stress p (positive in compression)
        double p = sig_geo.meanStress();

        // Deviatoric stress q = sqrt(3*J2)
        double q = sqrt(3.0 * sig_geo.getJ2());

        // Lode angle
        double theta = sig_geo.lodeAngle();

        // Scaling function f(θ) for Mohr-Coulomb shape in deviatoric plane
        double sin_phi = sin(phi);
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);

        double f_theta_denom = sqrt(3.0) * cos_theta - sin_theta * sin_phi;

        double f_theta;
        if (abs(f_theta_denom) > 1e-12) {
            f_theta = (3.0 - sin_phi) / (2.0 * f_theta_denom);
        } else {
            f_theta = 1.0;  // fallback for numerical stability
        }

        // Scaled deviatoric stress q*
        double q_star = q / f_theta;

        // Cap plastic potential (same as yield function for associated flow)
        // G_c = (q*/α)² + p² - p_c²
        double Gc = (q_star / alpha) * (q_star / alpha) + p * p - pc.value() * pc.value();

        return Gc;
    }

    PLASTIC_FLOW_DIRECTION
    {
        using namespace std;

        double ds = GET_PARAMETER_VALUE(MC_ds);

        double sigma_norm = sigma.norm();

        // Perturbation scaled by stress magnitude
        ds = std::max(ds, ds * sigma_norm);

        // Helper lambda for numerical differentiation of plastic potential
        auto computeNumericalFlowDirection = [this, &internal_variables_storage, &parameters_storage]
            (const VoigtVector& sig, double perturbation) -> VoigtVector {
            VoigtVector result;
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sig;
                VoigtVector SIG2 = sig;
                
                SIG1(i) += perturbation;
                SIG2(i) -= perturbation;

                double g1 = plasticPotential(SIG1, internal_variables_storage, parameters_storage);
                double g2 = plasticPotential(SIG2, internal_variables_storage, parameters_storage);

                result(i) =  (g1 - g2) / (2.0 * perturbation); //- is due to geotech convention
            }
            return result;
        };

        vv_out = computeNumericalFlowDirection(sigma, ds);

        return vv_out;
    }

    using internal_variables_t = std::tuple<CapPressureType>;

    using parameters_t = std::tuple<MC_phi, MC_ds, SS_alpha>;

private:

    static VoigtVector vv_out;
};

template<class CapPressureType>
VoigtVector StiffSoilCap_PF<CapPressureType>::vv_out;

#endif
