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

#ifndef StiffSoilShear_PF_H
#define StiffSoilShear_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
#include "../AllASDModelParameterTypes.h"
#include <cmath>

/**
 * @brief StiffSoilShear_PF - Plastic flow direction for the deviatoric (shear) mechanism
 * 
 * The deviatoric mechanism uses a non-associated flow rule. The plastic potential
 * has the same form as the yield function but uses the dilatancy angle ψ instead
 * of the friction angle φ for the volumetric component.
 * 
 * Flow rule (PLAXIS Eq. 15.5):
 *   ε̇_v^p = sin(ψ_m) * ε̇_q^p
 * 
 * The plastic potential is defined similarly to Mohr-Coulomb but with dilatancy:
 *   G_s = q / (E_i * (1 - q/q_a)) - q / E_ur - (deviatoric plastic strain)
 * 
 * For the flow direction, we use numerical differentiation of a Mohr-Coulomb
 * type potential with the dilatancy angle.
 * 
 * Options for mobilized dilatancy (Rowe's theory - Eq. 15.6-15.8):
 *   sin(ψ_m) = (sin(φ_m) - sin(φ_cv)) / (1 - sin(φ_m)*sin(φ_cv))
 * 
 * @tparam EpsQpShearType Internal variable type for plastic deviatoric strain
 */
template<class EpsQpShearType>
class StiffSoilShear_PF : public PlasticFlowBase<StiffSoilShear_PF<EpsQpShearType>>
{
public:

    static constexpr const char* NAME = "StiffSoilShear_PF";

    StiffSoilShear_PF():
        PlasticFlowBase<StiffSoilShear_PF<EpsQpShearType>>::PlasticFlowBase()
        {}

    /**
     * @brief Plastic potential function for deviatoric mechanism
     * 
     * Uses Mohr-Coulomb type potential with dilatancy angle ψ
     * G = (cos(θ) - sin(θ)*sin(ψ)/√3) * √J2 + I1*sin(ψ)/3 - c*cos(ψ)
     */
    double plasticPotential(const VoigtVector& sigma, double psi, double c) const
    {
        using namespace std;

        // Convert to geotechnical convention (compression positive)
        VoigtVector sig_geo = -sigma;

        double I1 = sig_geo.getI1();
        double J2 = sig_geo.getJ2();
        double lode_angle = sig_geo.lodeAngle();

        double sin_psi = sin(psi);
        double cos_psi = cos(psi);

        // Mohr-Coulomb type potential with dilatancy angle
        double rEquivalentStress = (cos(lode_angle) - sin(lode_angle) * sin_psi / sqrt(3.0)) * sqrt(J2) +
                                   I1 * sin_psi / 3.0;
        
        double rThresh = c * cos_psi;

        return rEquivalentStress - rThresh;
    }

    /**
     * @brief Compute mobilized dilatancy angle using Rowe's stress-dilatancy theory
     * 
     * sin(ψ_m) = (sin(φ_m) - sin(φ_cv)) / (1 - sin(φ_m)*sin(φ_cv))
     * sin(φ_m) = (σ1 - σ3) / (σ1 + σ3 - 2*c*cot(φ))
     * sin(φ_cv) = (sin(φ) - sin(ψ)) / (1 - sin(φ)*sin(ψ))
     */
    double computeMobilizedDilatancy(const VoigtVector& sigma, double phi, double psi, double c) const
    {
        using namespace std;

        // If ultimate dilatancy is zero, return zero
        if (abs(psi) < 1e-12) {
            return 0.0;
        }

        // Convert to geotechnical convention
        VoigtVector sig_geo = -sigma;
        auto [sigma3, sigma2, sigma1] = sig_geo.principalStresses();

        double sin_phi = sin(phi);
        double sin_psi = sin(psi);

        // Critical state friction angle (Eq. 15.8)
        double denom_cv = 1.0 - sin_phi * sin_psi;
        double sin_phi_cv = (abs(denom_cv) > 1e-12) ? (sin_phi - sin_psi) / denom_cv : 0.0;

        // Mobilized friction angle (Eq. 15.7)
        double cot_phi = (abs(sin_phi) > 1e-12) ? cos(phi) / sin_phi : 1e12;
        double denom_m = sigma1 + sigma3 - 2.0 * c * cot_phi;
        double sin_phi_m = (abs(denom_m) > 1e-12) ? (sigma1 - sigma3) / denom_m : 0.0;
        
        // Clamp to valid range
        sin_phi_m = std::max(-1.0, std::min(1.0, sin_phi_m));

        // Mobilized dilatancy (Eq. 15.6)
        double denom_psi = 1.0 - sin_phi_m * sin_phi_cv;
        double sin_psi_m = (abs(denom_psi) > 1e-12) ? (sin_phi_m - sin_phi_cv) / denom_psi : 0.0;

        // Clamp to valid range and ensure non-negative dilatancy
        sin_psi_m = std::max(0.0, std::min(1.0, sin_psi_m));

        return asin(sin_psi_m);
    }

    PLASTIC_FLOW_DIRECTION
    {
        using namespace std;

        double phi = GET_PARAMETER_VALUE(MC_phi) * M_PI / 180.0;
        double psi = GET_PARAMETER_VALUE(MC_psi) * M_PI / 180.0;
        double c = GET_PARAMETER_VALUE(MC_c);
        double ds = GET_PARAMETER_VALUE(MC_ds);

        double sigma_norm = sigma.norm();

        // Perturbation scaled by stress magnitude
        ds = ds * std::max(1.0, sigma_norm);

        // Option: Use mobilized dilatancy (Rowe's theory) or constant dilatancy
        // For now, use constant dilatancy angle (simpler approach)
        // To enable Rowe's theory, uncomment the following line:
        double psi_m = computeMobilizedDilatancy(sigma, phi, psi, c);
        // double psi_m = psi;  // Constant dilatancy

        // Helper lambda for numerical differentiation of plastic potential
        auto computeNumericalFlowDirection = [this, psi_m, c](const VoigtVector& sig, double perturbation) -> VoigtVector {
            VoigtVector result;
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sig;
                VoigtVector SIG2 = sig;
                
                SIG1(i) += perturbation;
                SIG2(i) -= perturbation;

                double g1 = plasticPotential(SIG1, psi_m, c);
                double g2 = plasticPotential(SIG2, psi_m, c);

                result(i) = -(g1 - g2) / (2.0 * perturbation); // - due to geotechnical convention
            }
            return result;
        };

        vv_out = computeNumericalFlowDirection(sigma, ds);

        return vv_out;
    }

    using internal_variables_t = std::tuple<EpsQpShearType>;

    using parameters_t = std::tuple<MC_phi, MC_psi, MC_c, MC_ds>;

private:

    static VoigtVector vv_out;
};

template<class EpsQpShearType>
VoigtVector StiffSoilShear_PF<EpsQpShearType>::vv_out;

#endif
