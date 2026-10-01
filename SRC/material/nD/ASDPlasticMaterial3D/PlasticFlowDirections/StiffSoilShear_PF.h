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

        // VoigtVector sig_geo = -sigma;
        VoigtVector sig_geo = sigma;

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
// PLASTIC_FLOW_DIRECTION
// {
//     using namespace std;
//
//     double phi = GET_PARAMETER_VALUE(MC_phi) * M_PI / 180.0;
//     double psi = GET_PARAMETER_VALUE(MC_psi) * M_PI / 180.0;
//     double c = GET_PARAMETER_VALUE(MC_c);
//
//         // Compute mobilized dilatancy
//         double psi_m = computeMobilizedDilatancy(sigma, phi, psi, c);
//         double sin_psi_m = sin(psi_m);
//
//         // Get deviatoric stress direction
//         VoigtVector s = sigma.deviator();
//         double J2 = sigma.getJ2();
//         double sqrt_J2 = sqrt(std::max(J2, 1e-20));
//
//         // Deviatoric flow direction (normalized deviatoric stress)
//         // m_dev = s / (2 * sqrt(J2)) gives ||m_dev|| = sqrt(J2/2)
//         // We want the deviatoric part to give ε_q^p when contracted
//         VoigtVector m_dev = s / (2.0 * sqrt_J2);
//
//         // Volumetric flow direction
//         // From Eq. 15.5: ε̇_v^p = sin(ψ_m) * ε̇_q^p
//         // The volumetric part is (1/3) * sin(ψ_m) * I
//         // But we need to scale properly with the deviatoric part
//
//         // For triaxial: if m_dev gives ε_q = m_1 - m_3, then
//         // m_vol should give ε_v = sin(ψ_m) * ε_q
//         // 
//         // Scale factor: ||s|| / (2*sqrt(J2)) for deviatoric part
//         // gives proper ε_q scaling
//
//         double scale = sqrt(3.0 / 2.0);  // Relates ||s|| to q
//
//         // Flow direction: deviatoric + volumetric
//         vv_out = m_dev;
//
//         // Add volumetric component (compression negative in mechanics convention)
//         // ε_v = ε_11 + ε_22 + ε_33, so add sin(ψ_m)/3 to each diagonal
//         double vol_contrib = -sin_psi_m / 3.0 * scale;  // negative for dilation
//         vv_out(0) += vol_contrib;
//         vv_out(1) += vol_contrib;
//         vv_out(2) += vol_contrib;
//
//         return vv_out;
//     }
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

                result(i) = (g1 - g2) / (2.0 * perturbation); 
            }
            return result;
        };

        vv_out = computeNumericalFlowDirection(sigma, ds);
        double norm = std::sqrt(tensor_dot_stress_like(vv_out, vv_out));
        vv_out /= norm;
        return vv_out;

        // return vv_out;
    }
    //
    // PLASTIC_FLOW_DIRECTION
    // {
    //     using namespace std;
    //
    //     double phi = GET_PARAMETER_VALUE(MC_phi) * M_PI / 180.0;
    //     double psi = GET_PARAMETER_VALUE(MC_psi) * M_PI / 180.0;
    //     double c = GET_PARAMETER_VALUE(MC_c);
    //
    //     // Compute mobilized dilatancy (Rowe's theory)
    //     double psi_m = computeMobilizedDilatancy(sigma, phi, psi, c);
    //     // Or use constant dilatancy: double psi_m = psi;
    //
    //     double sin_psi_m = sin(psi_m);
    //
    //     // Get principal stresses and directions in geotechnical convention
    //     VoigtVector sig_geo = -sigma;  // Convert to compression positive
    //
    //     auto [stresses, directions] = sig_geo.principalStressesAndDirections();
    //     auto [s1, s2, s3] = stresses;      // s1 <= s2 <= s3 (ascending)
    //     auto [n1, n2, n3] = directions;    // Corresponding unit eigenvectors
    //
    //     // In geotechnical convention with ascending order:
    //     //   s1 = minor principal (least compressive / most tensile)
    //     //   s3 = major principal (most compressive)
    //     // For triaxial compression: s3 = σ_axial, s1 = s2 = σ_lateral
    //
    //     // Principal plastic strain rates (PLAXIS-consistent):
    //     // These satisfy:
    //     //   m3 - m1 = 1  (gives ε_q = dλ for triaxial definition ε_q = ε_1 - ε_3)
    //     //   m1 + m2 + m3 = -sin(ψ_m)  (volumetric strain, negative = dilation)
    //     //
    //     // For triaxial compression (s3 > s1 = s2):
    //     //   m3 = major principal strain rate (most compressive, negative)
    //     //   m1 = minor principal strain rate (least compressive, positive for dilation)
    //
    //     double m3_p = -(2.0 + sin_psi_m) / 3.0;  // Major principal (compression, negative)
    //     double m1_p = -(sin_psi_m - 1.0) / 3.0;  // Minor principal (extension, positive for dilation)
    //
    //     // Intermediate principal - smooth interpolation
    //     double m2_p = -sin_psi_m / 3.0;  // = (m1_p + m3_p) / 2
    //
    //     // Alternatively, for Lode-angle dependent intermediate value:
    //     // double lode = sig_geo.lodeAngle();  // -π/6 to +π/6
    //     // double t = (lode + M_PI/6.0) / (M_PI/3.0);  // 0 at triax comp, 1 at triax ext
    //     // double m2_p = (1.0 - t) * m1_p + t * m3_p;
    //
    //     // Transform back to Cartesian coordinates
    //     // m_ij = sum_k( m_k * n_k_i * n_k_j )
    //     vv_out(0) = m1_p * n1(0) * n1(0) + m2_p * n2(0) * n2(0) + m3_p * n3(0) * n3(0);  // 11
    //     vv_out(1) = m1_p * n1(1) * n1(1) + m2_p * n2(1) * n2(1) + m3_p * n3(1) * n3(1);  // 22
    //     vv_out(2) = m1_p * n1(2) * n1(2) + m2_p * n2(2) * n2(2) + m3_p * n3(2) * n3(2);  // 33
    //     vv_out(3) = m1_p * n1(0) * n1(1) + m2_p * n2(0) * n2(1) + m3_p * n3(0) * n3(1);  // 12
    //     vv_out(4) = m1_p * n1(1) * n1(2) + m2_p * n2(1) * n2(2) + m3_p * n3(1) * n3(2);  // 23
    //     vv_out(5) = m1_p * n1(0) * n1(2) + m2_p * n2(0) * n2(2) + m3_p * n3(0) * n3(2);  // 13
    //
    //     // Note: vv_out is already in mechanics convention (tension positive) because
    //     // we applied the negative signs to m1_p, m2_p, m3_p
    //
    //     return vv_out;
    // }
    // PLASTIC_FLOW_DIRECTION
    // {
    //     using namespace std;
    //
    //     double phi = GET_PARAMETER_VALUE(MC_phi) * M_PI / 180.0;
    //     double psi = GET_PARAMETER_VALUE(MC_psi) * M_PI / 180.0;
    //     double c   = GET_PARAMETER_VALUE(MC_c);
    //
    //     // mobilized dilatancy (Rowe) or constant:
    //     double psi_m = computeMobilizedDilatancy(sigma, phi, psi, c);
    //     double sin_psi_m = sin(psi_m);
    //
    //     // principal stresses of sigma (mechanics convention)
    //     auto [s3, s2, s1] = sigma.principalStresses(); // ascending
    //     double denom = s1 - s3;
    //
    //     // deviatoric direction
    //     VoigtVector s = sigma.deviator();
    //
    //     vv_out.setZero();
    //
    //     if (std::abs(denom) > 1e-14)
    //     {
    //         // Normalize so that |m1 - m3| = 1  -> dεq^p = dλ
    //         vv_out = s / denom;
    //
    //         // Add volumetric part so that tr(m) = sin(psi_m)
    //         double a = sin_psi_m / 3.0;
    //         vv_out(0) += a;
    //         vv_out(1) += a;
    //         vv_out(2) += a;
    //     }
    //
    //     return vv_out;
    // }
    //


    using internal_variables_t = std::tuple<EpsQpShearType>;

    using parameters_t = std::tuple<MC_phi, MC_psi, MC_c, MC_ds>;

private:

    static VoigtVector vv_out;
};

template<class EpsQpShearType>
VoigtVector StiffSoilShear_PF<EpsQpShearType>::vv_out;

#endif
