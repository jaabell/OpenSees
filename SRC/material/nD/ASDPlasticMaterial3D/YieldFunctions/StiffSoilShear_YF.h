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

#ifndef StiffSoilShear_YF_H
#define StiffSoilShear_YF_H

#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>
#include "../AllASDModelParameterTypes.h"


double cot(double th){ return std::cos(th) / std::sin(th);};


template<class EpsQpShearType>
class StiffSoilShear_YF : public YieldFunctionBase<StiffSoilShear_YF<EpsQpShearType>> // CRTP
{
public:

    static constexpr const char* NAME = "StiffSoilShear_YF";


    StiffSoilShear_YF( ):
        YieldFunctionBase<StiffSoilShear_YF<EpsQpShearType>>::YieldFunctionBase() 
        {}
    // YIELD_FUNCTION 
    // {
    //     using namespace std;
    //
    //     double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
    //     double c = GET_PARAMETER_VALUE(MC_c);
    //     double E50_ref = GET_PARAMETER_VALUE(SS_E50_ref);
    //     double pref = GET_PARAMETER_VALUE(SS_pref);
    //     double Rf = GET_PARAMETER_VALUE(SS_Rf);
    //     double Eur_ref = GET_PARAMETER_VALUE(SS_Eur_ref);
    //     double m = GET_PARAMETER_VALUE(SS_m);
    //     double p_lim = 0.1 * pref;
    //
    //     // std::cout << "YF " << endl;
    //     // std::cout << "   phi = " << phi << endl;
    //     // std::cout << "   c = " << c << endl;
    //     // std::cout << "   E50_ref = " << E50_ref << endl;
    //     // std::cout << "   pref = " << pref << endl;
    //     // std::cout << "   Rf = " << Rf << endl;
    //     // std::cout << "   Eur_ref = " << Eur_ref << endl;
    //     // std::cout << "   m = " << m << endl;
    //
    //
    //     VoigtVector sig_geo = -sigma;   // Geotech stress 
    //     double q = sig_geo.stressDeviatorQ();
    //     auto [sigma3,sigma2,sigma1] = sig_geo.principalStresses(); // ascending order in geotech conventon
    //     sigma3 = std::max(sigma3, p_lim);
    //     // std::cout << "   sigma = " << sigma.transpose() << endl;
    //     // std::cout << "   sigma3 = " << sig_geo.transpose() << endl;
    //     // std::cout << "   sigma1 = " << sigma1 << endl;
    //     // std::cout << "   sigma2 = " << sigma2 << endl;
    //     // std::cout << "   sigma3 = " << sigma3 << endl;
    //     // std::cout << "   q = " << q << endl;
    //
    //     auto eps_qp_shear = GET_TRIAL_INTERNAL_VARIABLE(EpsQpShearType);
    //
    //     // cout << "  YF - eps_qp_shear = " << eps_qp_shear << endl;
    //
    //     double qf = (c * cot(phi) + sigma3) * 2 * sin(phi) / (1 - sin(phi));
    //     double qa = qf / Rf;
    //     double denom = (c*cos(phi) + pref*sin(phi));
    //     double E_ur = Eur_ref * pow((c*cos(phi) + sigma3*sin(phi)) / denom, m);
    //     double E50  = E50_ref * pow((c*cos(phi) + sigma3*sin(phi)) / denom, m);
    //
    //     double Ei = 2 * E50 / ( 2 - Rf);
    //
    //     // =====================================================================
    //     // PROPER MOHR-COULOMB FAILURE CAP
    //     //
    //     // The hyperbolic yield function is:
    //     //   F = γ_p(q) - eps_qp = 0
    //     // where γ_p(q) = q/(Ei*(1-q/qa)) - q/Eur
    //     //
    //     // Problem: As eps_qp increases (hardening), q can exceed qf without bound.
    //     //
    //     // Solution: Cap eps_qp at the failure value (eps_qp at q = qf).
    //     // This makes the material perfectly plastic at q = qf.
    //     //
    //     // At failure (q = qf, noting qf = Rf*qa so 1 - qf/qa = 1 - Rf):
    //     //   eps_qp_failure = qf/(Ei*(1-Rf)) - qf/Eur
    //     //
    //     // Modified yield function:
    //     //   F = γ_p(q) - min(eps_qp, eps_qp_failure)
    //     //
    //     // This ensures q cannot exceed qf regardless of plastic strain.
    //     // =====================================================================
    //
    //     // Compute plastic strain at failure
    //     double one_minus_Rf = 1.0 - Rf;
    //     double eps_qp_failure = qf / (Ei * one_minus_Rf) - qf / E_ur;
    //
    //     // Cap the effective plastic strain at failure value
    //     double eps_qp_effective = std::min(eps_qp_shear.value(), eps_qp_failure);
    //
    //     // Hyperbolic yield function with capped hardening
    //     // For numerical safety, also cap q to avoid singularity at qa
    //     double q_safe = std::min(q, 0.9999 * qa);
    //     double gamma_p = q_safe / (Ei * (1.0 - q_safe / qa)) - q_safe / E_ur;
    //
    //     double Fs = gamma_p - eps_qp_effective;
    //
    //     // std::cout << "    qf = " << qf << std::endl;
    //     // std::cout << "    qa = " << qa << std::endl;
    //     // std::cout << "    q = " << q << std::endl;
    //     // std::cout << "    eps_qp = " << eps_qp_shear.value() << std::endl;
    //     // std::cout << "    eps_qp_failure = " << eps_qp_failure << std::endl;
    //     // std::cout << "    eps_qp_effective = " << eps_qp_effective << std::endl;
    //     // std::cout << "--> Fs = " << Fs << std::endl;
    //
    //     //Handle tensile part
    //     double T = 0.0;
    //     double Ft = T - sigma3;
    //
    //     if (sigma3 > T + 1e-10 && sigma1 > 1e-10) {
    //       //pass
    //     } else {
    //       Fs = -1e10;
    //     }
    //     double F = std::max(Fs, Ft);
    //
    //     return F;
    // }

    YIELD_FUNCTION 
    {
        using namespace std;

        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double E50_ref = GET_PARAMETER_VALUE(SS_E50_ref);
        double pref = GET_PARAMETER_VALUE(SS_pref);
        double Rf = GET_PARAMETER_VALUE(SS_Rf);
        double Eur_ref = GET_PARAMETER_VALUE(SS_Eur_ref);
        double m = GET_PARAMETER_VALUE(SS_m);
        double p_lim = 0.1 * pref;

        // std::cout << "YF " << endl;
        // std::cout << "   phi = " << phi << endl;
        // std::cout << "   c = " << c << endl;
        // std::cout << "   E50_ref = " << E50_ref << endl;
        // std::cout << "   pref = " << pref << endl;
        // std::cout << "   Rf = " << Rf << endl;
        // std::cout << "   Eur_ref = " << Eur_ref << endl;
        // std::cout << "   m = " << m << endl;


        VoigtVector sig_geo = -sigma;   // Geotech stress 
        double q = sig_geo.stressDeviatorQ();
        auto [sigma3,sigma2,sigma1] = sig_geo.principalStresses(); // ascending order in geotech conventon
        sigma3 = std::max(sigma3, p_lim);
        // double q = (sigma1 - sigma3);
        // std::cout << "   sigma = " << sigma.transpose() << endl;
        // std::cout << "   sigma3 = " << sig_geo.transpose() << endl;
        // std::cout << "   sigma1 = " << sigma1 << endl;
        // std::cout << "   sigma2 = " << sigma2 << endl;
        // std::cout << "   sigma3 = " << sigma3 << endl;
        // std::cout << "   q = " << q << endl;

        auto eps_qp_shear = GET_TRIAL_INTERNAL_VARIABLE(EpsQpShearType);


        

        
        // cout << "  YF - eps_qp_shear = " << eps_qp_shear << endl;

        double qf = (c * cot(phi) + sigma3) * 2 * sin(phi) / (1 - sin(phi));
        // std::cout << "    qf = " << qf << std::endl;
        double qa = qf / Rf;
        double denom = (c*cos(phi) + pref*sin(phi));
        double E_ur = Eur_ref * pow((c*cos(phi) + sigma3*sin(phi)) / denom, m);
        double E50  = E50_ref * pow((c*cos(phi) + sigma3*sin(phi)) / denom, m);

        double Ei = 2 * E50 / ( 2 - Rf);
        // Cap the effective plastic strain at failur
        double one_minus_Rf = 1.0 - Rf;
        double eps_qp_failure = qf / (Ei * one_minus_Rf) - qf / E_ur;
        double eps_qp_effective = std::min(eps_qp_shear.value(), eps_qp_failure);
        // =====================================================================
        // SMOOTH REGULARIZATION: Avoid singularity without discontinuous switch
        //
        // The hyperbolic yield function has a singularity at q = qa (where qa > qf).
        // Strategy: Cap q at a safe value (q_max = 0.95*qf) to stay well below 
        // the asymptote, keeping the yield function and its gradients smooth.
        //
        // When q exceeds q_max, the yield function still uses the hyperbolic 
        // formula evaluated at q_max, but adds a linear penalty term for the 
        // excess (q - q_max). This creates a continuous function with continuous
        // first derivatives.
        // =====================================================================

        // Safe upper bound for q in the hyperbolic formula
        // Stay well below qf to avoid issues near the asymptote
        double q_max = 0.99 * qa;

        double Fs;

        if (q <= q_max) {
            // Normal hyperbolic yield function
            Fs = q / (Ei * (1.0 - q / qa)) - q / E_ur - eps_qp_effective;
        } else {
            // Above q_max: use hyperbolic formula at q_max + linear extension
            // This keeps gradients continuous
            double F_at_qmax = q_max / (Ei * (1.0 - q_max / qa)) - q_max / E_ur - eps_qp_effective;

            // Gradient of F w.r.t. q at q_max (for smooth continuation)
            // dF/dq = d/dq[q/(Ei*(1-q/qa))] - 1/Eur
            //       = [1*(1-q/qa) + q/qa] / [Ei*(1-q/qa)^2] - 1/Eur
            //       = 1 / [Ei*(1-q/qa)^2] - 1/Eur
            double one_minus_ratio = 1.0 - q_max / qa;
            double dF_dq_at_qmax = 1.0 / (Ei * one_minus_ratio * one_minus_ratio) - 1.0 / E_ur;

            // Linear extension from q_max
            Fs = F_at_qmax + dF_dq_at_qmax * (q - q_max);
        }

        // std::cout << "    qf = " << qf << std::endl;
        // std::cout << "    qa = " << qa << std::endl;
        // std::cout << "    q = " << q << std::endl;
        // std::cout << "    q_max = " << q_max << std::endl;
        // std::cout << "--> Fs = " << Fs << std::endl;

        //Handle tensile part
        double T = 0.0;
        double Ft = T - sigma3;

        if (sigma3 > T + 1e-10 && sigma1 > 1e-10) {
          //pass
        } else {
          Fs = -1e10;
        }
        double F = std::max(Fs, Ft);

        return F;
    }
    // YIELD_FUNCTION 
    // {
    //     using namespace std;
    //
    //     double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
    //     double c = GET_PARAMETER_VALUE(MC_c);
    //     double E50_ref = GET_PARAMETER_VALUE(SS_E50_ref);
    //     double pref = GET_PARAMETER_VALUE(SS_pref);
    //     double Rf = GET_PARAMETER_VALUE(SS_Rf);
    //     double Eur_ref = GET_PARAMETER_VALUE(SS_Eur_ref);
    //     double m = GET_PARAMETER_VALUE(SS_m);
    //     double p_lim = 0.1 * pref;
    //
    //     // std::cout << "YF " << endl;
    //     // std::cout << "   phi = " << phi << endl;
    //     // std::cout << "   c = " << c << endl;
    //     // std::cout << "   E50_ref = " << E50_ref << endl;
    //     // std::cout << "   pref = " << pref << endl;
    //     // std::cout << "   Rf = " << Rf << endl;
    //     // std::cout << "   Eur_ref = " << Eur_ref << endl;
    //     // std::cout << "   m = " << m << endl;
    //
    //
    //     VoigtVector sig_geo = -sigma;   // Geotech stress 
    //     double q = sig_geo.stressDeviatorQ();
    //     auto [sigma3,sigma2,sigma1] = sig_geo.principalStresses(); // ascending order in geotech conventon
    //     sigma3 = std::max(sigma3, p_lim);
    //     // std::cout << "   sigma = " << sigma.transpose() << endl;
    //     // std::cout << "   sigma3 = " << sig_geo.transpose() << endl;
    //     // std::cout << "   sigma1 = " << sigma1 << endl;
    //     // std::cout << "   sigma2 = " << sigma2 << endl;
    //     // std::cout << "   sigma3 = " << sigma3 << endl;
    //     // std::cout << "   q = " << q << endl;
    //
    //     auto eps_qp_shear = GET_TRIAL_INTERNAL_VARIABLE(EpsQpShearType);
    //
    //     // cout << "  YF - eps_qp_shear = " << eps_qp_shear << endl;
    //
    //     double qf = (c * cot(phi) + sigma3) * 2 * sin(phi) / (1 - sin(phi));
    //     double qa = qf / Rf;
    //     // double qf = (2.0 * c * cos(phi) + 2.0 * sigma3 * sin(phi)) / (1.0 - sin(phi));
    //     // double qa = qf / Rf;
    //     double denom = (c*cos(phi) + pref*sin(phi));
    //     double E_ur = Eur_ref * pow((c*cos(phi) + sigma3*sin(phi)) / denom, m);
    //     double E50  = E50_ref * pow((c*cos(phi) + sigma3*sin(phi)) / denom, m);
    //
    //     double Ei = 2 * E50 / ( 2 - Rf);
    //     if (q / qa > 1) {
    //      q = 0.99*qa; 
    //     }
    //     double Fs = q / (Ei * (1 - q / qa)) - q / E_ur - eps_qp_shear.value();
    //     // std::cout << "    qf = " << qf << std::endl;
    //     // std::cout << "    qa = " << qa << std::endl;
    //     // std::cout << "    denom = " << denom << std::endl;
    //     // std::cout << "    E_ur = " << E_ur << std::endl;
    //     // std::cout << "    E50 = " << E50 << std::endl;
    //     // std::cout << "    Ei = " << Ei << std::endl;
    //     // std::cout << "--> Fs = " << Fs << std::endl;
    //    // Fs *= E50_ref; 
    //
    //     //Handle tensile part
    //     double T = 0.0;
    //     double Ft = T - sigma3;
    //
    //     if (sigma3 > T + 1e-10 && sigma1 > 1e-10) {
    //       //pass
    //     } else {
    //       Fs = -1e10;
    //     }
    //     double F = std::max(Fs, Ft);
    //
    //     return F;
    // }

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

                result(i) = (yf1 - yf2) / (2*perturbation);
            }
            return result;
        };

        vv_out = computeNumericalDerivative(sigma, ds);
        double norm = std::sqrt(tensor_dot_stress_like(vv_out, vv_out));
        // vv_out /= norm;
        return vv_out;
    }

    YIELD_FUNCTION_HARDENING
    {
      auto eps_qp_shear = GET_TRIAL_INTERNAL_VARIABLE(EpsQpShearType);
      double HH = GET_INTERNAL_VARIABLE_HARDENING(EpsQpShearType).value();

      double E50_ref = GET_PARAMETER_VALUE(SS_E50_ref);

      using namespace std;
      // cout << "  YF HARD  eps_qp_shear = " << eps_qp_shear.value() << endl;
      // cout << "  YF HARD  HH = " << HH << endl;
      return -HH; 
      // return -HH*E50_ref;
    }

    using internal_variables_t = std::tuple<EpsQpShearType>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds,SS_E50_ref, SS_Eur_ref, SS_Rf, SS_m, SS_pref>;

private:


    static VoigtVector vv_out; //For returning VoigtVector's
};

template <class EpsQpShearType>
VoigtVector StiffSoilShear_YF<EpsQpShearType>::vv_out;


#endif
