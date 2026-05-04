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
        auto [sigma3,sigma2,sigma1] = sig_geo.principalStresses(); // ascending ordir in geotech conveniton
        // std::cout << "   sigma1 = " << sigma1 << endl;
        // std::cout << "   sigma2 = " << sigma2 << endl;
        // std::cout << "   sigma3 = " << sigma3 << endl;
        // std::cout << "   q = " << q << endl;

        auto eps_qp_shear = GET_TRIAL_INTERNAL_VARIABLE(EpsQpShearType);


        double qf = (c * cot(phi) + sigma1) * 2 * sin(phi) / (1 - sin(phi));
        // std::cout << "    qf = " << qf << std::endl;
        double qa = qf / Rf;
        // std::cout << "    qa = " << qa << std::endl;
        double denom = (c*cos(phi) + pref*sin(phi));
        // std::cout << "    denom = " << denom << std::endl;
        double E_ur = Eur_ref * pow((c*cos(phi) + sigma3*sin(phi)) / denom, m);
        // std::cout << "    E_ur = " << E_ur << std::endl;
        double E50  = E50_ref * pow((c*cos(phi) + sigma1*sin(phi)) / denom, m);
        // std::cout << "    E50 = " << E50 << std::endl;
        
        double Ei = 2 * E50 / ( 2 - Rf);
        // std::cout << "    Ei = " << Ei << std::endl;
        double Fs = q / (Ei * (1 - q / qa)) - q / E_ur - eps_qp_shear.value();
        // std::cout << "--> Fs = " << Fs << std::endl;
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

                result(i) = -(yf1 - yf2) / (2*perturbation); // - is due to geotech convention
            }
            return result;
        };

        vv_out = computeNumericalDerivative(sigma, ds);

        return vv_out;
    }

    YIELD_FUNCTION_HARDENING
    {
        return -1;
    }

    using internal_variables_t = std::tuple<EpsQpShearType>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds,SS_E50_ref, SS_Eur_ref, SS_Rf, SS_m, SS_pref>;

private:


    static VoigtVector vv_out; //For returning VoigtVector's
};

template <class EpsQpShearType>
VoigtVector StiffSoilShear_YF<EpsQpShearType>::vv_out;


#endif
