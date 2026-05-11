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
// ASDPlasticMaterial3D - Hoek-Brown Yield Function
//
// Implements the Hoek-Brown failure criterion for rock masses
// Based on: Hoek, E., Brown, E.T. (2018). The Hoek-Brown failure criterion 
// and GSI – 2018 edition. Journal of Rock Mechanics and Geotechnical Engineering.
//
// The generalized Hoek-Brown criterion is:
//   σ₁ = σ₃ + σci * (mb * σ₃/σci + s)^a
//
// Where:
//   σ₁, σ₃ = major and minor principal stresses (compression positive)
//   σci = unconfined compressive strength of intact rock
//   mb = reduced material constant for rock mass = mi * exp((GSI-100)/(28-14D))
//   s = rock mass parameter = exp((GSI-100)/(9-3D))
//   a = rock mass parameter = 0.5 + (1/6)*(exp(-GSI/15) - exp(-20/3))
//   GSI = Geological Strength Index (0-100)
//   D = Disturbance factor (0-1)
//   mi = material constant for intact rock
//
// Sign convention: Compression is POSITIVE (geomechanics convention)

#ifndef HoekBrown_YF_H
#define HoekBrown_YF_H

#include "../YieldFunctionBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
#include "../AllASDModelParameterTypes.h"
#include <cmath>
#include <iostream>
#include <algorithm>

template<class NO_HARDENING>
class HoekBrown_YF : public YieldFunctionBase<HoekBrown_YF<NO_HARDENING>> // CRTP
{
public:

    static constexpr const char* NAME = "HoekBrown_YF";

    HoekBrown_YF():
        YieldFunctionBase<HoekBrown_YF<NO_HARDENING>>::YieldFunctionBase()
    {}

    YIELD_FUNCTION
    {
        using namespace std;

        // Get Hoek-Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);  // Unconfined compressive strength
        double mb = GET_PARAMETER_VALUE(HB_mb);           // Reduced material constant
        double s = GET_PARAMETER_VALUE(HB_s);             // Rock mass parameter s
        double a = GET_PARAMETER_VALUE(HB_a);             // Rock mass parameter a

        // Get principal stresses (sorted: sigma1 >= sigma2 >= sigma3)
        // Note: In this code, compression is positive
        VoigtVector sigma_geo = -sigma;
        auto [sigma3, sigma2, sigma1] = sigma_geo.principalStresses();
        
        // The Hoek-Brown criterion in compression is:
        // f = σ₁ - σ₃ - σci * (mb * σ₃/σci + s)^a
        // 
        // We need to handle the case when mb * σ₃/σci + s < 0
        // which can happen in tension
        
        double arg = mb * sigma3 / sigma_ci + s;
        
        double yf;
        if (arg > 0) {
            // Standard Hoek-Brown
            yf = sigma1 - sigma3 - sigma_ci * pow(arg, a);
        } else {
            // In tension region - use tension cutoff behavior
            // The tensile strength is: σt = -s * σci / mb (for a = 0.5)
            // For general a, we use the limit as arg -> 0
            // When arg <= 0, we're beyond the tensile cutoff
            // Return a large positive value to indicate yielding
            double sigma_t = -s * sigma_ci / mb;  // Tensile strength (negative)
            yf = sigma1 - sigma3 - sigma_ci * s;  // Essentially tension cutoff
            // Alternative: linear extrapolation beyond tension cutoff
            // yf = sigma1 - sigma3;  // Pure shear criterion
        }

        return yf;
    }

    YIELD_FUNCTION_STRESS_DERIVATIVE
    {
        using namespace std;

        // Get Hoek-Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double mb = GET_PARAMETER_VALUE(HB_mb);
        double s = GET_PARAMETER_VALUE(HB_s);
        double a = GET_PARAMETER_VALUE(HB_a);
        double ds = GET_PARAMETER_VALUE(HB_ds);  // Perturbation for numerical derivative

        double sigma_norm = sigma.norm();
        ds = std::max(ds, ds * sigma_norm);

        // Use numerical differentiation for robustness
        // The analytical derivative involves principal stress directions which 
        // can be numerically sensitive
        auto computeNumericalDerivative = [this, &internal_variables_storage, &parameters_storage]
            (const VoigtVector& sig, double perturbation) -> VoigtVector {
            VoigtVector result;
            for (int i = 0; i < 6; ++i) {
                VoigtVector SIG1 = sig;
                VoigtVector SIG2 = sig;
                
                SIG1(i) += perturbation;
                SIG2(i) -= perturbation;

                double yf1 = YF(SIG1);
                double yf2 = YF(SIG2);

                result(i) = (yf1 - yf2) / (2.0 * perturbation);
            }
            return result;
        };

        if (ds > 0) {
            vv_out = computeNumericalDerivative(sigma, ds);
        } else {
            // Default to numerical with small perturbation
            vv_out = computeNumericalDerivative(sigma, 1e-8 * std::max(1.0, sigma_norm));
        }

        return vv_out;
    }

    YIELD_FUNCTION_HARDENING
    {
        // This model does not support hardening (perfectly plastic)
        return 0.0;
    }

    CHECK_APEX_REGION
    {
        using namespace std;

        // Get Hoek-Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double mb = GET_PARAMETER_VALUE(HB_mb);
        double s = GET_PARAMETER_VALUE(HB_s);
        double a = GET_PARAMETER_VALUE(HB_a);

        // The apex of the Hoek-Brown criterion in principal stress space
        // occurs at the tensile strength point
        // σt = -s * σci / mb (for a = 0.5)
        // For general a, it's approximately the same
        
        double sigma_t = -s * sigma_ci / mb;  // Tensile strength (negative value)
        
        // Get the mean stress
        double I1 = sigma.getI1();
        double p = I1 / 3.0;  // Mean stress (positive in compression)
        
        // Check if we're near the tensile apex
        // The apex is at p = sigma_t (all principal stresses equal to sigma_t)
        if (p < sigma_t * 1.01) {  // Small tolerance
            return true;
        }

        return false;
    }

    APEX_STRESS
    {
        using namespace std;

        // Get Hoek-Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double mb = GET_PARAMETER_VALUE(HB_mb);
        double s = GET_PARAMETER_VALUE(HB_s);

        // Tensile strength (apex stress state - hydrostatic tension)
        double sigma_t = -s * sigma_ci / mb;

        // Return hydrostatic stress state at tensile strength
        // (all normal stresses equal, no shear)
        vv_out = VoigtVector(sigma_t, sigma_t, sigma_t, 0, 0, 0);

        return vv_out;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<HB_sigci, HB_mb, HB_s, HB_a, HB_ds>;

private:
    static VoigtVector vv_out;  // For returning VoigtVectors
};

template <class NO_HARDENING>
VoigtVector HoekBrown_YF<NO_HARDENING>::vv_out;

// Declares this YF as featuring an apex
template<class NO_HARDENING>
struct yf_has_apex<HoekBrown_YF<NO_HARDENING>> : std::true_type {};

#endif
