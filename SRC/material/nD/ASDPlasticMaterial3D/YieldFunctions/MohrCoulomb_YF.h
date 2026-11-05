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

#ifndef MohrCoulomb_YF_H
#define MohrCoulomb_YF_H


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

#include "../YieldFunctionBase.h"
#include "cmath"
#include <iostream>
#include "../AllASDModelParameterTypes.h"
#include "../MohrCoulomb_Utils.h"


template<class NO_HARDENING>
class MohrCoulomb_YF : public YieldFunctionBase<MohrCoulomb_YF<NO_HARDENING>> // CRTP
{
public:

    static constexpr const char* NAME = "MohrCoulomb_YF";


    MohrCoulomb_YF( ):
        YieldFunctionBase<MohrCoulomb_YF<NO_HARDENING>>::YieldFunctionBase() 
        {}

    YIELD_FUNCTION 
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double theta_T = GET_PARAMETER_VALUE(MC_theta_T);
        return MohrCoulombUtils::function(sigma, phi, c, theta_T);
    }

    // central differences of the yield function above (corners rounded, see
    // MohrCoulomb_Utils.h), relative perturbation MC_ds
    YIELD_FUNCTION_STRESS_DERIVATIVE
    {  
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double ds = GET_PARAMETER_VALUE(MC_ds);
        double theta_T = GET_PARAMETER_VALUE(MC_theta_T);
        vv_out = MohrCoulombUtils::gradient(sigma, ds,
            [&](const VoigtVector& s) { return MohrCoulombUtils::function(s, phi, c, theta_T); });
        return vv_out;
    }

    YIELD_FUNCTION_HARDENING
    {
        // This model does not support hardening 
        return 0.0;
    }


    CHECK_APEX_REGION
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double I1 = sigma.getI1();
        double p = -I1 / 3;

        if (p < - c / tan(phi))
        {
            return true;
        }

        return false;
    }

    APEX_STRESS
    {
        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double p_apex = c / tan(phi);

        vv_out = VoigtVector(p_apex,p_apex,p_apex,0,0,0);

        return vv_out;
    }

  
    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds,MC_theta_T>;

private:


    static VoigtVector vv_out; //For returning VoigtVector's
};

template <class NO_HARDENING>
VoigtVector MohrCoulomb_YF<NO_HARDENING>::vv_out;

//Declares this YF as featuring an apex
template<class NO_HARDENING>
struct yf_has_apex<MohrCoulomb_YF<NO_HARDENING>> : std::true_type {};

#endif
