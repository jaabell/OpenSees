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

#ifndef VonMises_PF_H
#define VonMises_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"

#include <cmath>
#include <typeinfo>


using namespace std;

template<class AlphaHardeningType>
class VonMises_PF : public PlasticFlowBase<VonMises_PF<AlphaHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "VonMises_PF";

    VonMises_PF():
        PlasticFlowBase<VonMises_PF<AlphaHardeningType >>::PlasticFlowBase()  // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                { }

    PLASTIC_FLOW_DIRECTION
    {
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);

        auto s = sigma.deviator();
        vv_out = s - alpha;

        double den = sqrt(tensor_dot_stress_like(vv_out, vv_out));
        if (abs(den) > sqrt(tensor_dot_stress_like(s, s))*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            vv_out = vv_out / den;

        return vv_out;
    }

    // m = r / |r|, r = dev(sigma) - alpha:  dm/dsigma = (I - m strain_like(m)^T) P / |r|,
    // with P = d dev(sigma) / dsigma
    static constexpr bool has_dm_dsigma = true;
    PLASTIC_FLOW_DIRECTION_STRESS_DERIVATIVE
    {
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        VoigtVector r = sigma.deviator() - alpha;
        double rho = sqrt(tensor_dot_stress_like(r, r));
        VoigtMatrix dm = VoigtMatrix::Zero();
        if (!(rho > sqrt(tensor_dot_stress_like(sigma.deviator(), sigma.deviator()))*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON))
            return dm;
        VoigtVector mm = r / rho;
        VoigtMatrix P = VoigtMatrix::Identity();
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                P(a, b) -= 1.0 / 3.0;
        dm = (VoigtMatrix::Identity() - mm * strain_like(mm).transpose()) * P / rho;
        return dm;
    }

    using internal_variables_t = std::tuple<AlphaHardeningType>;
    using parameters_t = std::tuple<>;

private:

    static VoigtVector vv_out; 
};

template<class AlphaHardeningType>
VoigtVector VonMises_PF<AlphaHardeningType  >::vv_out;


#endif