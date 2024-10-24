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

#ifndef MohrCoulomb_PF_H
#define MohrCoulomb_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"

#include <cmath>
#include <typeinfo>


template<class NO_HARDENING>
class MohrCoulomb_PF : public PlasticFlowBase<MohrCoulomb_PF<NO_HARDENING>>// CRTP
{
public:

    static constexpr const char* NAME = "MohrCoulomb_PF";


    MohrCoulomb_PF( ):
        PlasticFlowBase<MohrCoulomb_PF<NO_HARDENING>>::PlasticFlowBase()  // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                { }

    double g(const VoigtVector& sigma, double phi) const
    {
        double I1 = sigma.getI1();
        double J2 = sigma.getJ2();
        double lode_angle = sigma.lodeAngle();

        double rEquivalentStress = (std::cos(lode_angle) - std::sin(lode_angle) * std::sin(phi) / std::sqrt(3.0))  * std::sqrt(J2) +
            I1 * std::sin(phi) / 3.0;

        return rEquivalentStress;
    }

    PLASTIC_FLOW_DIRECTION
    {

        double phi = GET_PARAMETER_VALUE(MC_phi)*M_PI/180;
        // double c = GET_PARAMETER_VALUE(MC_c);
        double ds = GET_PARAMETER_VALUE(MC_ds);

        using namespace std;

        double sigma_norm = sigma.norm();

        // ds = ds == 0 ? 1e-8 : ds;

        ds = std::max(ds, ds*sigma_norm);

        VoigtVector dg;

        for (int i = 0; i < 6; ++i) {
            VoigtVector SIG1 = sigma;
            VoigtVector SIG2 = sigma;
            
            // Increment SIG at index i by a small amount for numerical differentiation
            SIG1(i) += ds;
            SIG2(i) -= ds;

            // Compute the yield function at the perturbed state
            double g1 = g(SIG1, phi);
            double g2 = g(SIG2, phi);

            // Calculate the derivative
            dg(i) = (g1 - g2) / (2*ds);
        }

        double D = GET_PARAMETER_VALUE(Dilatancy);

        vv_out = dg.deviator();
        vv_out -= D * kronecker_delta() / 3;


        return vv_out;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;
    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds,Dilatancy>;

private:

    static VoigtVector vv_out; 

};


template<class NO_HARDENING>
VoigtVector MohrCoulomb_PF<NO_HARDENING>::vv_out;

#endif
