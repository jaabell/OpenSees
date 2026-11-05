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
#include "../MohrCoulomb_Utils.h"
#include <cmath>
#include <typeinfo>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

template<class NO_HARDENING>
class MohrCoulomb_PF : public PlasticFlowBase<MohrCoulomb_PF<NO_HARDENING>>// CRTP
{
public:

    static constexpr const char* NAME = "MohrCoulomb_PF";

    MohrCoulomb_PF( ):
        PlasticFlowBase<MohrCoulomb_PF<NO_HARDENING>>::PlasticFlowBase()
                { }
    // Mohr-Coulomb potential: the yield function's form with the dilatancy angle
    // psi in place of phi (corners rounded the same way, see MohrCoulomb_Utils.h)
    PLASTIC_FLOW_DIRECTION
    {
        double psi = GET_PARAMETER_VALUE(MC_psi)*M_PI/180;
        double c = GET_PARAMETER_VALUE(MC_c);
        double ds = GET_PARAMETER_VALUE(MC_ds);
        double theta_T = GET_PARAMETER_VALUE(MC_theta_T);
        vv_out = MohrCoulombUtils::gradient(sigma, ds,
            [&](const VoigtVector& s) { return MohrCoulombUtils::function(s, psi, c, theta_T); });
        return vv_out;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<MC_phi,MC_c,MC_ds, MC_psi, MC_theta_T>;

private:

    static VoigtVector vv_out; 
};

template<class NO_HARDENING>
VoigtVector MohrCoulomb_PF<NO_HARDENING>::vv_out;

#endif
