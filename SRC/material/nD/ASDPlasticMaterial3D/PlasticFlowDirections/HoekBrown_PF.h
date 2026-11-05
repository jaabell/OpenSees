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
// ASDPlasticMaterial3D - Hoek-Brown Plastic Flow Direction
//
// Implements a non-associated plastic flow rule for the Hoek-Brown criterion.
// The plastic potential uses a similar form to the yield function but with
// a different (typically lower) dilation parameter to control volumetric
// plastic strain.
//
// The plastic potential is:
//   g = σ₁ - σ₃ - σci * (mb_psi * σ₃/σci + s)^a
//
// Where mb_psi = mi_psi * exp((GSI-100)/(28-14D)) with mi_psi typically
// corresponding to a dilation angle ψ instead of friction angle φ.
//
// For rock masses, non-associated flow is essential because associated
// flow would predict unrealistically large dilation.
//
// Sign convention: Compression is POSITIVE (geomechanics convention)

#ifndef HoekBrown_PF_H
#define HoekBrown_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
#include "../AllASDModelParameterTypes.h"
#include "../HoekBrown_Utils.h"
#include "../HoekBrown_Invariants.h"
#include <cmath>
#include <typeinfo>

template<class NO_HARDENING>
class HoekBrown_PF : public PlasticFlowBase<HoekBrown_PF<NO_HARDENING>> // CRTP
{
public:

    static constexpr const char* NAME = "HoekBrown_PF";

    HoekBrown_PF():
        PlasticFlowBase<HoekBrown_PF<NO_HARDENING>>::PlasticFlowBase()
    {}

    PLASTIC_FLOW_DIRECTION
    {
        using namespace std;

        // Get Hoek-Brown parameters
        double sigma_ci = GET_PARAMETER_VALUE(HB_sigci);
        double mb_psi = GET_PARAMETER_VALUE(HB_mb_psi);  // Dilation parameter (reduced mb)
        double mb_unused, s, a;
        HoekBrownUtils::effective_rock_mass_parameters(
            GET_PARAMETER_VALUE(HB_mi), GET_PARAMETER_VALUE(HB_GSI), GET_PARAMETER_VALUE(HB_D),
            0.0, GET_PARAMETER_VALUE(HB_s), GET_PARAMETER_VALUE(HB_a),
            mb_unused, s, a);
        double ds = GET_PARAMETER_VALUE(HB_ds);
        double theta_T = GET_PARAMETER_VALUE(HB_theta_T);

        // gradient of the yield function's shear part with mb_psi in place of mb
        // (compression positive like the yield function, edges rounded; see
        // HoekBrown_Invariants.h)
        vv_out = HoekBrownInvariants::gradient(sigma, ds,
            [&](const VoigtVector& sig) { return HoekBrownInvariants::shear_function(sig, sigma_ci, mb_psi, s, a, theta_T); });

        return vv_out;
    }

    using internal_variables_t = std::tuple<NO_HARDENING>;

    using parameters_t = std::tuple<HB_sigci, HB_mb_psi, HB_GSI, HB_mi, HB_D, HB_s, HB_a, HB_ds, HB_theta_T>;

private:
    static VoigtVector vv_out;
};

template<class NO_HARDENING>
VoigtVector HoekBrown_PF<NO_HARDENING>::vv_out;

#endif
