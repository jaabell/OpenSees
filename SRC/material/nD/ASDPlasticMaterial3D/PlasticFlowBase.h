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

#ifndef PlasticFlowBase_H
#define PlasticFlowBase_H

#include <Channel.h>

#include "EigenAPI.h"


// Helper template to check if a class has a parameters_t type alias
template <typename T, typename = void>
struct pf_has_parameters_t : std::false_type {};

template <typename T>
struct pf_has_parameters_t<T, typename std::enable_if<!std::is_same<typename T::parameters_t, void>::value>::type> : std::true_type {};

// Helper template to check if a class has a internal_variables_t type alias
template <typename T, typename = void>
struct pf_has_internal_variables_t : std::false_type {};

template <typename T>
struct pf_has_internal_variables_t<T, typename std::enable_if<!std::is_same<typename T::internal_variables_t, void>::value>::type> : std::true_type {};


// A flow direction that also provides its closed-form derivative declares
//     static constexpr bool has_dm_dsigma = true;
// and a PLASTIC_FLOW_DIRECTION_STRESS_DERIVATIVE member returning
// dm_i/dsigma_j (m in tensor components, sigma_j the Voigt storage entry).
// Full_Backward_Euler uses it with derivatives Closed_Form.
template <typename T, typename = void>
struct pf_has_dm_dsigma : std::false_type {};
template <typename T>
struct pf_has_dm_dsigma<T, typename std::enable_if<T::has_dm_dsigma>::type> : std::true_type {};

#define PLASTIC_FLOW_DIRECTION_STRESS_DERIVATIVE template <typename StorageType, typename ParameterStorageType> \
    VoigtMatrix dm_dsigma_ij( \
        const VoigtVector &depsilon, \
        const VoigtVector& sigma, \
        const StorageType& internal_variables_storage,  \
        const ParameterStorageType& parameters_storage) const

#define PLASTIC_FLOW_DIRECTION template <typename StorageType, typename ParameterStorageType> \
    const VoigtVector& operator()( \
        const VoigtVector &depsilon, \
        const VoigtVector& sigma, \
        const StorageType& internal_variables_storage,  \
        const ParameterStorageType& parameters_storage) const

#define GET_TRIAL_INTERNAL_VARIABLE(type) \
    ((internal_variables_storage).template get<type>().trial_value)

#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value


template <class T>
class PlasticFlowBase
{
public:
    PlasticFlowBase() { 
        static_assert(pf_has_parameters_t<T>::value, "Derived class must have a 'parameters_t' type alias.");
        static_assert(pf_has_internal_variables_t<T>::value, "Derived class must have a 'internal_variables_t' type alias.");
    }

    PLASTIC_FLOW_DIRECTION
    {
        return static_cast<T*>(this)->operator()( depsilon,  sigma, internal_variables_storage, parameters_storage);
    }

    inline const char* getName() const { return static_cast<T*>(this)->NAME; }
};

#endif
