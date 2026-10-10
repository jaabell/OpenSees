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
// Hardening functions for the Stiff Soil (Hardening Soil) model

#ifndef StiffSoil_HardeningFunctions_H
#define StiffSoil_HardeningFunctions_H

#include "ASDPlasticMaterial3DGlobals.h"
#include "AllASDModelParameterTypes.h"
#include "HardeningFunction.h"

#include <cmath>
#include <tuple>

using namespace std;

/**
 * @brief Hardening policy for deviatoric plastic strain (shear mechanism)
 * 
 * The deviatoric yield function (PLAXIS Eq. 15.1) includes the plastic deviatoric
 * strain directly:
 * 
 *   F_s = q / (E_i * (1 - q/q_a)) - q / E_ur - ε_q^{p-shear} = 0
 * 
 * The evolution of the plastic deviatoric strain is:
 * 
 *   dε_q^{p-shear} = ||dev(dε^p)||_eq = √(2/3) * ||dev(m)|| * dλ
 * 
 * So the hardening function h = dε_q^{p-shear}/dλ = √(2/3) * ||dev(m)||
 * 
 * This is the equivalent deviatoric strain measure consistent with the
 * triaxial definition q = σ1 - σ3.
 */
struct StiffSoilShearHardeningPolicy {
    static constexpr const char* NAME = "StiffSoilShearHardening";
    
    HARDENING_FUNCTION_DEFINITION
    {
        // Extract deviatoric part of plastic flow direction
        VoigtVector m_dev = m.deviator();

        // Compute equivalent deviatoric strain increment
        // ||dev(m)||_eq = √(2/3) * √(dev(m) : dev(m))
        // Using stress-like inner product for the deviatoric tensor
        double m_dev_norm_sq = m_dev.v11() * m_dev.v11() + 
                               m_dev.v22() * m_dev.v22() + 
                               m_dev.v33() * m_dev.v33() +
                               2.0 * (m_dev.v12() * m_dev.v12() + 
                                      m_dev.v23() * m_dev.v23() + 
                                      m_dev.v13() * m_dev.v13());

        // h = √(2/3) * ||dev(m)||
        double h = sqrt(2.0 / 3.0) * sqrt(m_dev_norm_sq);

        return h;
        // return 0;
    }
    
    using parameters_t = tuple<>;  // No additional parameters needed
};


/**
 * @brief Hardening policy for cap pressure p_c (volumetric mechanism)
 * 
 * The cap hardening law relates volumetric plastic strain to p_c (PLAXIS Eq. 15.13):
 * 
 *   ε_v^{p-cap} = (β / (1-m)) * (p_c / p_ref)^{1-m}
 * 
 * Differentiating with respect to ε_v^{p-cap}:
 * 
 *   dε_v^{p-cap} = β * (p_c / p_ref)^{-m} * (dp_c / p_ref)
 * 
 * Rearranging for dp_c:
 * 
 *   dp_c = (p_ref / β) * (p_c / p_ref)^m * dε_v^{p-cap}
 * 
 * Since dε_v^{p-cap} = tr(m) * dλ, the hardening function is:
 * 
 *   h = dp_c/dλ = (p_ref / β) * (p_c / p_ref)^m * tr(m)
 * 
 * Note: tr(m) > 0 means volumetric compression (in mechanics convention),
 * which increases p_c (cap expands).
 */
struct StiffSoilCapHardeningPolicy {
    static constexpr const char* NAME = "StiffSoilCapHardening";
    
    HARDENING_FUNCTION_DEFINITION
    {
        // Get parameters
        double pref = GET_PARAMETER_VALUE(SS_pref);
        double beta = GET_PARAMETER_VALUE(SS_beta);
        double m_exp = GET_PARAMETER_VALUE(SS_m);  // Power exponent (note: 'm' is already used for flow direction)
        
        // Current cap pressure from the internal variable
        double pc = current_value.value();
        
        // Ensure p_c is positive and has a minimum value
        double pc_min = 0.1 * pref;
        pc = std::max(pc, pc_min);
        
        // Volumetric part of plastic flow direction
        // tr(m) = m_11 + m_22 + m_33
        double tr_m = m.trace();  
        
        // For the cap, we want compression (tr(m) > 0 in geotechnical convention)
        // to increase p_c. The flow direction from StiffSoilCap_PF gives
        // dG/dσ which in mechanics convention has tr(m) < 0 for compression.
        // We need to account for sign convention.
        // In geotechnical convention (compression positive), volumetric plastic
        // strain is positive for compression.
        
        // h = (p_ref / β) * (p_c / p_ref)^m * |tr(m)|
        // We use |tr(m)| because the cap only hardens under compression
        double ratio = pc / pref;
        double h = (pref / beta) * pow(ratio, m_exp) * std::abs(tr_m);
       
        // Cap should only expand (p_c increase) under compression
        // If tr(m) indicates tension/expansion, no hardening occurs
        // In mechanics convention (tension positive), tr(m) < 0 means compression
        // So we want h > 0 when tr(m) < 0
        if (tr_m > 0) {
            // Tension/expansion in mechanics convention - no cap hardening
            h = 0.0;
        }

        return h;
    }
    
    using parameters_t = tuple<SS_pref, SS_beta, SS_m>;
};


/**
 * @brief Alternative cap hardening with explicit volumetric strain tracking
 * 
 * This version allows for more control over the cap evolution by using
 * a simpler linear relationship between volumetric plastic strain and p_c.
 * 
 *   dp_c = H_cap * dε_v^p
 * 
 * where H_cap is a hardening modulus that can be stress-dependent.
 */
struct StiffSoilCapLinearHardeningPolicy {
    static constexpr const char* NAME = "StiffSoilCapLinearHardening";
    
    HARDENING_FUNCTION_DEFINITION
    {
        // Get parameters
        double pref = GET_PARAMETER_VALUE(SS_pref);
        double beta = GET_PARAMETER_VALUE(SS_beta);
        
        // Current cap pressure
        double pc = current_value.value();
        pc = std::max(pc, 0.1 * pref);
        
        // Volumetric part of plastic flow
        double tr_m = m.trace();
        
        // Simple linear hardening: H_cap = p_ref / beta
        double H_cap = pref / beta;
        
        // h = H_cap * |tr(m)|, only for compression
        double h = (tr_m < 0) ? H_cap * std::abs(tr_m) : 0.0;
        
        return h;
    }
    
    using parameters_t = tuple<SS_pref, SS_beta>;
};


// ============================================================================
// Aliases for HardeningFunction with specific policies
// ============================================================================

// Deviatoric (shear) mechanism hardening
using StiffSoilShearHardeningFunction = HardeningFunction<VoigtScalar, StiffSoilShearHardeningPolicy>;

// Cap (volumetric) mechanism hardening - nonlinear (PLAXIS formulation)
using StiffSoilCapHardeningFunction = HardeningFunction<VoigtScalar, StiffSoilCapHardeningPolicy>;

// Cap mechanism hardening - linear (simplified)
using StiffSoilCapLinearHardeningFunction = HardeningFunction<VoigtScalar, StiffSoilCapLinearHardeningPolicy>;


#endif // StiffSoil_HardeningFunctions_H


// ============================================================================
// IMPORTANT: Add the following to AllASDInternalVariableTypes.h
// ============================================================================
// 
// // Include the hardening functions
// #include "StiffSoil_HardeningFunctions.h"
// 
// // Plastic deviatoric strain for shear mechanism
// struct EpsQpShearName { static constexpr const char* name = "EpsQpShear"; };
// using EpsQpShear = InternalVariableType<VoigtScalar, StiffSoilShearHardeningFunction, EpsQpShearName>;
// 
// // Cap pressure for volumetric mechanism  
// struct CapPressureName { static constexpr const char* name = "CapPressure"; };
// using CapPressure = InternalVariableType<VoigtScalar, StiffSoilCapHardeningFunction, CapPressureName>;
// 
// // Cap pressure with linear hardening (alternative)
// using CapPressureLinear = InternalVariableType<VoigtScalar, StiffSoilCapLinearHardeningFunction, CapPressureName>;
// ============================================================================
