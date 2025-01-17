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



#ifndef ASDPlasticMaterial3D_H
#define ASDPlasticMaterial3D_H


#include "NDMaterial.h"
#include <G3Globals.h>
#include <iostream>
#include <Channel.h>
#include <Information.h>
#include <MaterialResponse.h>

#include "ASDPlasticMaterial3DTraits.h"

#include "ASDPlasticMaterial3DGlobals.h"


// #include "MaterialInternalVariables.h"
#include "YieldFunctions/AllYieldFunctions.h"
#include "PlasticFlowDirections/AllPlasticFlowDirections.h"
#include "ElasticityModels/AllElasticityModels.h"
#include "AllASDModelParameterTypes.h"
#include "AllASDInternalVariableTypes.h"
#include "AllASDHardeningFunctions.h"

#include "utuple_storage.h"

// for print p, q, theta
#include <tuple>
#include <utility> // For std::pair
#include <map> // For std::pair
#include <limits>
#include <type_traits>

#include "std_tuple_concat.h"

// for debugging printing
#include <fstream>

#define ASDPlasticMaterial3D_MAXITER_BRENT 50


using namespace ASDPlasticMaterial3DGlobals;

#define ASDP_TAG this->getTag()

template <
    class ElasticityType,
    class YieldFunctionType,
    class PlasticFlowType,
    int thisClassTag >
class ASDPlasticMaterial3D : public NDMaterial
{

public:

    // Concatenate the internal varibles into the storage
    using iv_concat_types = utuple_concat_unique_type <
                            typename YieldFunctionType::internal_variables_t,
                            typename PlasticFlowType::internal_variables_t >;
    using iv_storage_t = utuple_storage<iv_concat_types>;

    // Concatenate the model parameters into the parameters storage
    using extracted_parameters_t = utuple_concat_unique_type <
                                   ExtractNestedParameterTypes_t<typename YieldFunctionType::internal_variables_t>,
                                   ExtractNestedParameterTypes_t<typename PlasticFlowType::internal_variables_t>
                                   >;

    using parameters_concat_types = utuple_concat_type <
                                    typename YieldFunctionType::parameters_t,
                                    typename PlasticFlowType::parameters_t,
                                    typename ElasticityType::parameters_t,
                                    extracted_parameters_t,
                                    std::tuple<MassDensity>,
                                    std::tuple<InitialP0>
                                    >;
    using parameters_storage_t = utuple_storage<parameters_concat_types>;


    //==================================================================================================
    //  Constructors
    //==================================================================================================

    ASDPlasticMaterial3D( )
        : NDMaterial(0, thisClassTag)
    {

    }


    ASDPlasticMaterial3D(int tag)
        : NDMaterial(tag, thisClassTag)
    {


        TrialStress *= 0;
        CommitStress *= 0;
        TrialStrain *= 0;
        CommitStrain *= 0;
        TrialPlastic_Strain *= 0;
        CommitPlastic_Strain *= 0;

        first_step = true;
    }


    ~ASDPlasticMaterial3D(void)
    {

    }

    //==================================================================================================
    // To set internal variables values for the model
    //==================================================================================================
    auto getInternalVariablesNames() const
    {
        return iv_storage.getParameterNames();
    }

    int getInternalVariableSizeByName(const char * iv_name) const
    {
        return iv_storage.getInternalVariableSizeByName(iv_name);
    }

    int getInternalVariableIndexByName(const char * iv_name) const
    {
        return iv_storage.getInternalVariableIndexByName(iv_name);
    }

    auto setInternalVariableByName(const char * iv_name, int iv_size, double* iv_values)
    {
        cout << "  --->  Setting " << iv_name << " = ";
        for (int i = 0; i < iv_size; ++i)
        {
            cout << iv_values[i] << " ";
        }
        cout << endl;
        return iv_storage.setInternalVariableByName(iv_name, iv_size, iv_values);
    }



    //==================================================================================================
    // To set parameter values for the model
    //==================================================================================================
    auto getParameterNames() const
    {
        return parameters_storage.getParameterNames();
    }

    auto setParameterByName(const char * param_name, double param_value)
    {
        cout << "  --->  Setting " << param_name << " = " << param_value << endl;
        return parameters_storage.setParameterByName(param_name, param_value);
    }


    //==================================================================================================
    //  Class type function
    //==================================================================================================
    const char *getClassType(void) const
    {
        std::string name("ASDPlasticMaterial3D");

        return name.c_str();
    };

    double getRho(void)
    {
        return parameters_storage.template get<MassDensity>().value;
    }

    double getPressure(void)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        return -CommitStress.trace() / 3;
    }

    std::string getYFName() const {return YieldFunctionType::NAME;}
    std::string getPFName() const {return PlasticFlowType::NAME;}
    std::string getELName() const {return ElasticityType::NAME;}
    std::string getIVName() const {return iv_storage.getVariableNamesAndHardeningLaws();}


    //==================================================================================================
    //  Set Trial strain and trial strain increment
    //==================================================================================================
    // For total strain-based elements.
    // Receives the current total strain at a GP.
    // This function then computes the incremental strain (subtracting from the committed one)
    // and sets the increment.
    // Returns a success flag from the call to setTrialStrainIncr

    int setTrialStrain(const Vector &v)
    {

        if (first_step)
        {
            double p0 = parameters_storage.template get<InitialP0>().value;
            TrialStress(0) = p0;
            TrialStress(1) = p0;
            TrialStress(2) = p0;
            CommitStress(0) = p0;
            CommitStress(1) = p0;
            CommitStress(2) = p0;
        }


        TrialStrain = VoigtVector::fromStrain(v);
        return setTrialStrainIncr( TrialStrain - CommitStrain );
    }

    int setTrialStrainIncr( const VoigtVector &strain_increment )
    {

        int exitflag = -1;

        switch (INT_OPT_constitutive_integration_method[ASDP_TAG])
        {
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Not_Set :
            exitflag = -1;
            cerr << "CEP::setTrialStrainIncr - Integration method not set!\n" ;
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler :
            exitflag = this->Forward_Euler(strain_increment);
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Subincrement :
            exitflag = this->Forward_Euler_Subincrement(strain_increment);
            break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler :
        //     exitflag = this->Backward_Euler(strain_increment);;
        //     break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_ddlambda :
        //     exitflag = this->Backward_Euler_ddlambda(strain_increment);;
        //     break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_ddlambda_Subincrement :
        //     exitflag = this->Backward_Euler_ddlambda_Subincrement(strain_increment);;
        //     break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Crisfield :
        //     exitflag = this->Forward_Euler(strain_increment, true);
        //     break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Multistep_Forward_Euler_Crisfield :
        //     exitflag = this->Multistep_Forward_Euler(strain_increment, true);
        //     break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Modified_Euler_Error_Control :
        //     exitflag = this->Modified_Euler_Error_Control(strain_increment);
        //     break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control :
            exitflag = this->Runge_Kutta_45_Error_Control(strain_increment);;
            break;
        default:
            cerr << "ASDPlasticMaterial3D::setTrialStrainIncr - Integration method not available!\n" ;
            exitflag = -1;
        }

        return exitflag;
    }

//==================================================================================================
//  Getters
//==================================================================================================

    const VoigtMatrix &getTangentTensor( void )
    {

        return Stiffness;
    }

    const Vector &getStress(void)
    {
        static Vector result(6);
        TrialStress.toStress(result);
        return result;
    }



    const Vector &getStrain(void)
    {
        static Vector result(6);
        TrialStrain.toStrain(result);
        return result;
    }

    const Vector &getPstrain(void)
    {
        static Vector result(6);
        TrialPlastic_Strain.toStress(result);
        return result;
    }

    const Vector &getEQPstrain(void)
    {
        static Vector result(1);
        result(0) = std::sqrt(2./3.*TrialPlastic_Strain.squaredNorm());
        // opserr << "JOSE: getEQPstrain(void) result = " << result << endln;
        return result;
    }

    const Vector &getPStress(void)
    {
        static Vector result(1);
        result(0) = TrialStress.meanStress();
        // opserr << "JOSE: getPStress(void) result = " << result << endln;
        return result;
    }

    const Vector &getJ2Stress(void)
    {
        static Vector result(1);
        result(0) = TrialStress.getJ2();
        // opserr << "JOSE: getJ2Stress(void) result = " << result << endln;
        return result;
    }

    const Vector &getVolStrain(void)
    {
        static Vector result(1);
        result(0) = TrialStrain.getI1();
        // opserr << "JOSE: getVolStrain(void) result = " << result << endln;
        return result;
    }

    const Vector &getJ2Strain(void)
    {
        static Vector result(1);
        result(0) = TrialStrain.getJ2();
        // opserr << "JOSE: getJ2Strain(void) result = " << result << endln;
        return result;
    }

    const Vector &getInternalVariableByPos(int pos)
    {
        static Vector return_vector(6);
        int find_pos = 0;

        // Note by J. Abell on Wed 06 Dec 2023 11:44:24
        //
        // The lambda capture of return_vector which is declared static above triggers
        // the following warning in g++ 11.4.0
        //
        // warning: capture of variable ‘return_vector’ with non-automatic storage duration
        //
        // Explanation:
        //   Because usually the lambda functions are used for threaded applications
        //   this warning is there to help with race conditions on the return_vector.
        //   In this case the warning is benign because opensees is not threaded
        //   at this level.
        //   We want to keep the static allocation of return_vector to avoid
        //   many calls to malloc (new) every time this function is called
        //   for performance reasons, so we have to live with the warning.
        iv_storage.apply([&pos, &find_pos, this](auto & internal_variable)
        {
            if (pos == find_pos)
            {
                auto &iv = internal_variable.trial_value;
                int iv_size = iv.size();
                return_vector.resize(iv_size);
                for (int i = 0; i < iv_size; ++i)
                {
                    return_vector(i) = iv(i);
                }
            }
            find_pos ++;
        });

        return return_vector;
    }

    const VoigtVector &getStressTensor( void )
    {
        return TrialStress;
    }

    const VoigtVector &getStrainTensor( void )
    {
        return TrialStrain;
    }

    const VoigtVector &getPlasticStrainTensor( void )
    {
        return TrialPlastic_Strain;
    }

    const VoigtVector  &getCommittedStressTensor(void)
    {
        return CommitStress;
    }

    const VoigtVector &getCommittedStrainTensor(void)
    {
        return CommitStrain;
    }

    const VoigtVector &getCommittedPlasticStrainTensor(void)
    {
        return CommitPlastic_Strain;
    }

    void ComputeTangentStiffness()
    {
        if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Elastic)
        {
            VoigtMatrix Eelastic = et(CommitStress, parameters_storage);
            Stiffness = Eelastic;
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_FirstOrder)
        {
            compute_numerical_tangent_firstorder(TrialStrain-CommitStrain, Stiffness);
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_SecondOrder)
        {
            compute_numerical_tangent_secondorder(TrialStrain-CommitStrain, Stiffness);
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::
            Continuum)
        {

            VoigtMatrix Eelastic = et(TrialStress, parameters_storage);
            const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);

            double den_after_corrector = n.transpose() * Eelastic * m - hardening;

            VoigtMatrix Econtinuum = Eelastic - Eelastic * m * (n.transpose() * Eelastic) / den_after_corrector;
            Stiffness = Econtinuum;
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Secant)
        {

            VoigtMatrix Eelastic = et(TrialStress, parameters_storage);
            const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);

            double den_after_corrector = n.transpose() * Eelastic * m - hardening;

            VoigtMatrix Econtinuum = Eelastic - Eelastic * m * (n.transpose() * Eelastic) / den_after_corrector;
            Stiffness = (Econtinuum + Eelastic)/2;
        }
    }

    int compute_local_stress(
        const VoigtVector& local_stress, const VoigtVector& local_strain,
        const VoigtVector& strain_incr, VoigtVector& stress_incr) const
    {
        using namespace ASDPlasticMaterial3DGlobals;

        // Initialize local variables for stress and strain increments
        VoigtVector depsilon = strain_incr;  // Strain increment (perturbation)
        VoigtVector dsigma = VoigtVector();  // Stress increment
        VoigtVector trial_stress = VoigtVector();  // Trial stress

        // Compute elastic stiffness matrix based on the local stress and parameters
        VoigtMatrix Eelastic = et(local_stress, parameters_storage);  // Elasticity tensor
        
        // Compute the elastic stress increment: dsigma = E * depsilon
        dsigma = Eelastic * depsilon;
        
        // Compute the trial stress
        trial_stress = local_stress + dsigma;

        // // Evaluate the yield function for the current stress and trial stress
        double yf_val_start = yf(local_stress, iv_storage, parameters_storage);
        double yf_val_end = yf(trial_stress, iv_storage, parameters_storage);

        // Check if the material response is elastic or plastic
        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) {
            // Elastic response: no plastic correction
            stress_incr = dsigma;  // Stress increment is purely elastic
        } else {
            // Plastic response: need to apply plastic correction
            // Compute plastic correction by finding the intersection of the yield surface
            VoigtVector intersection_stress = local_stress;
            VoigtVector intersection_strain = local_strain;
            depsilon_elpl = strain_incr;

            if (yf_val_start < 0) {
                // Find the intersection of the yield surface between the start and trial stress
                double tol_yf = DOUBLE_OPT_f_absolute_tol[ASDP_TAG];
                double intersection_factor = compute_yf_crossing(
                    local_stress, trial_stress, 0.0, 1.0, tol_yf);
                
                // Ensure the intersection factor is within valid bounds [0, 1]
                intersection_factor = std::max(0.0, std::min(1.0, intersection_factor));

                // Compute the intersection stress and strain
                intersection_stress = local_stress * (1 - intersection_factor) +
                                      trial_stress * intersection_factor;
                intersection_strain = local_strain + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * strain_incr;
            }

            // The trial stress is updated based on the intersection
            // stress_incr = intersection_stress - local_stress;
            trial_stress = intersection_stress;

            Eelastic = et(intersection_stress, parameters_storage);
            trial_stress  += Eelastic * depsilon_elpl;

            //Compute normal to YF (n) and Plastic Flow direction (m)
            const VoigtVector& n = yf.df_dsigma_ij(intersection_stress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, intersection_stress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  intersection_stress, iv_storage, parameters_storage);
            double den = n.transpose() * Eelastic * m - hardening;

            //Compute the plastic multiplier
            if (abs(den) < MACHINE_EPSILON)
            {
                cout << "CEP - den = 0\n";
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor1("m", m);
                printTensor1("n", n);
                cout << "hardening = " << hardening << endl;
                cout << "den = " << den << endl;
                printTensor1("depsilon_elpl", depsilon_elpl);
                return -1;
            }

            double dLambda =  n.transpose() * Eelastic * depsilon_elpl;
            dLambda /= den;

            if (dLambda <= 0)
            {
                // cout << "CEP - dLambda = " << dLambda << " <= 0\n";
                // printTensor1("m", m);
                // printTensor1("n", n);
                // cout << "hardening = " << hardening << endl;
                // cout << "den = " << den << endl;
                // printTensor1("depsilon_elpl", depsilon_elpl);
                dLambda = 0;
            }

            trial_stress = trial_stress - dLambda * Eelastic * m;
        }

        stress_incr = trial_stress - CommitStress;

        return 0;  // Return success
    }



    int compute_numerical_tangent_firstorder(
        const VoigtVector& strain_incr, VoigtMatrix& tangent_matrix, double epsilon_ref = 1e-8, double delta_min = 1e-12)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        // cout << "strain_incr = " << strain_incr.transpose() << endl;
        // cout << "epsilon_ref = " << epsilon_ref << endl;
        // cout << "delta_min = " << delta_min << endl;

        // Number of strain and stress components (Voigt notation in 3D: 6 components)
        const int n = 6;

        // Initialize the local copies of stress and strain
        VoigtVector local_stress = CommitStress;
        VoigtVector local_strain = CommitStrain;

        // Allocate memory for perturbed stress and strain
        VoigtVector perturbed_stress;
        VoigtVector perturbed_strain;

        // Compute initial stress increment for the given strain increment (unperturbed)
        VoigtVector initial_stress_incr;
        compute_local_stress(local_stress, local_strain, strain_incr, initial_stress_incr);

        // cout << "local_stress = " << local_stress.transpose() << endl;
        // cout << "initial_stress_incr = " << initial_stress_incr.transpose() << endl;
        double delta = std::max(epsilon_ref * strain_incr.norm(), delta_min);
        // cout << "delta = " << delta << endl;

        // Loop over each strain component to compute the tangent matrix via finite differences
        for (int i = 0; i < n; ++i) {
            // Compute adaptive delta based on the current strain component

            // Perturb the i-th strain component by the adaptive delta
            VoigtVector strain_incr_perturbed = strain_incr;
            strain_incr_perturbed(i) += delta;

            // Compute the local stress for the perturbed strain increment
            compute_local_stress(local_stress, local_strain, strain_incr_perturbed, perturbed_stress);

            // cout << "        strain_incr_perturbed = " << strain_incr_perturbed.transpose() << endl;
            // cout << "        perturbed_stress = " << perturbed_stress.transpose() << endl;
            // Finite difference approximation of the tangent matrix (column i)
            for (int j = 0; j < n; ++j) {
                tangent_matrix(j, i) = (perturbed_stress(j) - initial_stress_incr(j)) / delta;
            }
        }

        // cout << "tangent_matrix = \n" << tangent_matrix << endl;

        return 0; // Return success
    }



    int compute_numerical_tangent_secondorder(
        const VoigtVector& strain_incr, VoigtMatrix& tangent_matrix, double epsilon_ref = 1e-8, double delta_min = 1e-12)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        // cout << "strain_incr = " << strain_incr.transpose() << endl;
        // cout << "epsilon_ref = " << epsilon_ref << endl;
        // cout << "delta_min = " << delta_min << endl;

        // Number of strain and stress components (Voigt notation in 3D: 6 components)
        const int n = 6;

        // Initialize the local copies of stress and strain
        VoigtVector local_stress = CommitStress;
        VoigtVector local_strain = CommitStrain;

        // Allocate memory for perturbed stress and strain
        VoigtVector perturbed_stress1;
        VoigtVector perturbed_stress2;
        VoigtVector perturbed_strain;

        // Compute initial stress increment for the given strain increment (unperturbed)
        // VoigtVector initial_stress_incr;
        // compute_local_stress(local_stress, local_strain, strain_incr, initial_stress_incr);

        // cout << "local_stress = " << local_stress.transpose() << endl;
        // cout << "initial_stress_incr = " << initial_stress_incr.transpose() << endl;
        double delta = std::max(epsilon_ref * strain_incr.norm(), delta_min);
        // cout << "delta = " << delta << endl;

        // Loop over each strain component to compute the tangent matrix via finite differences
        for (int i = 0; i < n; ++i) {
            // Compute adaptive delta based on the current strain component

            // Perturb the i-th strain component by the adaptive delta
            VoigtVector strain_incr_perturbed1 = strain_incr;
            VoigtVector strain_incr_perturbed2 = strain_incr;
            strain_incr_perturbed1(i) += delta;
            strain_incr_perturbed2(i) -= delta;

            // Compute the local stress for the perturbed strain increment
            compute_local_stress(local_stress, local_strain, strain_incr_perturbed1, perturbed_stress1);
            compute_local_stress(local_stress, local_strain, strain_incr_perturbed2, perturbed_stress2);

            // cout << "        strain_incr_perturbed = " << strain_incr_perturbed.transpose() << endl;
            // cout << "        perturbed_stress = " << perturbed_stress.transpose() << endl;
            // Finite difference approximation of the tangent matrix (column i)
            for (int j = 0; j < n; ++j) {
                tangent_matrix(j, i) = (perturbed_stress1(j) - perturbed_stress2(j)) / (2*delta);
            }
        }

        // cout << "tangent_matrix = \n" << tangent_matrix << endl;

        return 0; // Return success
    }

    const Matrix& getTangent()
    {
        static Matrix return_matrix(6, 6);

        copyToMatrixReference(Stiffness, return_matrix);

        return return_matrix;
    }


    const Matrix& getInitialTangent()
    {
        static Matrix return_matrix(6, 6);

        VoigtMatrix Eelastic = et(CommitStress, parameters_storage);
        Stiffness = Eelastic;

        copyToMatrixReference(this->Stiffness, return_matrix);

        return return_matrix;
    }



//==================================================================================================
//  State commiting and reversion
//==================================================================================================


    int commitState(void)
    {

        CommitStress = TrialStress;
        CommitStrain = TrialStrain;
        CommitPlastic_Strain = TrialPlastic_Strain;

        iv_storage.commit_all();

        if (first_step)
        {
            first_step = false;
        }

        if (GLOBAL_INT_max_iter[ASDP_TAG] > 0 || GLOBAL_DBL_max_error[ASDP_TAG] > 0.)
        {
            cout << "  () ASDP Integration Info. Tag = " << ASDP_TAG << " max_iter = " << GLOBAL_INT_max_iter[ASDP_TAG] << " max_error = " << GLOBAL_DBL_max_error[ASDP_TAG] << endl;
            GLOBAL_INT_max_iter[ASDP_TAG] = 0;
            GLOBAL_DBL_max_error[ASDP_TAG] = 0.;
        }

        return 0;
    }

    //Reverts the commited variables to the trials and calls revert on BET Classes.
    int revertToLastCommit(void)
    {

        TrialStress = CommitStress;
        TrialStrain = CommitStrain;
        TrialPlastic_Strain = CommitPlastic_Strain;

        iv_storage.revert_all();

        return 0;
    }

    int revertToStart(void)
    {
        cerr << "ASDPlasticMaterial3D::revertToStart - not implemented!!!\n" ;
        return -1;
    }

    NDMaterial *getCopy(void)
    {

        ASDPlasticMaterial3D <ElasticityType,
                           YieldFunctionType,
                           PlasticFlowType,
                           thisClassTag> *newmaterial = new ASDPlasticMaterial3D<ElasticityType,
        YieldFunctionType,
        PlasticFlowType,
        thisClassTag>(ASDP_TAG);
        newmaterial->TrialStrain = this->TrialStrain;
        newmaterial->TrialStress = this->TrialStress;
        newmaterial->TrialPlastic_Strain = this->TrialPlastic_Strain;
        newmaterial->CommitStress = this->CommitStress;
        newmaterial->CommitStrain = this->CommitStrain;
        newmaterial->CommitPlastic_Strain = this->CommitPlastic_Strain;
        newmaterial->iv_storage = this->iv_storage;
        newmaterial->parameters_storage = this->parameters_storage;

        return newmaterial;
    }


    NDMaterial *getCopy(const char *type) {
        if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
            ASDPlasticMaterial3D <ElasticityType,
                               YieldFunctionType,
                               PlasticFlowType,
                               thisClassTag> *newmaterial = new ASDPlasticMaterial3D<ElasticityType,
            YieldFunctionType,
            PlasticFlowType,
            thisClassTag>(ASDP_TAG);
            newmaterial->TrialStrain = this->TrialStrain;
            newmaterial->TrialStress = this->TrialStress;
            newmaterial->TrialPlastic_Strain = this->TrialPlastic_Strain;
            newmaterial->CommitStress = this->CommitStress;
            newmaterial->CommitStrain = this->CommitStrain;
            newmaterial->CommitPlastic_Strain = this->CommitPlastic_Strain;
            newmaterial->iv_storage = this->iv_storage;
            newmaterial->parameters_storage = this->parameters_storage;

            return newmaterial;
        } else
        {
            cout << "ASDPlasticMaterial3D::getCopy(const char *type) - Only 3D is currently supported. Use 3D elements!" << endl;
        }
        return 0;
    }


    Response *setResponse (const char **argv, int argc,
                           OPS_Stream & s)
    {

        static Vector return_vector(6);

        if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
            return new MaterialResponse(this, 1, this->getStress());
        else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
            return new MaterialResponse(this, 2, this->getStrain());
        else if (strcmp(argv[0], "pstrain") == 0 || strcmp(argv[0], "pstrains") == 0)
            return new MaterialResponse(this, 3, this->getPstrain());
        else if (strcmp(argv[0], "eqpstrain") == 0 )
            return new MaterialResponse(this, 4, this->getEQPstrain());        
        else if (strcmp(argv[0], "PStress") == 0 )
            return new MaterialResponse(this, 5, this->getPStress());        
        else if (strcmp(argv[0], "J2Stress") == 0 )
            return new MaterialResponse(this, 6, this->getJ2Stress());        
        else if (strcmp(argv[0], "VolStrain") == 0 )
            return new MaterialResponse(this, 7, this->getVolStrain());        
        else if (strcmp(argv[0], "J2Strain") == 0 )
            return new MaterialResponse(this, 8, this->getJ2Strain());
        else
        {
            const char *iv_name = argv[0];

            int iv_size = this->getInternalVariableSizeByName(iv_name);
            int pos = this->getInternalVariableIndexByName(iv_name);

            return new MaterialResponse(this, 1000 + pos, Vector(iv_size));
        }

        return 0;
    }


    int getResponse (int responseID, Information & matInformation)
    {


        if (matInformation.theVector == 0)
            return 0;

        if (responseID == -1)
        {
            return -1;
        }
        else if (responseID == 1)
        {
            *(matInformation.theVector) = getStress();
        }
        else if (responseID == 2)
            *(matInformation.theVector) = getStrain();
        else if (responseID == 3)
            *(matInformation.theVector) = getPstrain();
        else if (responseID == 4)
            *(matInformation.theVector) = getEQPstrain();        
        else if (responseID == 5)
            *(matInformation.theVector) = getPStress();        
        else if (responseID == 6)
            *(matInformation.theVector) = getJ2Stress();        
        else if (responseID == 7)
            *(matInformation.theVector) = getVolStrain();        
        else if (responseID == 8)
            *(matInformation.theVector) = getJ2Strain();
        else if (responseID >= 1000)
        {
            int pos = responseID - 1000;
            *(matInformation.theVector) = getInternalVariableByPos(pos);
        }

        return 0;
    }

    const char *getType(void) const {return "ThreeDimensional";}

    int sendSelf(int commitTag, Channel & theChannel)
    {
        cerr << "ASDPlasticMaterial3D::sendSelf - not implemented!!!\n" ;


        return 0;
    }

    int recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
    {
        cerr << "ASDPlasticMaterial3D::recvSelf - not implemented!!!\n" ;


        return 0;
    }

    void Print(OPS_Stream & s, int flag = 0) {
        s << "ASDPlasticMaterial3D" << endln;
        s << "  Yield Function          : " << yf.NAME << endln;
        s << "  Plastic flow direction  : " << pf.NAME << endln;
        s << "  Elasticity              : " << et.NAME << endln;
        s << "  # of Internal variables : " << (int) iv_storage.size() <<  endln;
        iv_storage.print_components();
        s << "  # of Parameters         : " << (int) parameters_storage.size() <<  endln;
        parameters_storage.print_components();
    }

    void Print(ostream & s, int flag = 0)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        s << "ASDPlasticMaterial3D" << endl;
        s << "  Yield Function          : " << yf.NAME << endl;
        s << "  Plastic flow direction  : " << pf.NAME << endl;
        s << "  Elasticity              : " << et.NAME << endl;
        s << "  # of Internal variables : " << iv_storage.size() <<  endl;
        iv_storage.print_components();
        s << "  # of Parameters         : " << parameters_storage.size() <<  endl;
        parameters_storage.print_components();

    }

    int getObjectSize()
    {
        int size = 0;

        // 6 3x3 VoigtVectors and 1 VoigtMatrix (3x3x3x3)
        size += (3 * 3 * 6 + 3 * 3 * 3 * 3) * sizeof(double);

        //Four pointers
        size += 4 * sizeof(YieldFunctionType*);

        //Whatever the base components size is
        size += sizeof(yf);//yf->getObjectSize();
        size += sizeof(et);//et->getObjectSize();
        size += sizeof(pf);//pf->getObjectSize();
        // size += sizeof(internal_variables);//internal_variables->getObjectSize();
        size += sizeof(NDMaterial);

        // size += static_cast<T*>(this)->getObjectSize();

        return size;
    }

    void setTrialStress(const VoigtVector & stress)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        TrialStress = stress;
    }

    bool set_constitutive_integration_method(int method, int tangent, double f_absolute_tol, double stress_absolute_tol, int n_max_iterations, int return_to_yield_surface, int rk45_niter_max, double rk45_dT_min)
    {
        if ( method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Not_Set
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Crisfield
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Multistep_Forward_Euler
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Multistep_Forward_Euler_Crisfield
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Modified_Euler_Error_Control
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Subincrement
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_ddlambda
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_ddlambda_Subincrement
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Full_Backward_Euler)
        {
            INT_OPT_constitutive_integration_method[ASDP_TAG] = (ASDPlasticMaterial3D_Constitutive_Integration_Method) method ;
            INT_OPT_tangent_operator_type[ASDP_TAG] = (ASDPlasticMaterial3D_Tangent_Operator_Type) tangent ;
            DOUBLE_OPT_f_absolute_tol[ASDP_TAG] = f_absolute_tol ;
            DBL_OPT_stress_absolute_tol[ASDP_TAG] = stress_absolute_tol ;
            INT_OPT_n_max_iterations[ASDP_TAG] = n_max_iterations ;
            INT_OPT_return_to_yield_surface[ASDP_TAG] = return_to_yield_surface ;
            DBL_OPT_RK45_dT_min[ASDP_TAG] = rk45_dT_min ;
            INT_OPT_RK45_niter_max[ASDP_TAG] = rk45_niter_max ;

            GLOBAL_INT_max_iter[ASDP_TAG] = 0;
            GLOBAL_DBL_max_error[ASDP_TAG] = 0.;

            cout << "set_constitutive_integration_method tag = " << ASDP_TAG << " ::: " << endl;
            cout << "   method = " << method << endl;
            cout << "   tanget_type = " << tangent << endl;
            cout << "   f_absolute_tol = " << f_absolute_tol << endl;
            cout << "   stress_absolute_tol = " << stress_absolute_tol << endl;
            cout << "   n_max_iterations = " << n_max_iterations << endl;
            cout << "   return_to_yield_surface = " << return_to_yield_surface << endl;
            cout << "   rk45_niter_max = " << rk45_niter_max << endl;
            cout << "   rk45_dT_min = " << rk45_dT_min << endl;

            return true;
        }
        else
        {
            cerr << "ASDPlasticMaterial3D::set_constitutive_integration_method - Unknown constitutive_integration_method\n";
            return false;
        }
    }

protected:

    void setTrialPlastic_Strain(const VoigtVector & strain)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        TrialPlastic_Strain = strain;
    }

    void setCommitStress(const VoigtVector & stress)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        CommitStress = stress;
    }

    void setCommitStrain(const VoigtVector & strain)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        CommitStrain = strain;
    }

    void setCommitPlastic_Strain(const VoigtVector & strain)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        CommitPlastic_Strain = strain;
    }

    void setStiffness(const VoigtMatrix & stiff)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        Stiffness = stiff;
    }


    void setStressTensor(VoigtVector & stress)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        CommitStress = stress;
        TrialStress = stress;
        return;
    }


private:


    int Forward_Euler(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;



        int errorcode = -1;

        static VoigtVector depsilon;
        depsilon *= 0;
        depsilon = strain_incr;

        const VoigtVector& sigma = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma *= 0;
        intersection_stress *= 0;
        intersection_strain *= 0;

        VoigtMatrix Eelastic = et(sigma, parameters_storage);

        dsigma = Eelastic * depsilon;

        TrialStress = sigma + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(sigma, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
            return 0;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = DOUBLE_OPT_f_absolute_tol[ASDP_TAG];
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = epsilon  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;

            Eelastic = et(intersection_stress, parameters_storage);
            TrialStress  += Eelastic * depsilon_elpl;

            //Compute normal to YF (n) and Plastic Flow direction (m)
            const VoigtVector& n = yf.df_dsigma_ij(intersection_stress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, intersection_stress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  intersection_stress, iv_storage, parameters_storage);
            double den = n.transpose() * Eelastic * m - hardening;

            //Compute the plastic multiplier
            if (abs(den) < MACHINE_EPSILON)
            {
                cout << "CEP - den = 0\n";
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor1("m", m);
                printTensor1("n", n);
                cout << "hardening = " << hardening << endl;
                cout << "den = " << den << endl;
                printTensor1("depsilon_elpl", depsilon_elpl);
                return -1;
            }

            double dLambda =  n.transpose() * Eelastic * depsilon_elpl;
            dLambda /= den;

            if (dLambda <= 0)
            {
                // cout << "CEP - dLambda = " << dLambda << " <= 0\n";
                // printTensor1("m", m);
                // printTensor1("n", n);
                // cout << "hardening = " << hardening << endl;
                // cout << "den = " << den << endl;
                // printTensor1("depsilon_elpl", depsilon_elpl);
                dLambda = 0;
            }

            // Update the trial plastic strain.
            TrialPlastic_Strain += dLambda * m;

            // This code iterates internal variables and updates the trial values
            iv_storage.apply([&m, &dLambda, this](auto & internal_variable)
            {
                auto h = internal_variable.hardening_function(depsilon_elpl, m, intersection_stress, parameters_storage);
                internal_variable.trial_value += dLambda * h;
            });


            // Deal with the APEX if needed
            // if constexpr (yf_has_apex<YieldFunctionType>::value) {
            // {
            //    // TODO
            // }

            //Correct the trial stress
            TrialStress = TrialStress - dLambda * Eelastic * m;


            // Returning to Yield surface as recommended by Crisfield...
            if (INT_OPT_return_to_yield_surface[ASDP_TAG])
            {
                double yf_val_corr = yf(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& n_corr = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_corr = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

                double hardening_corr = yf.hardening( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);
                double dLambda_corr = yf_val_corr / (
                                                     n_corr.transpose() * Eelastic * m_corr - hardening_corr
                                                 );
                TrialStress = TrialStress - dLambda_corr * Eelastic * m_corr;
                TrialPlastic_Strain += dLambda_corr * m_corr;
            }
            else
            {
                // Do nothing
            }
            // ============================================================================================
            // ============================================================================================

            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                // cout << "Numeric error!\n";
                // printTensor1("TrialStress = " , TrialStress);
                // printTensor1("CommitStress = " , CommitStress);
                // printTensor1("depsilon = " , depsilon);
                // printTensor1("dsigma   = " , dsigma);
                // printTensor1("intersection_stress = " , intersection_stress);
                // printTensor2("Eelastic = " , Eelastic);
                // printTensor2("Stiffness = " , Stiffness);
                // cout << "yf_val_start = " << yf_val_start << endl;
                // cout << "yf_val_end = " << yf_val_end << endl;
                // printTensor1("n = " , n );
                // printTensor1("m = " , m );
                // cout << "hardening  = " << hardening << endl;
                // cout << "den = " << den << endl;
                // cout << "dLambda = " << dLambda << endl;

                return -1;
            }
            else
            {
                return 0;
            }

            ComputeTangentStiffness();

        }

        return errorcode;
    }



    // int Forward_Euler_Subincrement(const VoigtVector &strain_incr, bool const& with_return2yield_surface)
    int Forward_Euler_Subincrement(const VoigtVector & strain_incr)
    {
       using namespace ASDPlasticMaterial3DGlobals;



        int errorcode = -1;

        static VoigtVector depsilon;
        depsilon *= 0;
        depsilon = strain_incr;

        const VoigtVector& sigma = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma *= 0;
        intersection_stress *= 0;
        intersection_strain *= 0;

        VoigtMatrix Eelastic = et(sigma, parameters_storage);

        dsigma = Eelastic * depsilon;

        TrialStress = sigma + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(sigma, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
            return 0;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = DOUBLE_OPT_f_absolute_tol[ASDP_TAG];
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = epsilon  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;

            int Nsubsteps = INT_OPT_n_max_iterations[ASDP_TAG];
            depsilon_elpl = depsilon_elpl/Nsubsteps;

            for (int substep = 0; substep < Nsubsteps; ++substep)
            {
          
                Eelastic = et(intersection_stress, parameters_storage);
                TrialStress  += Eelastic * depsilon_elpl;

                //Compute normal to YF (n) and Plastic Flow direction (m)
                const VoigtVector& n = yf.df_dsigma_ij(intersection_stress, iv_storage, parameters_storage);
                const VoigtVector& m = pf(depsilon_elpl, intersection_stress, iv_storage, parameters_storage);

                double hardening = yf.hardening( depsilon_elpl, m,  intersection_stress, iv_storage, parameters_storage);
                double den = n.transpose() * Eelastic * m - hardening;

                //Compute the plastic multiplier
                if (abs(den) < MACHINE_EPSILON)
                {
                    cout << "CEP - den = 0\n";
                    cout << "yf_val_start = " << yf_val_start << endl;
                    cout << "yf_val_end = " << yf_val_end << endl;
                    printTensor1("m", m);
                    printTensor1("n", n);
                    cout << "hardening = " << hardening << endl;
                    cout << "den = " << den << endl;
                    printTensor1("depsilon_elpl", depsilon_elpl);
                    return -1;
                }

                double dLambda =  n.transpose() * Eelastic * depsilon_elpl;
                dLambda /= den;

                if (dLambda <= 0)
                {
                    dLambda = 0;
                }

                // Update the trial plastic strain.
                TrialPlastic_Strain += dLambda * m;

                // This code iterates internal variables and updates the trial values
                iv_storage.apply([&m, &dLambda, this](auto & internal_variable)
                {
                    auto h = internal_variable.hardening_function(depsilon_elpl, m, intersection_stress, parameters_storage);
                    internal_variable.trial_value += dLambda * h;
                });

                //Correct the trial stress
                TrialStress = TrialStress - dLambda * Eelastic * m;
            }

            // Returning to Yield surface as recommended by Crisfield...
            if (INT_OPT_return_to_yield_surface[ASDP_TAG])
            {
                double yf_val_corr = yf(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& n_corr = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_corr = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

                double hardening_corr = yf.hardening( depsilon_elpl, m_corr,  TrialStress, iv_storage, parameters_storage);
                double dLambda_corr = yf_val_corr / (
                                                     n_corr.transpose() * Eelastic * m_corr - hardening_corr
                                                 );
                TrialStress = TrialStress - dLambda_corr * Eelastic * m_corr;
                TrialPlastic_Strain += dLambda_corr * m_corr;
            }
            else
            {
                // Do nothing
            }
            // ============================================================================================
            // ============================================================================================

            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                // cout << "Numeric error!\n";
                // printTensor1("TrialStress = " , TrialStress);
                // printTensor1("CommitStress = " , CommitStress);
                // printTensor1("depsilon = " , depsilon);
                // printTensor1("dsigma   = " , dsigma);
                // printTensor1("intersection_stress = " , intersection_stress);
                // printTensor2("Eelastic = " , Eelastic);
                // printTensor2("Stiffness = " , Stiffness);
                // cout << "yf_val_start = " << yf_val_start << endl;
                // cout << "yf_val_end = " << yf_val_end << endl;
                // printTensor1("n = " , n );
                // printTensor1("m = " , m );13
                // cout << "hardening  = " << hardening << endl;
                // cout << "den = " << den << endl;
                // cout << "dLambda = " << dLambda << endl;

                return -1;
            }
            else
            {
                return 0;
            }

            ComputeTangentStiffness();

        }

        return errorcode;
    }

    std::pair<double, VoigtVector> CalculateLambdaM(
        const VoigtVector & thisSigma,
        const VoigtVector & depsilon_elpl,
        const parameters_storage_t &this_parameters_storage,
        const iv_storage_t &this_iv_storage)
    {
        VoigtMatrix Eelastic = et(thisSigma, this_parameters_storage);
        const VoigtVector& n = yf.df_dsigma_ij(thisSigma, this_iv_storage, this_parameters_storage);
        VoigtVector m = pf(depsilon_elpl, thisSigma, this_iv_storage, this_parameters_storage);
        double hardening = yf.hardening(depsilon_elpl, m, thisSigma, this_iv_storage, this_parameters_storage);

        double den = n.transpose() * Eelastic * m - hardening;
        double dLambda = (den != 0) ? (n.transpose() * Eelastic * depsilon_elpl).value() / den : 0;

        if (dLambda != dLambda)
        {
            cout << "CalculateLambdaM error" << endl;
            cout << "yf = " << yf(thisSigma, this_iv_storage, this_parameters_storage) << endl;
            cout << "thisSigma = " << thisSigma.transpose() << endl;
            cout << "depsilon_elpl = " << depsilon_elpl.transpose() << endl;
            cout << "n = " << n.transpose() << endl;
            cout << "m = " << m.transpose() << endl;
            cout << "hardening = " << hardening << endl;
            cout << "den = " << den << endl;
            cout << "Eelastic = " << Eelastic << endl;

            cout << "IVSTORAGE" << endln;
            this_iv_storage.print_components();

            cout << "PARAMSTORAGE" << endln;
            this_parameters_storage.print_components();
        }


        return std::make_pair(dLambda, m);
    }

    int Runge_Kutta_45_Error_Control(const VoigtVector & strain_incr)
    {
        // cout << "Runge_Kutta_45_Error_Control" << endl;

        int errorcode = -1;

        static VoigtVector depsilon;
        depsilon *= 0;
        depsilon = strain_incr;


        iv_storage.revert_all();


        dsigma *= 0;
        intersection_stress *= 0;
        intersection_strain *= 0;

        VoigtMatrix Eelastic = et(CommitStress, parameters_storage);

        dsigma = Eelastic * depsilon;


        TrialStress = CommitStress + dsigma;
        TrialStrain = CommitStrain + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        double yf_val_start = yf(CommitStress, iv_storage, parameters_storage);
        double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

        VoigtVector start_stress = CommitStress;
        VoigtVector end_stress = TrialStress;

        intersection_stress = start_stress;

        if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
        {
            Stiffness = Eelastic;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = DOUBLE_OPT_f_absolute_tol[ASDP_TAG];
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = CommitStrain  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;
            double T = 0.0, dT = 1.0, dT_min = this->DBL_OPT_RK45_dT_min[ASDP_TAG], TolE = this->DBL_OPT_stress_absolute_tol[ASDP_TAG];


            VoigtVector next_Sigma = TrialStress;
            iv_storage_t next_iv_storage = iv_storage;
            VoigtVector next_EpsilonPl = CommitPlastic_Strain;

            VoigtVector prev_Sigma = TrialStress;
            VoigtVector prev_EpsilonPl = CommitPlastic_Strain;
            iv_storage_t prev_iv_storage = iv_storage;

            VoigtVector dSigma, dSigma1, dSigma2, dSigma3, dSigma4, dSigma5, dSigma6;
            VoigtVector dEpsilonPl, dEpsilonPl1, dEpsilonPl2, dEpsilonPl3, dEpsilonPl4, dEpsilonPl5, dEpsilonPl6;
            iv_storage_t iv_storage1 = iv_storage;
            iv_storage_t iv_storage2 = iv_storage;
            iv_storage_t iv_storage3 = iv_storage;
            iv_storage_t iv_storage4 = iv_storage;
            iv_storage_t iv_storage5 = iv_storage;


            int niter = 0;
            double maxStepError = 0;
            while (T < 1.0)
            {
                niter ++;

                VoigtVector dEPS = dT * depsilon_elpl;
                VoigtVector m;
                double dLambda;

                //Delta 1
                std::tie(dLambda, m) = CalculateLambdaM(prev_Sigma, dEPS, parameters_storage, prev_iv_storage);
                dSigma1  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl1 = dLambda * m;
                iv_storage1.apply([&m, &dLambda, &prev_Sigma, &dEPS, this](auto & iv1)
                {
                    auto h = iv1.hardening_function(dEPS, m, prev_Sigma, parameters_storage);
                    iv1.trial_value = iv1.committed_value + 0.5 * dLambda * h;
                });
                next_Sigma =  prev_Sigma + 0.5 * dSigma1;

                //Delta 2
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage1);
                dSigma2  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl2 = dLambda * m;
                iv_storage2.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, this](auto & iv2)
                {
                    using VT1 = std::decay_t<decltype(iv2)>;
                    const VT1 &iv1 = iv_storage1.template get<VT1>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;

                    auto h2 = iv2.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH2 = dLambda * h2;

                    iv2.trial_value = iv2.committed_value + 0.25 * (DH1 + DH2);
                });
                next_Sigma =  prev_Sigma + 0.25 * (dSigma1 + dSigma2);

                //Delta 3
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage2);
                dSigma3  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl3 = dLambda * m;
                iv_storage3.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, this](auto & iv3)
                {
                    using VT = std::decay_t<decltype(iv3)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;

                    auto h3 = iv3.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH3 = dLambda * h3;

                    iv3.trial_value = iv3.committed_value +  -DH2 + 2 * DH3;
                });
                next_Sigma =  prev_Sigma - dSigma2 + 2 * dSigma3;

                //Delta 4
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage3);
                dSigma4  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl4 = dLambda * m;
                iv_storage4.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, &iv_storage3, this](auto & iv4)
                {
                    using VT = std::decay_t<decltype(iv4)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;
                    const VT &iv3 = iv_storage3.template get<VT>();
                    auto DH3 = iv3.trial_value - iv3.committed_value;

                    auto h4 = iv4.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH4 = dLambda * h4;

                    iv4.trial_value = iv4.committed_value +  (7 * DH1 + 10 * DH2 + DH4) / 27;
                });
                next_Sigma =  prev_Sigma + (7 * dSigma1 + 10 * dSigma2 + dSigma4) / 27;

                //Delta 5
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage4);
                dSigma5  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl5 = dLambda * m;
                iv_storage5.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, &iv_storage3, &iv_storage4, this](auto & iv5)
                {
                    using VT = std::decay_t<decltype(iv5)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;
                    const VT &iv3 = iv_storage3.template get<VT>();
                    auto DH3 = iv3.trial_value - iv3.committed_value;
                    const VT &iv4 = iv_storage4.template get<VT>();
                    auto DH4 = iv4.trial_value - iv4.committed_value;

                    auto h5 = iv5.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH5 = dLambda * h5;

                    iv5.trial_value = iv5.committed_value +  (28 * DH1 - 125 * DH2 + 546 * DH3 + 54 * DH4 - 378 * DH5) / 625;
                });
                next_Sigma =  prev_Sigma + (28 * dSigma1 - 125 * dSigma2 + 546 * dSigma3 + 54 * dSigma4 - 378 * dSigma5) / 625;

                //Delta 6 - The final predictor
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage5);
                dSigma6  = Eelastic * (dEPS - dLambda * m);
                dEpsilonPl6 = dLambda * m;
                next_iv_storage.apply([&m, &dLambda, &next_Sigma, &dEPS, &iv_storage1, &iv_storage2, &iv_storage3, &iv_storage4, &iv_storage5, this](auto & niv)
                {
                    using VT = std::decay_t<decltype(niv)>;
                    const VT &iv1 = iv_storage1.template get<VT>();
                    auto DH1 = iv1.trial_value - iv1.committed_value;
                    const VT &iv2 = iv_storage2.template get<VT>();
                    auto DH2 = iv2.trial_value - iv2.committed_value;
                    const VT &iv3 = iv_storage3.template get<VT>();
                    auto DH3 = iv3.trial_value - iv3.committed_value;
                    const VT &iv4 = iv_storage4.template get<VT>();
                    auto DH4 = iv4.trial_value - iv4.committed_value;
                    const VT &iv5 = iv_storage5.template get<VT>();
                    auto DH5 = iv5.trial_value - iv5.committed_value;

                    auto h6 = niv.hardening_function(dEPS, m, next_Sigma, parameters_storage);
                    auto DH6 = dLambda * h6;

                    niv.trial_value = niv.committed_value +  ( DH1 +  4 * DH3 + DH4) / 6;
                });
                dSigma =  (dSigma1 + 4 * dSigma3 + dSigma4) / 6;
                dEpsilonPl =  (dEpsilonPl1 + 4 * dEpsilonPl3 + dEpsilonPl4) / 6;
                next_Sigma =  prev_Sigma + dSigma;
                next_EpsilonPl = prev_EpsilonPl + dEpsilonPl;


                if (dSigma != dSigma)
                {
                    cout << "ASDPlasticMaterial3D::Runge_Kutta_45_Error_Control Integration error" << endl;
                    cout << "T = " << T << endl;
                    cout << "dT = " << dT << endl;
                    cout << "dSigma1 = " << dSigma1.transpose() << endl;
                    cout << "dSigma2 = " << dSigma2.transpose() << endl;
                    cout << "dSigma3 = " << dSigma3.transpose() << endl;
                    cout << "dSigma4 = " << dSigma4.transpose() << endl;
                    cout << "dSigma5 = " << dSigma5.transpose() << endl;
                    cout << "dSigma6 = " << dSigma6.transpose() << endl;
                    cout << "m = " << m.transpose() << endl;
                    cout << "dEpsilonPl = " << dEpsilonPl.transpose() << endl;
                    cout << "TrialStress = " << TrialStress.transpose() << endl;
                    cout << "dEPS = " << dEPS.transpose() << endl;
                    exit(-1);
                }




                //Stress norm and stress error
                double stressNorm = next_Sigma.norm();
                double curStepError1 = (-42 * dSigma1 - 224 * dSigma3 - 21 * dSigma4 + 162 * dSigma5 + 125 * dSigma6).norm() / 336;
                if (stressNorm >= 0.5) { curStepError1 /= (2 * stressNorm); }
                //Internal variables norm and internal variables error (TODO)

                // double curStepError = curStepError1; //fmax(curStepError1, curStepError2);
                // double curStepError = curStepError1/stressNorm; //fmax(curStepError1, curStepError2);
                double curStepError = curStepError1; //fmax(curStepError1, curStepError2);

                //Check convergence and adjust integration timestep
                if (curStepError > TolE)
                {
                    double q = fmax(0.8 * pow(TolE / curStepError, 0.2), 0.1);

                    if (dT == dT_min) {

                        prev_Sigma = next_Sigma;
                        prev_EpsilonPl = next_EpsilonPl;
                        prev_iv_storage = next_iv_storage;

                        T += dT;
                    }
                    dT = fmax(q * dT, dT_min);
                }
                else {

                    prev_Sigma = next_Sigma;
                    prev_EpsilonPl = next_EpsilonPl;
                    prev_iv_storage = next_iv_storage;

                    double q = fmin(0.8 * pow(TolE / curStepError, 0.2), 2.0);
                    T += dT;
                    dT = fmax(q * dT, dT_min);
                    dT = fmin(dT, 1 - T);

                    maxStepError = max(maxStepError, curStepError);
                }

                if (niter > this->INT_OPT_RK45_niter_max[ASDP_TAG])
                {
                    cout << "ASDPlasticMaterial3D - tag = " << ASDP_TAG << " exceeded number of iterations. niter = " << niter << " niter_max = " <<this->INT_OPT_RK45_niter_max[ASDP_TAG] << " T= " << T << " dT = " << dT << endl;
                    // throw std::runtime_error("ASDPLasticMaterial3D - Unable to find a valid bracket in compute_yf_crossing");
                    return -1;
                }

            }

            GLOBAL_INT_max_iter[ASDP_TAG] = std::max(GLOBAL_INT_max_iter[ASDP_TAG], niter);
            GLOBAL_DBL_max_error[ASDP_TAG] = std::max(GLOBAL_DBL_max_error[ASDP_TAG], maxStepError);

            TrialStress = next_Sigma;
            TrialPlastic_Strain = next_EpsilonPl;
            iv_storage = next_iv_storage;

           //Return to Yield
            if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 1)  // Return to yield in one step
            {
                // In the evolve function, only dLambda and m are used. Other arguments are not used at all.
                // Make surface the internal variables are already updated. And then, return to the yield surface.
                double yf_val_after_corrector;
                int iter =0;
                // double TOL = 10*this->DBL_OPT_stress_absolute_tol[ASDP_TAG];
                // double NITER = this->INT_OPT_n_max_iterations[ASDP_TAG];
                // do
                {
                    yf_val_after_corrector = yf(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    // In the function below, depsilon_elpl is actually not used at all in hardening
                    // double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double dLambda_after_corrector = yf_val_after_corrector / (
                                                         n_after_corrector.transpose() * Eelastic * m_after_corrector - hardening_after_corrector
                                                     );
                    TrialStress = TrialStress - dLambda_after_corrector * Eelastic * m_after_corrector;
                    TrialPlastic_Strain += dLambda_after_corrector * m_after_corrector;

                    // iter++;
                }
                // while(yf_val_after_corrector > TOL && iter < NITER);
            }


           //Return to Yield with bisection
           else if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 2)  // Return to yield with iterations
           {
                // In the evolve function, only dLambda and m are used. Other arguments are not used at all.
                // Make surface the internal variables are already updated. And then, return to the yield surface.
                double y0  = yf(TrialStress, iv_storage, parameters_storage) ;
                int iter = 0;
                double TOL = this->DOUBLE_OPT_f_absolute_tol[ASDP_TAG];
                double NITER = this->INT_OPT_n_max_iterations[ASDP_TAG];
                // do
                if(y0 > 0 && iter < NITER)
                {
                    y0 = yf(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    // In the function below, depsilon_elpl is actually not used at all in hardening
                    // double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double dL = y0 / (
                                                         n_after_corrector.transpose() * Eelastic * m_after_corrector - hardening_after_corrector
                                                     );

                    
                    VoigtVector TS = TrialStress - dL * Eelastic * m_after_corrector;

                    double y1 = yf(TS, iv_storage, parameters_storage);

                    //Try to bracket solution
                    // cout << "   y0 = " << y0 << "  dL =" << dL << "   y1 =" << y1 << endl;
                    while( y1 > 0 && iter < NITER)
                    {
                        dL = dL * 1.1;
                        TS = TrialStress - dL * Eelastic * m_after_corrector;
                        y1 = yf(TS, iv_storage, parameters_storage);
                        iter ++;

                        // cout << "   iter = " << iter << "  dL =" << dL << "   y1 =" << y1 << endl;
                    }

                    iter = 0;

                    //Once solution is bracketed, use bisection to get to YS
                    if (y1 < 0)
                    {
                        double dL_min = 0;
                        double dL_max = dL;
                        double dL_mid = dL / 2;

                        VoigtVector TS2 = TrialStress - dL_mid * Eelastic * m_after_corrector;
                        double y_mid = yf(TS2, iv_storage, parameters_storage);

                        // cout << "   y_mid = " << y_mid << "  dL_min =" << dL_min << "   dL_mid =" << dL_mid << "   dL_max =" << dL_max << endl;
                        
                        while(abs(y_mid) > TOL && iter < NITER)
                        {
                            // cout << "   iter = " << iter << "   y_mid = " << y_mid << "  dL_min =" << dL_min << "   dL_mid =" << dL_mid << "   dL_max =" << dL_max << endl;
                            if (y_mid  > 0)
                            {
                                dL_min = dL_mid;
                            } else
                            {
                                dL_max = dL_mid;
                            }
                            dL_mid = 0.5*(dL_min + dL_max);
                            TS2 = TrialStress - dL_mid * Eelastic * m_after_corrector;
                            y_mid = yf(TS2, iv_storage, parameters_storage);
                            iter++;
                        }
                        dL = dL_mid;
                    }
                    

                    TrialStress = TrialStress - dL * Eelastic * m_after_corrector;
                    TrialPlastic_Strain += dL * m_after_corrector;

                    // iter++;
                }
            }

            // ============================================================================================
            // ============================================================================================

            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                cout << "ASDPlasticMaterial3D::Runge_Kutta_45_Error_Control  Numeric error!\n";
                printTensor1("TrialStress = " , TrialStress);
                printTensor1("CommitStress = " , CommitStress);
                printTensor1("depsilon = " , depsilon);
                printTensor1("dsigma   = " , dsigma);
                printTensor1("intersection_stress = " , intersection_stress);
                printTensor2("Eelastic = " , Eelastic);
                printTensor2("Stiffness = " , Stiffness);
                cout << "yf_val_start = " << yf_val_start << endl;
                cout << "yf_val_end = " << yf_val_end << endl;
                printTensor1("n = " , yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage) );
                printTensor1("m = " , pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage) );
                cout << "hardening  = " << yf.hardening( depsilon_elpl, pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage),  TrialStress, iv_storage, parameters_storage) << endl;

                errorcode = -1;

                exit(-1);
            }

            ComputeTangentStiffness();
        }

        return 0;
    }


    //Robust Brent algorithm
    double compute_yf_crossing(const VoigtVector & start_stress, const VoigtVector & end_stress, double x1, double x2, double tol) const
    {
        using namespace ASDPlasticMaterial3DGlobals;

        // Constants
        const double MAX_BRACKET_EXPANSION = 100.0; // Maximum factor to expand the bracket
        const double BRACKET_EXPANSION_FACTOR = 1.6; // Factor to expand the bracket each iteration
        const double EPS_MULTIPLE = 2.0;
        const double TOL_MULTIPLE = 0.5;
        const double EPS = std::numeric_limits<double>::epsilon();

        // Lambda functions
        auto calculateSigma = [](const VoigtVector & start, const VoigtVector & end, double multiplier) {
            return start * (1 - multiplier) + end * multiplier;
        };

        auto calculateYf = [this](const VoigtVector & sigma) {
            return yf(sigma, iv_storage, parameters_storage);
        };

        // Bracket Adjustment: Expand the bracket if necessary
        double a = x1;
        double b = x2;
        double c = x2;
        double d = 0;
        double e = 0.0;
        double fa = calculateYf(calculateSigma(start_stress, end_stress, a));
        double fb = calculateYf(calculateSigma(start_stress, end_stress, b));
        double fc = fb;


        for (double factor = 1.0; factor <= MAX_BRACKET_EXPANSION; factor *= BRACKET_EXPANSION_FACTOR) {
            if ((fb * fa) <= 0.0) {
                break; // Valid bracket found
            }
            if (std::abs(fa) < std::abs(fb)) {
                a -= (b - a) * factor;
                fa = calculateYf(calculateSigma(start_stress, end_stress, a));
            } else {
                b += (b - a) * factor;
                fb = calculateYf(calculateSigma(start_stress, end_stress, b));
            }
        }

        // Brent's Method Implementation
        if ((fb * fa) > 0.0) {
            throw std::runtime_error("ASDPLasticMaterial3D - Unable to find a valid bracket in compute_yf_crossing");
        }

        for (int iter = 1; iter <= ASDPlasticMaterial3D_MAXITER_BRENT; iter++) {
            if ((fb * fc) > 0.0) {
                c = a;   // Rename a, b, c and adjust bounding interval d
                fc = fa;
                e = d = b - a;
            }

            if (std::abs(fc) < std::abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            double tol1 = EPS_MULTIPLE * EPS * std::abs(b) + TOL_MULTIPLE * tol;
            double xm = 0.5 * (c - b);

            if (std::abs(xm) <= tol1 || fb == 0.0) {
                return b;  // Convergence
            }

            if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
                // Attempt inverse quadratic interpolation
                double s = fb / fa;
                double p, q;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    double r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (p > 0.0) q = -q;
                p = std::abs(p);
                double min1 = 3.0 * xm * q - std::abs(tol1 * q);
                double min2 = std::abs(e * q);

                if (2.0 * p < std::min(min1, min2)) {
                    e = d;
                    d = p / q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;  // Bounds decreasing too slowly, use bisection
                e = d;
            }

            a = b;
            fa = fb;

            if (std::abs(d) > tol1) {
                b += d;
            } else {
                b += (xm > 0.0 ? std::abs(tol1) : -std::abs(tol1));
            }

            fb = calculateYf(calculateSigma(start_stress, end_stress, b));

            // Check for NaN (Not a Number) values
            if (std::isnan(fb)) {
                throw std::runtime_error("compute_yf_crossing: NaN encountered in function evaluation.");
            }
        }

        throw std::runtime_error("Maximum iterations reached without convergence in compute_yf_crossing");
    }


protected:

    VoigtVector TrialStrain;
    VoigtVector TrialStress;
    VoigtVector TrialPlastic_Strain;

    VoigtVector CommitStress;
    VoigtVector CommitStrain;
    VoigtVector CommitPlastic_Strain;


    YieldFunctionType yf;
    ElasticityType    et;
    PlasticFlowType   pf;

    iv_storage_t iv_storage;
    parameters_storage_t parameters_storage;


protected:

    static std::map<int, ASDPlasticMaterial3D_Constitutive_Integration_Method> INT_OPT_constitutive_integration_method;     //
    static std::map<int, ASDPlasticMaterial3D_Tangent_Operator_Type> INT_OPT_tangent_operator_type;     //
    static std::map<int, double> DOUBLE_OPT_f_absolute_tol;
    static std::map<int, double> DBL_OPT_stress_absolute_tol;
    static std::map<int, int> INT_OPT_n_max_iterations;
    static std::map<int, int> INT_OPT_return_to_yield_surface;
    static std::map<int, double> DBL_OPT_RK45_dT_min;
    static std::map<int, int> INT_OPT_RK45_niter_max;

    static std::map<int, int> GLOBAL_INT_max_iter; 
    static std::map<int, double> GLOBAL_DBL_max_error; 

    bool first_step;

    static VoigtVector dsigma;
    static VoigtVector depsilon_elpl;    //Elastoplastic strain increment : For a strain increment that causes first yield, the step is divided into an elastic one (until yield) and an elastoplastic one.
    static VoigtVector intersection_stress;
    static VoigtVector intersection_strain;
    static VoigtMatrix Stiffness;


};

template < class E, class Y, class P, int tag>
std::map<int, ASDPlasticMaterial3D_Constitutive_Integration_Method> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_constitutive_integration_method;
template < class E, class Y, class P, int tag>
std::map<int, ASDPlasticMaterial3D_Tangent_Operator_Type> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_tangent_operator_type;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::DOUBLE_OPT_f_absolute_tol;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::DBL_OPT_stress_absolute_tol;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_n_max_iterations;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_return_to_yield_surface;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::DBL_OPT_RK45_dT_min;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_RK45_niter_max;

template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::GLOBAL_DBL_max_error;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::GLOBAL_INT_max_iter; 

template < class E, class Y, class P, int tag>
VoigtVector ASDPlasticMaterial3D< E,  Y,  P,  tag>::dsigma;

template < class E, class Y, class P, int tag>
VoigtVector ASDPlasticMaterial3D< E,  Y,  P,  tag>::depsilon_elpl;  //Used to compute the yield surface intersection.

template < class E, class Y, class P, int tag>
VoigtVector ASDPlasticMaterial3D< E,  Y,  P,  tag >::intersection_stress;  //Used to compute the yield surface intersection.

template < class E, class Y, class P, int tag>
VoigtVector ASDPlasticMaterial3D< E,  Y,  P,  tag>::intersection_strain;  //Used to compute the yield surface intersection.

template < class E, class Y, class P, int tag>
VoigtMatrix ASDPlasticMaterial3D< E,  Y,  P,  tag>::Stiffness;  //Used to compute the yield surface intersection.


#endif
