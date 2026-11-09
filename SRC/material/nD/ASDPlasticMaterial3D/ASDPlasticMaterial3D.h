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
#include <cstdio>
#include <string>
#include <Channel.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include "ASDPlasticMaterial3DGlobals.h"


// #include "MaterialInternalVariables.h"
#include "YieldFunctions/AllYieldFunctions.h"
#include "PlasticFlowDirections/AllPlasticFlowDirections.h"
#include "ElasticityModels/AllElasticityModels.h"
#include "AllASDModelParameterTypes.h"
#include "AllASDInternalVariableTypes.h"

#include "utuple_storage.h"

// for print p, q, theta
#include <tuple>
#include <utility> // For std::pair
#include <map> // For std::pair
#include <limits>
#include <type_traits>

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
        stress_set_externally = false;
        Stiffness.setZero();
    }


    ASDPlasticMaterial3D(int tag)
        : NDMaterial(tag, thisClassTag)
    {


        TrialStress.setZero();
        CommitStress.setZero();
        TrialStrain.setZero();
        CommitStrain.setZero();
        TrialPlastic_Strain.setZero();
        CommitPlastic_Strain.setZero();
        Stiffness.setZero();

        first_step = true;
        stress_set_externally = false;
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
        return "ASDPlasticMaterial3D";
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

        if (first_step && !stress_set_externally)
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
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler :
            exitflag = this->Backward_Euler(strain_increment);;
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_LineSearch :
            exitflag = this->Backward_Euler_LineSearch(strain_increment);;
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Full_Backward_Euler :
            exitflag = this->Full_Backward_Euler(strain_increment);
            break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_ddlambda_Subincrement :
        //     exitflag = this->Backward_Euler_ddlambda_Subincrement(strain_increment);;
        //     break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Crisfield :
        //     exitflag = this->Forward_Euler(strain_increment, true);
        //     break;
        // case ASDPlasticMaterial3D_Constitutive_Integration_Method::Multistep_Forward_Euler_Crisfield :
        //     exitflag = this->Multistep_Forward_Euler(strain_increment, true);
        //     break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Modified_Euler_Error_Control :
            exitflag = this->Modified_Euler_Error_Control(strain_increment);
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control :
            exitflag = this->Runge_Kutta_45_Error_Control(strain_increment);;
            break;
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control_old :
            exitflag = this->Runge_Kutta_45_Error_Control_old(strain_increment);;
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
        // plastic strain holds engineering shear: use its tensor components
        const VoigtVector epsp = tensor_like(TrialPlastic_Strain);
        result(0) = std::sqrt(2./3.*tensor_dot_stress_like(epsp, epsp));
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
        result(0) = tensor_like(TrialStrain).getJ2();  // strain holds engineering shear
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
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_FirstOrder
                 || INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic)
        {
            // Algorithmic is computed exactly by Full_Backward_Euler; other integrators difference it
            compute_numerical_tangent_firstorder(TrialStrain-CommitStrain, Stiffness);
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_SecondOrder)
        {
            compute_numerical_tangent_secondorder(TrialStrain-CommitStrain, Stiffness);
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Continuum)
        {

            VoigtMatrix Eelastic = et(TrialStress, parameters_storage);
            const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);

            double den_after_corrector = strain_like(n).dot(Eelastic * strain_like(m)) - hardening;

            VoigtMatrix Econtinuum = Eelastic - (Eelastic * strain_like(m)) * (strain_like(n).transpose() * Eelastic) / den_after_corrector;
            Stiffness = Econtinuum;
        }
        else if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Secant)
        {

            VoigtMatrix Eelastic = et(TrialStress, parameters_storage);
            const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
            const VoigtVector& m = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

            double hardening = yf.hardening( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);

            double den_after_corrector = strain_like(n).dot(Eelastic * strain_like(m)) - hardening;

            VoigtMatrix Econtinuum = Eelastic - (Eelastic * strain_like(m)) * (strain_like(n).transpose() * Eelastic) / den_after_corrector;
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
                double tol_yf = DBL_OPT_f_absolute_tol[ASDP_TAG];
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
            double den = strain_like(n).dot(Eelastic * strain_like(m)) - hardening;

            //Compute the plastic multiplier
            // n:E:m - H = 0 where the flow direction vanishes, e.g. at the Mohr-Coulomb
            // apex with psi = 0 (m = 0): no plastic flow can be computed and the stress
            // stays on the yield surface. The numerical tangents difference this
            // function and ignore its return value, so stress_incr is always set.
            if (abs(den) < MACHINE_EPSILON)
            {
                stress_incr = intersection_stress - local_stress;
                return -1;
            }

            double dLambda =  strain_like(n).dot(Eelastic * depsilon_elpl);
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

            trial_stress = trial_stress - dLambda * Eelastic * strain_like(m);
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

        // no strain set yet: the elastic tangent
        if (Stiffness.isZero(0.0))
            Stiffness = et(CommitStress, parameters_storage);

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

        // cerr << "ASDPlasticMaterial3D::revertToLastCommit !!!\n" ;


        // TrialStress = CommitStress;
        // TrialStrain = CommitStrain;
        // TrialPlastic_Strain = CommitPlastic_Strain;

        // iv_storage.revert_all();

        // if (GLOBAL_INT_max_iter[ASDP_TAG] > 0 || GLOBAL_DBL_max_error[ASDP_TAG] > 0.)
        // {
        //     // cout << "  () ASDP Integration Info. Tag = " << ASDP_TAG << " max_iter = " << GLOBAL_INT_max_iter[ASDP_TAG] << " max_error = " << GLOBAL_DBL_max_error[ASDP_TAG] << endl;
        //     GLOBAL_INT_max_iter[ASDP_TAG] = 0;
        //     GLOBAL_DBL_max_error[ASDP_TAG] = 0.;
        // }

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
        newmaterial->stress_set_externally = this->stress_set_externally;
        newmaterial->Stiffness = this->Stiffness;

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
            newmaterial->stress_set_externally = this->stress_set_externally;
            newmaterial->Stiffness = this->Stiffness;

            return newmaterial;
        } else
        {
            cout << "ASDPlasticMaterial3D::getCopy(const char *type) - Only 3D is currently supported. Use 3D elements!" << endl;
        }
        return 0;
    }


    int setParameter(const char **argv, int argc, Parameter &param)
    {


        
        // if (theMaterialTag == this->getTag()) {
        if (true) {
            
            // State variables (use specific response IDs)
            if (strcmp(argv[0], "stress") == 0) {
                return param.addObject(1, this);
            }
            else if (strcmp(argv[0], "strain") == 0) {
                return param.addObject(2, this);
            }
            else if (strcmp(argv[0], "plasticStrain") == 0) {
                return param.addObject(3, this);
            }
            else if (strcmp(argv[0], "trialStress") == 0) {
                return param.addObject(4, this);
            }
            else if (strcmp(argv[0], "trialStrain") == 0) {
                return param.addObject(5, this);
            }
            else if (strcmp(argv[0], "trialPlasticStrain") == 0) {
                return param.addObject(6, this);
            }
            else if (strcmp(argv[0], "K02D") == 0) {
                cout << "       ---->  K02D" << endl;
                return param.addObject(7, this);
            }
            else if (strcmp(argv[0], "K03D") == 0) {
                cout << "       ---->  K03D" << endl;
                return param.addObject(8, this);
            }
            else if (strcmp(argv[0], "trialStressIncrement") == 0) {
                return param.addObject(9, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementXX") == 0) {
                return param.addObject(10, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementYY") == 0) {
                return param.addObject(11, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementZZ") == 0) {
                return param.addObject(12, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementXY") == 0) {
                return param.addObject(13, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementYZ") == 0) {
                return param.addObject(14, this);
            }
            else if (strcmp(argv[0], "trialStressIncrementXZ") == 0) {
                return param.addObject(15, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementXX") == 0) {
                return param.addObject(16, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementYY") == 0) {
                return param.addObject(17, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementZZ") == 0) {
                return param.addObject(18, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementXY") == 0) {
                return param.addObject(19, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementYZ") == 0) {
                return param.addObject(20, this);
            }
            else if (strcmp(argv[0], "commitStressIncrementXZ") == 0) {
                return param.addObject(21, this);
            }
            else {
                // For all other parameter names, use the parameter system to pass the name
                // Store the parameter name in the Parameter object (if supported)
                // or use a generic response ID
                current_parameter_name = argv[0]; // Store for use in updateParameter
                return param.addObject(1000, this);
            }
        }
        
        return -1;
    }

    int updateParameter(int responseID, Information &info)
    {



        // State variables (committed values)
        if (responseID == 1) { // stress
            if (info.theType == VectorType) {
                const Vector& newStress = *(info.theVector);
                CommitStress = VoigtVector::fromStress(newStress);
                TrialStress = CommitStress;
                stress_set_externally = true;
            }
            return 0;
        }
        else if (responseID == 2) { // strain
            if (info.theType == VectorType) {
                const Vector& newStrain = *(info.theVector);
                CommitStrain = VoigtVector::fromStrain(newStrain);
                TrialStrain = CommitStrain;
            }
            return 0;
        }
        else if (responseID == 3) { // plasticStrain
            if (info.theType == VectorType) {
                const Vector& newPlasticStrain = *(info.theVector);
                CommitPlastic_Strain = VoigtVector::fromStrain(newPlasticStrain);
                TrialPlastic_Strain = CommitPlastic_Strain;
            }
            return 0;
        }
        // Trial state variables
        else if (responseID == 4) { // trialStress
            if (info.theType == VectorType) {
                const Vector& newTrialStress = *(info.theVector);
                TrialStress = VoigtVector::fromStress(newTrialStress);
                stress_set_externally = true;
            }
            return 0;
        }
        else if (responseID == 5) { // trialStrain
            if (info.theType == VectorType) {
                const Vector& newTrialStrain = *(info.theVector);
                TrialStrain = VoigtVector::fromStrain(newTrialStrain);
            }
            return 0;
        }
        else if (responseID == 6) { // trialPlasticStrain
            if (info.theType == VectorType) {
                const Vector& newTrialPlasticStrain = *(info.theVector);
                TrialPlastic_Strain = VoigtVector::fromStrain(newTrialPlasticStrain);
            }
            return 0;
        }
        else if (responseID == 7) { // K02D - Use Sigma_Y
            // cout << "responseID == 7 !! info.theType = " << info.theType << " DoubleType = " << DoubleType << endl;
            // if (info.theType == DoubleType) {
                const double& K02D = info.theDouble;
                CommitStress(0) = K02D * CommitStress(1);
                CommitStress(2) = K02D * CommitStress(1);
                stress_set_externally = true;
            // }
            return 0;
        }
        else if (responseID == 8) { // K03D - Use Sigma_Z
            // if (info.theType == DoubleType) {
                const double& K03D = info.theDouble;
                CommitStress(0) = K03D * CommitStress(2);
                CommitStress(1) = K03D * CommitStress(2);
                stress_set_externally = true;
            // }
            return 0;
        }
        else if (responseID == 9) { // trialStressIncrement
            if (info.theType == VectorType) {
                const Vector& newTrialStress = *(info.theVector);
                TrialStress += VoigtVector::fromStress(newTrialStress);
                stress_set_externally = true;
            }
            return 0;
        }
        else if (responseID == 10) { // trialStressIncrementXX
            TrialStress(0) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 11) { // trialStressIncrementYY
            TrialStress(1) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 12) { // trialStressIncrementZZ
            TrialStress(2) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 13) { // trialStressIncrementXY
            TrialStress(3) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 14) { // trialStressIncrementYZ
            TrialStress(4) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 15) { // trialStressIncrementXZ
            TrialStress(5) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 16) { // commitStressIncrementXX
            CommitStress(0) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 17) { // commitStressIncrementYY
            CommitStress(1) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 18) { // commitStressIncrementZZ
            CommitStress(2) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 19) { // commitStressIncrementXY
            CommitStress(3) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 20) { // commitStressIncrementYZ
            CommitStress(4) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        else if (responseID == 21) { // commitStressIncrementXZ
            CommitStress(5) += info.theDouble;
            stress_set_externally = true;
            return 0;
        }
        // Generic parameter update (model parameters and internal variables)
        else if (responseID == 1000) {
            // Use the stored parameter name from the most recent setParameter call
            const char* param_name = current_parameter_name.c_str();
            double param_value = info.theDouble;
            
            // Try to set as model parameter first
            // utuple_storage::setParameterByName handles non-existent parameters gracefully (does nothing)
            parameters_storage.setParameterByName(param_name, param_value);
            
            // Also try to set as internal variable (for scalar internal variables)
            // utuple_storage::setInternalVariableByName also handles non-existent variables gracefully
            iv_storage.setInternalVariableByName(param_name, 1, &param_value);
            
            return 0;
        }
        
        return -1;
    }



    Response *setResponse (const char **argv, int argc,
                           OPS_Stream & s)
    {

        static Vector return_vector(6);

        // the output tags give MPCO (and the xml recorders) the component names,
        // as in NDMaterial::setResponse; without them MPCO records no components
        auto tagged = [&](int id, const Vector& v, const char* const* names) {
            s.tag("NdMaterialOutput");
            s.attr("matType", this->getClassType());
            s.attr("matTag", this->getTag());
            for (int i = 0; i < 6; ++i)
                s.tag("ResponseType", names[i]);
            Response* resp = new MaterialResponse(this, id, v);
            s.endTag();
            return resp;
        };
        static const char* const stress_names[6]  = {"sigma11", "sigma22", "sigma33", "sigma12", "sigma23", "sigma13"};
        static const char* const strain_names[6]  = {"eps11", "eps22", "eps33", "eps12", "eps23", "eps13"};
        static const char* const pstrain_names[6] = {"epsP11", "epsP22", "epsP33", "epsP12", "epsP23", "epsP13"};

        if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
            return tagged(1, this->getStress(), stress_names);
        else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
            return tagged(2, this->getStrain(), strain_names);
        else if (strcmp(argv[0], "pstrain") == 0 || strcmp(argv[0], "pstrains") == 0
                 || strcmp(argv[0], "plasticStrain") == 0 || strcmp(argv[0], "plasticStrains") == 0)
            return tagged(3, this->getPstrain(), pstrain_names);
        else if (strcmp(argv[0], "eqpstrain") == 0 )
            return new MaterialResponse(this, 4, this->getEQPstrain());        
        else if (strcmp(argv[0], "equivalentPlasticStrain") == 0 || strcmp(argv[0], "EquivalentPlasticStrain") == 0)
        {
            // same as eqpstrain, under the name other ASD materials (and STKO's
            // MPCO key material.equivalentPlasticStrain) use; the tags give MPCO
            // its component name
            s.tag("NdMaterialOutput");
            s.attr("matType", this->getClassType());
            s.attr("matTag", this->getTag());
            s.tag("ResponseType", "PLE");
            Response* resp = new MaterialResponse(this, 4, this->getEQPstrain());
            s.endTag();
            return resp;
        }
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

            // unknown name: no response (Vector(-1) used to crash the recorder)
            if (pos < 0 || iv_size <= 0)
                return 0;

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
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control_old
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Subincrement
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_LineSearch
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Full_Backward_Euler)
        {
            INT_OPT_constitutive_integration_method[ASDP_TAG] = (ASDPlasticMaterial3D_Constitutive_Integration_Method) method ;
            INT_OPT_tangent_operator_type[ASDP_TAG] = (ASDPlasticMaterial3D_Tangent_Operator_Type) tangent ;
            DBL_OPT_f_absolute_tol[ASDP_TAG] = f_absolute_tol ;
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

    void set_derivative_options(int derivative_type, double fd_perturbation)
    {
        INT_OPT_derivative_type[ASDP_TAG] = (ASDPlasticMaterial3D_Derivative_Type) derivative_type;
        DBL_OPT_fd_perturbation[ASDP_TAG] = fd_perturbation > 0 ? fd_perturbation : 1e-6;
        FBE_fallbacks[ASDP_TAG] = 0;
        FBE_failures[ASDP_TAG] = 0;
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
        depsilon.setZero();
        depsilon = strain_incr;

        const VoigtVector& sigma = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();

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
                double tol_yf = DBL_OPT_f_absolute_tol[ASDP_TAG];
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
            double den = strain_like(n).dot(Eelastic * strain_like(m)) - hardening;

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

            double dLambda =  strain_like(n).dot(Eelastic * depsilon_elpl);
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
            TrialPlastic_Strain += dLambda * strain_like(m);

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
            TrialStress = TrialStress - dLambda * Eelastic * strain_like(m);


            // Returning to Yield surface as recommended by Crisfield...
            if (INT_OPT_return_to_yield_surface[ASDP_TAG])
            {
                double yf_val_corr = yf(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& n_corr = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_corr = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

                double hardening_corr = yf.hardening( depsilon_elpl, m,  TrialStress, iv_storage, parameters_storage);
                double dLambda_corr = yf_val_corr / (
                                                     strain_like(n_corr).dot(Eelastic * strain_like(m_corr)) - hardening_corr
                                                 );
                TrialStress = TrialStress - dLambda_corr * Eelastic * strain_like(m_corr);
                TrialPlastic_Strain += dLambda_corr * strain_like(m_corr);
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
        depsilon.setZero();
        depsilon = strain_incr;

        const VoigtVector& sigma = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();

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
                double tol_yf = DBL_OPT_f_absolute_tol[ASDP_TAG];
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
                double den = strain_like(n).dot(Eelastic * strain_like(m)) - hardening;

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

                double dLambda =  strain_like(n).dot(Eelastic * depsilon_elpl);
                dLambda /= den;

                if (dLambda <= 0)
                {
                    dLambda = 0;
                }

                // Update the trial plastic strain.
                TrialPlastic_Strain += dLambda * strain_like(m);

                // This code iterates internal variables and updates the trial values
                iv_storage.apply([&m, &dLambda, this](auto & internal_variable)
                {
                    auto h = internal_variable.hardening_function(depsilon_elpl, m, intersection_stress, parameters_storage);
                    internal_variable.trial_value += dLambda * h;
                });

                //Correct the trial stress
                TrialStress = TrialStress - dLambda * Eelastic * strain_like(m);
            }

            // Returning to Yield surface as recommended by Crisfield...
            if (INT_OPT_return_to_yield_surface[ASDP_TAG])
            {
                double yf_val_corr = yf(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& n_corr = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_corr = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);

                double hardening_corr = yf.hardening( depsilon_elpl, m_corr,  TrialStress, iv_storage, parameters_storage);
                double dLambda_corr = yf_val_corr / (
                                                     strain_like(n_corr).dot(Eelastic * strain_like(m_corr)) - hardening_corr
                                                 );
                TrialStress = TrialStress - dLambda_corr * Eelastic * strain_like(m_corr);
                TrialPlastic_Strain += dLambda_corr * strain_like(m_corr);
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



    // int Backward_Euler(const VoigtVector & strain_incr)
    // {
    //     using namespace ASDPlasticMaterial3DGlobals;

    //     int errorcode = -1;

    //     static VoigtVector depsilon;
    //     depsilon *= 0;
    //     depsilon = strain_incr;

    //     const VoigtVector& sigma = CommitStress;
    //     const VoigtVector& epsilon = CommitStrain;

    //     iv_storage.revert_all();

    //     dsigma *= 0;
    //     intersection_stress *= 0;
    //     intersection_strain *= 0;

    //     VoigtMatrix Eelastic = et(sigma, parameters_storage);

    //     // Initial elastic predictor
    //     dsigma = Eelastic * depsilon;
    //     TrialStress = sigma + dsigma;
    //     TrialStrain = CommitStrain + depsilon;
    //     TrialPlastic_Strain = CommitPlastic_Strain;

    //     double yf_val_start = yf(sigma, iv_storage, parameters_storage);
    //     double yf_val_end = yf(TrialStress, iv_storage, parameters_storage);

    //     VoigtVector start_stress = CommitStress;
    //     VoigtVector end_stress = TrialStress;

    //     intersection_stress = start_stress;

    //     if ((yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end) //Elasticity
    //     {
    //         Stiffness = Eelastic;
    //         return 0;
    //     }
    //     else  //Plasticity - Backward Euler Implementation
    //     {
    //         depsilon_elpl = depsilon;
            
    //         // Initialize backward Euler iteration variables
    //         VoigtVector sigma_trial = TrialStress;
    //         VoigtVector sigma_n_plus_1 = sigma_trial + Eelastic * depsilon_elpl; // Initial guess
    //         double dLambda = 0.0;
            
    //         // Newton-Raphson iteration for backward Euler
    //         int max_iterations = INT_OPT_n_max_iterations[ASDP_TAG];
    //         double f_tol = DBL_OPT_f_absolute_tol[ASDP_TAG];
    //         double stress_tol = DBL_OPT_stress_absolute_tol[ASDP_TAG];
            
    //         for (int iter = 0; iter < max_iterations; iter++)
    //         {
    //             // Evaluate yield function at current stress state
    //             double f_val = yf(sigma_n_plus_1, iv_storage, parameters_storage);
                
    //             // Check convergence on yield function
    //             if (abs(f_val) <= f_tol)
    //             {
    //                 TrialStress = sigma_n_plus_1;
    //                 TrialPlastic_Strain = CommitPlastic_Strain + dLambda * pf(depsilon_elpl, sigma_n_plus_1, iv_storage, parameters_storage);
                    
    //                 // Check for NaN
    //                 double norm_trial_stress = TrialStress.transpose() * TrialStress;
    //                 if (norm_trial_stress != norm_trial_stress)
    //                 {
    //                     return -1;
    //                 }
                    
    //                 ComputeTangentStiffness();
    //                 return 0;
    //             }
                
    //             // Compute derivatives for Newton-Raphson
    //             const VoigtVector& n = yf.df_dsigma_ij(sigma_n_plus_1, iv_storage, parameters_storage);
    //             const VoigtVector& m = pf(depsilon_elpl, sigma_n_plus_1, iv_storage, parameters_storage);
    //             double hardening = yf.hardening(depsilon_elpl, m, sigma_n_plus_1, iv_storage, parameters_storage);
                
    //             // Build residual vector and Jacobian for backward Euler
    //             // Residual: R1 = sigma_n+1 - sigma_trial - E * (depsilon - dLambda * m)
    //             //          R2 = f(sigma_n+1, alpha_n+1)
                
    //             VoigtVector R1 = sigma_n_plus_1 - sigma_trial - Eelastic * (depsilon_elpl - dLambda * m);
    //             double R2 = f_val;
                
    //             // Check stress residual convergence
    //             double stress_residual_norm = R1.norm();
    //             if (stress_residual_norm <= stress_tol && abs(R2) <= f_tol)
    //             {
    //                 TrialStress = sigma_n_plus_1;
    //                 TrialPlastic_Strain = CommitPlastic_Strain + dLambda * m;
                    
    //                 // Check for NaN
    //                 double norm_trial_stress = TrialStress.transpose() * TrialStress;
    //                 if (norm_trial_stress != norm_trial_stress)
    //                 {
    //                     return -1;
    //                 }
                    
    //                 ComputeTangentStiffness();
    //                 return 0;
    //             }
                
    //             // Build Jacobian matrix for Newton-Raphson
    //             // J11 = I + dLambda * E * dm/dsigma
    //             // J12 = E * m
    //             // J21 = df/dsigma
    //             // J22 = -hardening
                
    //             VoigtMatrix I = VoigtMatrix::Identity();
    //             VoigtMatrix J11 = I; // Simplified - could add dLambda * E * dm/dsigma for better convergence
    //             VoigtVector J12 = Eelastic * m;
    //             VoigtVector J21 = n;
    //             double J22 = -hardening;
                
    //             // Solve Newton system using block elimination
    //             // [J11  J12] [Δσ   ]   [R1]
    //             // [J21  J22] [ΔdL  ] = [R2]
                
    //             // Eliminate to get: (J22 - J21 * J11^-1 * J12) * ΔdL = R2 - J21 * J11^-1 * R1
    //             double den = J22 - J21.transpose() * J12;
                
    //             if (abs(den) < MACHINE_EPSILON)
    //             {
    //                 cout << "Backward Euler - Singular Jacobian, den = " << den << endl;
    //                 return -1;
    //             }
                
    //             double delta_dLambda = (R2 - J21.transpose() * R1) / den;
    //             VoigtVector delta_sigma = -R1 - delta_dLambda * J12;
                
    //             // Update variables
    //             sigma_n_plus_1 += delta_sigma;
    //             dLambda += delta_dLambda;
                
    //             // Ensure plastic multiplier is non-negative
    //             if (dLambda < 0)
    //             {
    //                 dLambda = 0;
    //             }
                
    //             // Update internal variables based on current plastic multiplier
    //             iv_storage.apply([&m, &dLambda, &sigma_n_plus_1,  this](auto & iv)
    //             {
    //                 auto h = iv.hardening_function(depsilon_elpl, m, sigma_n_plus_1, parameters_storage);
    //                 iv.trial_value = iv.committed_value + dLambda * h;
    //             });
    //         }
            
    //         // If we reach here, Newton-Raphson did not converge
    //         cout << "Backward Euler - Newton-Raphson did not converge after " << max_iterations << " iterations" << endl;
    //         return -1;
    //     }

    //     return errorcode;
    // }


// int Backward_Euler(const VoigtVector & strain_incr)
// {
//     using namespace ASDPlasticMaterial3DGlobals;

//     int errorcode = -1;

//     static VoigtVector depsilon;
//     depsilon *= 0;
//     depsilon = strain_incr;

//     const VoigtVector& sigma   = CommitStress;
//     const VoigtVector& epsilon = CommitStrain;

//     iv_storage.revert_all();

//     dsigma *= 0;
//     intersection_stress *= 0;
//     intersection_strain *= 0;

//     // Elastic stiffness at committed state
//     VoigtMatrix Eelastic = et(sigma, parameters_storage);

//     // Elastic predictor
//     dsigma      = Eelastic * depsilon;
//     TrialStress = sigma + dsigma;
//     TrialStrain = CommitStrain + depsilon;
//     TrialPlastic_Strain = CommitPlastic_Strain;

//     // Simple elastic check at committed internal vars
//     const double f_tol     = DBL_OPT_f_absolute_tol[ASDP_TAG];
//     const double stress_tol= DBL_OPT_stress_absolute_tol[ASDP_TAG];
//     int max_iterations     = INT_OPT_n_max_iterations[ASDP_TAG];

//     double f_trial = yf(TrialStress, iv_storage, parameters_storage);

//     // Keep some of your existing bookkeeping
//     intersection_stress = CommitStress;

//     if (f_trial <= f_tol) {
//         // Elastic step
//         Stiffness = Eelastic;
//         return 0;
//     }

//     // --- Plastic corrector (Backward Euler / return mapping) ---
//     // Keep your depsilon_elpl usage to preserve pf()/hardening() signatures
//     depsilon_elpl = depsilon;

//     // Unknowns
//     VoigtVector sigma_n1 = TrialStress; // start from trial
//     double dLambda = 0.0;

//     for (int iter = 0; iter < max_iterations; ++iter)
//     {
//         // Gradients at (sigma_n1, committed internal vars)
//         const VoigtVector n = yf.df_dsigma_ij(sigma_n1, iv_storage, parameters_storage);               // ∂f/∂σ
//         const VoigtVector m = pf(depsilon_elpl, sigma_n1, iv_storage, parameters_storage);             // ∂g/∂σ (flow)
//         const double      H = yf.hardening(depsilon_elpl, m, sigma_n1, iv_storage, parameters_storage);// effective hardening slope (scalar)

//         // Residuals:
//         // R1 = σ_{n+1} - σ_trial + E * (dΛ m) = 0
//         // R2 = f(σ_{n+1}, α_n) + H dΛ = 0
//         const VoigtVector R1 = sigma_n1 - TrialStress + Eelastic * (dLambda * m);
//         const double      f_sigma = yf(sigma_n1, iv_storage, parameters_storage);
//         const double      R2 = f_sigma + H * dLambda;

//         // Convergence check
//         if (R1.norm() <= stress_tol && std::fabs(R2) <= f_tol) {
//             // Finalize state
//             const VoigtVector m_fin = pf(depsilon_elpl, sigma_n1, iv_storage, parameters_storage);
//             const VoigtVector n_fin = yf.df_dsigma_ij(sigma_n1, iv_storage, parameters_storage);
//             const double      H_fin = yf.hardening(depsilon_elpl, m_fin, sigma_n1, iv_storage, parameters_storage);

//             TrialStress = sigma_n1;
//             TrialStrain = CommitStrain + depsilon;
//             TrialPlastic_Strain = CommitPlastic_Strain + dLambda * m_fin;

//             // NaN guard (keep same style as your code)
//             double norm_trial_stress = TrialStress.transpose() * TrialStress;
//             if (norm_trial_stress != norm_trial_stress) {
//                 return -1;
//             }

//             // Consistent algorithmic tangent (non-associative allowed)
//             const VoigtVector Em = Eelastic * m_fin;
//             const VoigtVector En = Eelastic * n_fin;
//             double denom_tan = H_fin + n_fin.transpose() * Em;
//             if (std::fabs(denom_tan) < MACHINE_EPSILON) {
//                 Stiffness = Eelastic; // fallback
//             } else {
//                 // rank-1 update: C_ep = E - (E m) ⊗ (E n) / (H + n^T E m)
//                 Stiffness = Eelastic - (Em * En.transpose()) / denom_tan;
//             }

//             // Update internal vars to trial values: α_{n+1} = α_n + dΛ h
//             iv_storage.apply([&](auto & iv)
//             {
//                 auto h = iv.hardening_function(depsilon_elpl, m_fin, sigma_n1, parameters_storage);
//                 iv.trial_value = iv.committed_value + dLambda * h;
//             });

//             // Maintain dsigma for downstream use
//             dsigma = TrialStress - CommitStress;
//             return 0;
//         }

//         // Newton step with J11 = I (neglecting dm/dσ for robustness)
//         const VoigtVector Em = Eelastic * m;
//         const double denom = H + n.transpose() * Em;  // Schur complement denominator

//         if (std::fabs(denom) < MACHINE_EPSILON) {
//             // singular / ill-conditioned -> fail gracefully
//             return -1;
//         }

//         // ΔdΛ = (n^T R1 - R2) / (H + n^T E m)
//         const double delta_dLambda = (n.transpose() * R1 - R2) / denom;

//         // Project to dΛ >= 0
//         double d_dLambda = delta_dLambda;
//         if (dLambda + d_dLambda < 0.0) d_dLambda = -dLambda;

//         // Δσ = -(R1 + E m ΔdΛ)
//         const VoigtVector delta_sigma = -(R1 + Em * d_dLambda);

//         // Update unknowns
//         sigma_n1 += delta_sigma;
//         dLambda  += d_dLambda;

//         // Basic sanity
//         double norm_s = sigma_n1.transpose() * sigma_n1;
//         if (norm_s != norm_s) {
//             return -1;
//         }
//     }

//     // If we reach here, Newton failed to converge
//     std::cout << "Backward Euler - Newton-Raphson did not converge after " << max_iterations << " iterations\n";
//     return -1;
// }



    // Return mapping to the apex of a cone-shaped yield surface (yf_has_apex_return).
    // At the apex f is not differentiable and the cutting-plane iteration stalls
    // (with psi = 0 the flow has no volumetric part and cannot reduce p at all).
    // The trial stress is in the apex region when removing its deviatoric part
    // along the flow direction at the trial stress, dsigma = -dLambda E W m,
    // leaves the mean stress at or beyond the apex:
    //     dLambda_q = (s_tr : s_tr) / (s_tr : dev(E W m)),
    //     p_tr - dLambda_q tr(E W m) / 3 >= p_apex          (p = I1 / 3, tension positive)
    // For psi = 0 this is p_tr >= p_apex. The stress is then the apex, the
    // plastic strain increment is C (sigma_tr - sigma_apex) exactly, and dLambda
    // is its projection on W m (used only to update the internal variables).
    // For a perfectly plastic cone the consistent tangent at the apex is zero.
    // It is replaced by E / 2 (the Secant tangent) for every tangent type but
    // Elastic (E): a near-zero tangent (1e-6 E was used first) leaves nodes
    // whose Gauss points are all at the apex, a pocket of soil in tension,
    // almost without stiffness, and the global solve then returns
    // displacement corrections of 1e10 and more without any solver warning.
    bool apex_return_applied = false;

    bool Apex_Return(const VoigtMatrix& Eelastic, const VoigtVector& depsilon)
    {
        if constexpr (!yf_has_apex_return<YieldFunctionType>::value)
        {
            return false;
        }
        else
        {
            using namespace ASDPlasticMaterial3DGlobals;

            const VoigtVector sigma_apex = yf.apex_stress(iv_storage, parameters_storage);
            if (!sigma_apex.allFinite())   // e.g. phi = 0: no apex
                return false;

            const double p_apex = sigma_apex.getI1() / 3.0;
            const double p_tr   = TrialStress.getI1() / 3.0;
            const VoigtVector s_tr = TrialStress.deviator();
            const double s_norm = std::sqrt(std::max(tensor_dot_stress_like(s_tr, s_tr), 0.0));
            const double scale  = std::max({std::abs(p_apex), std::abs(p_tr), s_norm});

            VoigtVector m = pf(depsilon, TrialStress, iv_storage, parameters_storage);
            if (!m.allFinite())
                m.setZero();
            const VoigtVector v = Eelastic * strain_like(m);
            const double sv = tensor_dot_stress_like(s_tr, v.deviator());

            bool in_apex_region;
            if (s_norm <= 1e-12 * scale || !(sv > 0.0))
                in_apex_region = p_tr >= p_apex;
            else
                in_apex_region = p_tr - (s_norm * s_norm / sv) * v.getI1() / 3.0 >= p_apex;
            if (!in_apex_region)
                return false;

            // plastic strain increment C (sigma_tr - sigma_apex), in strain-like (engineering) components
            Eigen::Matrix<double, 6, 6> E6;
            Eigen::Matrix<double, 6, 1> r6;
            for (int i = 0; i < 6; ++i) {
                r6(i) = TrialStress(i) - sigma_apex(i);
                for (int j = 0; j < 6; ++j) E6(i, j) = Eelastic(i, j);
            }
            const Eigen::Matrix<double, 6, 1> e6 = E6.ldlt().solve(r6);
            if (!e6.allFinite())
                return false;
            VoigtVector depsp;
            for (int i = 0; i < 6; ++i) depsp(i) = e6(i);

            const VoigtVector Wm = strain_like(m);
            const double WmWm = Wm.dot(Wm);
            const double dLambda = WmWm > MACHINE_EPSILON ? std::max(depsp.dot(Wm) / WmWm, 0.0) : 0.0;

            iv_storage.apply([&](auto & internal_variable)
            {
                auto h = internal_variable.hardening_function(depsilon, m, sigma_apex, parameters_storage);
                internal_variable.trial_value += dLambda * h;
            });

            // with hardening the apex moves: use this return only if the apex is still on the surface
            if (std::abs(yf(sigma_apex, iv_storage, parameters_storage)) > DBL_OPT_f_absolute_tol[ASDP_TAG])
            {
                iv_storage.revert_all();
                return false;
            }

            TrialStress = sigma_apex;
            TrialPlastic_Strain = CommitPlastic_Strain + depsp;

            if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Elastic)
                Stiffness = Eelastic;
            else
                Stiffness = 0.5 * Eelastic;

            apex_return_applied = true;
            return true;
        }
    }

    int Backward_Euler(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        int errorcode = -1;

        int    max_iter = INT_OPT_n_max_iterations[ASDP_TAG];
        double tol_yf   = DBL_OPT_f_absolute_tol[ASDP_TAG]; 

        apex_return_applied = false;
        be_last_iter = 0;

        // -------- setup
        static VoigtVector depsilon;
        depsilon = strain_incr;

        const VoigtVector& sigma   = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();

        VoigtMatrix Eelastic = et(sigma, parameters_storage);

        // -------- elastic predictor
        dsigma      = Eelastic * depsilon;
        TrialStress = sigma + dsigma;
        TrialStrain = epsilon + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        // cout << "BE - TrialStress = " << TrialStress.transpose() <<  endl;
        // cout << "BE - CommitStress = " << CommitStress.transpose() <<  endl;
        
        const double yf_val_start = yf(sigma,        iv_storage, parameters_storage);
        // cout << "----->  yf_val_start = " << yf_val_start <<  endl;
        
        const double yf_val_end   = yf(TrialStress,  iv_storage, parameters_storage);
        // cout << "----->  yf_val_end = " << yf_val_end << " tol_yf = " << tol_yf <<  endl;

        // purely elastic or moving deeper inside the surface
        if ( (yf_val_start <= 0.0 && yf_val_end <= 0.0) || (yf_val_start - yf_val_end > tol_yf) ) {
            // cout << "BE - ELASTIC!" << endl << endl;
            Stiffness = Eelastic;
            return 0;
        }


        // -------- trial stress beyond the apex of a cone (e.g. Mohr-Coulomb in tension)
        if (Apex_Return(Eelastic, depsilon))
            return 0;

        // -------- plastic correction (Backward Euler)

        double dLambda = 0.0;


        // cout << "BE - Plastic! Begin iterations----------" << endl << endl;

        be_last_iter = 0;
        for (int iter = 0; iter < max_iter; ++iter)
        {
            be_last_iter = iter + 1;
            // directions at current (trial) end state
            const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);

            // use total step depsilon as proxy argument; no derivatives of m required
            const VoigtVector& m = pf(depsilon, TrialStress, iv_storage, parameters_storage);

            // hardening at end state
            const double H = yf.hardening(depsilon, m, TrialStress, iv_storage, parameters_storage);

            // consistency residual Φ(σ_{n+1}) = 0
            const double Phi = yf(TrialStress, iv_storage, parameters_storage);

            if (std::abs(Phi) < tol_yf) {
                GLOBAL_INT_max_iter[ASDP_TAG] = std::max(GLOBAL_INT_max_iter[ASDP_TAG], iter);
                GLOBAL_DBL_max_error[ASDP_TAG] = std::max(GLOBAL_DBL_max_error[ASDP_TAG], std::abs(Phi));

                // cout << "  =>  n = " << n.transpose() << endl;
                // cout << "  =>  m = " << m.transpose() << endl;
                // cout << "  =>  H = " << H << endl;
                break; // converged
            }

            // dΦ/dλ ≈ - n^T E m + H
            // const double dPhi_dLambda = - (n.transpose() * Eelastic * m) + H;
            const double nEm = tensor_dot_stress_like(n, Eelastic * strain_like(m));
            const double dPhi_dLambda = H - nEm;   // == - n^T E m + H

            if (std::abs(dPhi_dLambda) < MACHINE_EPSILON) {
                // singular local tangent
                opserr << "ASDPlasticMaterial3D (tag " << ASDP_TAG << ") Backward_Euler: singular local tangent (n:E:m - H = 0)\n";
                return -1;
            }

            const double deltaLambda = - Phi / dPhi_dLambda;

            // cout << " ---- iter = " << iter << " / " << max_iter << endl;
            // cout << "  =>  n = " << n.transpose() << endl;
            // cout << "  =>  m = " << m.transpose() << endl;
            // cout << "  =>  H = " << H << endl;
            // cout << "  =>  nEm = " << nEm << endl;
            // cout << "  =>  Phi = " << Phi << endl;
            // cout << "  =>  dPhi_dLambda = H - nEm = " << dPhi_dLambda << endl;
            // cout << "  =>  deltaLambda = -Phi / dPhi_dLambda = " << deltaLambda << endl;
            // cout << "  =>  dLambda = " << dLambda << endl ;

            // keep λ >= 0
            if (dLambda + deltaLambda < 0.0) {
                // step cannot be plastic; fall back to elastic in this rare case
                Stiffness = Eelastic;
                return 0;
            }

            // incremental updates (use delta to avoid re-summing from commit each iter)
            dLambda += deltaLambda;
            TrialStress          = TrialStress - deltaLambda * (Eelastic * strain_like(m));
            TrialPlastic_Strain  = TrialPlastic_Strain + deltaLambda * strain_like(m);
            // cout << "  =>  dLambda + deltaLambda = " << dLambda << endl;
            // cout << "  =>  CommitStress = " << CommitStress.transpose() << endl;
            // cout << "  =>  TrialStress = " << TrialStress.transpose() << endl;
            // cout << "  =>  TrialPlastic_Strain = " << TrialPlastic_Strain.transpose() << endl;
            iv_storage.apply([&](auto & internal_variable)
            {
                auto h = internal_variable.hardening_function(depsilon, m, TrialStress, parameters_storage);
                internal_variable.trial_value += deltaLambda * h;
                // cout << "  => " <<  internal_variable << endl;
            });

            // NaN guard
            const double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (!(norm_trial_stress == norm_trial_stress)) { // NaN chec
               opserr << "ASDPlasticMaterial3D (tag " << ASDP_TAG << ") Backward_Euler: NaN stress\n";
                return -1;
            }
        }
        // cout << "BE - END iterations----------" << endl << endl;


        if (!skip_BE_tangent)   // Full_Backward_Euler uses this state only as its starting point
            ComputeTangentStiffness();

        return 0;
    }

    // Fully implicit Backward Euler (closest-point projection). Unlike
    // Backward_Euler, which iterates on dLambda alone and never differentiates
    // m or the hardening, this solves the whole return-mapping system with
    // Newton. Unknowns x = [sigma (6), q (every internal variable, flattened),
    // dLambda] and residual
    //     R_sigma = sigma - sigma_trial + dLambda E strain_like(m(sigma, q))
    //     R_q     = q - q_n - dLambda h(sigma, q)
    //     R_f     = f(sigma, q)
    // E is the elastic stiffness at the committed stress, as in Backward_Euler.
    // The Jacobian's sigma and q columns are central differences of R (option
    // derivatives Numerical). With derivatives Closed_Form, the closed forms
    // the components provide replace their blocks: df/dsigma from the yield
    // function, and dm/dsigma from flow directions that declare has_dm_dsigma.
    // The dLambda column is exact. Newton starts from the Backward_Euler
    // (cutting-plane) state; if it does not converge, as at the edges of
    // Mohr-Coulomb or Hoek-Brown where m jumps, that state is kept.
    // tangent_type Algorithmic gives the consistent tangent (J^-1)_sigma,sigma E.
    int Full_Backward_Euler(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        const int    max_iter = INT_OPT_n_max_iterations[ASDP_TAG];
        const double tol_yf   = DBL_OPT_f_absolute_tol[ASDP_TAG];
        const double tol_sig  = DBL_OPT_stress_absolute_tol[ASDP_TAG];
        const bool   closed_form = INT_OPT_derivative_type[ASDP_TAG] == ASDPlasticMaterial3D_Derivative_Type::Closed_Form;
        const double h_rel    = DBL_OPT_fd_perturbation[ASDP_TAG];

        static VoigtVector depsilon;
        depsilon = strain_incr;

        const VoigtVector& sigma_n   = CommitStress;
        const VoigtVector& epsilon_n = CommitStrain;

        iv_storage.revert_all();

        dsigma.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();

        const VoigtMatrix Eelastic = et(sigma_n, parameters_storage);

        // -------- elastic predictor (same test as Backward_Euler)
        const VoigtVector sigma_trial = sigma_n + Eelastic * depsilon;
        TrialStress = sigma_trial;
        TrialStrain = epsilon_n + depsilon;
        TrialPlastic_Strain = CommitPlastic_Strain;

        const double yf_val_start = yf(sigma_n,     iv_storage, parameters_storage);
        const double yf_val_end   = yf(sigma_trial, iv_storage, parameters_storage);
        if ( (yf_val_start <= 0.0 && yf_val_end <= 0.0) || (yf_val_start - yf_val_end > tol_yf) ) {
            Stiffness = Eelastic;
            return 0;
        }

        // -------- flatten the internal variables
        int nq = 0;
        iv_storage.apply([&](auto & iv) { nq += iv.size(); });
        const int N = 7 + nq;
        const int iL = 6 + nq;   // position of dLambda

        Eigen::VectorXd q_n(nq);
        std::vector<int> iv_of(nq);   // which internal variable each entry belongs to
        int n_iv = 0;
        {
            int k = 0;
            iv_storage.apply([&](auto & iv) {
                for (int j = 0; j < iv.size(); ++j) { iv_of[k] = n_iv; q_n(k++) = iv.committed_value(j); }
                ++n_iv;
            });
        }

        auto set_q = [&](const Eigen::VectorXd & q) {
            int k = 0;
            iv_storage.apply([&](auto & iv) { for (int j = 0; j < iv.size(); ++j) iv.trial_value(j) = q(k++); });
        };

        // residual at (sigma, q, dLambda); also returns m and h there
        VoigtVector m_last;
        Eigen::VectorXd h_last(nq);
        auto residual = [&](const VoigtVector & sig, const Eigen::VectorXd & q, double dL, Eigen::VectorXd & R) {
            set_q(q);
            const VoigtVector m = pf(depsilon, sig, iv_storage, parameters_storage);
            Eigen::VectorXd h(nq);
            int k = 0;
            iv_storage.apply([&](auto & iv) {
                auto hk = iv.hardening_function(depsilon, m, sig, parameters_storage);
                for (int j = 0; j < iv.size(); ++j) h(k++) = hk(j);
            });
            const VoigtVector Rsig = sig - sigma_trial + dL * (Eelastic * strain_like(m));
            for (int i = 0; i < 6; ++i) R(i) = Rsig(i);
            for (int j = 0; j < nq; ++j) R(6 + j) = q(j) - q_n(j) - dL * h(j);
            R(iL) = yf(sig, iv_storage, parameters_storage);
            m_last = m;
            h_last = h;
        };

        // scales for the finite differences and for the line search
        const double sig_scale = std::max({sigma_trial.cwiseAbs().maxCoeff(), sigma_n.cwiseAbs().maxCoeff(), tol_sig});
        // f changes by n . dsigma, so |n| sig_scale is the yield-function value of a stress-scale error
        const double n_scale   = strain_like(yf.df_dsigma_ij(sigma_trial, iv_storage, parameters_storage)).cwiseAbs().maxCoeff();
        const double f_scale   = std::max({std::abs(yf_val_end), n_scale * sig_scale, tol_yf});

        VoigtVector sig = sigma_trial;
        Eigen::VectorXd q = q_n;
        double dL = 0.0;

        // start from the cutting-plane Backward_Euler state: Newton from the
        // elastic trial overshoots on strongly curved surfaces (e.g. HoekBrown)
        // Backward_Euler returns 0 even when its loop ran out of iterations: its
        // state is used (as starting point, and as fallback) only if f is within
        // tolerance there
        bool cp_start = false;
        skip_BE_tangent = true;
        const int be_flag = Backward_Euler(strain_incr);
        skip_BE_tangent = false;
        if (apex_return_applied)   // the apex is the exact solution; there is no Newton system there
            return be_flag;
        const VoigtVector cp_sig = TrialStress;
        const VoigtVector cp_epsp = TrialPlastic_Strain;
        Eigen::VectorXd cp_q(nq);
        {
            int k = 0;
            iv_storage.apply([&](auto & iv) { for (int j = 0; j < iv.size(); ++j) cp_q(k++) = iv.trial_value(j); });
        }
        const bool cp_converged = be_flag == 0 && cp_sig.allFinite()
            && std::abs(yf(cp_sig, iv_storage, parameters_storage)) <= tol_yf;
        if (cp_converged && (cp_epsp - CommitPlastic_Strain).cwiseAbs().maxCoeff() > 0)
        {
            sig = cp_sig;
            q = cp_q;
            cp_start = true;
        }

        Eigen::VectorXd R(N), Rp(N), Rm(N);
        Eigen::MatrixXd J(N, N);
        Eigen::VectorXd q_step(nq);   // expected change of each q, for its perturbation size

        auto jacobian = [&]() {
            // R, m_last, h_last must be current at (sig, q, dL)
            const VoigtVector m0 = m_last;
            const Eigen::VectorXd h0 = h_last;
            for (int i = 0; i < 6; ++i) {
                const double d = h_rel * std::max(std::abs(sig(i)), sig_scale);
                VoigtVector sp = sig, sm = sig;
                sp(i) += d;
                sm(i) -= d;
                residual(sp, q, dL, Rp);
                residual(sm, q, dL, Rm);
                J.col(i) = (Rp - Rm) / (2 * d);
            }
            for (int j = 0; j < nq; ++j) {
                double d = h_rel * std::max({std::abs(q(j)), std::abs(q_n(j)), std::abs(q_step(j))});
                if (!(d > 0)) d = h_rel * sig_scale;
                Eigen::VectorXd qp = q, qm = q;
                qp(j) += d;
                qm(j) -= d;
                residual(sig, qp, dL, Rp);
                residual(sig, qm, dL, Rm);
                J.col(6 + j) = (Rp - Rm) / (2 * d);
            }
            // dLambda column: exact
            J.col(iL).setZero();
            const VoigtVector Em = Eelastic * strain_like(m0);
            for (int i = 0; i < 6; ++i) J(i, iL) = Em(i);
            for (int j = 0; j < nq; ++j) J(6 + j, iL) = -h0(j);

            if (closed_form) {
                set_q(q);
                const VoigtVector n = strain_like(yf.df_dsigma_ij(sig, iv_storage, parameters_storage));
                for (int i = 0; i < 6; ++i) J(iL, i) = n(i);
                if constexpr (pf_has_dm_dsigma<PlasticFlowType>::value) {
                    const VoigtMatrix dm = pf.dm_dsigma_ij(depsilon, sig, iv_storage, parameters_storage);
                    VoigtMatrix W = VoigtMatrix::Identity();
                    W(3, 3) = W(4, 4) = W(5, 5) = 2.0;
                    const VoigtMatrix Jss = VoigtMatrix::Identity() + dL * (Eelastic * W * dm);
                    J.block(0, 0, 6, 6) = Jss;
                }
            }
            // leave the state at (sig, q, dL)
            residual(sig, q, dL, R);
        };

        // scale of each entry: its own size or expected change, but at least 1e-8
        // of the largest entry of the same internal variable (its components share
        // units; round-off-sized entries must not set their own scale). An internal
        // variable with no value and no evolution so far (scale 0) stays at q_n
        // and is kept out of the merit and the convergence test.
        std::vector<double> iv_max(n_iv, 0.0);
        auto update_iv_max = [&]() {
            std::fill(iv_max.begin(), iv_max.end(), 0.0);
            for (int j = 0; j < nq; ++j)
                iv_max[iv_of[j]] = std::max({iv_max[iv_of[j]], std::abs(q_n(j)), std::abs(q_step(j))});
        };
        auto q_scale = [&](int j) {
            return std::max({std::abs(q_n(j)), std::abs(q_step(j)), 1e-8 * iv_max[iv_of[j]]});
        };

        auto merit = [&](const Eigen::VectorXd & Rv) {
            double r = 0.0;
            for (int i = 0; i < 6; ++i) r = std::max(r, std::abs(Rv(i)) / sig_scale);
            for (int j = 0; j < nq; ++j)
                if (q_scale(j) > 0) r = std::max(r, std::abs(Rv(6 + j)) / q_scale(j));
            return std::max(r, std::abs(Rv(iL)) / f_scale);
        };

        // converged: the absolute tolerances are met, or the scaled residual is at
        // round-off level (below what components that difference internally resolve)
        auto converged = [&](const Eigen::VectorXd & Rv) {
            if (merit(Rv) < 1e-10) return true;
            if (std::abs(Rv(iL)) > tol_yf) return false;
            for (int i = 0; i < 6; ++i) if (std::abs(Rv(i)) > tol_sig) return false;
            for (int j = 0; j < nq; ++j)
                if (q_scale(j) > 0 && std::abs(Rv(6 + j)) > 1e-10 * std::max({std::abs(q(j)), iv_max[iv_of[j]]})) return false;
            return true;
        };

        residual(sig, q, dL, R);
        if (cp_start)
        {
            const VoigtVector mW = strain_like(m_last);
            const double mm = mW.dot(mW);
            if (mm > 0)
                dL = std::max(0.0, (TrialPlastic_Strain - CommitPlastic_Strain).dot(mW) / mm);
            residual(sig, q, dL, R);
        }
        // first estimate of each q's change, from a linearized dLambda
        {
            const VoigtVector n = strain_like(yf.df_dsigma_ij(sig, iv_storage, parameters_storage));
            const double den = n.dot(Eelastic * strain_like(m_last));
            const double dL0 = den > 0 ? std::abs(R(iL)) / den : 0.0;
            for (int j = 0; j < nq; ++j)
                q_step(j) = std::max({std::abs(dL0 * h_last(j)), std::abs(dL * h_last(j)), std::abs(q(j) - q_n(j))});
            update_iv_max();
        }

        bool ok = false;
        bool newton_nonfinite_step = false;
        int iter = 0;
        for (; iter < max_iter; ++iter)
        {
            if (converged(R)) { ok = true; break; }

            jacobian();
            Eigen::PartialPivLU<Eigen::MatrixXd> lu(J);
            const Eigen::VectorXd dx = lu.solve(-R);
            if (!dx.allFinite()) { newton_nonfinite_step = true; break; }

            // backtracking on the scaled residual
            const double phi0 = merit(R);
            double alpha = 1.0;
            VoigtVector sig_new;
            Eigen::VectorXd q_new(nq);
            double dL_new = dL;
            bool decreased = false;
            for (int ls = 0; ls < 10; ++ls) {
                for (int i = 0; i < 6; ++i) sig_new(i) = sig(i) + alpha * dx(i);
                q_new = q + alpha * dx.segment(6, nq);
                for (int j = 0; j < nq; ++j) if (!(q_scale(j) > 0)) q_new(j) = q_n(j);
                dL_new = dL + alpha * dx(iL);
                residual(sig_new, q_new, dL_new, Rp);
                if (Rp.allFinite() && merit(Rp) < (1.0 - 1e-4 * alpha) * phi0) { decreased = true; break; }
                alpha *= 0.5;
            }
            // no descent left with f within tolerance and the scaled residual tiny:
            // Newton is at the noise floor of the components' own differences
            if (!decreased && phi0 < 1e-8 && std::abs(R(iL)) <= tol_yf) {
                residual(sig, q, dL, R);   // back to the current state
                ok = true;
                break;
            }
            sig = sig_new;
            q = q_new;
            dL = dL_new;
            R = Rp;
            for (int j = 0; j < nq; ++j) q_step(j) = std::max(std::abs(q_step(j)), std::abs(q(j) - q_n(j)));
            update_iv_max();
        }

        if (!ok || dL < 0.0 || !sig.allFinite()) {
            // keep the converged cutting-plane state for this step, if there is one.
            // Fallbacks (Backward_Euler state used) are reported 5 times per
            // combination and process; failures (-1 returned, which fails the
            // element and the domain update) 20 times, each with the full state
            // so the step can be replayed on one material point.
            const bool report = cp_converged ? (++FBE_fallbacks[ASDP_TAG] <= 5)
                                             : (++FBE_failures[ASDP_TAG] <= 20);
            if (report)
            {
                // stress invariants of each state, to tell a corner (Lode angle near
                // +-30 deg) from the apex (p near the tensile apex, q near 0)
                auto invariants = [](const VoigtVector& st) {
                    char buf[128];
                    std::snprintf(buf, sizeof(buf), "p = %.6g, q = %.6g, Lode = %.3f deg",
                                  st.getI1() / 3.0, std::sqrt(3.0 * std::max(st.getJ2(), 0.0)),
                                  st.lodeAngle() * 180.0 / M_PI);
                    return std::string(buf);
                };
                set_q(cp_q);
                const double f_cp = yf(cp_sig, iv_storage, parameters_storage);
                const int nrep = cp_converged ? FBE_fallbacks[ASDP_TAG] : FBE_failures[ASDP_TAG];
                const int cap = cp_converged ? 5 : 20;
                opserr << "ASDPlasticMaterial3D (tag " << ASDP_TAG << ") Full_Backward_Euler did not converge ("
                       << nrep << (nrep == cap ? ", further ones are counted, not reported" : "") << ")"
                       << (cp_converged ? " -- using the converged Backward_Euler state\n" : " -- Backward_Euler did not converge either: returning -1\n")
                       << "    p = I1/3 (tension positive), q = sqrt(3 J2)\n"
                       << "    committed: " << invariants(sigma_n).c_str() << ", f = " << yf_val_start << "\n"
                       << "    trial:     " << invariants(sigma_trial).c_str() << ", f = " << yf_val_end << "\n"
                       << "    Newton:    " << invariants(sig).c_str() << ", f = " << R(iL) << ", dLambda = " << dL << "\n"
                       << "    Backward_Euler: " << invariants(cp_sig).c_str() << ", f = " << f_cp << "\n";
                if (!cp_converged)
                {
                    // why each solver stopped
                    char buf[256];
                    if (newton_nonfinite_step)
                        std::snprintf(buf, sizeof(buf), "non-finite Newton step at iteration %d", iter);
                    else if (!ok)
                        std::snprintf(buf, sizeof(buf), "no convergence in %d iterations (merit %.3g)", iter, merit(R));
                    else if (!sig.allFinite())
                        std::snprintf(buf, sizeof(buf), "non-finite stress");
                    else
                        std::snprintf(buf, sizeof(buf), "converged to dLambda = %.6g < 0 after %d iterations", dL, iter);
                    opserr << "    Newton stopped: " << buf << "\n";
                    if (be_flag != 0)
                        std::snprintf(buf, sizeof(buf), "returned %d", be_flag);
                    else if (!cp_sig.allFinite())
                        std::snprintf(buf, sizeof(buf), "non-finite stress");
                    else
                        std::snprintf(buf, sizeof(buf), "|f| = %.6g > f_absolute_tol = %.6g after %d of %d iterations",
                                      std::abs(f_cp), tol_yf, be_last_iter, max_iter);
                    opserr << "    Backward_Euler stopped: " << buf << "\n";
                    // full Voigt vectors (11 22 33 12 23 13; strains with engineering shear)
                    auto row = [](const char* name, const VoigtVector& v) {
                        char b[512];
                        std::snprintf(b, sizeof(b), "    %-26s %.17g %.17g %.17g %.17g %.17g %.17g\n",
                                      name, v(0), v(1), v(2), v(3), v(4), v(5));
                        return std::string(b);
                    };
                    opserr << row("committed stress:", sigma_n).c_str()
                           << row("committed strain:", epsilon_n).c_str()
                           << row("committed plastic strain:", CommitPlastic_Strain).c_str()
                           << row("strain increment:", depsilon).c_str()
                           << row("trial stress:", sigma_trial).c_str()
                           << row("Newton stress:", sig).c_str()
                           << row("Backward_Euler stress:", cp_sig).c_str();
                    if (nq > 0)
                    {
                        opserr << "    committed internal variables:";
                        for (int j = 0; j < nq; ++j) { std::snprintf(buf, sizeof(buf), " %.17g", q_n(j)); opserr << buf; }
                        opserr << "\n";
                    }
                }
            }
            if (!cp_converged)
                return -1;
            TrialStress = cp_sig;
            TrialPlastic_Strain = cp_epsp;
            set_q(cp_q);
            depsilon_elpl = depsilon;
            ComputeTangentStiffness();
            return 0;
        }

        GLOBAL_INT_max_iter[ASDP_TAG] = std::max(GLOBAL_INT_max_iter[ASDP_TAG], iter);
        GLOBAL_DBL_max_error[ASDP_TAG] = std::max(GLOBAL_DBL_max_error[ASDP_TAG], std::abs(R(iL)));

        // -------- accept the converged state (residual() left q and m_last there)
        TrialStress = sig;
        TrialPlastic_Strain = CommitPlastic_Strain + dL * strain_like(m_last);
        set_q(q);

        if (INT_OPT_tangent_operator_type[ASDP_TAG] == ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic)
        {
            // dR/d(strain) = [-E; 0; 0], so dsigma/dstrain = (J^-1)_{sigma,sigma} E
            jacobian();
            Eigen::PartialPivLU<Eigen::MatrixXd> lu(J);
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(N, 6);
            B.block(0, 0, 6, 6) = Eelastic;
            const Eigen::MatrixXd X = lu.solve(B);
            Stiffness = X.block(0, 0, 6, 6);
            set_q(q);
        }
        else
        {
            depsilon_elpl = depsilon;
            ComputeTangentStiffness();
        }

        return 0;
    }

    int Backward_Euler_LineSearch(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        int errorcode = -1;

        // ------------------ setup ------------------
        static VoigtVector depsilon;
        depsilon.setZero();
        depsilon = strain_incr;

        const VoigtVector& sigma   = CommitStress;
        const VoigtVector& epsilon = CommitStrain;

        iv_storage.revert_all();

        dsigma.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();

        VoigtMatrix Eelastic = et(sigma, parameters_storage);

        // tolerancias
        const double tol_abs = (DBL_OPT_f_absolute_tol[ASDP_TAG] > 0.0) ? DBL_OPT_f_absolute_tol[ASDP_TAG] : 1e-10;
        const double tol_rel = 1e-8;
        const int    max_iter = 30;

        auto converged = [&](double Phi, double Phi_scale)->bool {
            double tol = std::max(tol_abs, tol_rel * std::max(1.0, Phi_scale));
            return std::abs(Phi) <= tol;
        };

        // ------------------ solver (un sub-paso Δε) ------------------
        auto solve_increment = [&](const VoigtVector& dEps)->bool
        {
            // Predictor elástico desde el estado commit
            TrialStrain = epsilon + dEps;
            TrialPlastic_Strain = CommitPlastic_Strain;
            TrialStress = sigma + Eelastic * dEps;

            // Estado interno trial = commit
            iv_storage.revert_all();

            const double yf_start = yf(sigma,       iv_storage, parameters_storage);
            const double yf_end   = yf(TrialStress, iv_storage, parameters_storage);

            // Misma lógica que usabas: puramente elástico o moviéndose "hacia adentro"
            if ( (yf_start <= 0.0 && yf_end <= 0.0) || (yf_start > yf_end) ) {
                Stiffness = Eelastic;
                return true;
            }

            // ---------- Newton escalar con backtracking ----------
            double dLambda = 0.0;

            // Arranque seguro (si el trial está fuera)
            {
                const VoigtVector& n0 = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m0 = pf(dEps,        TrialStress, iv_storage, parameters_storage);
                const double       H0 = yf.hardening(dEps, m0,       TrialStress, iv_storage, parameters_storage);
                const double       Phi0 = yf(TrialStress, iv_storage, parameters_storage);
                const double       nEm0 = strain_like(n0).dot(Eelastic * strain_like(m0));
                const double       denom0 = nEm0 - H0; // == -(H0 - nEm0) == -dPhi/dλ

                if (denom0 > MACHINE_EPSILON && Phi0 > 0.0) {
                    // dλ inicial ≥ 0
                    double dL0 = Phi0 / denom0;
                    // límite físico simple
                    double dLmax = dEps.norm() / std::max(m0.norm(), 1e-16);
                    dLambda = std::min(dL0, dLmax);

                    // aplicar arranque
                    TrialStress         = TrialStress - dLambda * (Eelastic * strain_like(m0));
                    TrialPlastic_Strain = TrialPlastic_Strain + dLambda * strain_like(m0);
                    iv_storage.apply([&](auto & iv){
                        auto h = iv.hardening_function(dEps, m0, TrialStress, parameters_storage);
                        iv.trial_value += dLambda * h;
                    });
                }
            }

            bool newton_ok = false;

            for (int iter = 0; iter < max_iter; ++iter)
            {
                // (opcional) si E depende fuerte de σ, descomenta:
                // Eelastic = et(TrialStress, parameters_storage);

                const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m = pf(dEps,        TrialStress, iv_storage, parameters_storage);
                const double       H = yf.hardening(dEps, m,        TrialStress, iv_storage, parameters_storage);

                const double Phi = yf(TrialStress, iv_storage, parameters_storage);
                const double Phi_scale = std::abs(strain_like(n).dot(TrialStress)) + 1.0;

                if (converged(Phi, Phi_scale)) { newton_ok = true; break; }

                const double nEm = strain_like(n).dot(Eelastic * strain_like(m));
                const double dPhi_dLambda = H - nEm;  // total deriv. aprox

                if (!std::isfinite(dPhi_dLambda) || std::abs(dPhi_dLambda) < 1e-20) {
                    // Singular → deja que el substepping maneje
                    newton_ok = false;
                    break;
                }

                // Paso de Newton propuesto
                double deltaLambda = - Phi / dPhi_dLambda;

                // Limitar tamaño (previene runaway)
                const double dLmax = dEps.norm() / std::max(m.norm(), 1e-16);
                if (std::abs(deltaLambda) > dLmax)
                    deltaLambda = (deltaLambda > 0.0 ? 1.0 : -1.0) * dLmax;

                // Enforce λ ≥ 0
                if (dLambda + deltaLambda < 0.0)
                    deltaLambda = -dLambda * 0.5; // reduce para no cruzar a negativo

                // -------- line search (backtracking) con Φ linealizado --------
                double alpha = 1.0;
                const double c = 1e-4;
                double dl_accepted = 0.0;

                for (int ls = 0; ls < 8; ++ls) {
                    const double dl = alpha * deltaLambda;

                    // Predicción lineal de Φ
                    const double Phi_pred = Phi + dPhi_dLambda * dl;

                    if (std::abs(Phi_pred) <= (1.0 - c*alpha) * std::abs(Phi)) {
                        dl_accepted = dl;
                        break;
                    }
                    alpha *= 0.5;
                }

                if (dl_accepted == 0.0) {
                    // No se pudo aceptar un paso útil → dejar a substepping
                    newton_ok = false;
                    break;
                }

                // Aplicar actualización aceptada
                const VoigtVector Em = Eelastic * strain_like(m);

                TrialStress         = TrialStress - dl_accepted * Em;
                TrialPlastic_Strain = TrialPlastic_Strain + dl_accepted * strain_like(m);

                iv_storage.apply([&](auto & iv){
                    auto h = iv.hardening_function(dEps, m, TrialStress, parameters_storage);
                    iv.trial_value += dl_accepted * h;
                });

                dLambda += dl_accepted;

                // Guard NaN
                if (!std::isfinite(TrialStress.squaredNorm())) {
                    newton_ok = false;
                    break;
                }
            }

            if (!newton_ok) return false;

            // -------- Optional one-shot polish: return-to-surface --------
            if (INT_OPT_return_to_yield_surface[ASDP_TAG]) {
                const VoigtVector& n_corr = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m_corr = pf(dEps,        TrialStress, iv_storage, parameters_storage);
                const double       H_corr = yf.hardening(dEps, m_corr,   TrialStress, iv_storage, parameters_storage);
                const double       Phi_corr = yf(TrialStress, iv_storage, parameters_storage);

                const VoigtVector Em_corr = Eelastic * strain_like(m_corr);
                const double nEm_corr = strain_like(n_corr).dot(Em_corr);
                const double denom_corr = nEm_corr - H_corr;

                if (std::abs(Phi_corr) > tol_abs && std::abs(denom_corr) > MACHINE_EPSILON) {
                    const double dLambda_corr =  Phi_corr / denom_corr;
                    TrialStress         = TrialStress - dLambda_corr * Em_corr;
                    TrialPlastic_Strain = TrialPlastic_Strain + dLambda_corr * strain_like(m_corr);

                    iv_storage.apply([&](auto & iv){
                        auto h = iv.hardening_function(dEps, m_corr, TrialStress, parameters_storage);
                        iv.trial_value += dLambda_corr * h;
                    });
                    dLambda += dLambda_corr;
                }
            }

            // -------- Tangente algorítmica consistente --------
            {
                const VoigtVector& n = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                const VoigtVector& m = pf(dEps,        TrialStress, iv_storage, parameters_storage);
                const double       H = yf.hardening(dEps, m,        TrialStress, iv_storage, parameters_storage);

                const VoigtVector Em = Eelastic * strain_like(m);
                const double nEm = strain_like(n).dot(Em);
                const double denom = nEm - H; // ojo: este es el de la fórmula de C_alg

                if (std::abs(denom) > MACHINE_EPSILON) {
                    // (6x6) = E - (6x1)*(1x6)/denom
                    const auto row_nE = (strain_like(n).transpose() * Eelastic); // 1x6
                    Stiffness = Eelastic - (Em * row_nE) / denom;
                } else {
                    Stiffness = Eelastic; // fallback
                }
            }

            return true;
        };

        // ------------------ substepping automático ------------------
        VoigtVector dEps = depsilon;

        // intenta paso completo; si falla, corta a la mitad repetidamente (hasta 1/32)
        bool ok = false;
        for (int split = 0; split <= 5; ++split)   // 0..5 → 1, 2, 4, 8, 16, 32 subpasos
        {
            // resetear a commit antes de intentar este tamaño de paso
            TrialStress = sigma + Eelastic * dEps;   // predictor para flags de arriba
            TrialStrain = epsilon + dEps;
            TrialPlastic_Strain = CommitPlastic_Strain;
            iv_storage.revert_all();

            if (solve_increment(dEps)) { ok = true; break; }

            // reducir paso y reintentar
            dEps *= 0.5;
        }

        if (!ok) return -1;

        return 0;
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

        double den = strain_like(n).dot(Eelastic * strain_like(m)) - hardening;
        double dLambda = (den != 0) ? strain_like(n).dot(Eelastic * depsilon_elpl) / den : 0;
        // the plastic multiplier is never negative: a substep heading back inside
        // the surface (unloading from a yielded state) is elastic
        if (dLambda < 0)
            dLambda = 0;

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

    int Runge_Kutta_45_Error_Control_old(const VoigtVector & strain_incr)
    {
        // cout << "Runge_Kutta_45_Error_Control" << endl;

        int errorcode = -1;

        static VoigtVector depsilon;
        depsilon.setZero();
        depsilon = strain_incr;


        iv_storage.revert_all();


        dsigma.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();

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
                double tol_yf = DBL_OPT_f_absolute_tol[ASDP_TAG];
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
                dSigma1  = Eelastic * (dEPS - dLambda * strain_like(m));
                dEpsilonPl1 = dLambda * strain_like(m);
                iv_storage1.apply([&m, &dLambda, &prev_Sigma, &dEPS, this](auto & iv1)
                {
                    auto h = iv1.hardening_function(dEPS, m, prev_Sigma, parameters_storage);
                    iv1.trial_value = iv1.committed_value + 0.5 * dLambda * h;
                });
                next_Sigma =  prev_Sigma + 0.5 * dSigma1;

                //Delta 2
                std::tie(dLambda, m) = CalculateLambdaM(next_Sigma, dEPS, parameters_storage, iv_storage1);
                dSigma2  = Eelastic * (dEPS - dLambda * strain_like(m));
                dEpsilonPl2 = dLambda * strain_like(m);
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
                dSigma3  = Eelastic * (dEPS - dLambda * strain_like(m));
                dEpsilonPl3 = dLambda * strain_like(m);
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
                dSigma4  = Eelastic * (dEPS - dLambda * strain_like(m));
                dEpsilonPl4 = dLambda * strain_like(m);
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
                dSigma5  = Eelastic * (dEPS - dLambda * strain_like(m));
                dEpsilonPl5 = dLambda * strain_like(m);
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
                dSigma6  = Eelastic * (dEPS - dLambda * strain_like(m));
                dEpsilonPl6 = dLambda * strain_like(m);
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
                                                         strain_like(n_after_corrector).dot(Eelastic * strain_like(m_after_corrector)) - hardening_after_corrector
                                                     );
                    TrialStress = TrialStress - dLambda_after_corrector * Eelastic * strain_like(m_after_corrector);
                    TrialPlastic_Strain += dLambda_after_corrector * strain_like(m_after_corrector);

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
                double TOL = this->DBL_OPT_f_absolute_tol[ASDP_TAG];
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
                                                        strain_like(n_after_corrector).dot(Eelastic * strain_like(m_after_corrector)) - hardening_after_corrector
                                                    );

                    
                    VoigtVector TS = TrialStress - dL * Eelastic * strain_like(m_after_corrector);

                    double y1 = yf(TS, iv_storage, parameters_storage);

                    //Try to bracket solution
                    // cout << "   y0 = " << y0 << "  dL =" << dL << "   y1 =" << y1 << endl;
                    while( y1 > 0 && iter < NITER)
                    {
                        dL = dL * 1.1;
                        TS = TrialStress - dL * Eelastic * strain_like(m_after_corrector);
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

                        VoigtVector TS2 = TrialStress - dL_mid * Eelastic * strain_like(m_after_corrector);
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
                            TS2 = TrialStress - dL_mid * Eelastic * strain_like(m_after_corrector);
                            y_mid = yf(TS2, iv_storage, parameters_storage);
                            iter++;
                        }
                        dL = dL_mid;
                    }
                    

                    TrialStress = TrialStress - dL * Eelastic * strain_like(m_after_corrector);
                    TrialPlastic_Strain += dL * strain_like(m_after_corrector);

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

    // Modified Euler with Error Control - Following RK45 coding style
    int Modified_Euler_Error_Control(const VoigtVector & strain_incr)
    {
        // Modified Euler coefficients (Heun's method)
        constexpr double a21 = 1.0;  // For predictor step
        constexpr double b1 = 0.5;   // Weight for corrector (average)
        constexpr double b2 = 0.5;   // Weight for corrector (average)

        int errorcode = -1;

        static VoigtVector depsilon;
        depsilon.setZero();
        depsilon = strain_incr;

        iv_storage.revert_all();

        dsigma.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();

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
            return 0;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = DBL_OPT_f_absolute_tol[ASDP_TAG];
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = CommitStrain  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;
            double T = 0.0, dT = 1.0, dT_min = this->DBL_OPT_RK45_dT_min[ASDP_TAG], TolE = this->DBL_OPT_stress_absolute_tol[ASDP_TAG];
            
            // Adaptive step control parameters
            const double safety_factor = 0.9;
            const double min_scale = 0.25;
            const double max_scale = 4.0;
            const double beta = 0.04;  // PI controller parameter
            double previous_error = TolE;

            VoigtVector current_Sigma = TrialStress;
            iv_storage_t current_iv_storage = iv_storage;
            VoigtVector current_EpsilonPl = CommitPlastic_Strain;

            // Storage for k values - stress derivatives
            VoigtVector k1_sigma, k2_sigma;
            // Storage for k values - plastic strain derivatives  
            VoigtVector k1_pstrain, k2_pstrain;
            // Storage for intermediate iv storages
            iv_storage_t iv_k1 = iv_storage, iv_k2 = iv_storage;

            // Storage for predictor and corrector solutions
            VoigtVector predictor_sigma, corrector_sigma;
            VoigtVector predictor_pstrain, corrector_pstrain;
            iv_storage_t predictor_iv = iv_storage, corrector_iv = iv_storage;

            int niter = 0;
            double maxStepError = 0;
            int max_iterations = this->INT_OPT_RK45_niter_max[ASDP_TAG];
            
            while (T < 1.0)
            {
                niter++;
                
                double effective_dT = std::min(dT, 1.0 - T);
                VoigtVector dEPS = effective_dT * depsilon_elpl;
                VoigtVector m;
                double dLambda;

                // Update elasticity matrix for current state
                Eelastic = et(current_Sigma, parameters_storage);

                // PREDICTOR STEP (Forward Euler)
                // k1 = f(t, y)
                std::tie(dLambda, m) = CalculateLambdaM(current_Sigma, dEPS, parameters_storage, current_iv_storage);
                k1_sigma = Eelastic * (dEPS - dLambda * strain_like(m));
                k1_pstrain = dLambda * strain_like(m);
                iv_k1 = current_iv_storage;
                iv_k1.apply([&m, &dLambda, &current_Sigma, &dEPS, this](auto & iv1)
                {
                    auto h = iv1.hardening_function(dEPS, m, current_Sigma, parameters_storage);
                    iv1.trial_value = iv1.committed_value + dLambda * h;
                });

                // Predictor solution: y_pred = y_n + h*k1
                predictor_sigma = current_Sigma + k1_sigma;
                predictor_pstrain = current_EpsilonPl + k1_pstrain;
                predictor_iv = current_iv_storage;
                predictor_iv.apply([&iv_k1, &current_iv_storage](auto & pred_var)
                {
                    using VT = std::decay_t<decltype(pred_var)>;
                    const VT &iv1_var = iv_k1.template get<VT>();
                    const VT &current_var = current_iv_storage.template get<VT>();
                    auto dk1 = iv1_var.trial_value - current_var.committed_value;
                    pred_var.trial_value = current_var.committed_value + dk1;
                });

                // CORRECTOR STEP (Modified Euler)
                // k2 = f(t + h, y_pred)
                Eelastic = et(predictor_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(predictor_sigma, dEPS, parameters_storage, predictor_iv);
                k2_sigma = Eelastic * (dEPS - dLambda * strain_like(m));
                k2_pstrain = dLambda * strain_like(m);
                iv_k2 = current_iv_storage;
                iv_k2.apply([&m, &dLambda, &predictor_sigma, &dEPS, this](auto & iv2)
                {
                    auto h = iv2.hardening_function(dEPS, m, predictor_sigma, parameters_storage);
                    iv2.trial_value = iv2.committed_value + dLambda * h;
                });

                // Corrector solution: y_corr = y_n + h/2*(k1 + k2)
                corrector_sigma = current_Sigma + b1 * k1_sigma + b2 * k2_sigma;
                corrector_pstrain = current_EpsilonPl + b1 * k1_pstrain + b2 * k2_pstrain;
                corrector_iv = current_iv_storage;
                corrector_iv.apply([&iv_k1, &iv_k2, &current_iv_storage, b1, b2](auto & corr_var)
                {
                    using VT = std::decay_t<decltype(corr_var)>;
                    const VT &iv1_var = iv_k1.template get<VT>();
                    const VT &iv2_var = iv_k2.template get<VT>();
                    const VT &current_var = current_iv_storage.template get<VT>();
                    auto dk1 = iv1_var.trial_value - current_var.committed_value;
                    auto dk2 = iv2_var.trial_value - current_var.committed_value;
                    corr_var.trial_value = current_var.committed_value + b1 * dk1 + b2 * dk2;
                });

                // Error estimation: difference between predictor and corrector
                VoigtVector sigma_error = corrector_sigma - predictor_sigma;
                VoigtVector pstrain_error = corrector_pstrain - predictor_pstrain;
                
                double step_error = sigma_error.norm() + pstrain_error.norm();
                
                // Normalize error by solution magnitude
                double solution_norm = corrector_sigma.norm() + corrector_pstrain.norm();
                if (solution_norm > 0.1) {
                    step_error /= solution_norm;
                }

                // Check for NaN
                if (std::isnan(step_error) || std::isnan(corrector_sigma.norm()) || std::isnan(corrector_pstrain.norm()))
                {
                    cout << "ASDPlasticMaterial3D::Modified_Euler_Error_Control - NaN encountered, reducing step size" << endl;
                    dT *= 0.5;
                    if (dT < dT_min) {
                        cout << "ASDPlasticMaterial3D::Modified_Euler_Error_Control - Minimum step size reached with NaN" << endl;
                        return -1;
                    }
                    continue;
                }

                // Step size control with PI controller
                double error_ratio = TolE / std::max(step_error, 1e-15);
                double scale_factor = safety_factor * std::pow(error_ratio, 0.5) * std::pow(previous_error / step_error, beta);
                scale_factor = std::max(min_scale, std::min(max_scale, scale_factor));

                // Accept or reject step
                if (step_error <= TolE || effective_dT <= dT_min) {
                    // Accept step - use corrector solution
                    current_Sigma = corrector_sigma;
                    current_EpsilonPl = corrector_pstrain;
                    current_iv_storage = corrector_iv;
                    
                    T += effective_dT;
                    maxStepError = std::max(maxStepError, step_error);
                    previous_error = step_error;
                    
                    // Validate yield function drift
                    double yf_val = yf(current_Sigma, current_iv_storage, parameters_storage);
                    if (yf_val > 10 * DBL_OPT_f_absolute_tol[ASDP_TAG]) {
                        // cout << "Warning: Yield function drift detected: f = " << yf_val << endl;
                    }
                }

                // Update step size for next iteration
                double new_dT = scale_factor * effective_dT;
                dT = std::max(dT_min, std::min(new_dT, 1.0 - T));

                if (niter > max_iterations)
                {
                    cout << "ASDPlasticMaterial3D - tag = " << ASDP_TAG << " Modified Euler exceeded number of iterations. niter = " << niter << " niter_max = " << max_iterations << " T= " << T << " dT = " << dT << endl;
                    return -1;
                }
            }

            GLOBAL_INT_max_iter[ASDP_TAG] = std::max(GLOBAL_INT_max_iter[ASDP_TAG], niter);
            GLOBAL_DBL_max_error[ASDP_TAG] = std::max(GLOBAL_DBL_max_error[ASDP_TAG], maxStepError);

            TrialStress = current_Sigma;
            TrialPlastic_Strain = current_EpsilonPl;
            iv_storage = current_iv_storage;

            //Return to Yield
            if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 1)  // Return to yield in one step
            {
                double yf_val_after_corrector = yf(TrialStress, iv_storage, parameters_storage);
                if (yf_val_after_corrector > DBL_OPT_f_absolute_tol[ASDP_TAG]) {
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double denominator = strain_like(n_after_corrector).dot(Eelastic * strain_like(m_after_corrector)) - hardening_after_corrector;
                    
                    if (std::abs(denominator) > MACHINE_EPSILON) {
                        double dLambda_after_corrector = yf_val_after_corrector / denominator;
                        TrialStress = TrialStress - dLambda_after_corrector * Eelastic * strain_like(m_after_corrector);
                        TrialPlastic_Strain += dLambda_after_corrector * strain_like(m_after_corrector);
                        
                        // Update internal variables
                        iv_storage.apply([&m_after_corrector, &dLambda_after_corrector, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dLambda_after_corrector * h;
                        });
                    }
                }
            }
            else if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 2)  // Return to yield with bisection
            {
                double y0 = yf(TrialStress, iv_storage, parameters_storage);
                int iter = 0;
                double TOL = this->DBL_OPT_f_absolute_tol[ASDP_TAG];
                int NITER = this->INT_OPT_n_max_iterations[ASDP_TAG];
                
                if(y0 > TOL && iter < NITER)
                {
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double denominator = strain_like(n_after_corrector).dot(Eelastic * strain_like(m_after_corrector)) - hardening_after_corrector;
                    
                    if (std::abs(denominator) > MACHINE_EPSILON) {
                        double dL = y0 / denominator;
                        VoigtVector TS = TrialStress - dL * Eelastic * strain_like(m_after_corrector);
                        double y1 = yf(TS, iv_storage, parameters_storage);

                        // Try to bracket solution
                        while( y1 > 0 && iter < NITER)
                        {
                            dL = dL * 1.1;
                            TS = TrialStress - dL * Eelastic * strain_like(m_after_corrector);
                            y1 = yf(TS, iv_storage, parameters_storage);
                            iter++;
                        }

                        iter = 0;

                        // Once solution is bracketed, use bisection to get to YS
                        if (y1 < 0)
                        {
                            double dL_min = 0;
                            double dL_max = dL;
                            double dL_mid = dL / 2;

                            VoigtVector TS2 = TrialStress - dL_mid * Eelastic * strain_like(m_after_corrector);
                            double y_mid = yf(TS2, iv_storage, parameters_storage);
                            
                            while(std::abs(y_mid) > TOL && iter < NITER)
                            {
                                if (y_mid > 0) {
                                    dL_min = dL_mid;
                                } else {
                                    dL_max = dL_mid;
                                }
                                dL_mid = 0.5*(dL_min + dL_max);
                                TS2 = TrialStress - dL_mid * Eelastic * strain_like(m_after_corrector);
                                y_mid = yf(TS2, iv_storage, parameters_storage);
                                iter++;
                            }
                            dL = dL_mid;
                        }
                        
                        TrialStress = TrialStress - dL * Eelastic * strain_like(m_after_corrector);
                        TrialPlastic_Strain += dL * strain_like(m_after_corrector);
                        
                        // Update internal variables
                        iv_storage.apply([&m_after_corrector, &dL, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dL * h;
                        });
                    }
                }
            }

            // Final validation
            double norm_trial_stress = TrialStress.transpose() * TrialStress;
            if (norm_trial_stress != norm_trial_stress) //check for nan
            {
                cout << "ASDPlasticMaterial3D::Modified_Euler_Error_Control  Numeric error!\n";
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

                return -1;
            }

            ComputeTangentStiffness();
        }

        return 0;
    }

    // Complete corrected RK45 implementation with proper internal variable handling
    int Runge_Kutta_45_Error_Control(const VoigtVector & strain_incr)
    {
        int RK45_EXIT_FLAG = 0;

        // Dormand-Prince RK45 coefficients
        constexpr double a21 = 1.0/5.0;
        constexpr double a31 = 3.0/40.0, a32 = 9.0/40.0;
        constexpr double a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0;
        constexpr double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0;
        constexpr double a61 = 9017.0/3168.0, a62 = -355.0/33.0, a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0;
        
        // 5th order solution coefficients
        constexpr double b1 = 35.0/384.0, b2 = 0.0, b3 = 500.0/1113.0, b4 = 125.0/192.0, b5 = -2187.0/6784.0, b6 = 11.0/84.0;
        
        // 4th order solution coefficients for error estimation
        constexpr double bhat1 = 5179.0/57600.0, bhat2 = 0.0, bhat3 = 7571.0/16695.0, bhat4 = 393.0/640.0;
        constexpr double bhat5 = -92097.0/339200.0, bhat6 = 187.0/2100.0, bhat7 = 1.0/40.0;

        int errorcode = -1;

        static VoigtVector depsilon;
        depsilon.setZero();
        depsilon = strain_incr;

        // CRITICAL FIX: Initialize trial values from committed values at start
        iv_storage.apply([](auto & iv) {
            iv.trial_value = iv.committed_value;
        });

        dsigma.setZero();
        intersection_stress.setZero();
        intersection_strain.setZero();

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
            return 0;
        }
        else  //Plasticity
        {
            depsilon_elpl = depsilon;
            if (yf_val_start < 0)
            {
                double tol_yf = DBL_OPT_f_absolute_tol[ASDP_TAG];
                double intersection_factor = compute_yf_crossing( start_stress, end_stress, 0.0, 1.0, tol_yf );

                intersection_factor = intersection_factor < 0 ? 0 : intersection_factor;
                intersection_factor = intersection_factor > 1 ? 1 : intersection_factor;

                intersection_stress = start_stress * (1 - intersection_factor) + end_stress * intersection_factor;
                intersection_strain = CommitStrain  + depsilon * intersection_factor;
                depsilon_elpl = (1 - intersection_factor) * depsilon;
            }

            TrialStress = intersection_stress;
            double T = 0.0, dT = 1.0, dT_min = this->DBL_OPT_RK45_dT_min[ASDP_TAG], TolE = this->DBL_OPT_stress_absolute_tol[ASDP_TAG];
            
            // Adaptive step control parameters
            const double safety_factor = 0.9;
            const double min_scale = 0.2;
            const double max_scale = 5.0;
            const double beta = 0.04;  // PI controller parameter
            double previous_error = TolE;

            VoigtVector current_Sigma = TrialStress;
            iv_storage_t current_iv_storage = iv_storage;
            VoigtVector current_EpsilonPl = CommitPlastic_Strain;

            // Storage for k values - stress derivatives
            VoigtVector k1_sigma, k2_sigma, k3_sigma, k4_sigma, k5_sigma, k6_sigma;
            // Storage for k values - plastic strain derivatives  
            VoigtVector k1_pstrain, k2_pstrain, k3_pstrain, k4_pstrain, k5_pstrain, k6_pstrain;
            
            // CORRECTED: Storage for internal variable derivatives (not absolute states)
            iv_storage_t iv_k1_derivatives, iv_k2_derivatives, iv_k3_derivatives;
            iv_storage_t iv_k4_derivatives, iv_k5_derivatives, iv_k6_derivatives;

            int niter = 0;
            double maxStepError = 0;
            
            while (T < 1.0)
            {
                niter++;
                
                double effective_dT = std::min(dT, 1.0 - T);
                VoigtVector dEPS = effective_dT * depsilon_elpl;
                VoigtVector m;
                double dLambda;

                // Update elasticity matrix for current state
                Eelastic = et(current_Sigma, parameters_storage);

                // k1 = f(t, y) - compute derivatives at current state
                std::tie(dLambda, m) = CalculateLambdaM(current_Sigma, dEPS, parameters_storage, current_iv_storage);
                k1_sigma = Eelastic * (dEPS - dLambda * strain_like(m));
                k1_pstrain = dLambda * strain_like(m);
                
                // CORRECTED: Store pure derivatives for internal variables
                iv_k1_derivatives = current_iv_storage;  // Copy structure
                iv_k1_derivatives.apply([&m, &dLambda, &current_Sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, current_Sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;  // Pure derivative
                });

                // k2 = f(t + c2*h, y + h*(a21*k1))
                VoigtVector y2_sigma = current_Sigma + a21 * k1_sigma;
                VoigtVector y2_pstrain = current_EpsilonPl + a21 * k1_pstrain;
                
                // CORRECTED: Properly propagate internal variable state using RK formula
                iv_storage_t iv2_state = current_iv_storage;
                iv2_state.apply([&iv_k1_derivatives, &current_iv_storage, a21](auto & iv2_var)
                {
                    using VT = std::decay_t<decltype(iv2_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    iv2_var.trial_value = current_var.trial_value + a21 * k1_deriv.trial_value;
                });
                
                Eelastic = et(y2_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y2_sigma, dEPS, parameters_storage, iv2_state);
                k2_sigma = Eelastic * (dEPS - dLambda * strain_like(m));
                k2_pstrain = dLambda * strain_like(m);
                
                iv_k2_derivatives = current_iv_storage;
                iv_k2_derivatives.apply([&m, &dLambda, &y2_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y2_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // k3 = f(t + c3*h, y + h*(a31*k1 + a32*k2))
                VoigtVector y3_sigma = current_Sigma + a31 * k1_sigma + a32 * k2_sigma;
                VoigtVector y3_pstrain = current_EpsilonPl + a31 * k1_pstrain + a32 * k2_pstrain;
                
                iv_storage_t iv3_state = current_iv_storage;
                iv3_state.apply([&iv_k1_derivatives, &iv_k2_derivatives, &current_iv_storage, a31, a32](auto & iv3_var)
                {
                    using VT = std::decay_t<decltype(iv3_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    iv3_var.trial_value = current_var.trial_value + a31 * k1_deriv.trial_value + a32 * k2_deriv.trial_value;
                });
                
                Eelastic = et(y3_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y3_sigma, dEPS, parameters_storage, iv3_state);
                k3_sigma = Eelastic * (dEPS - dLambda * strain_like(m));
                k3_pstrain = dLambda * strain_like(m);
                
                iv_k3_derivatives = current_iv_storage;
                iv_k3_derivatives.apply([&m, &dLambda, &y3_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y3_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // k4 = f(t + c4*h, y + h*(a41*k1 + a42*k2 + a43*k3))
                VoigtVector y4_sigma = current_Sigma + a41 * k1_sigma + a42 * k2_sigma + a43 * k3_sigma;
                VoigtVector y4_pstrain = current_EpsilonPl + a41 * k1_pstrain + a42 * k2_pstrain + a43 * k3_pstrain;
                
                iv_storage_t iv4_state = current_iv_storage;
                iv4_state.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &current_iv_storage, 
                               a41, a42, a43](auto & iv4_var)
                {
                    using VT = std::decay_t<decltype(iv4_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    iv4_var.trial_value = current_var.trial_value + a41 * k1_deriv.trial_value + 
                                         a42 * k2_deriv.trial_value + a43 * k3_deriv.trial_value;
                });
                
                Eelastic = et(y4_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y4_sigma, dEPS, parameters_storage, iv4_state);
                k4_sigma = Eelastic * (dEPS - dLambda * strain_like(m));
                k4_pstrain = dLambda * strain_like(m);
                
                iv_k4_derivatives = current_iv_storage;
                iv_k4_derivatives.apply([&m, &dLambda, &y4_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y4_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // k5 = f(t + c5*h, y + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4))
                VoigtVector y5_sigma = current_Sigma + a51 * k1_sigma + a52 * k2_sigma + a53 * k3_sigma + a54 * k4_sigma;
                VoigtVector y5_pstrain = current_EpsilonPl + a51 * k1_pstrain + a52 * k2_pstrain + 
                                        a53 * k3_pstrain + a54 * k4_pstrain;
                
                iv_storage_t iv5_state = current_iv_storage;
                iv5_state.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &iv_k4_derivatives, 
                               &current_iv_storage, a51, a52, a53, a54](auto & iv5_var)
                {
                    using VT = std::decay_t<decltype(iv5_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    const VT &k4_deriv = iv_k4_derivatives.template get<VT>();
                    iv5_var.trial_value = current_var.trial_value + a51 * k1_deriv.trial_value + 
                                         a52 * k2_deriv.trial_value + a53 * k3_deriv.trial_value + 
                                         a54 * k4_deriv.trial_value;
                });
                
                Eelastic = et(y5_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y5_sigma, dEPS, parameters_storage, iv5_state);
                k5_sigma = Eelastic * (dEPS - dLambda * strain_like(m));
                k5_pstrain = dLambda * strain_like(m);
                
                iv_k5_derivatives = current_iv_storage;
                iv_k5_derivatives.apply([&m, &dLambda, &y5_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y5_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // k6 = f(t + h, y + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5))
                VoigtVector y6_sigma = current_Sigma + a61 * k1_sigma + a62 * k2_sigma + a63 * k3_sigma + 
                                      a64 * k4_sigma + a65 * k5_sigma;
                VoigtVector y6_pstrain = current_EpsilonPl + a61 * k1_pstrain + a62 * k2_pstrain + 
                                        a63 * k3_pstrain + a64 * k4_pstrain + a65 * k5_pstrain;
                
                iv_storage_t iv6_state = current_iv_storage;
                iv6_state.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &iv_k4_derivatives, 
                               &iv_k5_derivatives, &current_iv_storage, a61, a62, a63, a64, a65](auto & iv6_var)
                {
                    using VT = std::decay_t<decltype(iv6_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    const VT &k4_deriv = iv_k4_derivatives.template get<VT>();
                    const VT &k5_deriv = iv_k5_derivatives.template get<VT>();
                    iv6_var.trial_value = current_var.trial_value + a61 * k1_deriv.trial_value + 
                                         a62 * k2_deriv.trial_value + a63 * k3_deriv.trial_value + 
                                         a64 * k4_deriv.trial_value + a65 * k5_deriv.trial_value;
                });
                
                Eelastic = et(y6_sigma, parameters_storage);
                std::tie(dLambda, m) = CalculateLambdaM(y6_sigma, dEPS, parameters_storage, iv6_state);
                k6_sigma = Eelastic * (dEPS - dLambda * strain_like(m));
                k6_pstrain = dLambda * strain_like(m);
                
                iv_k6_derivatives = current_iv_storage;
                iv_k6_derivatives.apply([&m, &dLambda, &y6_sigma, &dEPS, this](auto & iv_deriv)
                {
                    auto h = iv_deriv.hardening_function(dEPS, m, y6_sigma, parameters_storage);
                    iv_deriv.trial_value = dLambda * h;
                });

                // 5th order solution
                VoigtVector next_Sigma_5th = current_Sigma + (b1 * k1_sigma + b2 * k2_sigma + b3 * k3_sigma + 
                                                             b4 * k4_sigma + b5 * k5_sigma + b6 * k6_sigma);
                VoigtVector next_EpsilonPl_5th = current_EpsilonPl + (b1 * k1_pstrain + b2 * k2_pstrain + 
                                                                     b3 * k3_pstrain + b4 * k4_pstrain + 
                                                                     b5 * k5_pstrain + b6 * k6_pstrain);

                // CORRECTED: 5th order solution for internal variables including ALL terms
                iv_storage_t next_iv_5th = current_iv_storage;
                next_iv_5th.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &iv_k4_derivatives, 
                                 &iv_k5_derivatives, &iv_k6_derivatives, &current_iv_storage, 
                                 b1, b2, b3, b4, b5, b6](auto & next_var)
                {
                    using VT = std::decay_t<decltype(next_var)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    const VT &k4_deriv = iv_k4_derivatives.template get<VT>();
                    const VT &k5_deriv = iv_k5_derivatives.template get<VT>();
                    const VT &k6_deriv = iv_k6_derivatives.template get<VT>();
                    
                    next_var.trial_value = current_var.trial_value + 
                        b1 * k1_deriv.trial_value + b2 * k2_deriv.trial_value + b3 * k3_deriv.trial_value + 
                        b4 * k4_deriv.trial_value + b5 * k5_deriv.trial_value + b6 * k6_deriv.trial_value;
                });

                // 4th order solution for error estimation
                VoigtVector next_Sigma_4th = current_Sigma + (bhat1 * k1_sigma + bhat2 * k2_sigma + bhat3 * k3_sigma + 
                                                             bhat4 * k4_sigma + bhat5 * k5_sigma + bhat6 * k6_sigma);
                VoigtVector next_EpsilonPl_4th = current_EpsilonPl + (bhat1 * k1_pstrain + bhat2 * k2_pstrain + 
                                                                     bhat3 * k3_pstrain + bhat4 * k4_pstrain + 
                                                                     bhat5 * k5_pstrain + bhat6 * k6_pstrain);

                // 4th order solution for internal variables (for error estimation)
                iv_storage_t next_iv_4th = current_iv_storage;
                next_iv_4th.apply([&iv_k1_derivatives, &iv_k2_derivatives, &iv_k3_derivatives, &iv_k4_derivatives, 
                                 &iv_k5_derivatives, &iv_k6_derivatives, &current_iv_storage, 
                                 bhat1, bhat2, bhat3, bhat4, bhat5, bhat6](auto & next_var_4th)
                {
                    using VT = std::decay_t<decltype(next_var_4th)>;
                    const VT &current_var = current_iv_storage.template get<VT>();
                    const VT &k1_deriv = iv_k1_derivatives.template get<VT>();
                    const VT &k2_deriv = iv_k2_derivatives.template get<VT>();
                    const VT &k3_deriv = iv_k3_derivatives.template get<VT>();
                    const VT &k4_deriv = iv_k4_derivatives.template get<VT>();
                    const VT &k5_deriv = iv_k5_derivatives.template get<VT>();
                    const VT &k6_deriv = iv_k6_derivatives.template get<VT>();
                    
                    next_var_4th.trial_value = current_var.trial_value + 
                        bhat1 * k1_deriv.trial_value + bhat2 * k2_deriv.trial_value + bhat3 * k3_deriv.trial_value + 
                        bhat4 * k4_deriv.trial_value + bhat5 * k5_deriv.trial_value + bhat6 * k6_deriv.trial_value;
                });

                // Error estimation including internal variables
                VoigtVector sigma_error = next_Sigma_5th - next_Sigma_4th;
                VoigtVector pstrain_error = next_EpsilonPl_5th - next_EpsilonPl_4th;
                
                // ADDED: Include internal variable errors in step control
                double iv_error = 0.0;
                next_iv_5th.apply([&next_iv_4th, &iv_error](const auto & iv_5th)
                {
                    using VT = std::decay_t<decltype(iv_5th)>;
                    const VT &iv_4th = next_iv_4th.template get<VT>();
                    auto diff = iv_5th.trial_value - iv_4th.trial_value;
                    iv_error += diff.norm();
                });
                
                double step_error = sigma_error.norm() + pstrain_error.norm() + iv_error;
                
                // Normalize error by solution magnitude
                double solution_norm = next_Sigma_5th.norm() + next_EpsilonPl_5th.norm();
                double iv_norm = 0.0;
                next_iv_5th.apply([&iv_norm](const auto & iv)
                {
                    iv_norm += iv.trial_value.norm();
                });
                solution_norm += iv_norm;
                
                if (solution_norm > 0.1) {
                    step_error /= solution_norm;
                }

                // Check for NaN
                if (std::isnan(step_error) || std::isnan(next_Sigma_5th.norm()) || std::isnan(next_EpsilonPl_5th.norm()))
                {
                    cout << "ASDPlasticMaterial3D::RK45 - NaN encountered, reducing step size" << endl;
                    dT *= 0.5;
                    if (dT < dT_min) {
                        cout << "ASDPlasticMaterial3D::RK45 - Minimum step size reached with NaN" << endl;
                        return -1;
                    }
                    continue;
                }

                // Step size control with PI controller
                double error_ratio = TolE / std::max(step_error, 1e-15);
                double scale_factor = safety_factor * std::pow(error_ratio, 0.2) * std::pow(previous_error / step_error, beta);
                scale_factor = std::max(min_scale, std::min(max_scale, scale_factor));

                // Accept or reject step
                if (step_error <= TolE || effective_dT <= dT_min) {
                    // Accept step - use 5th order solution
                    current_Sigma = next_Sigma_5th;
                    current_EpsilonPl = next_EpsilonPl_5th;
                    current_iv_storage = next_iv_5th;
                    
                    T += effective_dT;
                    maxStepError = std::max(maxStepError, step_error);
                    previous_error = step_error;
                    
                    // Validate yield function drift
                    double yf_val = yf(current_Sigma, current_iv_storage, parameters_storage);
                    if (yf_val > 10 * TolE) {
                        // Optional: Add drift correction here
                        // cout << "Warning: Yield function drift detected: f = " << yf_val << endl;
                    }
                }

                // Update step size for next iteration
                double new_dT = scale_factor * effective_dT;
                dT = std::max(dT_min, std::min(new_dT, 1.0 - T));

                if (niter > this->INT_OPT_RK45_niter_max[ASDP_TAG])
                {
                    cout << "ASDPlasticMaterial3D - tag = " << ASDP_TAG << " RK45 exceeded number of iterations. niter = " << niter << " niter_max = " <<this->INT_OPT_RK45_niter_max[ASDP_TAG] << " T= " << T << " dT = " << dT << endl;
                    RK45_EXIT_FLAG = -1;
                    break;
                }
            }

            GLOBAL_INT_max_iter[ASDP_TAG] = std::max(GLOBAL_INT_max_iter[ASDP_TAG], niter);
            GLOBAL_DBL_max_error[ASDP_TAG] = std::max(GLOBAL_DBL_max_error[ASDP_TAG], maxStepError);

            TrialStress = current_Sigma;
            TrialPlastic_Strain = current_EpsilonPl;
            // iv_storage = current_iv_storage;
            iv_storage.updateTrialValueFromOther(current_iv_storage);


            //Return to Yield Surface
            if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 1)  // Return to yield in one step
            {
                double yf_val_after_corrector = yf(TrialStress, iv_storage, parameters_storage);
                if (yf_val_after_corrector > DBL_OPT_f_absolute_tol[ASDP_TAG]) {
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double denominator = strain_like(n_after_corrector).dot(Eelastic * strain_like(m_after_corrector)) - hardening_after_corrector;
                    
                    if (std::abs(denominator) > MACHINE_EPSILON) {
                        double dLambda_after_corrector = yf_val_after_corrector / denominator;
                        TrialStress = TrialStress - dLambda_after_corrector * Eelastic * strain_like(m_after_corrector);
                        TrialPlastic_Strain += dLambda_after_corrector * strain_like(m_after_corrector);
                        
                        // Update internal variables
                        iv_storage.apply([&m_after_corrector, &dLambda_after_corrector, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dLambda_after_corrector * h;
                        });
                    }
                }
            }
            else if (INT_OPT_return_to_yield_surface[ASDP_TAG] == 2)  // Return to yield with bisection
            {
                double y0 = yf(TrialStress, iv_storage, parameters_storage);
                int iter = 0;
                double TOL = this->DBL_OPT_f_absolute_tol[ASDP_TAG];
                int NITER = this->INT_OPT_n_max_iterations[ASDP_TAG];
                
                if(y0 > TOL && iter < NITER)
                {
                    const VoigtVector& n_after_corrector = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
                    const VoigtVector& m_after_corrector = pf(depsilon_elpl, TrialStress, iv_storage, parameters_storage);
                    double hardening_after_corrector = yf.hardening( depsilon_elpl, m_after_corrector,  TrialStress, iv_storage, parameters_storage);
                    double denominator = strain_like(n_after_corrector).dot(Eelastic * strain_like(m_after_corrector)) - hardening_after_corrector;
                    
                    if (std::abs(denominator) > MACHINE_EPSILON) {
                        double dL = y0 / denominator;
                        VoigtVector TS = TrialStress - dL * Eelastic * strain_like(m_after_corrector);
                        
                        // Create temporary iv storage for testing
                        iv_storage_t temp_iv = iv_storage;
                        temp_iv.apply([&m_after_corrector, &dL, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dL * h;
                        });
                        
                        double y1 = yf(TS, temp_iv, parameters_storage);

                        // Try to bracket solution
                        while( y1 > 0 && iter < NITER)
                        {
                            dL = dL * 1.1;
                            TS = TrialStress - dL * Eelastic * strain_like(m_after_corrector);
                            
                            temp_iv = iv_storage;
                            temp_iv.apply([&m_after_corrector, &dL, this](auto& iv) {
                                auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                                iv.trial_value += dL * h;
                            });
                            
                            y1 = yf(TS, temp_iv, parameters_storage);
                            iter++;
                        }

                        iter = 0;

                        // Once solution is bracketed, use bisection to get to YS
                        if (y1 < 0)
                        {
                            double dL_min = 0;
                            double dL_max = dL;
                            double dL_mid = dL / 2;

                            VoigtVector TS2 = TrialStress - dL_mid * Eelastic * strain_like(m_after_corrector);
                            
                            iv_storage_t mid_iv = iv_storage;
                            mid_iv.apply([&m_after_corrector, &dL_mid, this](auto& iv) {
                                auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                                iv.trial_value += dL_mid * h;
                            });
                            
                            double y_mid = yf(TS2, mid_iv, parameters_storage);
                            
                            while(std::abs(y_mid) > TOL && iter < NITER)
                            {
                                if (y_mid > 0) {
                                    dL_min = dL_mid;
                                } else {
                                    dL_max = dL_mid;
                                }
                                dL_mid = 0.5*(dL_min + dL_max);
                                TS2 = TrialStress - dL_mid * Eelastic * strain_like(m_after_corrector);
                                
                                mid_iv = iv_storage;
                                mid_iv.apply([&m_after_corrector, &dL_mid, this](auto& iv) {
                                    auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                                    iv.trial_value += dL_mid * h;
                                });
                                
                                y_mid = yf(TS2, mid_iv, parameters_storage);
                                iter++;
                            }
                            dL = dL_mid;
                        }
                        
                        TrialStress = TrialStress - dL * Eelastic * strain_like(m_after_corrector);
                        TrialPlastic_Strain += dL * strain_like(m_after_corrector);
                        
                        // Update internal variables with final correction
                        iv_storage.apply([&m_after_corrector, &dL, this](auto& iv) {
                            auto h = iv.hardening_function(depsilon_elpl, m_after_corrector, TrialStress, parameters_storage);
                            iv.trial_value += dL * h;
                        });
                    }
                }
            }

            // Final validation
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

                return -1;
            }

            // ADDED: Final consistency check for internal variables
            #ifdef DEBUG_INTERNAL_VARIABLES
            cout << "=== Final Internal Variable State ===" << endl;
            iv_storage.apply([](const auto & iv) {
                cout << iv.getName() << ": trial=" << iv.trial_value.transpose() 
                     << ", committed=" << iv.committed_value.transpose() << endl;
            });
            
            double final_yf = yf(TrialStress, iv_storage, parameters_storage);
            cout << "Final yield function value: " << final_yf << endl;
            #endif

            ComputeTangentStiffness();
        }

        return RK45_EXIT_FLAG;
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

    std::string current_parameter_name; // Stores the most recent parameter name from setParameter

protected:

    static std::map<int, ASDPlasticMaterial3D_Constitutive_Integration_Method> INT_OPT_constitutive_integration_method;     //
    static std::map<int, ASDPlasticMaterial3D_Tangent_Operator_Type> INT_OPT_tangent_operator_type;     //
    static std::map<int, double> DBL_OPT_f_absolute_tol;
    static std::map<int, double> DBL_OPT_stress_absolute_tol;
    static std::map<int, int> INT_OPT_n_max_iterations;
    static std::map<int, int> INT_OPT_return_to_yield_surface;
    static std::map<int, double> DBL_OPT_RK45_dT_min;
    static std::map<int, int> INT_OPT_RK45_niter_max;
    static std::map<int, ASDPlasticMaterial3D_Derivative_Type> INT_OPT_derivative_type;
    static std::map<int, double> DBL_OPT_fd_perturbation;
    static std::map<int, int> FBE_fallbacks;
    static std::map<int, int> FBE_failures;

    static std::map<int, int> GLOBAL_INT_max_iter; 
    static std::map<int, double> GLOBAL_DBL_max_error; 

    bool first_step;
    bool stress_set_externally;
    bool skip_BE_tangent = false;   // set by Full_Backward_Euler around its warm start
    int be_last_iter = 0;           // iterations of the last Backward_Euler loop, for the Full_Backward_Euler failure report

    static VoigtVector dsigma;
    static VoigtVector depsilon_elpl;    //Elastoplastic strain increment : For a strain increment that causes first yield, the step is divided into an elastic one (until yield) and an elastoplastic one.
    static VoigtVector intersection_stress;
    static VoigtVector intersection_strain;
    // per instance: each Gauss point has its own tangent (a static one was
    // shared by every instance of the same combination)
    VoigtMatrix Stiffness;


};

template < class E, class Y, class P, int tag>
std::map<int, ASDPlasticMaterial3D_Constitutive_Integration_Method> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_constitutive_integration_method;
template < class E, class Y, class P, int tag>
std::map<int, ASDPlasticMaterial3D_Tangent_Operator_Type> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_tangent_operator_type;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::DBL_OPT_f_absolute_tol;
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
std::map<int, ASDPlasticMaterial3D_Derivative_Type> ASDPlasticMaterial3D< E,  Y,  P,  tag>::INT_OPT_derivative_type;
template < class E, class Y, class P, int tag>
std::map<int, double> ASDPlasticMaterial3D< E,  Y,  P,  tag>::DBL_OPT_fd_perturbation;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::FBE_fallbacks;
template < class E, class Y, class P, int tag>
std::map<int, int> ASDPlasticMaterial3D< E,  Y,  P,  tag>::FBE_failures;

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



#endif
