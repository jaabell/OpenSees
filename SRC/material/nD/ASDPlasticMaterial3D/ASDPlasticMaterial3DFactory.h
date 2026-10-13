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


#ifndef ASDPlasticMaterial3DFactory_h
#define ASDPlasticMaterial3DFactory_h

// Shared by OPS_AllASDPlasticMaterial3Ds.cpp and the generated ASD_material_definitions_*.cpp
// compilation units: the templates that create and populate one ASDPlasticMaterial3D combination.

#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <list>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

#include "AllASDPlasticMaterial3Ds.h"

template <typename T>
void populate_ASDPlasticMaterial3D(T* instance);

typedef std::tuple<std::string, std::string, std::string, std::string> model_spec_t;







template<typename EL, typename YF, typename PF>
NDMaterial* createASDPlasticMaterial3D(int instance_tag, 
        const char* yf_type, const char* pf_type, const char* el_type, const char* iv_type, 
        std::list<NDMaterial*> &instance_pointers, std::list<model_spec_t> &available_models) 
{
    auto instance = new ASDPlasticMaterial3D<EL, YF, PF, ND_TAG_ASDPlasticMaterial3D>(instance_tag);


    available_models.push_back(std::make_tuple(
        instance->getYFName(),
        instance->getPFName(),
        instance->getELName(),
        instance->getIVName())
    );

    if(
        std::strcmp(yf_type,instance->getYFName().c_str())==0 &&
        std::strcmp(pf_type,instance->getPFName().c_str())==0 &&
        std::strcmp(el_type,instance->getELName().c_str())==0 &&
        std::strcmp(iv_type,instance->getIVName().c_str())==0 
        )
    {
        populate_ASDPlasticMaterial3D(instance);
        cout << "\n\nPrinting material info\n";
        instance->Print(opserr);
        cout << "\n\nDone creating ASDPlasticMaterial3D \n\n\n";
        instance_pointers.push_back(static_cast<NDMaterial*>(instance));
        return static_cast<NDMaterial*>(instance);
    }  
    else
    {
        delete instance;
        instance_pointers.push_back(nullptr);        
        return nullptr;
    }
}









template <typename T>
void populate_ASDPlasticMaterial3D(T* instance)
{

    int get_one_value = 1;

    cout << "\n\nDefined internal variables: \n";
    auto iv_names = instance->getInternalVariablesNames();
    for_each_in_tuple(iv_names, [instance](auto & name)
    {
        std::cout << "   "  << name << " size = " << instance->getInternalVariableSizeByName(name) << std::endl;
    });

    cout << "\n\nDefined model parameters\n";
    auto parameter_names = instance->getParameterNames();
    for_each_in_tuple(parameter_names, [](auto & name)
    {
        std::cout << "   "  <<  name << std::endl;
    });

    // Default integration options
    int method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler;
    int tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Secant;
    double f_absolute_tol = 1e-6; 
    double stress_absolute_tol = 1e-6; 
    int n_max_iterations = 100;
    int return_to_yield_surface = 1;
    double rk45_dT_min = 1e-2;
    int rk45_niter_max = 110;

    // Loop over input arguments
    while (OPS_GetNumRemainingInputArgs() > 0) {

        //Current command
        const char *cmd = OPS_GetString();


        // Specifying internal variables
        if (std::strcmp(cmd, "Begin_Internal_Variables") == 0)
        {
            cout << "\n\nReading internal variables from input\n";
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double iv_values[6];
                const char *iv_name = OPS_GetString();

                if (std::strcmp(iv_name, "End_Internal_Variables") == 0)
                {
                    cout << "\n\n Done reading internal variables from input\n";
                    break;
                }

                int iv_size = instance->getInternalVariableSizeByName(iv_name);
                OPS_GetDouble(&iv_size, iv_values);
                cout << iv_name << " = ";
                for (int i = 0; i < iv_size; ++i)
                {
                    cout << iv_values[i] << " ";
                }
                cout << endl;
                instance->setInternalVariableByName(iv_name, iv_size, &iv_values[0]);
            }

        }


        // Specifying model parameters
        if (std::strcmp(cmd, "Begin_Model_Parameters") == 0)
        {
            cout << "\n\nReading parameters from input\n";
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double param_value;
                const char *param_name = OPS_GetString();
                if (std::strcmp(param_name, "End_Model_Parameters") == 0)
                {
                    cout << "\n\n Done reading parameters from input\n";
                    break;
                }

                OPS_GetDouble(&get_one_value, &param_value);
                cout << param_name << " = " << param_value << endl;
                instance->setParameterByName(param_name, param_value);
            }
        }


        // set_constitutive_integration_method(int method, double f_absolute_tol, double stress_absolute_tol, int n_max_iterations)
        if (std::strcmp(cmd, "Begin_Integration_Options") == 0)
        {
            cout << "\n\nReading Integration Options\n";
            while (OPS_GetNumRemainingInputArgs() > 0) {
                const char *param_name = OPS_GetString();
                if (std::strcmp(param_name, "End_Integration_Options") == 0)
                {
                    cout << "\n\nDone reading Integration Options\n";
                    break;
                }

                if (std::strcmp(param_name, "f_absolute_tol") == 0)
                {
                    OPS_GetDouble(&get_one_value, &f_absolute_tol);
                    cout << "   Setting f_absolute_tol = " << f_absolute_tol << endl;
                }

                if (std::strcmp(param_name, "stress_absolute_tol") == 0)
                {
                    OPS_GetDouble(&get_one_value, &stress_absolute_tol);
                    cout << "   Setting stress_absolute_tol = " << stress_absolute_tol << endl;
                }

                if (std::strcmp(param_name, "n_max_iterations") == 0)
                {
                    OPS_GetInt(&get_one_value, &n_max_iterations);
                    cout << "   Setting n_max_iterations = " << n_max_iterations << endl;
                }
                if (std::strcmp(param_name, "rk45_dT_min") == 0)
                {
                    OPS_GetDouble(&get_one_value, &rk45_dT_min);
                    cout << "   Setting rk45_dT_min = " << rk45_dT_min << endl;
                }

                if (std::strcmp(param_name, "rk45_niter_max") == 0)
                {
                    OPS_GetInt(&get_one_value, &rk45_niter_max);
                    cout << "   Setting rk45_niter_max = " << rk45_niter_max << endl;
                }
                if (std::strcmp(param_name, "return_to_yield_surface") == 0)
                {
                    // OPS_GetInt(&get_one_value, &return_to_yield_surface);
                    const char *method_name = OPS_GetString();
                    if (std::strcmp(method_name, "Disabled") == 0)
                        return_to_yield_surface = 0;
                    else if (std::strcmp(method_name, "One_Step_Return") == 0)
                        return_to_yield_surface = 1;
                    else if (std::strcmp(method_name, "Iterative_Return") == 0)
                        return_to_yield_surface = 2;
                    else
                    	return_to_yield_surface = 1;

                    cout << "   Setting return_to_yield_surface = " << return_to_yield_surface << endl;
                }

                if (std::strcmp(param_name, "integration_method") == 0)
                {
                    const char *method_name = OPS_GetString();
                    if (std::strcmp(method_name, "Forward_Euler") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler;
                    else if (std::strcmp(method_name, "Forward_Euler_Subincrement") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Forward_Euler_Subincrement;
                    else if (std::strcmp(method_name, "Modified_Euler_Error_Control") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Modified_Euler_Error_Control;
                    else if (std::strcmp(method_name, "Runge_Kutta_45_Error_Control") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control;
                    else if (std::strcmp(method_name, "Runge_Kutta_45_Error_Control_old") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Runge_Kutta_45_Error_Control_old;
                    else if (std::strcmp(method_name, "Backward_Euler") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler;                    
                    else if (std::strcmp(method_name, "Backward_Euler_LineSearch") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_LineSearch;
                    else
                    {
                        cout << "\n\nWARNING! Unrecognised ASDPlasticMaterial3D_Constitutive_Integration_Method name " << method_name << endl;
                        cout << "Available methods: Backward_Euler, Backward_Euler_LineSearch, Forward_Euler, Forward_Euler_Subincrement, Modified_Euler_Error_Control, Runge_Kutta_45_Error_Control\n\n" << endl;
                        cout << "Defaulting to Modified_Euler_Error_Control\n\n" << endl;
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Modified_Euler_Error_Control;
                    }
                    cout << "   Setting integration method = " << method_name << " method_int = " << method << endl;
                }

                if (std::strcmp(param_name, "tangent_type") == 0)
                {
                    const char *tangent_type_name = OPS_GetString();
                    if (std::strcmp(tangent_type_name, "Elastic") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Elastic;
                    else if (std::strcmp(tangent_type_name, "Continuum") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Continuum;
                    else if (std::strcmp(tangent_type_name, "Secant") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Secant;
                    else if (std::strcmp(tangent_type_name, "Numerical_Algorithmic_FirstOrder") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_FirstOrder;
                    else if (std::strcmp(tangent_type_name, "Numerical_Algorithmic_SecondOrder") == 0)
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Numerical_Algorithmic_SecondOrder;
                    else
                    {
                        cout << "WARNING! Unrecognised ASDPlasticMaterial3D_Tangent_Operator_Type name " << tangent_type_name << endl;
                        cout << "Defaulting to Elastic" << endl;
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Elastic;
                    }
                    cout << "   Setting tangent type = " << tangent_type_name << " tangent int = " << tangent << endl;
                    
                }

            }
        }
    }

    instance->set_constitutive_integration_method(method, tangent, f_absolute_tol, stress_absolute_tol, n_max_iterations, return_to_yield_surface, rk45_niter_max, rk45_dT_min);
}

#endif // ASDPlasticMaterial3DFactory_h
