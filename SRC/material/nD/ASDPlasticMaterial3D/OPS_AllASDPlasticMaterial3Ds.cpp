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

#ifdef _EIGEN3

#include "ASDPlasticMaterial3DFactory.h"
#include "ASD_material_definitions.h"



NDMaterial*  ASDPlasticMaterial3DFactory(int tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type, std::list<model_spec_t> &available_models);



void print_usage(void)
{
    opserr <<
       "SYNTAX:\n"
       "nDMaterial ASDPlasticMaterial3D $tag \\ \n"
       "    YF_type \\ \n"
       "    PF_type \\ \n"
       "    EL_type \\ \n"
       "    IV_type \\ \n"
       "Begin_Internal_Variables \\ \n"
       "    IV1 (initial value(s)) \\ \n"
       "    IV2 (initial value(s)) \\ \n"
       "    ....\n"
       "End_Internal_Variables \\ \n"
       "Begin_Model_Parameters \\ \n"
       "    PARAM1 (value) \\\n"
       "    PARAM2 (value) \\\n"
       "    ....\n"
       "End_Model_Parameters \\ \n"
       "Begin_Integration_Options \\ \n"
       "    f_absolute_tol (double value)\\ \n"
       "    stress_absolute_tol (double value)\\ \n"
       "    n_max_iterations (int value)\\ \n"
       "    rk45_dT_min (double value)\\ \n"
       "    rk45_niter_max (int value)\\ \n"
       "    return_to_yield_surface (string) : Disabled | One_Step_Return (default) | Iterative_Return\\ \n"
       "    integration_method (string) : Backward_Euler (default) | Backward_Euler_LineSearch | Forward_Euler |\n"
       "        Forward_Euler_Subincrement | Modified_Euler_Error_Control | Runge_Kutta_45_Error_Control\\ \n"
       "    tangent_type (string) : Secant (default) | Elastic | Continuum |\n"
       "        Numerical_Algorithmic_FirstOrder | Numerical_Algorithmic_SecondOrder\\ \n"
       "End_Integration_Options \\ \n"
       "\n"
       "To list the available YF/PF/EL/IV combinations (optionally for one yield function):\n"
       "nDMaterial ASDPlasticMaterial3D $tag list <YF_type>\n"
       "\n";
}

void *OPS_AllASDPlasticMaterial3Ds(void)
{
    // some kudos
    static bool first_done = false;
    if (!first_done) {
        opserr << "\n\nUsing ASDPlasticMaterial3D - Developed by: Jose Abell (UANDES), Massimo Petracca and Guido Camata (ASDEA Software Technology)\n\n";
        first_done = true;
    }

    // check arguments
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 2) {
        "nDMaterial ASDPlasticMaterial3D Error: Few arguments \n";
        opserr << "    numArgs = " << numArgs << endln << endln;
        print_usage();
        return nullptr;
    }

    int numData;

    int tag = 0;
    numData = 1;
    if (OPS_GetInt(&numData, &tag) != 0)  {
        opserr << "nDMaterial ASDPlasticMaterial3D Error: invalid 'tag'.\n";
        return nullptr;
    }


    const char *yf_type = nullptr;
    const char *pf_type = nullptr;
    const char *el_type = nullptr;
    const char *iv_type = nullptr;

    // Now use conditional checks based on numArgs to assign values
    yf_type = numArgs >= 2 ? OPS_GetString() : " X ";
    pf_type = numArgs >= 3 ? OPS_GetString() : " X ";
    el_type = numArgs >= 4 ? OPS_GetString() : " X ";
    iv_type = numArgs >= 5 ? OPS_GetString() : " X ";
    
    opserr << "    numArgs = " << numArgs << endln << endln;


    cout << "Searching for instance with:\n";


    cout << "yf_type = " << yf_type << "\n";
    cout << "pf_type = " << pf_type << "\n";
    cout << "el_type = " << el_type << "\n";
    cout << "iv_type = " << iv_type << "\n";

    std::list<model_spec_t> available_models;

    NDMaterial* instance = ASDPlasticMaterial3DFactory(tag, yf_type, pf_type, el_type, iv_type, available_models);

    if(std::strcmp(yf_type, "list")==0)
    {
        cout << "Available models for YF = " << pf_type << endln;
        for(model_spec_t &model : available_models)
        {
            std::string model_yf_type = std::get<0>(model);
            std::string model_pf_type = std::get<1>(model);
            std::string model_el_type = std::get<2>(model);
            std::string model_iv_type = std::get<3>(model);
            
            if (std::strcmp(pf_type, model_yf_type.c_str())==0 || std::strcmp(pf_type, " X ")==0 )
            {
                cout << "  YF = " << model_yf_type << endl;
                cout << "  PF = " << model_pf_type << endl;
                cout << "  EL = " << model_el_type << endl;
                cout << "  IV = " << model_iv_type << endl << endln;
            }
        }
    }

    if(instance==nullptr)
    {

        bool matches_yf = false;
        bool matches_pf = false;
        bool matches_el = false;


        for(model_spec_t &model : available_models)
        {
            std::string model_yf_type = std::get<0>(model);
            std::string model_pf_type = std::get<1>(model);
            std::string model_el_type = std::get<2>(model);
            // std::string model_iv_type = std::get<3>(model);
            if (std::strcmp(yf_type, model_yf_type.c_str())==0)
                matches_yf = true;
            if (std::strcmp(pf_type, model_pf_type.c_str())==0)
                matches_pf = true;
            if (std::strcmp(el_type, model_el_type.c_str())==0)
                matches_el = true;
        }


        cout << "\n";
        cout << "ASDPlasticMaterial3D -- ERROR! Material not found for input specification:\n";
        cout << "yf_type = " << yf_type << (matches_yf ? " :) " : " ") <<"\n";
        cout << "pf_type = " << pf_type << (matches_pf ? " :) " : " ") <<"\n";
        cout << "el_type = " << el_type << (matches_el ? " :) " : " ") <<"\n";
        cout << "iv_type = " << iv_type  <<"\n";
        cout << "\n";

        cout << endl;
        print_usage();
        cout << endl;

    }


    return instance;

}



NDMaterial*  ASDPlasticMaterial3DFactory(int instance_tag, const char * yf_type, const char * pf_type, const char * el_type, const char * iv_type, std::list<model_spec_t> &available_models)
{

    std::list<NDMaterial*> instance_pointers;


    // The combinations are instantiated in the generated ASD_material_definitions_*.cpp units
    ASDPlasticMaterial3D_all_definitions(instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);

    //Search for the valid pointer and return that one
    for(auto instance : instance_pointers)
    {
        if(instance != nullptr)
        {
            return instance;
        }
    }

    return nullptr;
}

#endif // _EIGEN3
