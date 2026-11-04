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


#ifndef ASDPlasticMaterial3DGlobals_H
#define ASDPlasticMaterial3DGlobals_H

#include "EigenAPI.h"
#include <limits>
#include <tuple>
#include <string>
namespace ASDPlasticMaterial3DGlobals
{


constexpr double SQRT_2  = 1.4142135623730951;//sqrt(2.0 ); //in case unsupported by compiler
constexpr double SQRT_3  = 1.7320508075688772;//sqrt(3.0); //in case unsupported by compiler
constexpr double SQRT_2_over_3  = 0.816496580928;//sqrt(2.0 / 3.0); //in case unsupported by compiler
//constexpr double SQRT_2_over_3  = sqrt(2.0 / 3.0);
constexpr double SQRT_2_over_27 = 0.272165526976;//sqrt(2.0 / 27.0); //in case unsupported by compiler
//constexpr double SQRT_2_over_27 = sqrt(2.0 / 27.0);

constexpr double MACHINE_EPSILON = std::numeric_limits<double>::epsilon();


enum struct ASDPlasticMaterial3D_Constitutive_Integration_Method : int
{
    Not_Set,
    Forward_Euler,
    Forward_Euler_Crisfield,
    Multistep_Forward_Euler,
    Multistep_Forward_Euler_Crisfield,
    Modified_Euler_Error_Control,
    Runge_Kutta_45_Error_Control,
    Runge_Kutta_45_Error_Control_old,
    Backward_Euler,
    Full_Backward_Euler,
    Forward_Euler_Subincrement,
    Backward_Euler_LineSearch,
};

// Derivatives used by Full_Backward_Euler's Newton Jacobian
enum struct ASDPlasticMaterial3D_Derivative_Type : int
{
    Numerical,      // central differences of the return-mapping residual
    Closed_Form,    // closed forms the components provide, differences for the rest
};

enum struct ASDPlasticMaterial3D_Tangent_Operator_Type : int
{
    Elastic,
    Continuum,
    Secant,
    Algorithmic,
    Numerical_Algorithmic_FirstOrder,
    Numerical_Algorithmic_SecondOrder,

};

inline void printTensor1(std::string const& name, VoigtVector const& v)
{
    // This is in good format but take 3 lines.
    // stderr will print immediately, not like cout (may be reordered by CPU).
    fprintf(stderr, "%s = \n", name.c_str());
    fprintf(stderr, "[%16.8f \t %16.8f \t %16.8f \t %16.8f \t %16.8f \t %16.8f]\n",   v(0), v(1), v(2),   v(3), v(4), v(5));

}

inline void printTensor2(std::string const& name, VoigtMatrix const& v)
{
    std::cout << name << " = [ " ;
    for (int ii = 0; ii < 6; ii++)
    {
        for (int jj = 0; jj < 6; jj++)
        {
            std::cout << v(ii, jj) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "]" << std::endl;
}


// Macaulay Bracket < >  operator. (Integral of Heaviside function)
inline double macaulay_bracket(double x)
{
    return x > 0 ? x : 0;
}


}

// Voigt conventions. Stresses, yield-function gradients n = df/dsigma and flow
// directions m = dg/dsigma hold tensor components (shear entry = sigma_12).
// Strains (total and plastic) hold OpenSees' engineering shear (gamma_12 =
// 2 eps_12), and the elasticity matrices map engineering strain to stress.
// So a tensor-form m is the plastic strain direction strain_like(m) = W m,
// W = diag(1,1,1,2,2,2), and n:dsigma = strain_like(n).dot(dsigma).
inline VoigtVector strain_like(const VoigtVector& t)
{
    VoigtVector r = t;
    r(3) *= 2;
    r(4) *= 2;
    r(5) *= 2;
    return r;
}

// Inverse of strain_like: engineering-shear strain to tensor components.
inline VoigtVector tensor_like(const VoigtVector& e)
{
    VoigtVector r = e;
    r(3) *= 0.5;
    r(4) *= 0.5;
    r(5) *= 0.5;
    return r;
}

#endif