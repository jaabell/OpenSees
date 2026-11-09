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

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell (Ladruño project)
//
// References:
//   J. Chung, G.M. Hulbert (1993). A time integration algorithm for structural
//     dynamics with improved numerical dissipation: the generalized-alpha
//     method. Journal of Applied Mechanics 60(2):371-375.
//     doi:10.1115/1.2900803
//   M. Kleiber, H. Antunez, T.D. Hien, P. Kowalczyk (1997). Parameter
//     Sensitivity in Nonlinear Mechanics: Theory and Finite Element
//     Computations. Wiley, Chichester. ISBN 0-471-96854-4.
//   J.P. Conte, P.K. Vijalapura, M. Meghella (2003). Consistent finite-element
//     response sensitivity analysis. Journal of Engineering Mechanics
//     129(12):1380-1393. doi:10.1061/(ASCE)0733-9399(2003)129:12(1380)
//   T. Haukaas, A. Der Kiureghian (2005). Parameter sensitivity and importance
//     measures in nonlinear finite element reliability analysis. Journal of
//     Engineering Mechanics 131(10):1013-1026.
//     doi:10.1061/(ASCE)0733-9399(2005)131:10(1013)
//
// Description: LadrunoGeneralizedAlpha is the generalized-alpha integrator
// (Chung and Hulbert 1993) with Direct Differentiation Method (DDM) response
// sensitivities. The primal algorithm is inherited from GeneralizedAlpha
// unchanged. In the OpenSees convention alphaF = 1 - alpha_f and
// alphaM = 1 - alpha_m of Chung and Hulbert, and the residual is
//
//     R = F(t_n + alphaF*dt) - P(U_aF) - C*Udot_aF - M*Udotdot_aM
//     U_aF       = (1-alphaF)*U_n       + alphaF*U_{n+1}
//     Udot_aF    = (1-alphaF)*Udot_n    + alphaF*Udot_{n+1}
//     Udotdot_aM = (1-alphaM)*Udotdot_n + alphaM*Udotdot_{n+1}
//
// with U, Udot and Udotdot following the Newmark relations. The sensitivity
// propagation (saveSensitivity, commitSensitivity, formSensitivityRHS) is the
// Newmark one. The right-hand side of the sensitivity equation differs:
//   - the mass multiplier carries the alphaM weighting of dUdotdot_aM/dh,
//   - the damping multiplier carries the alphaF weighting of dUdot_aF/dh,
//   - dM/dh acts on Udotdot_aM and dC/dh on Udot_aF,
//   - an element stiffness term -(1-alphaF)*K*dU_n/dh appears, the
//     (1-alphaF) part of dU_aF/dh.
// The left-hand side is alphaF*K + alphaF*c2*C + alphaM*c3*M, the exact
// Jacobian of R with respect to U_{n+1}. At alphaM = alphaF = 1 the class
// reduces to the Newmark DDM.
//
// computeSensitivities() re-forms the current tangent once before the
// sensitivity solves. The factorization left by the primal algorithm can be
// stale (ModifiedNewton) or built from the initial stiffness (-initial), and
// either would bias the gradients.

#ifndef LadrunoGeneralizedAlpha_h
#define LadrunoGeneralizedAlpha_h

#include <GeneralizedAlpha.h>
#include <Vector.h>

class DOF_Group;
class FE_Element;

class LadrunoGeneralizedAlpha : public GeneralizedAlpha
{
public:
    // constructors (same signatures and defaults as GeneralizedAlpha)
    LadrunoGeneralizedAlpha();
    LadrunoGeneralizedAlpha(double alphaM, double alphaF);
    LadrunoGeneralizedAlpha(double alphaM, double alphaF, double beta, double gamma);

    // destructor
    ~LadrunoGeneralizedAlpha();

    // the residual hooks branch on sensitivityFlag; with the flag off they
    // call the TransientIntegrator path, so a run without sensitivities is
    // identical to GeneralizedAlpha.
    int formEleResidual(FE_Element *theEle);
    int formNodUnbalance(DOF_Group *theDof);

    // the tangent hooks branch on sensTangentFlag; with the flag off they
    // call GeneralizedAlpha. With the flag on they form the current tangent
    // used by the sensitivity solves.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);

    const char *getClassType(void) const;
    void Print(OPS_Stream &s, int flag = 0);

    // AddingSensitivity:BEGIN //////////////////////////////////
    int revertToStart();
    int formSensitivityRHS(int gradNum);
    int formIndependentSensitivityRHS();
    int saveSensitivity(const Vector &v, int gradNum, int numGrads);
    int commitSensitivity(int gradNum, int numGrads);
    int computeSensitivities();
    // AddingSensitivity:END ////////////////////////////////////

protected:

private:
    // DDM workspace, same roles as the Newmark sensitivity members
    int sensitivityFlag;
    int gradNumber;
    Vector *massMatrixMultiplicator;
    Vector *dampingMatrixMultiplicator;
    int assemblyFlag;
    Vector independentRHS;
    int sensTangentFlag;   // 1 while forming the tangent for the sensitivity solves
};

#endif
