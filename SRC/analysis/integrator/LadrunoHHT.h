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
//   H.M. Hilber, T.J.R. Hughes, R.L. Taylor (1977). Improved numerical
//     dissipation for time integration algorithms in structural dynamics.
//     Earthquake Engineering and Structural Dynamics 5(3):283-292.
//     doi:10.1002/eqe.4290050306
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
// Description: LadrunoHHT is the HHT-alpha integrator (Hilber et al. 1977)
// with Direct Differentiation Method (DDM) response sensitivities. The primal
// algorithm (newStep, update, commit, tangent) is inherited from HHT
// unchanged; the class adds the sensitivity interface that OpenSees provides
// on Newmark, so that reliability and gradient-based analyses can use the
// numerically damped HHT scheme.
//
// HHT evaluates the dynamic residual at the intermediate state
//
//     R = F(t_n + alpha*dt) - P(U_alpha) - C*Udot_alpha - M*Udotdot_{n+1}
//     U_alpha    = (1-alpha)*U_n    + alpha*U_{n+1}
//     Udot_alpha = (1-alpha)*Udot_n + alpha*Udot_{n+1}
//
// while U, Udot and Udotdot follow the Newmark relations. Differentiating R
// with respect to a parameter h at fixed U_{n+1} (Conte et al. 2003) gives the
// same sensitivity propagation as Newmark (saveSensitivity, commitSensitivity,
// formSensitivityRHS, computeSensitivities). Only the right-hand side differs:
//   - the damping multiplier carries the alpha weighting of dUdot_alpha/dh,
//   - dC/dh acts on Udot_alpha instead of Udot_{n+1},
//   - an element stiffness term -(1-alpha)*K*dU_n/dh appears, the (1-alpha)
//     part of dU_alpha/dh; it vanishes at alpha = 1, where the class reduces
//     to the Newmark DDM.
// The mass terms keep the Newmark form because HHT evaluates the inertia at
// t_{n+1}. The left-hand side is the HHT tangent alpha*K + alpha*c2*C + c3*M,
// which is the exact Jacobian of R with respect to U_{n+1}.

#ifndef LadrunoHHT_h
#define LadrunoHHT_h

#include <HHT.h>
#include <Vector.h>

class DOF_Group;
class FE_Element;

class LadrunoHHT : public HHT
{
public:
    // constructors (same signatures and defaults as HHT)
    LadrunoHHT();
    LadrunoHHT(double alpha);
    LadrunoHHT(double alpha, double beta, double gamma);

    // destructor
    ~LadrunoHHT();

    // the residual hooks branch on sensitivityFlag; with the flag off they
    // call the TransientIntegrator path, so a run without sensitivities is
    // identical to HHT.
    int formEleResidual(FE_Element *theEle);
    int formNodUnbalance(DOF_Group *theDof);

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
};

#endif
