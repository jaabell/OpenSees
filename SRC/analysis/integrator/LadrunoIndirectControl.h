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
//   Batoz, J.L. and Dhatt, G. (1979). Incremental displacement algorithms
//     for nonlinear problems. International Journal for Numerical Methods
//     in Engineering, 14(8), 1262-1267. doi:10.1002/nme.1620140811
//   de Borst, R. (1987). Computation of post-bifurcation and post-failure
//     behavior of strain-softening solids. Computers & Structures, 25(2),
//     211-224. doi:10.1016/0045-7949(87)90144-1
//   Ramm, E. (1981). Strategies for tracing the nonlinear response near
//     limit points. In W. Wunderlich, E. Stein and K.-J. Bathe (eds.),
//     Nonlinear Finite Element Analysis in Structural Mechanics, Springer,
//     Berlin, pp. 63-89. doi:10.1007/978-3-642-81589-8_5
//
// Description: LadrunoIndirectControl, indirect displacement control
// (de Borst 1987). DisplacementControl prescribes the increment of one
// nodal DOF; this integrator prescribes the increment of a weighted sum of
// DOFs
//     zeta = c . U = sum_k coef_k U(node_k, dof_k),
// for example a crack mouth opening u_A - u_B. A relative quantity of this
// kind can stay monotone through snap-back, where the load factor and the
// individual DOFs both reverse.
//
// Constraint per step:  c . DeltaU_step = Delta-zeta (the increment incr)
//   predictor:  K dUhat = phat ;  dLambda = Delta-zeta / (c . dUhat) ;
//               dU = dLambda dUhat
//   corrector:  K dUhat = phat ;  dLambda = -(c . dUbar) / (c . dUhat) ;
//               dU = dUbar + dLambda dUhat
// with phat the reference load and dUbar the residual displacement of the
// iteration (Batoz and Dhatt 1979). With one entry and coef = 1 this is
// DisplacementControl. -iter numIter dmin dmax scales the increment by
// numIter / Jlast, Jlast being the iterations of the last step (Ramm 1981),
// and clamps it to [dmin, dmax]. The DDM sensitivity terms are not
// implemented.

#ifndef LadrunoIndirectControl_h
#define LadrunoIndirectControl_h

#include <StaticIntegrator.h>
#include <Vector.h>
#include <ID.h>

class LinearSOE;
class AnalysisModel;
class Domain;
class FE_Element;

class LadrunoIndirectControl : public StaticIntegrator
{
  public:
    // nodes/dofs/coefs define the control quantity c.U (dofs are 1-based on input,
    // stored 0-based). incr = control increment Delta-zeta. numIncr>0 enables the
    // increment adaptation; min/max clamp the signed increment.
    LadrunoIndirectControl(const ID &nodes, const ID &dofs, const Vector &coefs,
                           double incr, int numIncr = 1,
                           double min = 0.0, double max = 0.0);
    LadrunoIndirectControl();   // for the FEM_ObjectBroker (must recvSelf)
    ~LadrunoIndirectControl();

    int newStep(void);
    int update(const Vector &deltaU);
    int domainChanged(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

  protected:

  private:
    double controlDot(const Vector &v) const;   // c.v over the controlled eqns

    // --- control definition (input) ---
    ID     ctrlNode;            // node tags
    ID     ctrlDof;             // 0-based dof indices
    Vector ctrlCoef;           // participation coefficients
    double theIncrement;        // control increment Delta-zeta (adapted if numIncr>0)
    int    specNumIncrStep;     // desired corrector iterations/step (Ramm target)
    int    numIncrLastStep;     // corrector iterations of the last step (Jlast)
    double minIncrement;        // increment-magnitude clamp
    double maxIncrement;

    // --- resolved at domainChanged ---
    ID     ctrlEqn;             // global equation number per control entry (-1 if fixed or constrained)

    // per-step state, as in DisplacementControl
    Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep;
    Vector *phat;               // reference load vector
    double  deltaLambdaStep, currentLambda;
};

#endif
