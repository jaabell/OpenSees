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
//   Riks, E. (1979). An incremental approach to the solution of snapping
//     and buckling problems. International Journal of Solids and
//     Structures, 15(7), 529-551. doi:10.1016/0020-7683(79)90081-7
//   Crisfield, M.A. (1981). A fast incremental/iterative solution procedure
//     that handles "snap-through". Computers & Structures, 13(1-3), 55-62.
//     doi:10.1016/0045-7949(81)90108-5
//   Ramm, E. (1981). Strategies for tracing the nonlinear response near
//     limit points. In W. Wunderlich, E. Stein and K.-J. Bathe (eds.),
//     Nonlinear Finite Element Analysis in Structural Mechanics, Springer,
//     Berlin, pp. 63-89. doi:10.1007/978-3-642-81589-8_5
//   Dassault Systemes Simulia. Abaqus Analysis User's Guide, "Automatic
//     stabilization of unstable problems" (STABILIZE parameter of *STATIC).
//
// Description: LadrunoArcLength, a static integrator that contains the
// cylindrical arc-length method of ArcLength and adds three optional
// features, all off by default. With no option the class reproduces
// ArcLength step for step.
//
// Constraint (as ArcLength, Crisfield 1981):
//   first iteration   dU^T dU + alpha^2 dLambda^2 = arcLength^2
//   later iterations  the ArcLength quadratic in dLambda, root chosen by
//                     the positive angle with the step increment.
//
// 1. Iteration-count adaptation of the radius (-adapt Jd ellMin ellMax
//    <-p exp>), the rule of Ramm (1981) that LoadControl,
//    DisplacementControl and MinUnbalDispNorm already apply to their
//    increments:
//      arcLength <- clamp(arcLength * (Jd / Jlast)^p, ellMin, ellMax)
//    with Jlast the corrector iterations of the last converged step. p = 1
//    is the LoadControl rule; p = 1/2 is the variant of Crisfield (1981).
//
// 2. Viscous stabilization (-stabilize <f>), the static stabilization of
//    Abaqus (STABILIZE). The arc-length constraint is not used in this mode:
//    the step is incremental load control with load increment arcLength,
//    and the equations solved are
//      (K + (c/dt) M*) dU = lambda p - f_int - (c/dt) M* DeltaU_step
//    with M* the diagonal fictitious mass of LadrunoFictitiousMass.h and
//    DeltaU_step the displacement increment of the step. The viscous term
//    keeps the iteration matrix nonsingular through limit points and local
//    instabilities. Unless -cVisc gives c, c is calibrated at the first
//    converged step so that the viscous work of that step is the fraction f
//    (default 2e-4) of the strain energy 0.5 |DeltaU_step . p| |lambda|;
//    -adaptStab rescales c after each step to keep the cumulative fraction
//    at f.
//
// 3. Cut-and-retry from a script. commit() stores the converged per-step
//    state {deltaLambdaStep, deltaUstep, currentLambda, sign}, and
//    revertToLastStep(), which StaticAnalysis::analyze calls after a failed
//    step, restores it. ArcLength itself overwrites deltaLambdaStep and
//    deltaUstep in newStep() before a failure can be detected, so without
//    the stored copy a retried step would start from a polluted state. The
//    radius is not restored: it belongs to the script, which changes it
//    through the `ladrunoArcLength` command (reduceStep, increaseStep,
//    setArcLength) between attempts.
//
// The DDM sensitivity terms of ArcLength are not implemented.

#ifndef LadrunoArcLength_h
#define LadrunoArcLength_h

#include <StaticIntegrator.h>
#include <Vector.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;

class LadrunoArcLength : public StaticIntegrator
{
  public:
    // all arguments after alpha are defaulted, so the broker can construct
    // the object with LadrunoArcLength(1.0) before recvSelf()
    LadrunoArcLength(double arcLength, double alpha = 1.0,
                     bool adapt = false, int Jd = 5,
                     double ellMin = 0.0, double ellMax = 0.0,
                     double pExp = 1.0,
                     bool stabilize = false, double fTarget = 2.0e-4,
                     bool adaptStab = false, int massMode = 0,
                     double massScale = 1.0, double cVisc = 0.0,
                     double dtPseudo = 1.0);
    ~LadrunoArcLength();

    int newStep(void);
    int update(const Vector &deltaU);
    int domainChanged(void);

    // -stabilize: K + (c/dt) M* on the diagonal, viscous force on the residual
    int formTangent(int statusFlag = CURRENT_TANGENT);
    int formUnbalance(void);
    int commit(void);
    int revertToLastStep(void);

    // radius mutators for script-driven cut-and-retry
    int    reduceStep(double f);       // arcLength *= f, floored at ellMin
    int    increaseStep(double f);     // arcLength *= f, capped at ellMax if set
    int    setArcLength(double arcLength);

    // read-only state for the `ladrunoArcLength` command
    double getArcLength(void) const;
    double getDeltaLambdaStep(void) const;
    double getCurrentLambda(void) const;
    int    getSignLastDeltaLambdaStep(void) const;
    double getDeltaUstepNorm(void) const;

    // In -stabilize mode the SOE residual is lambda p - f_int - f_v. These
    // return the residual before the viscous force f_v is subtracted, for
    // the LadrunoStabilizedUnbalance convergence test. The scalar form
    // returns -1 and the vector form 0 when the integrator is not
    // stabilizing.
    int getStabilizedTrueResidual(double &normOut) const;
    const Vector *getStabilizedTrueResidualVector(void) const;

    // stabilization energy: cumulative viscous work, the reference strain
    // energy of the calibration step, their ratio, and a runtime scale
    // factor on c (factor > 0)
    double getStabilizationDissipatedEnergy(void) const;
    double getReferenceStrainEnergy(void) const;
    double getStabilizationDissipationRatio(void) const;
    int    scaleCVisc(double factor);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

  protected:

  private:
    int buildArtificialMass(void);

    // ArcLength state
    double arcLength2;                 // arc radius squared
    double alpha2;                     // load-factor scaling squared
    double a, b, c, b24ac;             // quadratic coefficients
    Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep;
    Vector *phat;                      // reference load vector
    double deltaLambdaStep, currentLambda;
    int    signLastDeltaLambdaStep;

    // radius adaptation
    bool   adapt;
    int    Jd;                         // desired corrector iterations per step
    double ellMin, ellMax;             // radius (or load increment) clamp
    double pExp;                       // exponent of (Jd/Jlast)
    int    numIncrLastStep;            // corrector iterations of the last step

    // viscous stabilization
    bool   stabilize;
    bool   adaptStab;                  // rescale c after each step
    double fTarget;                    // dissipated energy fraction
    double cVisc;                      // viscous factor c
    double dtPseudo;                   // pseudo time step dt
    double cOverDt;                    // c / dt, the coefficient applied
    bool   cCalibrated;                // c fixed (given or calibrated)
    double Estrain0;                   // strain energy of the calibration step
    double dissipVisc;                 // cumulative viscous work
    double residualTrueNorm;           // ||lambda p - f_int||_2 before f_v
    Vector residualTrue;               // lambda p - f_int before f_v
    double dLambda0;                   // load increment of the stabilized step
    int    massMode;                   // 0 gershgorin, 1 lumped, 2 unity
    double massScale;                  // scale for the lumped mass
    int    nEqn;
    Vector *Mstar;                     // diagonal fictitious mass

    // converged per-step state restored by revertToLastStep()
    double  cArcLength2;
    double  cDeltaLambdaStep;
    double  cCurrentLambda;
    int     cSign;
    Vector *cDeltaUstep;
    bool    haveSnapshot;              // false until the first commit()
};

#endif
