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
//   Belytschko, T., Liu, W.K., Moran, B. and Elkhodary, K.I. (2014).
//     Nonlinear Finite Elements for Continua and Structures, 2nd ed.,
//     Wiley, Ch. 6 (central difference method, stability with damping).
//   Hughes, T.J.R. (1987). The Finite Element Method: Linear Static and
//     Dynamic Finite Element Analysis, Prentice-Hall, Ch. 9 (amplification
//     matrix, period error of the central difference method).
//   Bathe, K.J. (1996). Finite Element Procedures, Prentice-Hall, Sec. 9.2
//     (starting procedure of the central difference method).
//
// Description: explicit leap-frog central difference integrator with a
// starting procedure computed from the initial state, a critical time step
// estimate, and full-step velocity output.
//
// Scheme (mass alone on the left-hand side, damping lagged at the known
// half-step velocity):
//
//   a_n       = M^-1 ( P_n - C v_{n-1/2} - F_int(u_n) )
//   v_{n+1/2} = v_{n-1/2} + dt a_n
//   u_{n+1}   = u_n + dt v_{n+1/2}
//
// This is Newmark(beta = 0, gamma = 1/2): second order, no algorithmic
// dissipation, relative period error -Omega^2/24 (period shortening).
//
// Starting procedure. On the first newStep() the initial acceleration
//   a_0 = M^-1 ( P_0 - C v_0 - F_int(u_0) )
// is solved and the half-step velocity is seeded as v_{-1/2} = v_0 - dt/2 a_0.
// This is done in newStep() and not in domainChanged(), which has neither a
// time step nor a factored system; the existing CentralDifference classes
// assume u_{-1} = u_0 there instead.
//
// Velocities. getVel() returns the half-step velocity v_{n-1/2} known at
// solve time; it is the velocity of the modal damping force and of the
// element damping forces. Nodes and recorders receive the full-step velocity
// v_{n+1} = v_{n+1/2} + dt/2 a_{n+1}.
//
// Critical time step. computeCriticalTimeStep() (CriticalTimeStep.h) runs its
// own element eigensolves, so dt_cr is available after domainChanged(), before
// the first step. The central-difference limit is 2/omega_max; with damping
// the reduced limit (2/omega_max)(sqrt(1 + xi^2) - xi) is reported. A damped
// limit below 0.9 of the undamped one prints a note: stiffness-proportional
// Rayleigh damping (xi = betaK omega/2) reduces the step to about
// 2/(betaK omega_max^2), while mass-proportional damping barely moves it.
//
// Recipe: lumped (diagonal) element mass, system Diagonal, algorithm Linear
// (one solve per step; a second update() in a step is an error) and
// dt < dt_cr. For a coupled or consistent-mass central difference use
// integrator NewmarkExplicit 0.5.

#ifndef CentralDifferenceLadruno_h
#define CentralDifferenceLadruno_h

#include <TransientIntegrator.h>
#include <CriticalTimeStep.h>
#include <LadrunoProjectionConsumer.h>

class DOF_Group;
class FE_Element;
class Vector;
class LinearSOE;
class LadrunoConstraintProjector;

class CentralDifferenceLadruno : public TransientIntegrator,
                                 public LadrunoProjectionConsumer
{
public:
    CentralDifferenceLadruno();
    CentralDifferenceLadruno(int compute_critical_timestep,
                             bool verbose, bool cflAbort,
                             double divergenceFactor,
                             bool cflUseTangent, int cflRecomputeEvery,
                             CTSLumping lumping);

    ~CentralDifferenceLadruno();

    // Central-difference critical time step (the damped limit when damping
    // is active); <= 0 if not computed or not applicable.
    double getCriticalTimeStep(void) const override;

    // Mass alone on the left-hand side; no damping, no stiffness.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);

    // The assembled tangent is the mass matrix alone (the modal damping
    // matrix branch of TransientIntegrator::formTangent is scaled by
    // getCFactor() = 0), so it is constant between events that change the
    // mass or resize the system. Once formed it is kept, and the factored
    // diagonal of system Diagonal stays valid, as with
    // `algorithm Linear -factorOnce`. The cache is invalidated by
    // domainChanged() and revertToLastStep(). updateParameter on a density
    // changes the mass without a domain change: use -noMassCache for that.
    int formTangent(int statFlag);

    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    // Re-sync with the reverted Domain after a failed step. The analysis
    // calls Domain::revertToLastCommit() and then this method. The leap-frog
    // advance in newStep() has already changed the state vectors, so they
    // are re-seeded from the committed DOF state and the starting procedure
    // is re-armed: the retry rebuilds v_{-1/2} = v_n - dt/2 a_n for the dt it
    // uses. The restart is exact for undamped models; with damping it is a
    // second-order restart (the starter evaluates C v at v_n instead of the
    // lagged half-step velocity).
    int revertToLastStep(void);

    // Half-step velocity v_{n-1/2} (modal damping).
    const Vector &getVel(void);

    // Time step of the current step (> 0 once newStep() has run). Penalty
    // contact formulations size a Courant-stable penalty from it.
    double getCurrentDeltaT(void) const { return deltaT; }

    // -noMassCache: assemble the mass every step.
    void setMassCacheEnabled(bool on) { useMassCache = on; if (!on) massTangentValid = false; }

    // -commitSolveState (default off): skip the second element state
    // determination in update(). Element state is then committed as the
    // solve used it (u_{n+1}, v_{n+1/2}) instead of being recomputed with the
    // full-step velocity v_{n+1}. Identical for rate-independent materials;
    // it changes results for elements that pass a velocity to their material
    // (zeroLength family, TwoNodeLink, truss family, multiple-spring elements
    // with viscous materials).
    void setCommitSolveState(bool on) { commitSolveState = on; }

    // LadrunoProjectionConsumer: stores the (non-owning) acceleration
    // projector of a projection constraint handler and invalidates the mass
    // tangent cache, since a new projector comes with a new constraint set
    // and DOF numbering. With a projector, diag(M) is read in the starting
    // procedure between the mass-only formTangent() and the first solve, the
    // starting acceleration a_0 and every solved a_{n+1} are projected
    // before they form velocities, and domainChanged() checks (or, with
    // -projectICs, sets) the initial displacement and velocity.
    void setConstraintProjector(LadrunoConstraintProjector *theProjector) override;

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

private:
    // Seed Ut, Vhalf and Aprev from the committed DOF state (Vhalf receives
    // the full-step v_n; the starting procedure re-centres it).
    void seedFromCommitted(AnalysisModel *theModel);

    double deltaT;

    // leap-frog state
    Vector *Ut;        // displacement u_n (advanced to u_{n+1} in newStep)
    Vector *Vhalf;     // half-step velocity v_{n-1/2}
    Vector *Aprev;     // acceleration a_n carried to the next newStep
    Vector *Vfull;     // full-step velocity output v_{n+1}
    Vector *Azero;     // zero vector: setAccel(0) before each solve, so the
                       //   residual carries no inertia term

    int updateCount;   // one solve per step
    bool firstStep;    // starting procedure pending

    // critical time step
    int    compute_critical_timestep;     // 0 off, 1 (re)compute, 2 computed
    double damped_minimum_critical_timestep;
    double undamped_minimum_critical_timestep;
    int    damped_critical_element_tag;
    int    undamped_critical_element_tag;

    // diagnostics and run control
    bool   verbose;            // per-step report (default off)
    bool   cflAbort;           // abort the step if dt exceeds dt_cr
    double divergenceFactor;   // > 0: abort on kinetic-energy growth
    double prevKE;             // running maximum of the 0.5 v.v proxy
    double committedPrevKE;    // prevKE at newStep() entry (revert)

    // critical time step refresh
    bool   cflUseTangent;      // tangent instead of initial stiffness
    int    cflRecomputeEvery;  // refresh every N steps (0 = once)
    int    cflStepCount;
    int    committedCflStepCount;
    bool   cflFirstComputation;
    CTSLumping lumping;        // element mass lumping of the pencil (-lump)
    bool   betaKWarned;

    // mass tangent cache (transient, not serialized)
    bool   useMassCache;
    bool   massTangentValid;
    bool   massCacheNoted;

    // -commitSolveState (transient, not serialized)
    bool   commitSolveState;
    bool   commitSolveStateNoted;

    // constraint projector (not owned; null unless a projection handler
    // provides one)
    LadrunoConstraintProjector *theProjector;
    bool    massBuilt;         // projector has read diag(M) since the last change
    Vector *Aproj;             // projected acceleration a_{n+1} (scratch)
};

#endif
