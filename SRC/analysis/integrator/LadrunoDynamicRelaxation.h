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
//   Underwood, P. (1983). Dynamic relaxation. In T. Belytschko and
//     T.J.R. Hughes (eds.), Computational Methods for Transient Analysis,
//     North-Holland, Amsterdam, pp. 245-265.
//   Cundall, P.A. (1976). Explicit finite difference methods in
//     geomechanics. Proceedings of the 2nd International Conference on
//     Numerical Methods in Geomechanics, Blacksburg, VA, pp. 132-150.
//   Papadrakakis, M. (1981). A method for the automatic evaluation of the
//     dynamic relaxation parameters. Computer Methods in Applied Mechanics
//     and Engineering, 25(1), 35-48. doi:10.1016/0045-7825(81)90066-9
//
// Description: LadrunoDynamicRelaxation, a quasi-static solver by explicit
// dynamic relaxation (Underwood 1983). A fictitious-mass leap-frog march is
// damped until it comes to rest at the static equilibrium state:
//
//   M* a_n = f_ext - f_int(u_n)
//   v_{n+1/2} = v_{n-1/2} + dt a_n ,   u_{n+1} = u_n + dt v_{n+1/2}
//
// Kinetic damping (Cundall 1976, the default): all velocities are set to
// zero at each peak of the kinetic energy 1/2 v^T M* v. -damping viscous
// <zeta> applies instead a mass-proportional damping C* = cVisc M* with
// ratio zeta of the frequency bound of M* (cVisc = 2 zeta omega), lagged at
// v_{n+1/2}.
//
// M* is a diagonal fictitious mass owned by the integrator
// (LadrunoFictitiousMass.h), by default the Gershgorin row sums of the
// element tangents scaled by dt^2 / (4 f^2), which bounds omega_max dt <= 2f
// (Papadrakakis 1981). It does not use Element::getMass(), which is zero
// for models without density. formTangent() assembles M* on the diagonal of
// the SOE, so the solve is a diagonal inverse and K is never factorized.
// The class runs under DirectIntegrationAnalysis with `algorithm Linear`
// and `system Diagonal`; analyze(n, dt) advances n relaxation steps and the
// script decides when the state has come to rest (residualNorm,
// kineticEnergy through the `ladrunoDR` command).
//
// Stability margin (-massSafety f, default 0.5): with f = 1 the Gershgorin
// bound is attained with equality, omega_max dt = 2, which is the
// central-difference stability limit, where the amplification matrix has a
// double eigenvalue -1 and round-off grows linearly instead of decaying. A
// relaxation takes many more steps than a physical explicit run, and a hold
// on an exact equilibrium can drift away from it without any error being
// reported. f < 1 keeps the march inside the limit at a cost of f^2 in
// relaxation per step.
//
// The bound holds for the tangent M* was built from. When the model
// stiffens (geometric stiffening, elastic unloading after yield) the
// margin is spent. M* is rebuilt at each kinetic-energy peak (where the
// velocities are zero, so no energy is injected) unless -noAutoRefresh is
// given, and every -recompute N steps. Independently of the rebuild policy,
// every -marginEvery N steps (default 500) the row sums of the current
// tangent are compared with M*; `ladrunoDR stabilityMargin` returns the
// largest (omega_max dt / 2)^2 seen since the last domainChanged(). A value
// above 1 means the march has crossed the stability limit; -1 means no
// Gershgorin mass (lumped or unity) and -2 that the probe is off.
//
// revertToLastStep() restores the march state stored at the last commit
// (u, v_{n-1/2}, a_n, M* and the damping bookkeeping), so a failed step
// that is retried continues exactly as an uninterrupted run.

#ifndef LadrunoDynamicRelaxation_h
#define LadrunoDynamicRelaxation_h

#include <TransientIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class DOF_Group;
class Vector;

// default safety factor f of the Gershgorin mass, omega_max dt <= 2f
#define LADRUNO_DR_DEFAULT_MASS_SAFETY 0.5

// default step cadence of the stability-margin probe (one element tangent
// pass every N steps, without replacing M*); 0 switches the probe off
#define LADRUNO_DR_DEFAULT_MARGIN_EVERY 500

class LadrunoDynamicRelaxation : public TransientIntegrator
{
  public:
    // massMode: 0 gershgorin (default), 1 lumped (row sums of the element
    // mass times massScale), 2 unity. All arguments are defaulted, so the
    // broker can construct the object before recvSelf().
    LadrunoDynamicRelaxation(int massMode = 0, double dtPseudo = 1.0,
                             double massScale = 1.0, int recomputeEvery = 0,
                             bool interp = false, double divergenceFactor = 0.0,
                             bool verbose = false,
                             int dampMode = 0, double zetaTarget = 1.0,
                             bool autoRefresh = true,
                             double massSafety = LADRUNO_DR_DEFAULT_MASS_SAFETY,
                             int marginEvery = LADRUNO_DR_DEFAULT_MARGIN_EVERY);
    ~LadrunoDynamicRelaxation();

    int newStep(double deltaT);
    int update(const Vector &deltaU);
    int commit(void);
    int domainChanged(void);
    // restore the march state of the last commit after a failed step
    int revertToLastStep(void);

    // assemble M* on the diagonal of the SOE (K is not assembled)
    int formTangent(int statusFlag = CURRENT_TANGENT);
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);

    const Vector &getVel(void);

    // ||f_ext - f_int||_inf of the last step, the force criterion for rest
    // (the physical kinetic energy is zero or meaningless with a fictitious
    // mass), and the fictitious kinetic energy 1/2 v^T M* v
    double getResidualNorm(void) const;
    double getKineticEnergy(void) const;

    // largest (omega_max dt / 2)^2 = max_i (dt^2/4) sum_j |K_ij| / M*_i
    // measured since the last domainChanged(), with K the tangent at the
    // time of the measurement and M* the mass in use. Equals massSafety^2
    // while the tangent does not change; above 1 the march has crossed the
    // explicit stability limit. Measured at each M* rebuild and every
    // marginEvery steps. Negative values mean not measured: -1 for a lumped
    // or unity mass, -2 with the probe off.
    double getStabilityMargin(void) const;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

  protected:

  private:
    int buildFictitiousMass(void);
    // compare the current tangent with M* without changing M*
    int measureStabilityMargin(void);
    // record one margin sample (current-tangent mass vs mass in use)
    void applyMargin(const Vector &live, const Vector &march);
    // store the march state at a commit or at domainChanged()
    void takeSnapshot(void);

    // options
    int    massMode;                 // 0 gershgorin, 1 lumped, 2 unity
    double dtPseudo;                 // pseudo time step used to scale M*
    double massScale;                // scale of the lumped mass
    int    recomputeEvery;           // rebuild M* every N steps (0 = never)
    bool   interp;                   // reserved
    double divergenceFactor;         // abort when KE grows by this factor (0 = off)
    bool   verbose;

    // fictitious mass and leap-frog state
    Vector *Mstar;                   // diagonal fictitious mass
    Vector *Mprobe;                  // row sums for the margin probe
    Vector *Ut, *Vhalf, *Aprev, *Vfull, *Azero;
    bool   firstStep;
    int    updateCount;              // one solve per step
    int    size;
    int    stepCount;                // -recompute cadence

    // damping bookkeeping
    double prevKE;                   // kinetic energy of the previous step
    double kineticEnergy;            // 1/2 v^T M* v
    double residualNorm;             // ||f_ext - f_int||_inf = ||M* a||_inf
    double deltaT;

    int    dampMode;                 // 0 kinetic, 1 viscous
    double zetaTarget;               // viscous damping ratio
    double cVisc;                    // C* = cVisc M*
    bool   autoRefresh;              // rebuild M* at kinetic-energy peaks

    // stability margin
    double massSafety;               // f in omega_max dt <= 2f (gershgorin)
    double stabMargin;               // largest measured margin, < 0 none
    bool   marginWarned;             // warning printed once
    int    marginEvery;              // probe cadence in steps (0 = off)
    int    marginCount;              // steps since the last probe
    bool   haveBuiltOnce;            // M* has been built

    // march state at the last commit. M* is included because -recompute
    // and the kinetic-energy-peak rebuild change it inside a step. The
    // margin diagnostics are not reverted.
    Vector *UtC, *VhalfC, *AprevC, *MstarC;
    bool   snapFirstStep;
    int    snapStepCount;
    double snapPrevKE, snapKE, snapRes, snapCVisc;
    bool   haveSnap;
};

#endif
