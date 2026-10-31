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
// Description: explicit leap-frog central difference integrator. See
// CentralDifferenceLadruno.h for the scheme and its references.

#include <CentralDifferenceLadruno.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Domain.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <classTags.h>
#include <cmath>
#include <limits>
#include <string.h>
#include <stdlib.h>

// A damped critical step below this fraction of the undamped one prints the
// stiffness-proportional damping note. Mass-proportional damping
// (xi = alphaM/(2 omega)) barely changes the step at omega_max; betaK
// (xi = betaK omega/2) reduces it to about 2/(betaK omega_max^2).
#define CDL_BETAK_WARN_RATIO 0.90

void *OPS_CentralDifferenceLadruno(void)
{
    // integrator CentralDifferenceLadruno <-cfl> <-cflAbort> <-tangent>
    //     <-recompute N> <-lump rowsum|diagonal|hrz> <-verbose> <-divergence f>
    //     <-noMassCache> <-commitSolveState>
    int compute_critical_timestep = 0;
    bool verbose = false;
    bool cflAbort = false;
    double divergenceFactor = 0.0;
    bool cflUseTangent = false;
    int cflRecomputeEvery = 0;
    bool massCache = true;
    bool commitSolve = false;
    // the consistent diagonal is positive for the rotational DOFs of beams
    // and shells, where the row sum can vanish
    CTSLumping lumping = CTSLumping::Diagonal;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *arg = OPS_GetString();
        if (strcmp(arg, "-verbose") == 0) {
            verbose = true;
        } else if (strcmp(arg, "-cflAbort") == 0) {
            cflAbort = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-cfl") == 0 ||
                   strcmp(arg, "-criticalTimestep") == 0) {
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-divergence") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int nd = 1;
                OPS_GetDoubleInput(&nd, &divergenceFactor);
            }
        } else if (strcmp(arg, "-tangent") == 0) {
            cflUseTangent = true;
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-recompute") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int nd = 1;
                OPS_GetIntInput(&nd, &cflRecomputeEvery);
            }
            if (cflRecomputeEvery <= 0)
                opserr << "WARNING CentralDifferenceLadruno -recompute expects a positive "
                          "integer N (steps between dt_cr refreshes); dt_cr computed once\n";
            cflUseTangent = true;   // refreshing from the initial stiffness gives the same value
            compute_critical_timestep = 1;
        } else if (strcmp(arg, "-noMassCache") == 0) {
            massCache = false;
        } else if (strcmp(arg, "-commitSolveState") == 0) {
            commitSolve = true;
        } else if (strcmp(arg, "-lump") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *m = OPS_GetString();
                if (strcmp(m, "diagonal") == 0)      lumping = CTSLumping::Diagonal;
                else if (strcmp(m, "rowsum") == 0)   lumping = CTSLumping::RowSum;
                else if (strcmp(m, "hrz") == 0)      lumping = CTSLumping::HRZ;
                else opserr << "WARNING CentralDifferenceLadruno - unknown -lump " << m
                            << " (use rowsum|diagonal|hrz; keeping diagonal)\n";
            }
        } else {
            opserr << "WARNING CentralDifferenceLadruno - unknown option " << arg
                   << " (ignored)\n";
        }
    }

    CentralDifferenceLadruno *theIntegrator =
        new CentralDifferenceLadruno(compute_critical_timestep, verbose, cflAbort,
                                     divergenceFactor, cflUseTangent,
                                     cflRecomputeEvery, lumping);
    if (theIntegrator == 0) {
        opserr << "WARNING - out of memory creating CentralDifferenceLadruno integrator\n";
        return theIntegrator;
    }
    theIntegrator->setMassCacheEnabled(massCache);
    theIntegrator->setCommitSolveState(commitSolve);
    return theIntegrator;
}

CentralDifferenceLadruno::CentralDifferenceLadruno()
    : TransientIntegrator(INTEGRATOR_TAGS_CentralDifferenceLadruno),
      deltaT(0.0),
      Ut(0), Vhalf(0), Aprev(0), Vfull(0), Azero(0),
      updateCount(0), firstStep(true),
      compute_critical_timestep(0),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0),
      verbose(false), cflAbort(false), divergenceFactor(0.0), prevKE(0.0),
      committedPrevKE(0.0),
      cflUseTangent(false), cflRecomputeEvery(0), cflStepCount(0),
      committedCflStepCount(0),
      cflFirstComputation(true), lumping(CTSLumping::Diagonal), betaKWarned(false),
      useMassCache(true), massTangentValid(false), massCacheNoted(false),
      commitSolveState(false), commitSolveStateNoted(false),
      theProjector(0)
{
}

CentralDifferenceLadruno::CentralDifferenceLadruno(
    int compute_critical_timestep_, bool verbose_, bool cflAbort_,
    double divergenceFactor_, bool cflUseTangent_, int cflRecomputeEvery_,
    CTSLumping lumping_)
    : TransientIntegrator(INTEGRATOR_TAGS_CentralDifferenceLadruno),
      deltaT(0.0),
      Ut(0), Vhalf(0), Aprev(0), Vfull(0), Azero(0),
      updateCount(0), firstStep(true),
      compute_critical_timestep(compute_critical_timestep_),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0),
      verbose(verbose_), cflAbort(cflAbort_), divergenceFactor(divergenceFactor_),
      prevKE(0.0), committedPrevKE(0.0), cflUseTangent(cflUseTangent_),
      cflRecomputeEvery(cflRecomputeEvery_), cflStepCount(0),
      committedCflStepCount(0),
      cflFirstComputation(true), lumping(lumping_), betaKWarned(false),
      useMassCache(true), massTangentValid(false), massCacheNoted(false),
      commitSolveState(false), commitSolveStateNoted(false),
      theProjector(0)
{
}

CentralDifferenceLadruno::~CentralDifferenceLadruno()
{
    if (Ut    != 0) delete Ut;
    if (Vhalf != 0) delete Vhalf;
    if (Aprev != 0) delete Aprev;
    if (Vfull != 0) delete Vfull;
    if (Azero != 0) delete Azero;
    // theProjector is owned by the constraint handler
}

void
CentralDifferenceLadruno::setConstraintProjector(LadrunoConstraintProjector *p)
{
    theProjector = p;
    massTangentValid = false;
}

int CentralDifferenceLadruno::formTangent(int statFlag)
{
    if (useMassCache && massTangentValid) {
        statusFlag = statFlag;
        return 0;
    }
    int res = this->TransientIntegrator::formTangent(statFlag);
    if (res == 0 && useMassCache) {
        massTangentValid = true;
        // updateParameter on a density changes M without a domain change,
        // the one event the cache cannot see; note it once when the model
        // declares parameters.
        if (!massCacheNoted) {
            massCacheNoted = true;
            AnalysisModel *theModel = this->getAnalysisModel();
            Domain *theDomain = (theModel != 0) ? theModel->getDomainPtr() : 0;
            if (theDomain != 0 && theDomain->getNumParameters() > 0)
                opserr << "CentralDifferenceLadruno: the mass tangent is assembled once and "
                          "the model declares parameters. An updateParameter that changes a "
                          "density is not picked up; use -noMassCache (or trigger a domain "
                          "change) for density updates during the analysis.\n";
        }
    }
    return res;
}

// Mass alone on the left-hand side: a diagonal M^-1 solve. Damping enters the
// residual at the lagged half-step velocity set on the nodes.
int CentralDifferenceLadruno::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    theEle->addMtoTang();
    return 0;
}

int CentralDifferenceLadruno::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    theDof->addMtoTang();
    return 0;
}

// Reported critical step: the damped limit when damping reduces it, otherwise
// 2/omega_max; -1 when no element contributed.
static inline double cdl_reported(double damped, double undamped)
{
    const double inf = std::numeric_limits<double>::infinity();
    if (undamped == inf || undamped <= 0.0)
        return -1.0;
    if (damped != inf && damped > 0.0 && damped < undamped)
        return damped;
    return undamped;
}

int CentralDifferenceLadruno::domainChanged()
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    if (theModel == 0 || theLinSOE == 0) {
        opserr << "CentralDifferenceLadruno::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    if (size == 0) {
        opserr << "CentralDifferenceLadruno::domainChanged - invalid size\n";
        return -1;
    }

    if (Ut == 0 || Ut->Size() != size) {
        if (Ut    != 0) delete Ut;
        if (Vhalf != 0) delete Vhalf;
        if (Aprev != 0) delete Aprev;
        if (Vfull != 0) delete Vfull;
        if (Azero != 0) delete Azero;

        Ut    = new Vector(size);
        Vhalf = new Vector(size);
        Aprev = new Vector(size);
        Vfull = new Vector(size);
        Azero = new Vector(size);

        if (Ut == 0 || Ut->Size() != size ||
            Vhalf == 0 || Vhalf->Size() != size ||
            Aprev == 0 || Aprev->Size() != size ||
            Vfull == 0 || Vfull->Size() != size ||
            Azero == 0 || Azero->Size() != size) {
            opserr << "CentralDifferenceLadruno::domainChanged - out of memory\n";
            return -1;
        }
    }
    Azero->Zero();

    // u_n, v (provisionally v_0) and a (provisionally the committed a_0) from
    // the committed DOF state; the starting procedure runs on the first
    // newStep().
    this->seedFromCommitted(theModel);

    firstStep = true;
    cflStepCount = 0;

    // the kinetic-energy baseline of the previous mesh does not apply
    prevKE = 0.0;
    committedPrevKE = 0.0;

    // any change of mesh, mass or system size passes through here
    massTangentValid = false;

    // dt_cr is computed here in every case: the element eigensolves need
    // neither dt nor the factored system, so criticalTimeStep is valid before
    // the first step. compute_critical_timestep controls only the per-step
    // report, -cflAbort and -recompute.
    CTSResult r = computeCriticalTimeStep(theModel, cflUseTangent, lumping);
    damped_minimum_critical_timestep   = r.damped_dt;
    undamped_minimum_critical_timestep = r.undamped_dt;
    damped_critical_element_tag        = r.damped_tag;
    undamped_critical_element_tag      = r.undamped_tag;
    if (compute_critical_timestep == 2)
        compute_critical_timestep = 1;

    const double inf = std::numeric_limits<double>::infinity();
    const double dt_cr = cdl_reported(damped_minimum_critical_timestep,
                                      undamped_minimum_critical_timestep);

    if (!betaKWarned &&
        damped_minimum_critical_timestep != inf &&
        undamped_minimum_critical_timestep != inf &&
        undamped_minimum_critical_timestep > 0.0 &&
        damped_minimum_critical_timestep <
            CDL_BETAK_WARN_RATIO * undamped_minimum_critical_timestep) {
        opserr << "CentralDifferenceLadruno - NOTE damping is significantly reducing the "
                  "critical time step (damped dt_cr = " << damped_minimum_critical_timestep
               << " vs undamped " << undamped_minimum_critical_timestep << ").\n"
                  "  In an explicit scheme stiffness-proportional (betaK) Rayleigh damping "
                  "reduces dt_cr to about 2/(betaK*omega_max^2); mass-proportional alphaM "
                  "(xi = alphaM/2omega) is benign at high frequency. Reporting the damped dt_cr.\n";
        betaKWarned = true;
    }

    if (compute_critical_timestep || verbose) {
        opserr << "CentralDifferenceLadruno: critical time step estimate"
               << (cflUseTangent ? " (tangent)" : "") << "\n"
               << "  central-difference limit (2/omega_max): "
               << undamped_minimum_critical_timestep
               << " @ element #" << undamped_critical_element_tag << "\n";
        if (damped_minimum_critical_timestep < undamped_minimum_critical_timestep)
            opserr << "  damped limit (reported): " << damped_minimum_critical_timestep
                   << " @ element #" << damped_critical_element_tag << "\n";
        if (dt_cr <= 0.0)
            opserr << "  note: dt_cr not applicable (no element produced a finite, "
                      "positive estimate, e.g. a model with nodal mass only).\n";
    }

    return 0;
}

void CentralDifferenceLadruno::seedFromCommitted(AnalysisModel *theModel)
{
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) (*Ut)(loc) = disp(i);
        }
        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) (*Vhalf)(loc) = vel(i);
        }
        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) (*Aprev)(loc) = accel(i);
        }
    }
}

// The analysis has already called Domain::revertToLastCommit() (trial state,
// time, loads); only the integrator state is restored here and the time is
// not touched. The committed DOF values were not changed by the failed step.
int CentralDifferenceLadruno::revertToLastStep(void)
{
    if (Ut == 0)                                  // domainChanged never ran
        return 0;
    massTangentValid = false;
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)
        return 0;

    this->seedFromCommitted(theModel);

    // the next newStep(dt_retry) rebuilds v_{-1/2} = v_n - dt_retry/2 a_n;
    // this also repairs a starting procedure that failed after changing Vhalf
    firstStep = true;
    updateCount = 0;
    prevKE = committedPrevKE;
    cflStepCount = committedCflStepCount;
    return 0;
}

int CentralDifferenceLadruno::newStep(double _deltaT)
{
    deltaT = _deltaT;
    if (deltaT <= 0.0) {
        opserr << "CentralDifferenceLadruno::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -1;
    }

    if (Ut == 0 || Vhalf == 0 || Aprev == 0) {
        opserr << "CentralDifferenceLadruno::newStep() - domainChanged() failed or not called\n";
        return -2;
    }

    // scalar state for revertToLastStep (the vector state restarts from the
    // committed DOF values)
    committedPrevKE       = prevKE;
    committedCflStepCount = cflStepCount;

    updateCount = 0;

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    // -recompute N: refresh dt_cr from the tangent stiffness
    if (cflRecomputeEvery > 0 && compute_critical_timestep == 2) {
        if (++cflStepCount >= cflRecomputeEvery) {
            compute_critical_timestep = 1;
            cflStepCount = 0;
        }
    }
    if (compute_critical_timestep == 1) {
        CTSResult r = computeCriticalTimeStep(theModel, cflUseTangent, lumping);
        damped_minimum_critical_timestep   = r.damped_dt;
        undamped_minimum_critical_timestep = r.undamped_dt;
        damped_critical_element_tag        = r.damped_tag;
        undamped_critical_element_tag      = r.undamped_tag;
        compute_critical_timestep = 2;

        const double dt_cr = cdl_reported(damped_minimum_critical_timestep,
                                          undamped_minimum_critical_timestep);
        if (cflFirstComputation || verbose) {
            opserr << "CentralDifferenceLadruno: critical time step (2/omega_max)"
                   << (cflUseTangent ? " (tangent)" : "")
                   << ": " << dt_cr
                   << " @ element #"
                   << ((dt_cr == damped_minimum_critical_timestep)
                           ? damped_critical_element_tag
                           : undamped_critical_element_tag)
                   << "\n";
            if (dt_cr > 0.0 && deltaT > dt_cr)
                opserr << "  WARNING dt = " << deltaT
                       << " exceeds the central-difference stability limit " << dt_cr
                       << " - expect INSTABILITY.\n";
        }
        cflFirstComputation = false;
    }

    if (cflAbort) {
        const double dt_cr = cdl_reported(damped_minimum_critical_timestep,
                                          undamped_minimum_critical_timestep);
        if (dt_cr > 0.0 && deltaT > dt_cr) {
            opserr << "CentralDifferenceLadruno::newStep() - ABORT: dt = " << deltaT
                   << " > central-difference stability limit " << dt_cr << "\n";
            return -2;
        }
    }

    // Starting procedure: initial acceleration at the committed state and
    // v_{-1/2} = v_0 - dt/2 a_0. domainChanged() left Vhalf = v_0 and
    // Ut = u_0. The system is formed and factored here.
    if (firstStep) {
        firstStep = false;

        // residual at (u_0, v_0, a = 0): r = P_0 - C v_0 - F_int(u_0)
        theModel->setResponse(*Ut, *Vhalf, *Azero);
        double t0 = theModel->getCurrentDomainTime();
        if (theModel->updateDomain(t0, deltaT) < 0) {
            opserr << "CentralDifferenceLadruno::newStep() - starter: failed to update domain\n";
            return -3;
        }
        if (this->formTangent(CURRENT_TANGENT) < 0) {
            opserr << "CentralDifferenceLadruno::newStep() - starter: formTangent failed\n";
            return -3;
        }
        if (this->formUnbalance() < 0) {
            opserr << "CentralDifferenceLadruno::newStep() - starter: formUnbalance failed\n";
            return -3;
        }
        if (theLinSOE->solve() < 0) {
            opserr << "CentralDifferenceLadruno::newStep() - starter: SOE solve failed\n";
            return -3;
        }
        *Aprev = theLinSOE->getX();                     // a_0
        Vhalf->addVector(1.0, *Aprev, -0.5 * deltaT);   // v_{-1/2} = v_0 - dt/2 a_0
    }

    // leap-frog: v_{n+1/2} = v_{n-1/2} + dt a_n ;  u_{n+1} = u_n + dt v_{n+1/2}
    Vhalf->addVector(1.0, *Aprev, deltaT);
    Ut->addVector(1.0, *Vhalf, deltaT);

    // trial state for the solve of a_{n+1}: zero acceleration (no inertia
    // term in the residual), damping at the lagged v_{n+1/2}
    theModel->setResponse(*Ut, *Vhalf, *Azero);

    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0) {
        opserr << "CentralDifferenceLadruno::newStep() - failed to update the domain\n";
        return -3;
    }

    return 0;
}

int CentralDifferenceLadruno::update(const Vector &U)
{
    updateCount++;
    if (updateCount > 1) {
        opserr << "WARNING CentralDifferenceLadruno::update() - called more than once -";
        opserr << " the central-difference scheme requires a LINEAR solution algorithm "
                  "(exactly one solve/step)\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING CentralDifferenceLadruno::update() - no AnalysisModel set\n";
        return -2;
    }
    if (Ut == 0) {
        opserr << "WARNING CentralDifferenceLadruno::update() - domainChanged() failed or not called\n";
        return -3;
    }
    if (U.Size() != Ut->Size()) {
        opserr << "WARNING CentralDifferenceLadruno::update() - Vectors of incompatible size ";
        opserr << " expecting " << Ut->Size() << " obtained " << U.Size() << endln;
        return -4;
    }

    // U is the solved acceleration a_{n+1}. A NaN or Inf entry stops the
    // step before it propagates (checked entry by entry: pNorm(0) skips NaN).
    const double A_max = U.pNorm(0);   // verbose report only
    if (!vectorIsFinite(U)) {
        opserr << "CentralDifferenceLadruno::update() - ABORT: non-finite acceleration "
                  "(NaN/Inf) - the integration has diverged (dt likely too large).\n";
        return -5;
    }

    // full-step velocity v_{n+1} = v_{n+1/2} + dt/2 a_{n+1} (the centred
    // difference (u_{n+2} - u_n)/(2 dt)); Vhalf holds v_{n+1/2}
    *Vfull = *Vhalf;
    Vfull->addVector(1.0, U, 0.5 * deltaT);

    // Full-step state (u_{n+1}, v_{n+1}, a_{n+1}) to the nodes. This is
    // needed also under -commitSolveState: committed nodal velocities and
    // recorders read it. Only the element state determination below may be
    // skipped.
    theModel->setResponse(*Ut, *Vfull, U);
    if (!commitSolveState) {
        // second element state determination at the same displacement (only
        // v and a changed); identical state for rate-independent materials
        if (theModel->updateDomain() < 0) {
            opserr << "CentralDifferenceLadruno::update() - failed to update the domain\n";
            return -6;
        }
    } else if (!commitSolveStateNoted) {
        commitSolveStateNoted = true;
        opserr << "CentralDifferenceLadruno: -commitSolveState - element state is committed "
                  "at the solve (half-step) velocity v_{n+1/2}. Identical for "
                  "rate-independent materials; rate-dependent materials (viscous materials "
                  "in zeroLength, TwoNodeLink, truss and multiple-spring elements) commit a "
                  "v_{n+1/2}-rate state.\n";
    }

    // Optional stop on kinetic-energy growth. The baseline is the running
    // maximum of the proxy, not the previous step: in free vibration the
    // kinetic energy passes near zero at every displacement extreme, and a
    // previous-step baseline would trip on the regrowth from that trough.
    if (divergenceFactor > 0.0) {
        double ke = 0.5 * ((*Vfull) ^ (*Vfull));
        if (prevKE > 0.0 && ke > divergenceFactor * prevKE) {
            opserr << "CentralDifferenceLadruno::update() - ABORT: kinetic-energy proxy grew to "
                   << (ke / prevKE) << "x its running maximum (> " << divergenceFactor
                   << ") - spurious energy growth / instability.\n";
            return -7;
        }
        if (ke > prevKE) prevKE = ke;
    }

    if (verbose)
        opserr << "CentralDifferenceLadruno::update() max|a| = " << A_max << endln;

    // a_{n+1} drives the velocity advance of the next newStep()
    *Aprev = U;

    return 0;
}

int CentralDifferenceLadruno::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING CentralDifferenceLadruno::commit() - no AnalysisModel set\n";
        return -1;
    }

    // time was advanced in newStep()
    return theModel->commitDomain();
}

const Vector &CentralDifferenceLadruno::getVel(void)
{
    if (Vhalf == 0) {
        // queried before domainChanged() allocated the state
        opserr << "WARNING CentralDifferenceLadruno::getVel() - called before "
                  "domainChanged(); returning an empty vector\n";
        static Vector emptyVel(0);
        return emptyVel;
    }
    return *Vhalf;
}

double CentralDifferenceLadruno::getCriticalTimeStep(void) const
{
    return cdl_reported(damped_minimum_critical_timestep,
                        undamped_minimum_critical_timestep);
}

int CentralDifferenceLadruno::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(7);
    data(0) = (double)compute_critical_timestep;
    data(1) = verbose ? 1.0 : 0.0;
    data(2) = cflAbort ? 1.0 : 0.0;
    data(3) = divergenceFactor;
    data(4) = cflUseTangent ? 1.0 : 0.0;
    data(5) = (double)cflRecomputeEvery;
    data(6) = (double)(int)lumping;   // CTSLumping code

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "WARNING CentralDifferenceLadruno::sendSelf() - could not send data\n";
        return -1;
    }
    return 0;
}

int CentralDifferenceLadruno::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(7);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "WARNING CentralDifferenceLadruno::recvSelf() - could not receive data\n";
        return -1;
    }
    compute_critical_timestep = (int)data(0);
    verbose          = (data(1) != 0.0);
    cflAbort         = (data(2) != 0.0);
    divergenceFactor = data(3);
    cflUseTangent    = (data(4) != 0.0);
    cflRecomputeEvery = (int)data(5);
    int lc = (int)data(6);
    lumping = (lc == 0) ? CTSLumping::RowSum
            : (lc == 2) ? CTSLumping::HRZ
                        : CTSLumping::Diagonal;
    if (lc < 0 || lc > 2)
        opserr << "WARNING CentralDifferenceLadruno::recvSelf - unknown lumping "
               << "code " << lc << " (using Diagonal)\n";

    return 0;
}

void CentralDifferenceLadruno::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    s << "CentralDifferenceLadruno - explicit leap-frog central difference\n";
    if (theModel != 0)
        s << "  currentTime: " << theModel->getCurrentDomainTime() << endln;
    s << "  recipe: lumped/diagonal mass, system Diagonal, algorithm Linear, dt < dt_cr\n";
    if (compute_critical_timestep > 0) {
        s << "  critical dt (2/omega_max): " << undamped_minimum_critical_timestep << "\n";
        if (damped_minimum_critical_timestep < undamped_minimum_critical_timestep)
            s << "  critical dt (damped, reported): " << damped_minimum_critical_timestep << "\n";
    }
}
