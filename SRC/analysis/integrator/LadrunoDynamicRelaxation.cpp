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
// Description: see LadrunoDynamicRelaxation.h.

#include <LadrunoDynamicRelaxation.h>
#include <LadrunoFictitiousMass.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Element.h>
#include <classTags.h>
#include <elementAPI.h>
#include <math.h>
#include <string.h>
#include <cmath>

// true when every entry is finite. Vector::pNorm(0) skips NaN entries in its
// max comparison, so it cannot detect a NaN.
static bool
ladrunoDRVectorIsFinite(const Vector &v)
{
    for (int i = 0; i < v.Size(); i++)
        if (!std::isfinite(v(i)))
            return false;
    return true;
}

// Read an optional trailing number. On failure the argument is left in
// place for the option loop (the Python reader advances past a token it
// cannot convert, the Tcl reader does not).
static bool
ladrunoDROptionalDouble(double &value)
{
    int before = OPS_GetNumRemainingInputArgs();
    if (before < 1)
        return false;
    int nd = 1;
    double v;
    if (OPS_GetDoubleInput(&nd, &v) == 0) {
        value = v;
        return true;
    }
    if (OPS_GetNumRemainingInputArgs() < before)
        OPS_ResetCurrentInputArg(-1);
    return false;
}

void *OPS_LadrunoDynamicRelaxation(void)
{
    // Usage: integrator LadrunoDynamicRelaxation
    //          <-mass gershgorin|lumped $scale|unity> <-massSafety $f>
    //          <-marginEvery $N> <-dt $dt>
    //          <-recompute $N> <-damping kinetic|viscous <$zeta>>
    //          <-noAutoRefresh> <-interp> <-divergence $f> <-verbose>
    int    massMode = 0;            // gershgorin
    double dtPseudo = 1.0, massScale = 1.0, divergenceFactor = 0.0;
    int    recomputeEvery = 0;
    bool   interp = false, verbose = false;
    int    dampMode = 0;            // 0 kinetic (default) | 1 viscous-critical
    double zetaTarget = 1.0;        // viscous damping ratio (1 = critical)
    bool   autoRefresh = true;      // rebuild M* at KE peaks (gershgorin)
    double massSafety = LADRUNO_DR_DEFAULT_MASS_SAFETY;   // omega_max*dt <= 2f
    bool   massSafetySet = false;   // explicit -massSafety (for the mode warning)
    int    marginEvery = LADRUNO_DR_DEFAULT_MARGIN_EVERY;  // margin-probe cadence

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-mass") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING LadrunoDynamicRelaxation -mass expects gershgorin|lumped|unity\n";
                return 0;
            }
            const char *m = OPS_GetString();
            if (strcmp(m, "gershgorin") == 0)      massMode = 0;
            else if (strcmp(m, "lumped") == 0) {
                massMode = 1;
                ladrunoDROptionalDouble(massScale);
            } else if (strcmp(m, "unity") == 0)    massMode = 2;
            else opserr << "WARNING LadrunoDynamicRelaxation - unknown -mass " << m
                        << " (keeping gershgorin)\n";
        } else if (strcmp(opt, "-massSafety") == 0) {
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &massSafety) < 0) {
                opserr << "WARNING LadrunoDynamicRelaxation -massSafety failed to read f\n";
                return 0;
            }
            massSafetySet = true;
        } else if (strcmp(opt, "-marginEvery") == 0) {
            int nd = 1;
            if (OPS_GetIntInput(&nd, &marginEvery) < 0) {
                opserr << "WARNING LadrunoDynamicRelaxation -marginEvery failed to read N\n";
                return 0;
            }
            if (marginEvery < 0) marginEvery = 0;
        } else if (strcmp(opt, "-dt") == 0) {
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &dtPseudo) < 0) {
                opserr << "WARNING LadrunoDynamicRelaxation -dt failed to read value\n";
                return 0;
            }
        } else if (strcmp(opt, "-recompute") == 0) {
            int nd = 1;
            if (OPS_GetIntInput(&nd, &recomputeEvery) < 0) {
                opserr << "WARNING LadrunoDynamicRelaxation -recompute failed to read N\n";
                return 0;
            }
        } else if (strcmp(opt, "-damping") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING LadrunoDynamicRelaxation -damping expects kinetic|viscous\n";
                return 0;
            }
            const char *d = OPS_GetString();
            if (strcmp(d, "kinetic") == 0)      dampMode = 0;
            else if (strcmp(d, "viscous") == 0) {
                dampMode = 1;
                ladrunoDROptionalDouble(zetaTarget);   // optional zeta
            } else opserr << "WARNING LadrunoDynamicRelaxation - unknown -damping " << d
                          << " (keeping kinetic)\n";
        } else if (strcmp(opt, "-noAutoRefresh") == 0) {
            autoRefresh = false;
        } else if (strcmp(opt, "-interp") == 0) {
            interp = true;
        } else if (strcmp(opt, "-divergence") == 0) {
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &divergenceFactor) < 0) divergenceFactor = 0.0;
        } else if (strcmp(opt, "-verbose") == 0) {
            verbose = true;
        } else {
            opserr << "WARNING LadrunoDynamicRelaxation - unknown option " << opt
                   << " (ignored)\n";
        }
    }

    // -massSafety is a safety factor in (0,1]. f = 1 is the central-difference
    // stability limit itself, accepted with a warning.
    if (!(massSafety > 0.0)) {
        opserr << "WARNING LadrunoDynamicRelaxation -massSafety " << massSafety
               << " is not positive; using the default "
               << LADRUNO_DR_DEFAULT_MASS_SAFETY << "\n";
        massSafety = LADRUNO_DR_DEFAULT_MASS_SAFETY;
    } else if (massSafety > 1.0) {
        opserr << "WARNING LadrunoDynamicRelaxation -massSafety " << massSafety
               << " > 1 is ABOVE the central-difference stability limit "
                  "(omega*dt > 2); clamping to 1\n";
        massSafety = 1.0;
    }
    // massSafety scales only the gershgorin mass
    if (massSafety >= 1.0 && massMode == 0)
        opserr << "WARNING LadrunoDynamicRelaxation -massSafety 1 marches EXACTLY "
                  "ON the explicit stability boundary (omega_max*dt = 2), where the "
                  "march amplifies round-off instead of damping it. A long relaxation "
                  "can walk away from an exact equilibrium and report NO error. "
                  "Prefer the default "
                  "-massSafety " << LADRUNO_DR_DEFAULT_MASS_SAFETY << ".\n";
    if (massSafetySet && massMode != 0)
        opserr << "WARNING LadrunoDynamicRelaxation -massSafety only scales the "
                  "gershgorin mass; with -mass lumped/unity it is IGNORED (use the "
                  "-mass lumped $scale factor instead)\n";

    return new LadrunoDynamicRelaxation(massMode, dtPseudo, massScale,
                                        recomputeEvery, interp, divergenceFactor,
                                        verbose, dampMode, zetaTarget, autoRefresh,
                                        massSafety, marginEvery);
}

LadrunoDynamicRelaxation::LadrunoDynamicRelaxation(int massMode_, double dtPseudo_,
                                                   double massScale_, int recomputeEvery_,
                                                   bool interp_, double divergenceFactor_,
                                                   bool verbose_, int dampMode_,
                                                   double zetaTarget_, bool autoRefresh_,
                                                   double massSafety_, int marginEvery_)
 :TransientIntegrator(INTEGRATOR_TAGS_LadrunoDynamicRelaxation),
  massMode(massMode_), dtPseudo(dtPseudo_), massScale(massScale_),
  recomputeEvery(recomputeEvery_), interp(interp_),
  divergenceFactor(divergenceFactor_), verbose(verbose_),
  Mstar(0), Mprobe(0), Ut(0), Vhalf(0), Aprev(0), Vfull(0), Azero(0),
  firstStep(true), updateCount(0), size(0), stepCount(0),
  prevKE(0.0), kineticEnergy(0.0), residualNorm(0.0), deltaT(0.0),
  dampMode(dampMode_), zetaTarget(zetaTarget_), cVisc(0.0), autoRefresh(autoRefresh_),
  massSafety(massSafety_), stabMargin(-1.0), marginWarned(false),
  marginEvery(marginEvery_ < 0 ? 0 : marginEvery_), marginCount(0),
  haveBuiltOnce(false),
  UtC(0), VhalfC(0), AprevC(0), MstarC(0),
  snapFirstStep(true), snapStepCount(0), snapPrevKE(0.0), snapKE(0.0),
  snapRes(0.0), snapCVisc(0.0), haveSnap(false)
{
    // guard the direct-construction path too (the parser validates separately)
    if (!(massSafety > 0.0) || massSafety > 1.0)
        massSafety = LADRUNO_DR_DEFAULT_MASS_SAFETY;
}

LadrunoDynamicRelaxation::~LadrunoDynamicRelaxation()
{
    if (Mstar != 0) delete Mstar;
    if (Mprobe != 0) delete Mprobe;
    if (Ut    != 0) delete Ut;
    if (Vhalf != 0) delete Vhalf;
    if (Aprev != 0) delete Aprev;
    if (Vfull != 0) delete Vfull;
    if (Azero != 0) delete Azero;
    if (UtC    != 0) delete UtC;
    if (VhalfC != 0) delete VhalfC;
    if (AprevC != 0) delete AprevC;
    if (MstarC != 0) delete MstarC;
}

// The SOE matrix is M* alone (zeroed first), so the solve is a diagonal
// inverse. The element and nodal tangents are not assembled.
int
LadrunoDynamicRelaxation::formTangent(int statFlag)
{
    LinearSOE *theSOE = this->getLinearSOE();
    if (theSOE == 0 || Mstar == 0) {
        opserr << "WARNING LadrunoDynamicRelaxation::formTangent() - no SOE or M* not built\n";
        return -1;
    }
    Ladruno::addDiagonalToSOE(theSOE, *Mstar, size, 1.0, /*zeroFirst=*/true);
    return 0;
}

int LadrunoDynamicRelaxation::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();   // unused: formTangent() assembles M* directly
    return 0;
}

int LadrunoDynamicRelaxation::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();   // unused: formTangent() assembles M* directly
    return 0;
}

// gershgorin: m_i = dt^2/(4 f^2) sum_j |K_ij|, omega_max dt <= 2f
// lumped:     row sums of the element mass times massScale
// unity:      m_i = 1
int
LadrunoDynamicRelaxation::buildFictitiousMass(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0 || Mstar == 0) return -1;

    // keep the mass in use, so that a rebuild is also a margin sample
    bool hadPrevMass = (massMode == 0 && haveBuiltOnce && Mprobe != 0
                        && Mprobe->Size() == size);
    if (hadPrevMass) *Mprobe = *Mstar;

    // gershgorin prefactor dt^2/(4 massSafety^2), so omega_max dt <= 2 massSafety.
    // With viscous damping of ratio z the explicit stability limit of the damped
    // step is reduced by s = sqrt(1+z^2) - z, so M* is divided by s^2 as well.
    // The lumped and unity masses ignore the prefactor (and massSafety).
    double f = 0.25 * dtPseudo * dtPseudo / (massSafety * massSafety);
    if (dampMode == 1) {
        double s = sqrt(1.0 + zetaTarget * zetaTarget) - zetaTarget;
        if (s > 0.0) f /= (s * s);
    }
    if (Ladruno::buildGershgorinDiagonal(theModel, *Mstar, massMode, massScale, f) < 0)
        return -1;

    // C* = cVisc M*. After the scaling above omega dt is about 2 s fSafe, and
    // the damping coefficient for ratio zeta of that frequency is
    // cVisc = 2 zeta omega = 4 zeta s fSafe / dt, which keeps dt cVisc < 2
    // (0.83 at zeta = 1, fSafe = 0.5). For the lumped and unity masses the
    // frequency is not known and fSafe = 1.
    if (dampMode == 1) {
        double s = sqrt(1.0 + zetaTarget * zetaTarget) - zetaTarget;
        double fSafe = (massMode == 0) ? massSafety : 1.0;
        cVisc = 4.0 * zetaTarget * s * fSafe / dtPseudo;
    }

    // on a rebuild, M* now describes the current tangent and Mprobe the mass in
    // use until now: compare them (one pass over a vector)
    if (massMode == 0 && haveBuiltOnce && hadPrevMass)
        this->applyMargin(*Mstar, *Mprobe);

    haveBuiltOnce = true;
    marginCount = 0;
    if (massMode == 0) {
        if (stabMargin < 0.0) stabMargin = massSafety * massSafety;
    } else {
        stabMargin = -1.0;    // no Gershgorin bound for lumped or unity
    }

    return 0;
}

// margin_i = live_i f^2 / march_i, with both vectors built with the same
// prefactor (see measureStabilityMargin()).
void
LadrunoDynamicRelaxation::applyMargin(const Vector &live, const Vector &march)
{
    double mmax = 0.0;
    for (int i = 0; i < size; i++)
        if (live(i) > mmax) mmax = live(i);
    // skip DOFs on the 1e-8 floor of the builder: they carry no stiffness, and
    // the floored value changes between samples
    double cut = mmax * 1.0e-6;
    double worst = 0.0;
    for (int i = 0; i < size; i++) {
        double m = march(i);
        if (!(m > 0.0) || !(live(i) > cut)) continue;
        double r = live(i) * massSafety * massSafety / m;
        if (r > worst) worst = r;
    }
    if (!(worst > 0.0)) return;

    // keep the largest value since domainChanged(): a rebuild resets the
    // current margin to massSafety^2, which would hide an earlier excursion
    if (worst > stabMargin) stabMargin = worst;

    if (worst > 1.0 + 1.0e-6 && !marginWarned) {
        marginWarned = true;
        opserr << "WARNING LadrunoDynamicRelaxation - the fictitious mass NO LONGER "
                  "BOUNDS the explicit step: the current tangent gives a Gershgorin "
                  "stability margin of " << worst << " (> 1 means omega*dt > 2). "
                  "The march is amplifying round-off, not damping it, and no error "
                  "is reported when it does: the relaxed state can be wrong. "
                  "Lower -massSafety (currently " << massSafety
               << ", try " << 0.5 * massSafety << ") and/or refresh M* more often "
                  "(-recompute N / drop -noAutoRefresh). Printed once per integrator.\n";
    }
}

// Stability-margin probe. The row sums s of the current tangent go into
// Mprobe with the prefactor M* was built with, m_probe = dt^2/(4 f^2) s, so
//   margin_i = s_i dt^2 / (4 M*_i) = m_probe_i f^2 / M*_i
// (the viscous 1/s^2 factor cancels with the reduced limit), which is
// (omega_max dt / 2)^2 for the mass in use: massSafety^2 while the tangent is
// unchanged, above 1 once the model has stiffened past M*. M* is not
// changed, and the probe runs on its own cadence so that it also measures
// under -noAutoRefresh without -recompute.
int
LadrunoDynamicRelaxation::measureStabilityMargin(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0 || Mstar == 0 || Mprobe == 0 || massMode != 0 || !haveBuiltOnce)
        return -1;

    double f = 0.25 * dtPseudo * dtPseudo / (massSafety * massSafety);
    if (dampMode == 1) {
        double s = sqrt(1.0 + zetaTarget * zetaTarget) - zetaTarget;
        if (s > 0.0) f /= (s * s);
    }
    if (Ladruno::buildGershgorinDiagonal(theModel, *Mprobe, 0, massScale, f) < 0)
        return -1;

    this->applyMargin(*Mprobe, *Mstar);
    return 0;
}

int
LadrunoDynamicRelaxation::domainChanged(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "LadrunoDynamicRelaxation::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    size = x.Size();
    if (size == 0) {
        opserr << "LadrunoDynamicRelaxation::domainChanged - invalid size\n";
        return -1;
    }

    if (Ut == 0 || Ut->Size() != size) {
        if (Ut    != 0) delete Ut;
        if (Vhalf != 0) delete Vhalf;
        if (Aprev != 0) delete Aprev;
        if (Vfull != 0) delete Vfull;
        if (Azero != 0) delete Azero;
        if (Mstar != 0) delete Mstar;
        if (Mprobe != 0) delete Mprobe;
        Ut    = new Vector(size);
        Vhalf = new Vector(size);
        Aprev = new Vector(size);
        Vfull = new Vector(size);
        Azero = new Vector(size);
        Mstar = new Vector(size);
        Mprobe = new Vector(size);
        if (Ut == 0 || Vhalf == 0 || Aprev == 0 || Vfull == 0 || Azero == 0 ||
            Mstar == 0 || Mprobe == 0 || Mstar->Size() != size) {
            opserr << "LadrunoDynamicRelaxation::domainChanged - out of memory\n";
            return -1;
        }
    }
    Azero->Zero();

    // seed u_n, v_{n-1/2} = v_0 and a_n = a_0 from the committed DOF state;
    // the starter a_0, v_{-1/2} is computed in the first newStep()
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();
        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; i++) { int loc = id(i); if (loc >= 0) (*Ut)(loc) = disp(i); }
        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; i++) { int loc = id(i); if (loc >= 0) (*Vhalf)(loc) = vel(i); }
        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; i++) { int loc = id(i); if (loc >= 0) (*Aprev)(loc) = accel(i); }
    }

    firstStep = true;
    stepCount = 0;
    prevKE = 0.0;
    // the equations may have been renumbered: restart the margin record
    stabMargin = -1.0;
    haveBuiltOnce = false;
    marginCount = 0;

    if (this->buildFictitiousMass() < 0) {
        opserr << "LadrunoDynamicRelaxation::domainChanged - failed to build fictitious mass\n";
        return -1;
    }
    // the seeded state is the one a failed first step returns to
    this->takeSnapshot();
    if (verbose) {
        double mn = (*Mstar)(0), mx = (*Mstar)(0);
        for (int i = 1; i < size; i++) { if ((*Mstar)(i) < mn) mn = (*Mstar)(i);
                                         if ((*Mstar)(i) > mx) mx = (*Mstar)(i); }
        opserr << "LadrunoDynamicRelaxation: fictitious M* (mode " << massMode
               << ") range [" << mn << ", " << mx << "], dt = " << dtPseudo
               << ", massSafety = " << massSafety << " (omega*dt <= "
               << 2.0 * massSafety << ")" << endln;
    }
    return 0;
}

int
LadrunoDynamicRelaxation::newStep(double _deltaT)
{
    deltaT = _deltaT;
    if (deltaT <= 0.0) {
        opserr << "LadrunoDynamicRelaxation::newStep() - dt = " << deltaT << " <= 0\n";
        return -1;
    }
    if (Ut == 0 || Vhalf == 0 || Aprev == 0 || Mstar == 0) {
        opserr << "LadrunoDynamicRelaxation::newStep() - domainChanged() not called\n";
        return -2;
    }
    updateCount = 0;

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    // stability-margin probe, before and independent of any M* rebuild
    if (marginEvery > 0 && massMode == 0 && haveBuiltOnce && !firstStep) {
        if (++marginCount >= marginEvery) {
            marginCount = 0;
            this->measureStabilityMargin();
        }
    }

    // -recompute N: rebuild M* from the current tangent and zero the velocities,
    // so that the rebuild injects no kinetic energy
    if (recomputeEvery > 0 && !firstStep) {
        if (++stepCount >= recomputeEvery) {
            stepCount = 0;
            this->buildFictitiousMass();
            Vhalf->Zero();
            prevKE = 0.0;
        }
    }

    // first step: a_0 = M*^-1 (P_0 - F_int(u_0)), v_{-1/2} = v_0 - dt/2 a_0
    if (firstStep) {
        firstStep = false;
        theModel->setResponse(*Ut, *Vhalf, *Azero);
        double t0 = theModel->getCurrentDomainTime();
        if (theModel->updateDomain(t0, deltaT) < 0) return -3;
        if (this->formTangent(CURRENT_TANGENT) < 0) return -3;
        if (this->formUnbalance() < 0) return -3;
        if (theLinSOE->solve() < 0) return -3;
        *Aprev = theLinSOE->getX();
        Vhalf->addVector(1.0, *Aprev, -0.5 * deltaT);
    }

    // leap-frog: v_{n+1/2} = v_{n-1/2} + dt a_n ; u_{n+1} = u_n + dt v_{n+1/2}
    Vhalf->addVector(1.0, *Aprev, deltaT);
    Ut->addVector(1.0, *Vhalf, deltaT);

    // trial state for the solve of a_{n+1}; with zero acceleration the residual
    // B = P - F_int is the static unbalance
    theModel->setResponse(*Ut, *Vhalf, *Azero);
    double time = theModel->getCurrentDomainTime() + deltaT;
    if (theModel->updateDomain(time, deltaT) < 0) {
        opserr << "LadrunoDynamicRelaxation::newStep() - failed to update domain\n";
        return -3;
    }
    return 0;
}

int
LadrunoDynamicRelaxation::update(const Vector &U)
{
    updateCount++;
    if (updateCount > 1) {
        opserr << "WARNING LadrunoDynamicRelaxation::update() - called more than once "
                  "in a step; use `algorithm Linear` (one solve per step)\n";
        return -1;
    }
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0 || Ut == 0) return -2;
    if (U.Size() != size) {
        opserr << "WARNING LadrunoDynamicRelaxation::update() - size mismatch\n";
        return -4;
    }

    // U is the solved acceleration a_{n+1}; stop on NaN or Inf
    const double A_max = U.pNorm(0);   // verbose output only
    if (!ladrunoDRVectorIsFinite(U)) {
        opserr << "LadrunoDynamicRelaxation::update() - ABORT: non-finite acceleration "
                  "(dt likely too large / M* too small)\n";
        return -5;
    }

    // static unbalance ||f_ext - f_int||_inf = ||M* a||_inf, before damping
    residualNorm = 0.0;
    for (int i = 0; i < size; i++) {
        double r = fabs((*Mstar)(i) * U(i));
        if (r > residualNorm) residualNorm = r;
    }

    // viscous damping: force -cVisc M* v_{n+1/2} at the known half-step
    // velocity, so the step stays explicit: a_work = a - cVisc v_{n+1/2}
    Vector aWork(U);
    if (dampMode == 1)
        for (int i = 0; i < size; i++) aWork(i) -= cVisc * (*Vhalf)(i);

    // full-step velocity v_{n+1} = v_{n+1/2} + dt/2 a_work
    *Vfull = *Vhalf;
    Vfull->addVector(1.0, aWork, 0.5 * deltaT);

    // kinetic energy 1/2 v^T M* v
    kineticEnergy = 0.0;
    for (int i = 0; i < size; i++) {
        double vi = (*Vfull)(i);
        kineticEnergy += (*Mstar)(i) * vi * vi;
    }
    kineticEnergy *= 0.5;

    if (dampMode == 0) {
        // kinetic damping (Cundall 1976): zero all velocities at each peak
        if (prevKE > 0.0 && kineticEnergy < prevKE) {
            Vhalf->Zero();
            Vfull->Zero();
            prevKE = 0.0;
            // rebuild M* from the current tangent while the velocities are
            // zero (no energy injected); follows the stiffening of the model,
            // e.g. on the far branch after snap-through
            if (autoRefresh && massMode == 0)
                this->buildFictitiousMass();
        } else {
            prevKE = kineticEnergy;
        }
    } else {
        prevKE = kineticEnergy;      // viscous: velocities are never zeroed
    }

    // push the (possibly damped) full-step snapshot to the nodes
    theModel->setResponse(*Ut, *Vfull, aWork);
    if (theModel->updateDomain() < 0) {
        opserr << "LadrunoDynamicRelaxation::update() - failed to update domain\n";
        return -6;
    }

    if (divergenceFactor > 0.0 && prevKE > 0.0) {
        if (kineticEnergy > divergenceFactor * prevKE) {
            opserr << "LadrunoDynamicRelaxation::update() - ABORT: KE grew "
                   << (kineticEnergy / prevKE) << "x in one step (> " << divergenceFactor << ")\n";
            return -7;
        }
    }
    if (verbose)
        opserr << "LadrunoDynamicRelaxation::update() max|a|=" << A_max
               << " RES=" << residualNorm << " KE=" << kineticEnergy << endln;

    // carry the (damped) accel so the next leap-frog half-step advance uses it
    *Aprev = aWork;
    return 0;
}

int
LadrunoDynamicRelaxation::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING LadrunoDynamicRelaxation::commit() - no AnalysisModel\n";
        return -1;
    }
    const int rc = theModel->commitDomain();
    if (rc == 0)
        this->takeSnapshot();
    return rc;
}

// The analysis has already called Domain::revertToLastCommit() (trial state
// and time); only the march state of the integrator is restored here. The
// stored state is used rather than a re-seed from the committed DOF state:
// a re-seed reruns the first-step starter, v_{-1/2} = v_0 - dt/2 a_0, which
// differs from the uninterrupted march after a kinetic-damping reset.
int
LadrunoDynamicRelaxation::revertToLastStep(void)
{
    if (Ut == 0)                                  // domainChanged never ran
        return 0;
    if (!haveSnap || UtC == 0 || UtC->Size() != size) {
        // no stored state for this numbering: seed from the committed DOF
        // state and rerun the starter
        AnalysisModel *theModel = this->getAnalysisModel();
        if (theModel == 0)
            return 0;
        DOF_GrpIter &theDOFs = theModel->getDOFs();
        DOF_Group *dofPtr;
        while ((dofPtr = theDOFs()) != 0) {
            const ID &id = dofPtr->getID();
            const Vector &d = dofPtr->getCommittedDisp();
            const Vector &v = dofPtr->getCommittedVel();
            const Vector &a = dofPtr->getCommittedAccel();
            for (int i = 0; i < id.Size(); i++) {
                int loc = id(i);
                if (loc >= 0) { (*Ut)(loc) = d(i); (*Vhalf)(loc) = v(i); (*Aprev)(loc) = a(i); }
            }
        }
        firstStep = true;
        updateCount = 0;
        prevKE = 0.0;
        return 0;
    }
    *Ut    = *UtC;
    *Vhalf = *VhalfC;
    *Aprev = *AprevC;
    *Mstar = *MstarC;
    firstStep     = snapFirstStep;
    stepCount     = snapStepCount;
    prevKE        = snapPrevKE;
    kineticEnergy = snapKE;
    residualNorm  = snapRes;
    cVisc         = snapCVisc;
    updateCount   = 0;
    return 0;
}

void
LadrunoDynamicRelaxation::takeSnapshot(void)
{
    if (Ut == 0 || Mstar == 0)
        return;
    if (UtC == 0 || UtC->Size() != size) {
        if (UtC    != 0) delete UtC;
        if (VhalfC != 0) delete VhalfC;
        if (AprevC != 0) delete AprevC;
        if (MstarC != 0) delete MstarC;
        UtC = new Vector(size); VhalfC = new Vector(size);
        AprevC = new Vector(size); MstarC = new Vector(size);
    }
    *UtC = *Ut; *VhalfC = *Vhalf; *AprevC = *Aprev; *MstarC = *Mstar;
    snapFirstStep = firstStep;
    snapStepCount = stepCount;
    snapPrevKE    = prevKE;
    snapKE        = kineticEnergy;
    snapRes       = residualNorm;
    snapCVisc     = cVisc;
    haveSnap      = true;
}


const Vector &
LadrunoDynamicRelaxation::getVel(void)
{
    return *Vhalf;
}

// max_i |f_ext - f_int|_i of the last step
double LadrunoDynamicRelaxation::getResidualNorm(void) const
{
    return residualNorm;
}

double LadrunoDynamicRelaxation::getKineticEnergy(void) const
{
    return kineticEnergy;
}

// see the header; -2 when the probe is off
double LadrunoDynamicRelaxation::getStabilityMargin(void) const
{
    if (massMode == 0 && marginEvery <= 0) return -2.0;
    return stabMargin;
}

int
LadrunoDynamicRelaxation::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(12);
    data(0) = (double)massMode;
    data(1) = dtPseudo;
    data(2) = massScale;
    data(3) = (double)recomputeEvery;
    data(4) = interp ? 1.0 : 0.0;
    data(5) = divergenceFactor;
    data(6) = verbose ? 1.0 : 0.0;
    data(7) = (double)dampMode;
    data(8) = zetaTarget;
    data(9) = autoRefresh ? 1.0 : 0.0;
    data(10) = massSafety;
    data(11) = (double)marginEvery;
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoDynamicRelaxation::sendSelf() - failed to send\n";
        return -1;
    }
    return 0;
}

int
LadrunoDynamicRelaxation::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(12);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoDynamicRelaxation::recvSelf() - failed to receive\n";
        return -1;
    }
    massMode        = (int)data(0);
    dtPseudo        = data(1);
    massScale       = data(2);
    recomputeEvery  = (int)data(3);
    interp          = data(4) != 0.0;
    divergenceFactor = data(5);
    verbose         = data(6) != 0.0;
    dampMode        = (int)data(7);
    zetaTarget      = data(8);
    autoRefresh     = data(9) != 0.0;
    massSafety      = data(10);
    marginEvery     = (int)data(11);
    if (marginEvery < 0) marginEvery = 0;
    // an out-of-range value means the two ends disagree on the data layout
    if (!(massSafety > 0.0) || massSafety > 1.0) {
        opserr << "WARNING LadrunoDynamicRelaxation::recvSelf() - received massSafety "
               << massSafety << " outside (0,1]; falling back to "
               << LADRUNO_DR_DEFAULT_MASS_SAFETY << ". This rank will NOT march with "
                  "the same fictitious mass as the sender -- suspect a mixed-build "
                  "actor pair (the data vector holds 12 doubles).\n";
        massSafety = LADRUNO_DR_DEFAULT_MASS_SAFETY;
    }
    stabMargin      = -1.0;      // per process; measured from this actor's first build
    marginWarned    = false;
    marginCount     = 0;
    haveBuiltOnce   = false;
    return 0;
}

void
LadrunoDynamicRelaxation::Print(OPS_Stream &s, int flag)
{
    s << "LadrunoDynamicRelaxation - quasi-static dynamic relaxation\n";
    s << "\t mass mode: " << (massMode == 0 ? "gershgorin" : massMode == 1 ? "lumped" : "unity")
      << "  dt: " << dtPseudo << "  recompute: " << recomputeEvery
      << "  autoRefresh: " << (autoRefresh ? "on" : "off") << endln;
    s << "\t massSafety: " << massSafety << " (omega*dt <= " << 2.0 * massSafety
      << ")  stability margin: " << this->getStabilityMargin()
      << "  probe every: " << marginEvery << " steps" << endln;
    s << "\t damping: " << (dampMode == 0 ? "kinetic" : "viscous")
      << "  zeta: " << zetaTarget << "  cVisc: " << cVisc << endln;
    s << "\t last RES: " << residualNorm << "  KE: " << kineticEnergy << endln;
}
