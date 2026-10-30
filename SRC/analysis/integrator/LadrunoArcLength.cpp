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
//   Crisfield, M.A. (1981). A fast incremental/iterative solution procedure
//     that handles "snap-through". Computers & Structures, 13(1-3), 55-62.
//     doi:10.1016/0045-7949(81)90108-5
//   Ramm, E. (1981). Strategies for tracing the nonlinear response near
//     limit points. In W. Wunderlich, E. Stein and K.-J. Bathe (eds.),
//     Nonlinear Finite Element Analysis in Structural Mechanics, Springer,
//     Berlin, pp. 63-89. doi:10.1007/978-3-642-81589-8_5
//   Dassault Systemes Simulia. Abaqus Analysis User's Guide, "Automatic
//     stabilization of unstable problems".
//
// Description: see LadrunoArcLength.h. The predictor and corrector are
// those of ArcLength (fmk, 07/98) without the DDM sensitivity terms.

#include <LadrunoArcLength.h>
#include <LadrunoFictitiousMass.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Element.h>
#include <math.h>
#include <string.h>
#include <elementAPI.h>
#include <classTags.h>

// Read an optional trailing number. On failure the argument is left in
// place for the option loop (the Python reader advances past a token it
// cannot convert, the Tcl reader does not).
static bool
ladrunoArcLengthOptionalDouble(double &value)
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

void* OPS_LadrunoArcLength()
{
    // integrator LadrunoArcLength arcLength alpha
    //     <-adapt Jd ellMin ellMax> <-p exp>
    //     <-stabilize <f>> <-adaptStab> <-cVisc c>
    //     <-mass gershgorin|lumped <scale>|unity>
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING integrator LadrunoArcLength arcLength alpha "
                  "<-adapt Jd ellMin ellMax> <-p exp> <-stabilize <f>> "
                  "<-adaptStab> <-cVisc c> <-mass gershgorin|lumped <scale>|unity>\n";
        return 0;
    }

    int numdata = 1;
    double arcLength = 0.0, alpha = 1.0;
    if (OPS_GetDoubleInput(&numdata, &arcLength) < 0) {
        opserr << "WARNING integrator LadrunoArcLength failed to read arcLength\n";
        return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &alpha) < 0) {
        opserr << "WARNING integrator LadrunoArcLength failed to read alpha\n";
        return 0;
    }

    bool   adapt = false;
    int    Jd = 5;
    double ellMin = 0.0, ellMax = 0.0, pExp = 1.0;
    bool   stabilize = false, adaptStab = false;
    double fTarget = 2.0e-4, cVisc = 0.0, dtPseudo = 1.0, massScale = 1.0;
    int    massMode = 0;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-adapt") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "WARNING integrator LadrunoArcLength -adapt expects "
                          "Jd ellMin ellMax\n";
                return 0;
            }
            int idata = 1;
            if (OPS_GetIntInput(&idata, &Jd) < 0) {
                opserr << "WARNING integrator LadrunoArcLength -adapt failed to read Jd\n";
                return 0;
            }
            int ndata = 2;
            double mm[2] = {0.0, 0.0};
            if (OPS_GetDoubleInput(&ndata, mm) < 0) {
                opserr << "WARNING integrator LadrunoArcLength -adapt failed to read "
                          "ellMin ellMax\n";
                return 0;
            }
            ellMin = mm[0];
            ellMax = mm[1];
            adapt = true;
        } else if (strcmp(opt, "-p") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING integrator LadrunoArcLength -p expects exp\n";
                return 0;
            }
            int idata = 1;
            if (OPS_GetDoubleInput(&idata, &pExp) < 0) {
                opserr << "WARNING integrator LadrunoArcLength -p failed to read exp\n";
                return 0;
            }
        } else if (strcmp(opt, "-stabilize") == 0) {
            stabilize = true;
            ladrunoArcLengthOptionalDouble(fTarget);   // optional energy fraction
        } else if (strcmp(opt, "-adaptStab") == 0) {
            adaptStab = true;
        } else if (strcmp(opt, "-cVisc") == 0) {
            int nd = 1;
            if (OPS_GetDoubleInput(&nd, &cVisc) < 0) {
                opserr << "WARNING integrator LadrunoArcLength -cVisc failed to read value\n";
                return 0;
            }
        } else if (strcmp(opt, "-mass") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING integrator LadrunoArcLength -mass expects gershgorin|lumped|unity\n";
                return 0;
            }
            const char *m = OPS_GetString();
            if (strcmp(m, "gershgorin") == 0)      massMode = 0;
            else if (strcmp(m, "lumped") == 0) {
                massMode = 1;
                ladrunoArcLengthOptionalDouble(massScale);
            } else if (strcmp(m, "unity") == 0)    massMode = 2;
            else opserr << "WARNING integrator LadrunoArcLength - unknown -mass " << m
                        << " (keeping gershgorin)\n";
        } else {
            opserr << "WARNING integrator LadrunoArcLength - unknown option " << opt
                   << " (ignored)\n";
        }
    }

    return new LadrunoArcLength(arcLength, alpha, adapt, Jd, ellMin, ellMax, pExp,
                                stabilize, fTarget, adaptStab, massMode, massScale,
                                cVisc, dtPseudo);
}

LadrunoArcLength::LadrunoArcLength(double arcLength, double alpha,
                                   bool adapt_, int Jd_,
                                   double ellMin_, double ellMax_, double pExp_,
                                   bool stabilize_, double fTarget_, bool adaptStab_,
                                   int massMode_, double massScale_, double cVisc_,
                                   double dtPseudo_)
 :StaticIntegrator(INTEGRATOR_TAGS_LadrunoArcLength),
  arcLength2(arcLength*arcLength), alpha2(alpha*alpha),
  a(0.0), b(0.0), c(0.0), b24ac(0.0),
  deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), phat(0),
  deltaLambdaStep(0.0), currentLambda(0.0), signLastDeltaLambdaStep(1),
  adapt(adapt_), Jd(Jd_ > 0 ? Jd_ : 1), ellMin(ellMin_), ellMax(ellMax_),
  pExp(pExp_),
  // Jlast = Jd makes the factor of the first step (Jd/Jlast)^p = 1, as
  // LoadControl does with numIncrLastStep = numIncr
  numIncrLastStep(Jd_ > 0 ? Jd_ : 1),
  stabilize(stabilize_), adaptStab(adaptStab_), fTarget(fTarget_), cVisc(cVisc_),
  dtPseudo(dtPseudo_ > 0.0 ? dtPseudo_ : 1.0),
  cOverDt(cVisc_ / (dtPseudo_ > 0.0 ? dtPseudo_ : 1.0)),
  // a given c (-cVisc) is used as is; otherwise c is calibrated from
  // fTarget at the first converged stabilized step
  cCalibrated(cVisc_ != 0.0), Estrain0(0.0), dissipVisc(0.0), residualTrueNorm(0.0),
  dLambda0(fabs(arcLength)), massMode(massMode_), massScale(massScale_),
  nEqn(0), Mstar(0),
  cArcLength2(arcLength*arcLength), cDeltaLambdaStep(0.0), cCurrentLambda(0.0),
  cSign(1), cDeltaUstep(0), haveSnapshot(false)
{

}

LadrunoArcLength::~LadrunoArcLength()
{
    if (deltaUhat  != 0) delete deltaUhat;
    if (deltaU     != 0) delete deltaU;
    if (deltaUstep != 0) delete deltaUstep;
    if (deltaUbar  != 0) delete deltaUbar;
    if (phat       != 0) delete phat;
    if (Mstar      != 0) delete Mstar;
    if (cDeltaUstep != 0) delete cDeltaUstep;
}

// diagonal fictitious mass without the dt^2/4 prefactor: the scale of the
// viscous term is carried by cOverDt
int
LadrunoArcLength::buildArtificialMass(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0 || Mstar == 0) return -1;
    return Ladruno::buildGershgorinDiagonal(theModel, *Mstar, massMode, massScale, 0.0);
}

// K, plus cOverDt*M* on the diagonal when stabilizing
int
LadrunoArcLength::formTangent(int statFlag)
{
    int res = this->IncrementalIntegrator::formTangent(statFlag);
    if (res < 0)
        return res;
    if (stabilize && Mstar != 0)
        Ladruno::addDiagonalToSOE(this->getLinearSOE(), *Mstar, nEqn, cOverDt,
                                  /*zeroFirst=*/false);
    return 0;
}

// lambda p - f_int, minus cOverDt*M*.DeltaU_step when stabilizing. The
// residual before the viscous force is kept for the convergence test.
int
LadrunoArcLength::formUnbalance(void)
{
    int res = this->IncrementalIntegrator::formUnbalance();
    if (res < 0)
        return res;
    if (stabilize && Mstar != 0 && deltaUstep != 0) {
        LinearSOE *soe = this->getLinearSOE();
        const Vector &B = soe->getB();
        residualTrueNorm = B.Norm();
        residualTrue = B;
        Vector R(B);
        for (int i = 0; i < nEqn; i++)
            R(i) -= cOverDt * (*Mstar)(i) * (*deltaUstep)(i);
        soe->setB(R);
    }
    return 0;
}

int
LadrunoArcLength::newStep(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING LadrunoArcLength::newStep() - "
                  "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }

    // get the current load factor
    currentLambda = theModel->getCurrentDomainTime();

    // base sign of the load change on what happened last step
    if (deltaLambdaStep < 0)
        signLastDeltaLambdaStep = -1;
    else
        signLastDeltaLambdaStep = +1;

    // -stabilize: incremental load control on the regularized system. The
    // load factor is fixed within the step; -adapt scales the load increment
    // dLambda0 instead of the arc radius.
    if (stabilize) {
        if (adapt) {
            int jlast = numIncrLastStep > 0 ? numIncrLastStep : 1;
            dLambda0 *= pow((double)Jd / (double)jlast, pExp);
            if (dLambda0 < ellMin) dLambda0 = ellMin;
            else if (ellMax > 0.0 && dLambda0 > ellMax) dLambda0 = ellMax;
        }
        numIncrLastStep = 0;
        deltaUstep->Zero();
        deltaLambdaStep = signLastDeltaLambdaStep * dLambda0;
        currentLambda += deltaLambdaStep;
        theModel->applyLoadDomain(currentLambda);
        if (theModel->updateDomain() < 0) {
            opserr << "LadrunoArcLength::newStep() - model failed to update\n";
            return -1;
        }
        return 0;
    }

    // radius adaptation (Ramm 1981):
    //   arcLength <- clamp(arcLength * (Jd / Jlast)^p, ellMin, ellMax)
    if (adapt) {
        int jlast = numIncrLastStep > 0 ? numIncrLastStep : 1;
        double factor = pow((double)Jd / (double)jlast, pExp);
        double ell = sqrt(arcLength2) * factor;
        if (ell < ellMin) ell = ellMin;
        else if (ell > ellMax) ell = ellMax;
        arcLength2 = ell * ell;
    }
    numIncrLastStep = 0;

    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
    if (theLinSOE->solve() < 0) {
        opserr << "LadrunoArcLength::newStep(void) - failed in solver\n";
        return -1;
    }
    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    // determine delta lambda(1) == dlambda
    double dLambda = sqrt(arcLength2 / ((dUhat ^ dUhat) + alpha2));
    dLambda *= signLastDeltaLambdaStep;
    deltaLambdaStep = dLambda;
    currentLambda += dLambda;

    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);
    theModel->applyLoadDomain(currentLambda);
    if (theModel->updateDomain() < 0) {
        opserr << "LadrunoArcLength::newStep() - model failed to update for new dU\n";
        return -1;
    }

    return 0;
}

int
LadrunoArcLength::update(const Vector &dU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING LadrunoArcLength::update() - "
                  "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }

    // -stabilize: load-control corrector; accumulate DeltaU_step for f_v
    if (stabilize) {
        (*deltaU) = dU;
        (*deltaUstep) += dU;
        theModel->incrDisp(dU);
        theModel->applyLoadDomain(currentLambda);
        if (theModel->updateDomain() < 0) {
            opserr << "LadrunoArcLength::update() - model failed to update\n";
            return -1;
        }
        theLinSOE->setX(dU);
        numIncrLastStep++;
        return 0;
    }

    (*deltaUbar) = dU;   // have to do this as the SOE is gonna change

    // determine dUhat
    theLinSOE->setB(*phat);
    if (theLinSOE->solve() < 0) {
        opserr << "LadrunoArcLength::update() - failed in solver\n";
        return -1;
    }
    (*deltaUhat) = theLinSOE->getX();

    // determine the coefficients of the quadratic equation
    a = alpha2 + ((*deltaUhat) ^ (*deltaUhat));
    b = alpha2 * deltaLambdaStep
        + ((*deltaUhat) ^ (*deltaUbar))
        + ((*deltaUstep) ^ (*deltaUhat));
    b *= 2.0;
    c = 2.0 * ((*deltaUstep) ^ (*deltaUbar)) + ((*deltaUbar) ^ (*deltaUbar));

    // check for a solution to quadratic
    b24ac = b * b - 4.0 * a * c;
    if (b24ac < 0) {
        opserr << "LadrunoArcLength::update() - imaginary roots due to multiple "
                  "instability directions - initial load increment was too large\n";
        opserr << "a: " << a << " b: " << b << " c: " << c << " b24ac: " << b24ac << endln;
        return -1;
    }
    double a2 = 2.0 * a;
    if (a2 == 0.0) {
        opserr << "LadrunoArcLength::update() - zero denominator,"
                  " alpha was set to 0.0 and zero reference load\n";
        return -2;
    }

    // determine the roots of the quadratic
    double sqrtb24ac = sqrt(b24ac);
    double dlambda1 = (-b + sqrtb24ac) / a2;
    double dlambda2 = (-b - sqrtb24ac) / a2;

    double val = (*deltaUhat) ^ (*deltaUstep);
    double theta1 = ((*deltaUstep) ^ (*deltaUstep)) + ((*deltaUbar) ^ (*deltaUstep));
    theta1 += dlambda1 * val;

    // choose dLambda based on angle between incremental displacement before
    // and after this step -- want positive
    double dLambda;
    if (theta1 > 0)
        dLambda = dlambda1;
    else
        dLambda = dlambda2;

    // determine delta U(i)
    (*deltaU) = (*deltaUbar);
    deltaU->addVector(1.0, *deltaUhat, dLambda);

    // update dU and dlambda
    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda += dLambda;

    // update the model
    theModel->incrDisp(*deltaU);
    theModel->applyLoadDomain(currentLambda);
    if (theModel->updateDomain() < 0) {
        opserr << "LadrunoArcLength::update() - model failed to update for new dU\n";
        return -1;
    }

    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    numIncrLastStep++;

    return 0;
}

int
LadrunoArcLength::domainChanged(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING LadrunoArcLength::domainChanged() - "
                  "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }
    int size = theModel->getNumEqn();   // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) {
        if (deltaUhat != 0) delete deltaUhat;
        deltaUhat = new Vector(size);
    }
    if (deltaUbar == 0 || deltaUbar->Size() != size) {
        if (deltaUbar != 0) delete deltaUbar;
        deltaUbar = new Vector(size);
    }
    if (deltaU == 0 || deltaU->Size() != size) {
        if (deltaU != 0) delete deltaU;
        deltaU = new Vector(size);
    }
    if (deltaUstep == 0 || deltaUstep->Size() != size) {
        if (deltaUstep != 0) delete deltaUstep;
        deltaUstep = new Vector(size);
    }
    if (phat == 0 || phat->Size() != size) {
        if (phat != 0) delete phat;
        phat = new Vector(size);
    }
    if (cDeltaUstep == 0 || cDeltaUstep->Size() != size) {
        if (cDeltaUstep != 0) delete cDeltaUstep;
        cDeltaUstep = new Vector(size);
    }
    if (deltaUhat == 0 || deltaUbar == 0 || deltaU == 0 ||
        deltaUstep == 0 || phat == 0 || cDeltaUstep == 0 ||
        deltaUhat->Size() != size || phat->Size() != size) {
        opserr << "FATAL LadrunoArcLength::domainChanged() - ran out of memory "
                  "for Vectors of size " << size << endln;
        return -1;
    }

    // determine phat:
    // do this by incrementing lambda by 1, applying load and getting phat
    // from the unbalance. Mstar is not yet built for a new size, and
    // deltaUstep is zero after a resize, so formUnbalance() adds no viscous
    // force here.
    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);
    if (Mstar != 0)
        Mstar->Zero();
    this->formUnbalance();   // NOTE: assumes unbalance at last was 0
    (*phat) = theLinSOE->getB();
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);

    // check there is a reference load
    int haveLoad = 0;
    for (int i = 0; i < size; i++)
        if ((*phat)(i) != 0.0) {
            haveLoad = 1;
            i = size;
        }
    if (haveLoad == 0) {
        opserr << "WARNING LadrunoArcLength::domainChanged() - zero reference load\n";
        return -1;
    }

    // -stabilize: build M* after phat is known
    nEqn = size;
    if (stabilize) {
        if (Mstar == 0 || Mstar->Size() != size) {
            if (Mstar != 0) delete Mstar;
            Mstar = new Vector(size);
        }
        this->buildArtificialMass();
    }

    // the state a failure before the first commit() returns to
    cArcLength2      = arcLength2;
    cDeltaLambdaStep = deltaLambdaStep;
    cCurrentLambda   = currentLambda;
    cSign            = signLastDeltaLambdaStep;
    cDeltaUstep->Zero();

    return 0;
}

int
LadrunoArcLength::commit(void)
{
    if (stabilize && Mstar != 0 && deltaUstep != 0) {
        // viscous work of the step, W_v = sum cOverDt M*_i DeltaU_i^2
        double Wv = 0.0;
        for (int i = 0; i < nEqn; i++)
            Wv += cOverDt * (*Mstar)(i) * (*deltaUstep)(i) * (*deltaUstep)(i);
        dissipVisc += Wv;

        // calibrate c once from the energy fraction f (unless -cVisc given)
        if (!cCalibrated) {
            Estrain0 = 0.5 * fabs((*deltaUstep) ^ (*phat)) * fabs(currentLambda);
            double Wunit = 0.0;
            for (int i = 0; i < nEqn; i++)
                Wunit += (*Mstar)(i) * (*deltaUstep)(i) * (*deltaUstep)(i) / dtPseudo;
            cVisc   = (Wunit > 0.0) ? fTarget * Estrain0 / Wunit : 0.0;
            cOverDt = cVisc / dtPseudo;
            cCalibrated = true;
        } else if (adaptStab && Estrain0 > 0.0) {
            double ratio = dissipVisc / Estrain0;
            if (ratio > 0.0) { cVisc *= (fTarget / ratio); cOverDt = cVisc / dtPseudo; }
        }
    }

    // store the converged per-step state for revertToLastStep()
    cArcLength2      = arcLength2;
    cDeltaLambdaStep = deltaLambdaStep;
    cCurrentLambda   = currentLambda;
    cSign            = signLastDeltaLambdaStep;
    if (cDeltaUstep != 0 && deltaUstep != 0 &&
        cDeltaUstep->Size() == deltaUstep->Size())
        *cDeltaUstep = *deltaUstep;
    haveSnapshot = true;

    return this->StaticIntegrator::commit();
}

int
LadrunoArcLength::revertToLastStep(void)
{
    // Restore the per-step state of the last commit(). The radius (arcLength2
    // and dLambda0) is left as the script set it: restoring it would undo each
    // reduceStep() and a retry loop could never shrink the step. With -adapt
    // the next newStep() derives the radius from numIncrLastStep anyway.
    // No-op before the first commit().
    if (!haveSnapshot)
        return 0;

    deltaLambdaStep         = cDeltaLambdaStep;
    currentLambda           = cCurrentLambda;
    signLastDeltaLambdaStep = cSign;
    if (deltaUstep != 0 && cDeltaUstep != 0 &&
        deltaUstep->Size() == cDeltaUstep->Size())
        *deltaUstep = *cDeltaUstep;
    return 0;
}

// Radius mutators. A script drives the retry loop, e.g.
//   ok = analyze(1); while ok != 0: ladrunoArcLength reduceStep 0.5; ok = analyze(1)
int
LadrunoArcLength::reduceStep(double f)
{
    if (f <= 0.0) {
        opserr << "WARNING LadrunoArcLength::reduceStep - factor must be > 0 "
                  "(got " << f << "); ignored\n";
        return -1;
    }
    double newAL2 = arcLength2 * f * f;
    // floor at ellMin when given, else at a tiny absolute radius
    double floor2 = (ellMin > 0.0) ? ellMin * ellMin : 1.0e-30;
    if (newAL2 < floor2) newAL2 = floor2;
    arcLength2 = newAL2;
    dLambda0  *= f;                       // load increment of -stabilize
    return 0;
}

int
LadrunoArcLength::increaseStep(double f)
{
    if (f <= 0.0) {
        opserr << "WARNING LadrunoArcLength::increaseStep - factor must be > 0 "
                  "(got " << f << "); ignored\n";
        return -1;
    }
    double newAL2 = arcLength2 * f * f;
    if (ellMax > 0.0) {
        double ceil2 = ellMax * ellMax;
        if (newAL2 > ceil2) newAL2 = ceil2;
    }
    arcLength2 = newAL2;
    dLambda0  *= f;
    return 0;
}

int
LadrunoArcLength::setArcLength(double arcLength)
{
    if (arcLength <= 0.0) {
        opserr << "WARNING LadrunoArcLength::setArcLength - arcLength must be > 0 "
                  "(got " << arcLength << "); ignored\n";
        return -1;
    }
    arcLength2 = arcLength * arcLength;
    dLambda0   = arcLength;
    return 0;
}

double LadrunoArcLength::getArcLength(void) const           { return sqrt(arcLength2); }
double LadrunoArcLength::getDeltaLambdaStep(void) const     { return deltaLambdaStep; }
double LadrunoArcLength::getCurrentLambda(void) const       { return currentLambda; }
int    LadrunoArcLength::getSignLastDeltaLambdaStep(void) const { return signLastDeltaLambdaStep; }
double LadrunoArcLength::getDeltaUstepNorm(void) const
{
    return (deltaUstep != 0) ? deltaUstep->Norm() : 0.0;
}

int
LadrunoArcLength::getStabilizedTrueResidual(double &normOut) const
{
    if (!stabilize) return -1;
    normOut = residualTrueNorm;
    return 0;
}

const Vector *
LadrunoArcLength::getStabilizedTrueResidualVector(void) const
{
    return stabilize ? &residualTrue : 0;
}

double LadrunoArcLength::getStabilizationDissipatedEnergy(void) const
{
    return dissipVisc;
}

double LadrunoArcLength::getReferenceStrainEnergy(void) const
{
    return Estrain0;
}

double LadrunoArcLength::getStabilizationDissipationRatio(void) const
{
    return (Estrain0 > 0.0) ? dissipVisc / Estrain0 : 0.0;
}

int LadrunoArcLength::scaleCVisc(double factor)
{
    if (factor <= 0.0) {
        opserr << "WARNING LadrunoArcLength::scaleCVisc - factor must be > 0 "
                  "(got " << factor << "); ignored\n";
        return -1;
    }
    cVisc   *= factor;
    cOverDt  = (dtPseudo > 0.0) ? cVisc / dtPseudo : cVisc;
    return 0;
}

int
LadrunoArcLength::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(27);
    data(0) = arcLength2;
    data(1) = alpha2;
    data(2) = deltaLambdaStep;
    data(3) = currentLambda;
    data(4) = signLastDeltaLambdaStep;
    data(5) = adapt ? 1.0 : 0.0;
    data(6) = (double)Jd;
    data(7) = ellMin;
    data(8) = ellMax;
    data(9) = pExp;
    data(10) = (double)numIncrLastStep;
    data(11) = stabilize ? 1.0 : 0.0;
    data(12) = adaptStab ? 1.0 : 0.0;
    data(13) = fTarget;
    data(14) = cVisc;
    data(15) = dtPseudo;
    data(16) = cOverDt;
    data(17) = cCalibrated ? 1.0 : 0.0;
    data(18) = (double)massMode;
    data(19) = massScale;
    // the per-step vectors are rebuilt by domainChanged(), as in ArcLength
    data(20) = dLambda0;
    data(21) = cArcLength2;
    data(22) = cDeltaLambdaStep;
    data(23) = cCurrentLambda;
    data(24) = (double)cSign;
    data(25) = haveSnapshot ? 1.0 : 0.0;
    data(26) = 0.0;   // reserved

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoArcLength::sendSelf() - failed to send the data\n";
        return -1;
    }
    return 0;
}

int
LadrunoArcLength::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(27);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "LadrunoArcLength::recvSelf() - failed to receive the data\n";
        return -1;
    }
    arcLength2              = data(0);
    alpha2                  = data(1);
    deltaLambdaStep         = data(2);
    currentLambda           = data(3);
    signLastDeltaLambdaStep = (int)data(4);
    adapt                   = data(5) != 0.0;
    Jd                      = (int)data(6);
    ellMin                  = data(7);
    ellMax                  = data(8);
    pExp                    = data(9);
    numIncrLastStep         = (int)data(10);
    stabilize               = data(11) != 0.0;
    adaptStab               = data(12) != 0.0;
    fTarget                 = data(13);
    cVisc                   = data(14);
    dtPseudo                = data(15);
    cOverDt                 = data(16);
    cCalibrated             = data(17) != 0.0;
    massMode                = (int)data(18);
    massScale               = data(19);
    dLambda0                = data(20);
    cArcLength2             = data(21);
    cDeltaLambdaStep        = data(22);
    cCurrentLambda          = data(23);
    cSign                   = (int)data(24);
    haveSnapshot            = data(25) != 0.0;
    return 0;
}

void
LadrunoArcLength::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double cLambda = theModel->getCurrentDomainTime();
        s << "\t LadrunoArcLength - currentLambda: " << cLambda;
        s << "  arcLength: " << sqrt(arcLength2) << "  alpha: " << sqrt(alpha2);
        if (adapt)
            s << "  [adapt Jd=" << Jd << " ell=[" << ellMin << "," << ellMax
              << "] p=" << pExp << "]";
        if (stabilize)
            s << "  [stabilize f=" << fTarget << " c=" << cVisc
              << " diss/E=" << (Estrain0 > 0.0 ? dissipVisc / Estrain0 : 0.0)
              << " trueRes=" << residualTrueNorm << "]";
        s << endln;
    } else
        s << "\t LadrunoArcLength - no associated AnalysisModel\n";
}
