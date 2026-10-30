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
// Description: see LadrunoStabilizedUnbalance.h. The iteration logic and
// print flags follow CTestNormUnbalance.

#include <LadrunoStabilizedUnbalance.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <IncrementalIntegrator.h>
#include <LadrunoArcLength.h>
#include <classTags.h>
#include <elementAPI.h>

void *OPS_LadrunoStabilizedUnbalance(void)
{
    // Usage: test LadrunoStabilizedUnbalance $tol $maxIter <$printFlag $normType>
    // Drop-in for `test NormUnbalance ...`; norms the TRUE static unbalance when
    // the active integrator is a stabilizing LadrunoArcLength.
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING test LadrunoStabilizedUnbalance $tol $maxIter "
                  "<$printFlag $normType>\n";
        return 0;
    }
    double tol = 1.0e-6;
    int    numData = 1;
    if (OPS_GetDoubleInput(&numData, &tol) < 0) {
        opserr << "WARNING test LadrunoStabilizedUnbalance failed to read tol\n";
        return 0;
    }
    numData = OPS_GetNumRemainingInputArgs();
    if (numData > 3) numData = 3;
    int data[3] = {25, 0, 2};   // maxIter, printFlag, normType
    if (OPS_GetIntInput(&numData, &data[0]) < 0) {
        opserr << "WARNING test LadrunoStabilizedUnbalance failed to read int values\n";
        return 0;
    }
    return new LadrunoStabilizedUnbalance(tol, data[0], data[1], data[2]);
}

LadrunoStabilizedUnbalance::LadrunoStabilizedUnbalance()
 :ConvergenceTest(CONVERGENCE_TEST_LadrunoStabilizedUnbalance),
  theSOE(0), theIntegrator(0), tol(0.0), maxTol(OPS_MAXTOL),
  maxNumIter(0), currentIter(0), printFlag(0), nType(2), norms(1)
{

}

LadrunoStabilizedUnbalance::LadrunoStabilizedUnbalance(double theTol, int maxIter,
                                                       int printIt, int normType,
                                                       double max)
 :ConvergenceTest(CONVERGENCE_TEST_LadrunoStabilizedUnbalance),
  theSOE(0), theIntegrator(0), tol(theTol), maxTol(max),
  maxNumIter(maxIter), currentIter(0), printFlag(printIt), nType(normType),
  norms(maxIter > 0 ? maxIter : 1)
{

}

LadrunoStabilizedUnbalance::~LadrunoStabilizedUnbalance()
{

}

ConvergenceTest *
LadrunoStabilizedUnbalance::getCopy(int iterations)
{
    LadrunoStabilizedUnbalance *theCopy =
        new LadrunoStabilizedUnbalance(tol, iterations, printFlag, nType, maxTol);
    theCopy->theSOE = this->theSOE;
    theCopy->theIntegrator = this->theIntegrator;
    return theCopy;
}

void
LadrunoStabilizedUnbalance::setTolerance(double newTol)
{
    tol = newTol;
}

int
LadrunoStabilizedUnbalance::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE        = theAlgo.getLinearSOEptr();
    theIntegrator = theAlgo.getIncrementalIntegratorPtr();   // for the true residual
    return 0;
}

// the residual norm to test against: the LadrunoArcLength true unbalance when it
// is stabilizing, else the stock SOE ||B||.
double
LadrunoStabilizedUnbalance::currentResidualNorm(bool &usedTrue) const
{
    usedTrue = false;
    if (theIntegrator != 0 &&
        theIntegrator->getClassTag() == INTEGRATOR_TAGS_LadrunoArcLength) {
        LadrunoArcLength *la = (LadrunoArcLength *)theIntegrator;
        const Vector *v = la->getStabilizedTrueResidualVector();
        if (v != 0) { usedTrue = true; return v->pNorm(nType); }   // honor the configured nType
    }
    return (theSOE != 0) ? theSOE->getB().pNorm(nType) : 0.0;
}

int
LadrunoStabilizedUnbalance::test(void)
{
    if (theSOE == 0) {
        opserr << "WARNING LadrunoStabilizedUnbalance::test() - no SOE set.\n";
        return -2;
    }
    if (currentIter == 0) {
        opserr << "WARNING LadrunoStabilizedUnbalance::test() - start() not invoked.\n";
        return -2;
    }

    bool usedTrue = false;
    double norm = this->currentResidualNorm(usedTrue);
    if (currentIter <= maxNumIter) norms(currentIter - 1) = norm;

    if (printFlag == 1) {
        opserr << "LadrunoStabilizedUnbalance::test() - iteration: " << currentIter
               << (usedTrue ? " TRUE |R|: " : " SOE |B|: ") << norm
               << " (max: " << tol << ")\n";
    }

    // converged
    if (norm <= tol) {
        if (printFlag == 2 || printFlag == 6) {
            opserr << "LadrunoStabilizedUnbalance::test() - iteration: " << currentIter
                   << (usedTrue ? " TRUE |R|: " : " SOE |B|: ") << norm
                   << " (max: " << tol << ")\n";
        }
        return currentIter;
    }
    // failed but go on
    else if ((printFlag == 5 || printFlag == 6) && currentIter >= maxNumIter) {
        opserr << "WARNING LadrunoStabilizedUnbalance::test() - failed to converge "
                  "but going on, current Norm: " << norm << " (max: " << tol << ")\n";
        return currentIter;
    }
    // failed
    else if (currentIter >= maxNumIter || norm > maxTol) {
        opserr << "WARNING LadrunoStabilizedUnbalance::test() - failed to converge after "
               << currentIter << " iterations, current Norm: " << norm
               << " (max: " << tol << ")\n";
        currentIter++;
        return -2;
    }
    // keep going
    else {
        currentIter++;
        return -1;
    }
}

int
LadrunoStabilizedUnbalance::start(void)
{
    if (theSOE == 0) {
        opserr << "WARNING LadrunoStabilizedUnbalance::start() - no SOE set.\n";
        return -1;
    }
    currentIter = 1;
    norms.Zero();
    return 0;
}

int    LadrunoStabilizedUnbalance::getNumTests(void)       { return currentIter; }
int    LadrunoStabilizedUnbalance::getMaxNumTests(void)    { return maxNumIter; }
double LadrunoStabilizedUnbalance::getRatioNumToMax(void)
{
    double div = (maxNumIter > 0) ? maxNumIter : 1;
    return currentIter / div;
}
const Vector &LadrunoStabilizedUnbalance::getNorms(void)   { return norms; }

int
LadrunoStabilizedUnbalance::sendSelf(int cTag, Channel &theChannel)
{
    static Vector x(5);
    x(0) = tol;
    x(1) = maxNumIter;
    x(2) = printFlag;
    x(3) = nType;
    x(4) = maxTol;
    if (theChannel.sendVector(this->getDbTag(), cTag, x) < 0) {
        opserr << "LadrunoStabilizedUnbalance::sendSelf() - failed to send data\n";
        return -1;
    }
    return 0;
}

int
LadrunoStabilizedUnbalance::recvSelf(int cTag, Channel &theChannel,
                                     FEM_ObjectBroker &theBroker)
{
    static Vector x(5);
    if (theChannel.recvVector(this->getDbTag(), cTag, x) < 0) {
        opserr << "LadrunoStabilizedUnbalance::recvSelf() - failed to receive data\n";
        tol = 1.0e-6; maxNumIter = 25; printFlag = 0; nType = 2; maxTol = OPS_MAXTOL;
        return -1;
    }
    tol        = x(0);
    maxNumIter = (int)x(1);
    printFlag  = (int)x(2);
    nType      = (int)x(3);
    maxTol     = x(4);
    norms.resize(maxNumIter > 0 ? maxNumIter : 1);
    return 0;
}
