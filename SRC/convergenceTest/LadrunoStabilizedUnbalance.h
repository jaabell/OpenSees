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
//   Dassault Systemes Simulia. Abaqus Analysis User's Guide, "Automatic
//     stabilization of unstable problems".
//
// Description: LadrunoStabilizedUnbalance, a NormUnbalance convergence
// test for the -stabilize mode of LadrunoArcLength.
//
// In -stabilize mode the integrator subtracts an artificial viscous force
// f_v from the residual, so the SOE vector B normed by CTestNormUnbalance is
// lambda p - f_int - f_v. Newton can drive that vector to zero while the
// static unbalance lambda p - f_int is still of the size of f_v. This test
// norms lambda p - f_int instead (with the configured norm type), so a
// converged step is an equilibrium state to the tolerance and the viscous
// force is below it.
//
// setEquiSolnAlgo() stores the LinearSOE and the IncrementalIntegrator of
// the algorithm. When the integrator is a LadrunoArcLength in -stabilize
// mode its residual before f_v is used; otherwise the test norms the SOE B
// and behaves as NormUnbalance (without the maxIncr option).

#ifndef LadrunoStabilizedUnbalance_h
#define LadrunoStabilizedUnbalance_h

#include <ConvergenceTest.h>
#include <Vector.h>

class EquiSolnAlgo;
class LinearSOE;
class IncrementalIntegrator;

class LadrunoStabilizedUnbalance : public ConvergenceTest
{
  public:
    LadrunoStabilizedUnbalance();
    LadrunoStabilizedUnbalance(double tol, int maxNumIter, int printFlag,
                               int normType = 2, double maxTol = OPS_MAXTOL);
    ~LadrunoStabilizedUnbalance();

    ConvergenceTest *getCopy(int iterations);

    void setTolerance(double newTol);
    int setEquiSolnAlgo(EquiSolnAlgo &theAlgo);

    int test(void);
    int start(void);

    int getNumTests(void);
    int getMaxNumTests(void);
    double getRatioNumToMax(void);
    const Vector &getNorms(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  protected:
    // the norm tested: the residual before f_v when a stabilizing
    // LadrunoArcLength is active (usedTrue = true), else the SOE ||B||
    double currentResidualNorm(bool &usedTrue) const;

  private:
    LinearSOE             *theSOE;
    IncrementalIntegrator *theIntegrator;
    double tol;
    double maxTol;
    int    maxNumIter;
    int    currentIter;
    int    printFlag;
    int    nType;
    Vector norms;
};

#endif
