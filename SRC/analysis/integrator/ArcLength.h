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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:00:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/ArcLength.h,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/ArcLength.h
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for ArcLength.
// ArcLength is an algorithmic class for performing a static analysis
// using the arc length scheme, that is within a load step the following
// constraint is enforced: 
//  i=1        delta U^T delta U + alpha^2 delta lambda^2 = delta s^2
//  i>1        dU^T delta U + alpha^2 dLambda delta lambda = 0
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and arcLength is a control parameter.
//
// What: "@(#) ArcLength.h, revA"

#ifndef ArcLength_h
#define ArcLength_h

#include <StaticIntegrator.h>
#include <Vector.h>
class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class Domain;
class ArcLength : public StaticIntegrator
{
  public:
    ArcLength(double arcLength, double alpha = 1.0);
 //   ArcLength(int node, int dof,Domain *domain);
    ~ArcLength();

    int newStep(void);
    int update(const Vector &deltaU);
    int domainChanged(void);

    int commit(void);
    int revertToLastStep(void);

    // Continuous-time (continuation) mode -- OFF by default; see
    // ContinuationLambda.h and DisplacementControl::setContinuationTime.
    //
    // Progress-normalised exactly as in DisplacementControl, with the ARC
    // LENGTH as the progress measure instead of the controlled displacement:
    //
    //   DisplacementControl   dt = D |theIncrement| / |U_target|
    //   ArcLength             dt = D  sqrt(arcLength2) / |S_target|
    //
    // The mapping is one to one because the arc length has both properties the
    // controlled displacement has and lambda does not: the per-step increment
    // is PRESCRIBED (it is the method's own parameter, not an unknown), and it
    // is MONOTONE (ds > 0 always, being the root of a sum of squares -- it
    // keeps increasing even where the method retraces the path). The corrector
    // here is the spherical Crisfield constraint (see ::update: the -s^2 term
    // is absent from c because the previous iterate already lies on the
    // sphere), so a committed step advances the arc by exactly sqrt(arcLength2)
    // and the sum over N steps is exactly N*sqrt(arcLength2). Declare
    // S_target = N*s and the stage lands on its nominal end time exactly.
    //
    // A positive dtFixed is still accepted as the cruder alternative, used when
    // the stage declares no total arc length.
    int setContinuationTime(int lambdaChannel, double duration,
			    double target, double dtFixed);

    // Take the sign of the predictor from the ANGLE with the previous converged
    // step instead of from the sign of the previous lambda increment.
    //
    //   legacy      sign(dLambda_1) = sign(deltaLambdaStep of the last step)
    //   this        sign(dLambda_1) = sign(dU_committed . Uhat)
    //
    // The legacy rule is why ArcLength cannot pass a limit point: past the
    // extremum lambda starts increasing again, the sign flips, and the method
    // walks back down the path it just climbed. The replacement is not a new
    // idea imported from outside -- it is the SAME criterion the corrector in
    // ::update already uses to pick between the two roots of the quadratic
    // (there: dUstep . dUstep_new > 0), simply applied to the predictor too.
    //
    // OPT-IN, because it changes the mechanics: with the flag off this class is
    // bit-for-bit what it was.
    int setSignFromAngle(bool flag);

    // Change the arc-length increment IN PLACE, mid-stage. Re-issuing
    // "integrator ArcLength" instead builds a fresh object, and this class
    // carries the direction of travel of the previous converged step in
    // committedUstep (and in deltaLambdaStep for the legacy rule): a new object
    // starts with both zeroed, so the predictor loses the way it was going.
    // That is what makes step-size recovery unsafe for an arc-length method,
    // independently of any timeline bookkeeping.
    //
    // With -target the timeline follows automatically: dt = D s / S_target uses
    // the NEW s from the next step on, so the stage still ends on T_end when the
    // arc consumed reaches S_target.
    int setArcLength(double arcLength);


    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    

   //////////////////////////Sensitivity Begin/////////////////
      void formTangDispSensitivity(int gradNumber);
      double formdLambdaDh(int gradNumber);
      double getLambdaSensitivity(int gradNumber);
      int formSensitivityRHS(int gradNum);// it's been modified to compute dLambdadh and dUdh
      int formIndependentSensitivityRHS();
      int saveSensitivity(const Vector &v, int gradNum, int numGrads);
      int saveLambdaSensitivity(double dlambdadh, int gradNum, int numGrads);
      int commitSensitivity(int gradNum, int numGrads);
      int computeSensitivities(void);// this function is modified to obtain both dLambdadh and dUdh 
      int formEleResidual(FE_Element *theEle);
  bool computeSensitivityAtEachIteration();// A key that return 1 for loadControl and 2 for DisplacementControl
 void formResidualDispSensitivity( int gradNumber);


  protected:
    
  private:
    double arcLength2;
    double alpha2;
    double a,b,c,b24ac;
        Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep;
	Vector *deltaUstep2; // will be in the sensitivity part
    Vector *phat; // the reference load vector
    Vector *dUhatdh, *dphatdh,*dLAMBDAdh, *dUIJdh,*dDeltaUstepdh,*sensU,*Residual;
    double deltaLambdaStep,dDeltaLambdaStepdh, currentLambda, dlambdaJdh;
    int signLastDeltaLambdaStep;
    
    double dLAMBDA, dLAMBDA2; // Need it to be called when deriving dLambda1dh.
    double dlambda1dh;
    int gradNumber;

    // ---- continuous-time (continuation) mode; see setContinuationTime ----
    bool useContinuationTime; // false => legacy behaviour, bit for bit
    int lambdaChannel;
    double stageDuration;     // D_stage on the global timeline
    double stageTarget;       // total ARC LENGTH of the stage
    double dtFixed;           // used when stageTarget == 0
    double timeStep;          // domain time frozen for the current step
    double committedLambda;

    double continuationDt(void) const;

    // ---- angle-based predictor sign; see setSignFromAngle ------------------
    bool signFromAngle;        // false => legacy sign rule, bit for bit
    Vector *committedUstep;    // dU of the last CONVERGED step, kept separately
			       // from deltaUstep so a failed step cannot
			       // pollute the direction of travel
   int sensitivityFlag;
};

#endif

