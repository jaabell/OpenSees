
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
                                                                        
// Written: M. Salehi opensees.net@gmail.com
// Created: 02/19
// Revision: A

//refs
//Structural Engineeringand Mechanics   Volume 48, Number 6, December25 2013, pages 849 - 878
//DOI: https://doi.org/10.12989/sem.2013.48.6.849	
//
//Comprehensive evaluation of structural geometrical nonlinear solution techniques Part I : Formulation and characteristics of the methods
//M.Rezaiee - Pajand, M.Ghalishooyan and M.Salehi - Ahmadabad
//FULLTEXT : https://www.researchgate.net/publication/264146397_Comprehensive_evaluation_of_structural_geometrical_nonlinear_solution_techniques_Part_I_Formulation_and_characteristics_of_the_methods


//Structural Engineeringand Mechanics   Volume 48, Number 6, December25 2013, pages 879 - 914
//DOI: https://doi.org/10.12989/sem.2013.48.6.879	
//
//Comprehensive evaluation of structural geometrical nonlinear solution techniques Part II : Comparing efficiencies of the methods
//M.Rezaiee - Pajand, M.Ghalishooyan and M.Salehi - Ahmadabad
//FULLTEXT : https://www.researchgate.net/publication/263361974_Comprehensive_evaluation_of_structural_geometrical_nonlinear_solution_techniques_Part_II_Comparing_efficiencies_of_the_methods


#ifndef EQPath_h
#define EQPath_h

#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;

#define SIGN_LAST_STEP      1
#define CHANGE_DETERMINANT  2

class EQPath : public StaticIntegrator
{
  public:
    EQPath(double arcLeng,int type);

    ~EQPath();

    int newStep(void);
    int update(const Vector &deltaU);
    int domainChanged(void);

    int commit(void);
    int revertToLastStep(void);

    // Continuous-time (continuation) mode -- OFF by default.
    // See ContinuationLambda.h. Prescribed-dt variant only, same stopgap and
    // Progress-normalised on the arc length, exactly as ArcLength and
    // DisplacementControl: -target is the total the stage covers in the
    // method's own progress measure. Here that is |du|, which the constraint
    // sets to arclen exactly (newStep: dLambda = sign*arclen/|uq0| and
    // du = dLambda*uq0, so |du| = arclen), so
    //     dt = D * arclen / |S_target|
    // A positive dtFixed is kept as the cruder alternative.
    int setContinuationTime(int lambdaChannel, double duration,
			    double target, double dtFixed);

    // Change the arc-length increment IN PLACE, mid-stage; see
    // ArcLength::setArcLength. This class keeps the previous step's direction in
    // du (newStep: sign from du^uq0), so a fresh object loses it just the same.
    int setArcLength(double arcLength);


    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:
    double arclen,dl,m;
    double sign;
    int type,changed,nitr;
    Vector *du,*du0, *uq, *uq0, *uqn, *ur;
    Vector *q;

    // ---- continuous-time (continuation) mode ----
    // NOTE: unlike the other continuation integrators, EQPath reads the load
    // factor from the domain as a LOCAL in both newStep() and update(), so it
    // needs its own running accumulator here.
    bool useContinuationTime;
    int lambdaChannel;
    double dtFixed;
    double timeStep;
    double contLambda;        // running lambda within the step
    double stageDuration;     // D_stage on the global timeline
    double stageTarget;       // total arc length of the stage

    double continuationDt(void) const;
    double committedLambda;
};

#endif


