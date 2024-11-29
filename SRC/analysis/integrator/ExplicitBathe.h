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

// $Revision$
// $Date$
// $URL$

#ifndef ExplicitBathe_h
#define ExplicitBathe_h

// Written: User
// Created: Date
// Revision: A
//
// Description: This file contains the class definition for ExplicitBathe.
// ExplicitBathe is an algorithmic class for performing a transient analysis
// using the explicit Bathe time integration scheme.

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class ExplicitBathe : public TransientIntegrator
{
public:
    // constructors
    ExplicitBathe();
    ExplicitBathe(double p);//, double q0, double q1, double q2, double s);
    
    // destructor
    ~ExplicitBathe();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    
    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    const Vector &getVel(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);

protected:

    int formRHS_tpdt();

private:
    double deltaT;

    // Explicit Bathe parameters
    double p;
    double q0;
    double q1;
    double q2;
    double s;

    // State variables
    Vector *U_t;                    // Acceleration at time t
    Vector *V_t;                    // Velocity at time t 
    Vector *A_t;                    // Acceleration at time t 
    Vector *R_t;                    // Forces at time t 
    Vector *U_tpdt;                 // Acceleration at time t + p Dt
    Vector *V_tpdt;                 // Velocity at time t + p Dt
    Vector *A_tpdt;                 // Acceleration at time t + p Dt
    Vector *Rhat_tpdt;              //    Equivalent Forces Hat at time t + p Dt
    Vector *Rfunnyhat_tpdt;         //    Equivalent Forces Funnyhat time t + p Dt
    Vector *U_tdt;                  // Acceleration at time t + Dt
    Vector *V_tdt;                  // Velocity at time t + Dt
    Vector *A_tdt;                  // Acceleration at time t + Dt
    Vector *R_tdt;                  // Forces at time t + Dt

    int updateCount;

    double a0;
    double a1;
    double a2;
    double a3;
    double a4;
    double a5;
    double a6;
    double a7;
};

#endif
