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


#ifndef PetscSolverMP_h
#define PetscSolverMP_h

#include <petscksp.h>
#include <LinearSOESolver.h>


// Uncomment or define -D_DEBUG_PetscSolverMP to enable debugging
// #define _DEBUG_PetscSolverMP

#ifdef _DEBUG_PetscSolverMP
#define PetscSolverMP_DEBUGOUT cout // or any other ostream
#else
#define PetscSolverMP_DEBUGOUT 0 && cout
#endif

// Uncomment to enable logger that will creat petsc_log.txt on each
// run, to profile PETSc. 
// #define _ENABLE_PETSC_LOGGER

class PetscSOEMP;

class PetscSolverMP : public LinearSOESolver
{
public:
    PetscSolverMP();
    PetscSolverMP(KSPType method, PCType preconditioner);
    PetscSolverMP(KSPType method, PCType preconditioner, double rTol, double aTol, double dTol, int maxIts, MatType mat=MATMPIAIJ);//Guanzhou
    ~PetscSolverMP();

    int solve(void);
    int setSize(void);
    virtual int setLinearSOE(PetscSOEMP &theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                    FEM_ObjectBroker &theBroker);

    friend class ActorPetscSOEMP;
    friend class ShadowPetscSOEMP;

protected:
    PetscSOEMP *theSOE;

private:
    KSP ksp;
    PC pc;
    int its;
    KSPType method;
    PCType preconditioner;
    double rTol;
    double aTol;
    double dTol;
    int maxIts;
    MatType matType;

    bool is_KSP_initialized;
};

#endif

