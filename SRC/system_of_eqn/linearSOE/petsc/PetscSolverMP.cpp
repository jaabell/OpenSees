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


#include <OPS_Globals.h>
#include <PetscSolverMP.h>
#include <PetscSOEMP.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Timer.h>
#include <ID.h>
#include <Vector.h>
#include <string.h>
#include "petscksp.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <petsctime.h>

using namespace std;
//-------------------------------------

PetscSolverMP::PetscSolverMP()
    : LinearSOESolver(SOLVER_TAGS_PetscSolverMP),
      rTol(PETSC_DEFAULT), aTol(PETSC_DEFAULT), dTol(PETSC_DEFAULT), maxIts(PETSC_DEFAULT), matType(MATMPIAIJ),
      is_KSP_initialized(false)
{

}

PetscSolverMP::PetscSolverMP(KSPType meth, PCType pre)
    : LinearSOESolver(SOLVER_TAGS_PetscSolverMP), method(meth), preconditioner(pre),
      rTol(PETSC_DEFAULT), aTol(PETSC_DEFAULT), dTol(PETSC_DEFAULT), maxIts(PETSC_DEFAULT), matType(MATMPIAIJ),
      is_KSP_initialized(false)
{

}

PetscSolverMP::PetscSolverMP(KSPType meth, PCType pre, double relTol, double absTol, double divTol, int maxIterations, MatType mat)
    : LinearSOESolver(SOLVER_TAGS_PetscSolverMP), method(meth), preconditioner(pre),
      rTol(relTol), aTol(absTol), dTol(divTol), maxIts(maxIterations), matType(mat),
      is_KSP_initialized(false)
{

}

PetscSolverMP::~PetscSolverMP()
{
    // KSPDestroy(&ksp);
    // CHKERRQ(ierr);
}


int
PetscSolverMP::solve(void)
{

    int size = theSOE->size;
    // int numProcesses = theSOE->numProcesses;
    int processID = theSOE->processID;
    int ierr;

    if (processID >= 0)
    {

        //MatSetType(theSOE->A, matType);
        PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") MatAssemblyBegin\n";
        ierr = MatAssemblyBegin(theSOE->A, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);

        ierr = MatAssemblyEnd(theSOE->A, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") MatAssemblyEnd\n";


        if (not is_KSP_initialized)
        {

            PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") KSPCreate Begin\n";
            KSPCreate(PETSC_COMM_WORLD, &ksp);
            KSPSetFromOptions(ksp);
            PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") KSPCreate End\n";


            is_KSP_initialized = true;
        }
        
        PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") KSPSetOperators Begin\n";
        KSPSetOperators(ksp, theSOE->A, theSOE->A);
        PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") KSPSetOperators End\n";

    }



    // Everyone sends the assembled vector B to the rank 0 processor. Which sends it back...
    static Vector recvVector(1);


    Vector *vectX = theSOE->vectX;
    Vector *vectB = theSOE->vectB;

    // zero X
    vectX->Zero();

    //
    // form B on each
    //
    PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") FormBBegin\n";
    int numChannels = theSOE->numChannels;
    Channel **theChannels = theSOE->theChannels;

    if (processID != 0)
    {
        Channel *theChannel = theChannels[0];

        theChannel->sendVector(0, 0, *vectB);
        theChannel->recvVector(0, 0, *vectB);

        double vectBNorm  = vectB->Norm();
        if (std::isnan(vectBNorm))
        {
            cout << "norm(vectB) = " << vectBNorm << endl;
            return -1;
        }

    }
    else //Master process assembles b vector. This is a bottleneck and can be improved with
        //a collective reduction
    {

        if (recvVector.Size() != size)
        {
            recvVector.resize(size);
        }

        //Better done with ALLREDUCE
        for (int j = 0; j < numChannels; j++)
        {
            Channel *theChannel = theChannels[j];
            theChannel->recvVector(0, 0, recvVector);
            *vectB += recvVector;

            double rvNorm  = recvVector.Norm() ;
            if (std::isnan(rvNorm))
            {
                cout << "norm(recvVector) = " << rvNorm << endl;
                return -1;
            }
        }

        //Better done with a BCast
        for (int j = 0; j < numChannels; j++)
        {
            Channel *theChannel = theChannels[j];
            theChannel->sendVector(0, 0, *vectB);
        }
    }
    PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") FormBEnd\n";


    //
    // solve and mark as having been solved
    //

    double t1 = 0;
    double t2 = 0;
    PetscTime(&t1);
    PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") SolveBegin\n";
    if (processID >= 0)
    {
        ierr = KSPSolve(ksp, theSOE->b, theSOE->x);
        CHKERRQ(ierr);
        theSOE->isFactored = 1;
    }
    PetscTime(&t2);
    double delta_time = t2 - t1;
    PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") SolveEnd dt = " << delta_time << "s \n";





    //
    // if parallel, we must form the total X: each processor has startRow through endRow-1
    //  
    vectX = theSOE->vectX;

    numChannels = theSOE->numChannels;
    theChannels = theSOE->theChannels;
    PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") SendXBegin\n";
    if (processID != 0)
    {
        Channel *theChannel = theChannels[0];

        theChannel->sendVector(0, 0, *vectX);
        theChannel->recvVector(0, 0, *vectX);

    }
    else //Again, the master process forms the global X vector which is then sent to all processors
    {

        if (recvVector.Size() != size)
        {
            recvVector.resize(size);
        }

        //Better done with ALLREDUCE
        for (int j = 0; j < numChannels; j++)
        {
            Channel *theChannel = theChannels[j];
            theChannel->recvVector(0, 0, recvVector);
            *vectX += recvVector;
        }

        //Better done with a BCast
        for (int j = 0; j < numChannels; j++)
        {
            Channel *theChannel = theChannels[j];
            theChannel->sendVector(0, 0, *vectX);
        }
    }

    PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") SendXEnd\n";
    // Destroy KSP and collect the error at P0
    KSPDestroy(&ksp);
    is_KSP_initialized = false;

    PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") Logging Begin\n";

#ifdef _ENABLE_PETSC_LOGGER
        PetscViewer    viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, "petsc_log.txt" , &viewer);
        PetscLogView(viewer);
#endif
        
    PetscSolverMP_DEBUGOUT << "PetscSolverMP::solve (" << processID << ") Logging End - Done solve\n";

    return ierr;
}






int
PetscSolverMP::setSize()
{
    /*
     * Create linear solver context
     */
    if (theSOE->processID >= 0)
    {
        KSPCreate(PETSC_COMM_WORLD, &ksp);
    }


    return 0;
}


int
PetscSolverMP::setLinearSOE(PetscSOEMP &theSys)
{
    theSOE = &theSys;
    return 0;
}


int
PetscSolverMP::sendSelf(int cTag, Channel &theChannel)
{
   
    return 0;
}

int
PetscSolverMP::recvSelf(int cTag, Channel &theChannel,
                         FEM_ObjectBroker &theBroker)
{
   

    return 0;
}



