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
                                                                        
                                                                       
#ifndef PARDISOGenLinSOE_h
#define PARDISOGenLinSOE_h

// Written: M. Salehi opensees.net@gmail.com
// website : http://opensees.net
// Created: 02/19
// Revision: A

#include <LinearSOE.h>
#include <Vector.h>

class PARDISOGenLinSolver;

// matType selects the storage and, through the solver, the PARDISO mtype:
//
//   0  unsymmetric                    full CSR             mtype 11
//   1  symmetric positive definite    upper-triangle CSR   mtype  2
//   2  symmetric indefinite           upper-triangle CSR   mtype -2
//
// The numbering matches `system Mumps -matrixType`. For the symmetric types
// PARDISO reads the upper triangle only, row by row, with ascending column
// indices and the diagonal always present (Intel oneMKL Developer Reference,
// "Sparse Matrix Storage Formats"); setSize() builds and checks that layout.
//
// With matType 1 or 2 only the col >= row half of each element matrix is
// assembled. A genuinely unsymmetric tangent (contact, non-associated flow,
// follower loads, corotational transformations) is then solved with its upper
// triangle reflected, which can converge to a wrong answer. addA() samples the
// element matrices for this and warns once; the default stays 0.

class PARDISOGenLinSOE : public LinearSOE
{
  public:
    PARDISOGenLinSOE(PARDISOGenLinSolver &theSolver);
    PARDISOGenLinSOE(PARDISOGenLinSolver &theSolver, int matType);

    ~PARDISOGenLinSOE();

    int getNumEqn(void) const;
    int getMatType(void) const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA(void);
    void zeroB(void);
    
    const Vector &getX(void);
    const Vector &getB(void);    
    double normRHS(void);

    void setX(int loc, double value);        
    void setX(const Vector &x);        
    int setPARDISOGenLinSolver(PARDISOGenLinSolver &newSolver);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
	friend class PARDISOGenLinSolver;
  protected:
    
  private:
    int size;            // order of A
    int nnz;             // number of non-zeros in A
    double *A, *B, *X;   // 1d arrays containing coefficients of A, B and X
    int *colA, *rowStartA; // int arrays containing info about coefficients in A
    Vector *vectX;
    Vector *vectB;    
    int Asize, Bsize;    // size of the 1d array holding A
    bool factored;
    int matType;         // 0 unsymmetric, 1 SPD, 2 symmetric indefinite
    int asymWarned;      // half-storage asymmetry already reported
    int asymBudget;      // element-matrix checks left in this sampling window
    int asymPass;        // tangent assemblies seen for the current pattern
    int missWarned;      // an addA entry without a CSR slot already reported
};


#endif

