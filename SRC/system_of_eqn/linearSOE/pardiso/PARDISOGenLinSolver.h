// File: ~/system_of_eqn/pardiso/PARDISOGenLinSolver.h
//
// Written: M. Salehi opensees.net@gmail.com
// website : http://opensees.net
// Created: 02/19
// Revision: A
//
// Description: This file contains the class definition for
// PARDISOLinSolver. It solves the Sparse General SOE by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) PARDISOGenLinSolver.h, revA"


#ifndef PARDISOGenLinSolver_h
#define PARDISOGenLinSolver_h

#include <LinearSOESolver.h>
#include <PARDISOGenLinSOE.h>

// Parameter meanings: Intel oneMKL Developer Reference, "pardiso iparm
// Parameter". Comments below use the zero-based C index iparm[i], which is
// iparm(i+1) in Intel's one-based table.

class PARDISOGenLinSolver : public LinearSOESolver
{
  public:
	  PARDISOGenLinSolver();
    ~PARDISOGenLinSolver();

    int solve(void);
    int setSize(void);

    int setLinearSOE(PARDISOGenLinSOE &theSOE);

    // Print PARDISO's memory, fill and flop counters after every numerical
    // factorization.
    void setStats(int on);

    // Use the retained factorization as a preconditioner for a changed matrix
    // (iparm[3] = 10*L + K, phase 23). digits is Intel's L: the CGS iteration
    // stops at ||dx_i||/||dx_0|| < 10^-L. 0 disables it (the default).
    void setKrylov(int digits);

    // Conditional numerical reproducibility (CNR). branch is an MKL_CBWR_*
    // code. mkl_cbwr_set() is called here, before PARDISO runs, and
    // iparm[33] is set to the MKL thread count at every symbolic phase.
    // The mode is process-wide; MKL refuses to change it once its BLAS/LAPACK
    // dispatch is initialized, and the MKL_CBWR environment variable is the
    // fallback. keepEnv = 1 keeps a branch already set through MKL_CBWR.
    // Returns 0 when CNR is in force afterwards, -1 otherwise.
    int setDeterministic(int branch, int keepEnv);
    // MKL_CBWR name ("AUTO", "COMPATIBLE", "AVX2,STRICT", ...) to its code,
    // case-insensitive; -1 for an unknown name.
    static int cbwrBranchFromName(const char *name);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag,
		 Channel &theChannel,
		 FEM_ObjectBroker &theBroker);
  protected:

  private:
	  PARDISOGenLinSOE *theSOE;

	  // The PARDISO handle persists across solves so the symbolic and numeric
	  // factorizations can be reused. pt must be zeroed before the first call
	  // and never copied.
	  void *pt[64];
	  int   iparm[64];
	  int   mtype;         // derived from the SOE matType at the symbolic phase
	  bool  init;          // phase 11 has run, so a phase -1 release is owed
	  bool  needsSymbolic; // the sparsity pattern changed
	  // The destructor must not read theSOE: ~LinearSOE() deletes the solver
	  // after the derived SOE has freed its arrays. The order is cached here.
	  int   cachedN;
	  int   reportStats;

	  // Preconditioned CGS state. haveFactors: a factorization exists in the
	  // handle. factorsCurrent: it is a factorization of the matrix now in the
	  // SOE. A CGS success leaves the previous factors in place, so a later
	  // phase-33-only solve would use an older matrix; factorsCurrent forbids it.
	  int   krylovL;
	  int   krylovK;       // 1 = CGS (mtype 11), 2 = CG (mtype 2), 0 = off
	  bool  haveFactors;
	  bool  factorsCurrent;
	  int   cgsCalls;
	  int   cgsWins;
	  bool  cgsAdviceDone;

	  int   cnrBranch;     // -1 = off (default), else the requested branch
	  bool  cnrInForce;
	  bool  cnrNoticeDone;
};

#endif

