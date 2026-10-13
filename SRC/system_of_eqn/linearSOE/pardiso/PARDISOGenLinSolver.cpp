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
// What: "@(#) PARDISOGenLinSolver.cpp, revA"
//
// The PARDISO phases are split so that work is reused across solves:
//   phase 11 (reordering, symbolic)  once per sparsity pattern (setSize)
//   phase 22 (numerical factorization) only when A has changed
//   phase 33 (solve, refinement)      every call
//   phase -1 (release)                once, in the destructor
// Algorithms that keep the tangent (ModifiedNewton, Initial, KrylovNewton)
// then pay one factorization per tangent instead of one per solve.
//
// Reference: Intel oneMKL Developer Reference, "pardiso" and "pardiso iparm
// Parameter".


#include <PARDISOGenLinSolver.h>
#include <PARDISOGenLinSOE.h>
#include <math.h>
#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <mkl_pardiso.h>
#include <mkl_types.h>
#include <mkl_service.h>

// MKL CBWR branches, spelled as for the MKL_CBWR environment variable
// (Intel oneMKL Developer Guide, "Obtaining Numerically Reproducible Results").
static const struct { const char *name; int code; } pardiso_cbwr_table[] = {
	{"OFF",           MKL_CBWR_OFF},
	{"BRANCH_OFF",    MKL_CBWR_BRANCH_OFF},
	{"AUTO",          MKL_CBWR_AUTO},
	{"COMPATIBLE",    MKL_CBWR_COMPATIBLE},
	{"SSE2",          MKL_CBWR_SSE2},
	{"SSE3",          MKL_CBWR_SSE3},
	{"SSSE3",         MKL_CBWR_SSSE3},
	{"SSE4_1",        MKL_CBWR_SSE4_1},
	{"SSE4_2",        MKL_CBWR_SSE4_2},
	{"AVX",           MKL_CBWR_AVX},
	{"AVX2",          MKL_CBWR_AVX2},
	{"AVX512_MIC",    MKL_CBWR_AVX512_MIC},
	{"AVX512",        MKL_CBWR_AVX512},
	{"AVX512_MIC_E1", MKL_CBWR_AVX512_MIC_E1},
	{"AVX512_E1",     MKL_CBWR_AVX512_E1},
#ifdef MKL_CBWR_AVX10
	{"AVX10",         MKL_CBWR_AVX10},   // oneMKL 2025.0 and later
#endif
};

static const char *
pardiso_cbwr_name(int code)
{
	for (const auto &e : pardiso_cbwr_table)
		if (e.code == code) return e.name;
	return "UNKNOWN";
}

// case-insensitive compare of exactly n characters
static bool
pardiso_cbwr_ieq(const char *a, const char *b, size_t n)
{
	for (size_t i = 0; i < n; i++) {
		char ca = a[i], cb = b[i];
		if (ca >= 'a' && ca <= 'z') ca = ca - 'a' + 'A';
		if (cb >= 'a' && cb <= 'z') cb = cb - 'a' + 'A';
		if (ca != cb) return false;
	}
	return true;
}

PARDISOGenLinSolver::PARDISOGenLinSolver()
:LinearSOESolver(SOLVER_TAGS_PARDISOGenLinSolver),
 theSOE(0), mtype(11), init(false), needsSymbolic(false), cachedN(0),
 reportStats(0),
 krylovL(0), krylovK(0), haveFactors(false), factorsCurrent(false),
 cgsCalls(0), cgsWins(0), cgsAdviceDone(false),
 cnrBranch(-1), cnrInForce(false), cnrNoticeDone(false)
{
	for (int i = 0; i < 64; i++) {
		pt[i] = 0;
		iparm[i] = 0;
	}
}


PARDISOGenLinSolver::~PARDISOGenLinSolver()
{
	// Release PARDISO's memory once. theSOE is not touched here (see the
	// header); the release phase does not read the matrix.
	if (init == true) {
		int maxfct = 1, mnum = 1, phase = -1, error = 0, msglvl = 0, nrhs = 1;
		int n = cachedN;
		double ddum = 0.0; int idum = 0;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
			&n, &ddum, &idum, &idum, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
		init = false;
	}
}


static void
pardiso_report(const char *whatPhase, int error, int mtype)
{
	opserr << "WARNING PARDISOGenLinSolver::solve() - error " << error
	       << " during " << whatPhase << " (mtype " << mtype << "): ";
	switch (error) {
	case  -1: opserr << "input inconsistent\n"; break;
	case  -2: opserr << "not enough memory\n"; break;
	case  -3: opserr << "reordering problem\n"; break;
	case  -4:
		opserr << "zero pivot, numerical factorization or iterative "
		          "refinement problem\n";
		if (mtype == 2)
			opserr << "     -matrixType 1 assumes a positive definite matrix; "
			          "use -matrixType 2 for an indefinite\n     (softening "
			          "or buckling) tangent.\n";
		break;
	case  -5: opserr << "unclassified internal error\n"; break;
	case  -6: opserr << "reordering failed\n"; break;
	case  -7: opserr << "diagonal matrix is singular\n"; break;
	case  -8: opserr << "32-bit integer overflow problem\n"; break;
	case -10: opserr << "error opening OOC files\n"; break;
	default:  opserr << "see the oneMKL PARDISO error table\n"; break;
	}
}


// iparm[9] makes PARDISO replace pivots smaller than 10^-iparm[9]*||A||
// instead of failing, so a nearly singular matrix returns error 0 and the
// solution of a perturbed matrix. Report it once per factorization.
static void
pardiso_perturbed(const int *iparm, int mtype)
{
	if (iparm[13] > 0)
		opserr << "WARNING PARDISOGenLinSolver: PARDISO perturbed "
		       << iparm[13] << " pivot(s) during factorization (mtype "
		       << mtype << ", threshold 1e-" << iparm[9]
		       << "); the matrix is nearly singular\n";
}


int
PARDISOGenLinSolver::solve(void)
{
	if (theSOE == 0) {
		opserr << "WARNING PARDISOLinSolver::solve(void)- ";
		opserr << " No LinearSOE object has been set\n";
		return -1;
	}

	int     n    = theSOE->size;
	int    *ia   = theSOE->rowStartA;   // one-based CSR row pointers
	int    *ja   = theSOE->colA;        // one-based CSR column indices
	double *a    = theSOE->A;
	double *Xptr = theSOE->X;
	double *Bptr = theSOE->B;

	if (n <= 0 || ia == 0 || ja == 0 || a == 0) {
		opserr << "WARNING PARDISOGenLinSolver::solve() - the SOE has no "
		          "equations (n=" << n << "); nothing to factor\n";
		return -1;
	}

	int maxfct = 1, mnum = 1, msglvl = 0, nrhs = 1, error = 0;
	double ddum = 0.0; int idum = 0;

	// ---- reordering and symbolic factorization, once per pattern ----------
	if (needsSymbolic == true || init == false) {

		if (init == true) {   // discard the previous pattern first
			int phase = -1;
			PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja,
				&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
			init = false;
			for (int i = 0; i < 64; i++) pt[i] = 0;
		}

		// mtype follows the storage the SOE built, so the two cannot disagree
		const int soeMatType = theSOE->getMatType();
		mtype = (soeMatType == 1) ? 2 : (soeMatType == 2 ? -2 : 11);
		const bool symmetric = (mtype != 11);

		// iparm is rebuilt for every new pattern; any option stored in it
		// must be re-applied here.
		for (int i = 0; i < 64; i++) iparm[i] = 0;
		iparm[0]  =  1;  /* use the values below, not the solver defaults */
		iparm[1]  =  2;  /* nested dissection from METIS */
		iparm[3]  =  0;  /* no iterative-direct algorithm */
		iparm[4]  =  0;  /* no user permutation */
		iparm[5]  =  0;  /* solution in x, b unchanged */
		iparm[7]  =  2;  /* maximum iterative refinement steps */
		/* Pivoting and scaling use Intel's documented defaults per mtype:
		   unsymmetric: perturbation 1e-13, scaling and weighted matching on;
		   symmetric:   perturbation 1e-8, Bunch-Kaufman 1x1/2x2 pivoting,
		                scaling and matching off. */
		iparm[9]  = symmetric ?  8 : 13;
		iparm[10] = symmetric ?  0 :  1;
		iparm[12] = symmetric ?  0 :  1;
		if (symmetric)
			iparm[20] = 1;  /* Bunch-Kaufman pivoting for indefinite matrices */
		/* nonzeros in the factors and factorization Mflops are reported only
		   when requested with a negative value; the flop count costs extra
		   analysis time, so only with -stats. */
		iparm[17] = reportStats ? -1 : 0;
		iparm[18] = reportStats ? -1 : 0;
		iparm[34] =  0;  /* one-based indexing */
		/* CNR: iparm[33] > 0 is the thread count PARDISO reproduces results
		   for (in-core mode). METIS (iparm[1] = 2) is compatible with it. */
		if (cnrBranch >= 0)
			iparm[33] = mkl_get_max_threads();

		int phase = 11;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			pardiso_report("symbolic factorization", error, mtype);
			return -1;
		}

		init = true;
		needsSymbolic = false;
		cachedN = n;

		// report once what MKL says is in force, not what was requested
		if (cnrBranch >= 0 && cnrNoticeDone == false) {
			cnrNoticeDone = true;
			const int inForce = mkl_cbwr_get(MKL_CBWR_BRANCH);
			const int all     = mkl_cbwr_get(MKL_CBWR_ALL);
			cnrInForce = (inForce > MKL_CBWR_BRANCH_OFF);
			opserr << "PARDISO deterministic mode: MKL CNR branch "
			       << pardiso_cbwr_name(inForce);
			if (inForce == MKL_CBWR_AUTO) {
				const int resolved = mkl_cbwr_get_auto_branch();
				if (resolved > MKL_CBWR_AUTO)
					opserr << " (-> " << pardiso_cbwr_name(resolved) << ")";
			}
			if (all & MKL_CBWR_STRICT)
				opserr << ",STRICT";
			opserr << ", iparm(34)=" << iparm[33] << " thread(s), CNR "
			       << (cnrInForce ? "ACTIVE" : "NOT ACTIVE; results are not "
			           "guaranteed reproducible, relaunch with MKL_CBWR=AUTO")
			       << "\n";
		}
		theSOE->factored = false;   // a new pattern always needs phase 22

		haveFactors = false;
		factorsCurrent = false;
		// Intel documents K=1 (CGS) for unsymmetric and K=2 (CG) for
		// symmetric positive definite matrices only; there is no mode for
		// mtype -2.
		if (krylovL > 0) {
			krylovK = (mtype == 11) ? 1 : (mtype == 2 ? 2 : 0);
			if (krylovK == 0)
				opserr << "WARNING PARDISOGenLinSolver: -krylov is not available "
				          "for -matrixType 2 (mtype -2); PARDISO documents the "
				          "preconditioned iteration for mtype 11 and 2 only. "
				          "Continuing with direct factorization.\n";
		}
	}

	// ---- preconditioned CGS with the retained factors (phase 23) ----------
	// Used when the stored factors are not a factorization of the current A:
	// either A was reassembled, or the last solve was a CGS success. Phase 23
	// is required: PARDISO refactors automatically when the iteration fails
	// only in phase 23.
	bool solvedByKrylov = false;
	bool didFactorNow = false;

	if (krylovK != 0 && haveFactors == true &&
	    (theSOE->factored == false || factorsCurrent == false)) {

		iparm[3] = 10 * krylovL + krylovK;
		int phase = 23;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, Bptr, Xptr, &error);
		iparm[3] = 0;

		if (error != 0) {
			factorsCurrent = false;
			pardiso_report("CGS solve (phase 23)", error, mtype);
			return -2;
		}

		cgsCalls++;
		solvedByKrylov = true;
		theSOE->factored = true;

		if (iparm[19] > 0) {
			// converged; the stored factors are still those of an older A
			cgsWins++;
			factorsCurrent = false;
			if (reportStats)
				opserr << "PARDISO -krylov: CGS converged in " << iparm[19]
				       << " iteration(s); factorization reused\n";
		} else {
			// PARDISO refactored; iparm[19] = -it_cgs*10 - cgs_error
			factorsCurrent = true;
			const int cgsErr = -iparm[19] % 10;
			const int cgsIts = -iparm[19] / 10;
			if (reportStats)
				opserr << "PARDISO -krylov: CGS stopped after " << cgsIts
				       << " iteration(s) (cgs_error " << cgsErr
				       << "); PARDISO refactored  [" << cgsWins << "/"
				       << cgsCalls << " CGS successes so far]\n";
			// cgs_error 5: factorization is fast enough that CGS does not pay
			if (cgsErr == 5 && cgsAdviceDone == false) {
				cgsAdviceDone = true;
				opserr << "WARNING PARDISOGenLinSolver: PARDISO reports "
				          "cgs_error 5 (the factorization is fast enough that "
				          "CGS is not worthwhile); consider dropping -krylov. "
				          "Reported once.\n";
			}
			pardiso_perturbed(iparm, mtype);
		}
	}

	// ---- numerical factorization, only when A changed ---------------------
	else if (theSOE->factored == false) {
		int phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			pardiso_report("numerical factorization", error, mtype);
			return -2;
		}
		theSOE->factored = true;
		haveFactors = true;
		factorsCurrent = true;
		didFactorNow = true;
		pardiso_perturbed(iparm, mtype);
	}

	// ---- back substitution and iterative refinement, every call -----------
	if (solvedByKrylov == false) {
		int phase = 33;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
			&idum, &nrhs, iparm, &msglvl, Bptr, Xptr, &error);
		if (error != 0) {
			pardiso_report("solution", error, mtype);
			return -3;
		}
	}

	// ---- -stats, after every numerical factorization ----------------------
	// Read after phase 33 because iparm(17) is the peak over factorization
	// and solution. Labels use Intel's one-based iparm() numbering.
	if (reportStats && didFactorNow) {
		opserr << "PARDISO stats: n=" << n << " nnz(A)=" << theSOE->nnz
		       << " matrixType=" << mtype
		       << " threads=" << mkl_get_max_threads() << "\n";
		opserr << "  factor entries iparm(18)  = " << iparm[17] << "\n";
		opserr << "  peak memory KB iparm(15)  = " << iparm[14] << "\n";
		opserr << "  perm memory KB iparm(16)  = " << iparm[15] << "\n";
		opserr << "  fact memory KB iparm(17)  = " << iparm[16] << "\n";
		opserr << "  factor Mflops  iparm(19)  = " << iparm[18] << "\n";
		opserr << "  perturbed pivots = " << iparm[13]
		       << "   refinement steps = " << iparm[6] << "\n";
	}

	return 0;
}


int
PARDISOGenLinSolver::setSize()
{
	// the pattern changed; phase 11 runs on the next solve()
	needsSymbolic = true;
	return 0;
}


int
PARDISOGenLinSolver::setLinearSOE(PARDISOGenLinSOE &theLinearSOE)
{
    theSOE = &theLinearSOE;
    // state derived from the previous SOE no longer applies
    needsSymbolic = true;
    haveFactors = false;
    factorsCurrent = false;
    return 0;
}


void
PARDISOGenLinSolver::setStats(int on)
{
    reportStats = on;
}


void
PARDISOGenLinSolver::setKrylov(int digits)
{
    // iparm[3] = 10*L + K must stay non-negative
    if (digits < 0) {
        opserr << "WARNING PARDISOGenLinSolver::setKrylov() - negative digits ("
               << digits << "); -krylov disabled\n";
        digits = 0;
    }
    // PARDISO caps the iteration at 150 steps; a tolerance tighter than the
    // iteration can reach only forces the refactorization fallback.
    if (digits > 9) {
        opserr << "WARNING PARDISOGenLinSolver::setKrylov() - digits (" << digits
               << ") clamped to 9\n";
        digits = 9;
    }
    krylovL = digits;
    // K is otherwise set only at the symbolic phase; clear it when disabling
    if (krylovL == 0)
        krylovK = 0;
}


int
PARDISOGenLinSolver::cbwrBranchFromName(const char *name)
{
	if (name == 0) return -1;
	// optional ",STRICT" suffix, as in MKL_CBWR=AVX2,STRICT
	const char *comma = strchr(name, ',');
	const size_t len = comma ? (size_t)(comma - name) : strlen(name);
	int strict = 0;
	if (comma) {
		if (strlen(comma + 1) != 6 || !pardiso_cbwr_ieq(comma + 1, "STRICT", 6))
			return -1;
		strict = MKL_CBWR_STRICT;
	}
	for (const auto &e : pardiso_cbwr_table) {
		// OFF and BRANCH_OFF can be read back but not requested
		if (e.code <= MKL_CBWR_BRANCH_OFF) continue;
		if (strlen(e.name) == len && pardiso_cbwr_ieq(name, e.name, len))
			return e.code | strict;
	}
	return -1;
}


int
PARDISOGenLinSolver::setDeterministic(int branch, int keepEnv)
{
	cnrBranch = branch;
	cnrNoticeDone = false;
	if (branch < 0) {          // off: MKL's process-wide mode is left alone
		cnrInForce = false;
		return 0;
	}

	const int current = mkl_cbwr_get(MKL_CBWR_ALL);
	const bool alreadyOn = (mkl_cbwr_get(MKL_CBWR_BRANCH) > MKL_CBWR_BRANCH_OFF);

	// Keep a branch fixed by MKL_CBWR, and do not re-set the value already in
	// force: a second set after MKL has computed can fail.
	if ((keepEnv && alreadyOn) || current == branch) {
		cnrInForce = alreadyOn;
		return 0;
	}

	const int rc = mkl_cbwr_set(branch);
	cnrInForce = (mkl_cbwr_get(MKL_CBWR_BRANCH) > MKL_CBWR_BRANCH_OFF);
	if (rc == MKL_CBWR_SUCCESS && cnrInForce)
		return 0;

	const char *want = pardiso_cbwr_name(branch & ~MKL_CBWR_STRICT);
	opserr << "WARNING system Pardiso -deterministic: mkl_cbwr_set(" << want
	       << ") failed (rc " << rc;
	switch (rc) {
	case MKL_CBWR_ERR_MODE_CHANGE_FAILURE:
		opserr << "): the CNR mode is process-wide and MKL refuses to change "
		          "it once its BLAS/LAPACK dispatch is initialized (for "
		          "example by an earlier eigen solve).\n     Set the "
		          "environment variable MKL_CBWR=" << want
		       << " before OpenSees starts.\n";
		break;
	case MKL_CBWR_ERR_UNSUPPORTED_BRANCH:
		opserr << "): this CPU cannot run that code branch; the instruction-"
		          "set branches require an Intel CPU. Use AUTO, or COMPATIBLE "
		          "for a branch every x86 CPU can run.\n";
		break;
	default:
		opserr << "): see the oneMKL CBWR error codes).\n";
		break;
	}
	opserr << "     MKL CNR branch in force: "
	       << pardiso_cbwr_name(mkl_cbwr_get(MKL_CBWR_BRANCH))
	       << "; results are not guaranteed reproducible.\n";
	return -1;
}


int
PARDISOGenLinSolver::sendSelf(int cTAg, Channel &theChannel)
{
    // doing nothing
    return 0;
}


int
PARDISOGenLinSolver::recvSelf(int cTag,
			     Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}
