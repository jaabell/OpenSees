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
// Description: per-element critical time step for explicit integrators.
// See CriticalTimeStep.h for the method and its references.

#include <CriticalTimeStep.h>
#include <LadrunoMassLumping.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <OPS_Globals.h>

#include <cmath>
#include <limits>

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#include <mpi.h>
#endif

// LAPACK: general (DGGEV) and symmetric-definite expert (DSYGVX) generalized
// eigensolvers. DSYGVX is used instead of DSYGV because the reference LAPACK
// in OTHER/LAPACK provides dsygvx.f but not dsygv.f (SymmGeneralizedEigenSolver
// links the same routine). With RANGE='I' and IL=IU=N it returns only the
// largest eigenvalue.
#ifdef _WIN32
extern "C" int DGGEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                     double *BETA, double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);
extern "C" int DSYGVX(int *ITYPE, char *JOBZ, char *RANGE, char *UPLO,
                      int *N, double *A, int *LDA, double *B, int *LDB,
                      double *VL, double *VU, int *IL, int *IU,
                      double *ABSTOL, int *M, double *W, double *Z,
                      int *LDZ, double *WORK, int *LWORK, int *IWORK,
                      int *IFAIL, int *INFO);
#define LAPACK_DGGEV DGGEV
#define LAPACK_DSYGVX DSYGVX
#else
extern "C" int dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                      double *BETA, double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);
extern "C" int dsygvx_(int *ITYPE, char *JOBZ, char *RANGE, char *UPLO,
                       int *N, double *A, int *LDA, double *B, int *LDB,
                       double *VL, double *VU, int *IL, int *IU,
                       double *ABSTOL, int *M, double *W, double *Z,
                       int *LDZ, double *WORK, int *LWORK, int *IWORK,
                       int *IFAIL, int *INFO);
#define LAPACK_DGGEV dggev_
#define LAPACK_DSYGVX dsygvx_
#endif

// Largest generalized eigenvalue of K v = lambda M v, with M diagonal and both
// matrices packed column-major. Returns < 0 on failure. DSYGVX is tried first
// when every diagonal mass is positive; otherwise, or if DSYGVX fails, DGGEV
// is used.
static double maxGeneralizedEigenvalue(int n, double *K_data, double *M_data)
{
    // attempt 1: DSYGVX (M positive definite). DSYGVX overwrites its
    // arguments, so it works on copies and K_data/M_data stay intact for the
    // DGGEV fallback.
    bool mPositive = true;
    for (int i = 0; i < n; ++i)
        if (M_data[i * n + i] <= 0.0) { mPositive = false; break; }

    if (mPositive) {
        double *A = new double[n * n];
        double *B = new double[n * n];
        double *w = new double[n];
        double *z = new double[n];        // eigenvectors not requested
        int    *iwork = new int[5 * n];
        int    *ifail = new int[n];
        for (int k = 0; k < n * n; ++k) { A[k] = K_data[k]; B[k] = M_data[k]; }

        int itype = 1;      // A x = lambda B x
        char jobz = 'N';    // eigenvalues only
        char range = 'I';   // selected by index
        char uplo = 'U';
        int lda = n, ldb = n, ldz = n;
        int il = n, iu = n; // largest eigenvalue only (1-based, ascending)
        double vl = 0.0, vu = 0.0;
        double abstol = 0.0;
        int m = 0, info = 0;

        int lwork = -1; double wkopt = 0.0;
        LAPACK_DSYGVX(&itype, &jobz, &range, &uplo, &n, A, &lda, B, &ldb,
                      &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz,
                      &wkopt, &lwork, iwork, ifail, &info);
        lwork = (info == 0) ? (int)wkopt : (8 * n);
        if (lwork < 8 * n) lwork = 8 * n;
        double *work = new double[lwork];
        LAPACK_DSYGVX(&itype, &jobz, &range, &uplo, &n, A, &lda, B, &ldb,
                      &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz,
                      work, &lwork, iwork, ifail, &info);

        double lambdaMax = -1.0;
        if (info == 0 && m >= 1)
            lambdaMax = w[0];

        delete[] A; delete[] B; delete[] w; delete[] z;
        delete[] iwork; delete[] ifail; delete[] work;
        if (info == 0 && m >= 1)
            return lambdaMax;
        // B was not positive definite after all: fall through to DGGEV
    }

    // attempt 2: DGGEV (general; tolerates an indefinite or singular M)
    double *A = new double[n * n];
    double *B = new double[n * n];
    for (int k = 0; k < n * n; ++k) { A[k] = K_data[k]; B[k] = M_data[k]; }

    double *alphar = new double[n];
    double *alphai = new double[n];
    double *beta   = new double[n];
    double *vl = 0, *vr = 0;
    char jobvl = 'N', jobvr = 'N';
    int lda = n, ldb = n, info = 0;

    int lwork = -1; double wkopt = 0.0;
    LAPACK_DGGEV(&jobvl, &jobvr, &n, A, &lda, B, &ldb,
                 alphar, alphai, beta, vl, &lda, vr, &ldb, &wkopt, &lwork, &info);
    lwork = (info == 0) ? (int)wkopt : (8 * n);
    if (lwork < 1) lwork = 1;
    double *work = new double[lwork];
    LAPACK_DGGEV(&jobvl, &jobvr, &n, A, &lda, B, &ldb,
                 alphar, alphai, beta, vl, &lda, vr, &ldb, work, &lwork, &info);

    double lambdaMax = -1.0;
    if (info == 0) {
        // A small |beta| relative to the largest one is a near-infinite
        // eigenvalue of a massless mode and is discarded.
        double betaMax = 0.0;
        for (int i = 0; i < n; ++i)
            if (std::fabs(beta[i]) > betaMax) betaMax = std::fabs(beta[i]);
        const double betaTol = 1.0e-12 * betaMax;

        for (int i = 0; i < n; ++i) {
            if (std::fabs(beta[i]) > betaTol) {
                // the symmetric pencil has real eigenvalues; a complex pair
                // is round-off of the general driver and is skipped
                if (std::fabs(alphai[i]) > 1.0e-8 * std::fabs(alphar[i]))
                    continue;
                double lambda = alphar[i] / beta[i];
                if (lambda > lambdaMax) lambdaMax = lambda;
            }
        }
    }

    delete[] A; delete[] B;
    delete[] alphar; delete[] alphai; delete[] beta; delete[] work;
    return lambdaMax;
}

void lumpElementMass(Element *ele, const Matrix &M, CTSLumping lumping,
                     double *mdiag)
{
    int n = M.noRows();
    if (n <= 0) return;

    if (lumping == CTSLumping::Diagonal) {
        for (int i = 0; i < n; ++i)
            mdiag[i] = M(i, i);
    } else if (lumping == CTSLumping::HRZ) {
        // DOF directions from the node layout: the first ndm DOFs of each node
        // (ndm from that node's coordinates) are translational, the rest
        // (rotations, pressure) are tagged -1
        ID dofDir(n);
        Node **nodes = ele->getNodePtrs();
        int numNodes = ele->getNumExternalNodes();
        int pos = 0;
        bool built = (nodes != 0);
        for (int a = 0; a < numNodes && built && pos < n; ++a) {
            Node *nd = nodes[a];
            if (nd == 0) { built = false; break; }
            int ndmA = nd->getCrds().Size();
            int ndfA = nd->getNumberDOF();
            for (int k = 0; k < ndfA && pos < n; ++k, ++pos)
                dofDir(pos) = (k < ndmA) ? k : -1;
        }
        if (built && pos == n) {
            int rc = Ladruno::hrzLump(M, dofDir, mdiag, n);
            if (rc != Ladruno::HRZ_OK)
                opserr << "WARNING CriticalTimeStep - element " << ele->getTag()
                       << ": HRZ lump fell back to diagonal-of-consistent\n";
        } else {
            // DOF layout not inferable (element DOFs other than the node DOFs)
            opserr << "WARNING CriticalTimeStep - element " << ele->getTag()
                   << ": -lump hrz could not tag DOFs, using diagonal\n";
            for (int i = 0; i < n; ++i) mdiag[i] = M(i, i);
        }
    } else {
        // row sum
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j)
                sum += M(i, j);
            mdiag[i] = sum;
        }
    }
}

double elementLambdaMax(Element *ele, bool useTangent,
                        const double *mdiag, int n)
{
    const Matrix &K = useTangent ? ele->getTangentStiff() : ele->getInitialStiff();
    if (K.noRows() != n) return -1.0;

    // column-major packing for LAPACK
    double *K_data = new double[n * n];
    double *M_data = new double[n * n];
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            K_data[j * n + i] = K(i, j);
            M_data[j * n + i] = (i == j) ? mdiag[i] : 0.0;
        }
    }

    double lambdaMax = maxGeneralizedEigenvalue(n, K_data, M_data);

    delete[] K_data;
    delete[] M_data;
    return lambdaMax;
}

CTSResult computeCriticalTimeStep(AnalysisModel *theModel,
                                  bool useTangent,
                                  CTSLumping lumping)
{
    CTSResult r;
    r.damped_dt      = std::numeric_limits<double>::infinity();
    r.undamped_dt    = std::numeric_limits<double>::infinity();
    r.damped_tag     = -1;
    r.undamped_tag   = -1;
    r.n_contributing = 0;
    r.n_scanned      = 0;

    if (theModel == 0) return r;
    Domain *theDomain = theModel->getDomainPtr();
    if (theDomain == 0) return r;

    Element *ele;
    ElementIter &elements = theDomain->getElements();

    while ((ele = elements()) != 0) {
        // An element may carry a stability bound that its own pencil cannot
        // express (a penalty coupling whose DOFs are massless in the element
        // problem gives lambda_max = 0). A positive self-reported step
        // governs that element and the eigensolve is skipped; such an element
        // carries no Rayleigh damping of its own, so damped and undamped
        // values coincide.
        double selfDt = ele->getExplicitCriticalTimeStep();
        if (selfDt > 0.0) {
            r.n_scanned++;
            r.n_contributing++;
            if (selfDt < r.undamped_dt) { r.undamped_dt = selfDt; r.undamped_tag = ele->getTag(); }
            if (selfDt < r.damped_dt)   { r.damped_dt   = selfDt; r.damped_tag   = ele->getTag(); }
            continue;
        }

        // mass first, copied to a diagonal, before the stiffness is fetched
        // (see the aliasing note in CriticalTimeStep.h)
        const Matrix &M = ele->getMass();
        int n = M.noRows();
        if (n == 0) { r.n_scanned++; continue; }

        double *mdiag = new double[n];
        lumpElementMass(ele, M, lumping, mdiag);
        double lambdaMax = elementLambdaMax(ele, useTangent, mdiag, n);
        delete[] mdiag;
        r.n_scanned++;

        if (lambdaMax > 0.0) {
            double w_max = std::sqrt(lambdaMax);

            Vector coefs = ele->getRayleighDampingFactors();
            double alphaM = coefs(0);
            double betaK  = stiffnessRayleighBeta(ele);
            double xi = 0.5 * (alphaM / w_max + betaK * w_max);

            double undamped_dt = 2.0 / w_max;
            double damped_dt   = 2.0 / w_max * (std::sqrt(1.0 + xi * xi) - xi);

            r.n_contributing++;
            if (damped_dt < r.damped_dt) {
                r.damped_dt  = damped_dt;
                r.damped_tag = ele->getTag();
            }
            if (undamped_dt < r.undamped_dt) {
                r.undamped_dt  = undamped_dt;
                r.undamped_tag = ele->getTag();
            }
        }
    }

    // Parallel: the governing step is the minimum over all processes. The
    // element tags are process-local; a process whose value is not the global
    // minimum reports an unknown tag (-1).
#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    {
        double local[2] = { r.damped_dt, r.undamped_dt };
        double global[2] = { local[0], local[1] };
        MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        if (global[0] != r.damped_dt)   r.damped_tag   = -1;
        if (global[1] != r.undamped_dt) r.undamped_tag = -1;
        r.damped_dt   = global[0];
        r.undamped_dt = global[1];
    }
#endif

    return r;
}

double stiffnessRayleighBeta(Element *ele)
{
    Vector coefs = ele->getRayleighDampingFactors();  // (aM, bK, bK0, bKc)
    double beta = 0.0;
    if (coefs(1) > 0.0) beta += coefs(1);
    if (coefs(2) > 0.0) beta += coefs(2);
    if (coefs(3) > 0.0) beta += coefs(3);
    return beta;
}

bool vectorIsFinite(const Vector &v)
{
    for (int i = 0; i < v.Size(); i++)
        if (!std::isfinite(v(i)))
            return false;
    return true;
}
