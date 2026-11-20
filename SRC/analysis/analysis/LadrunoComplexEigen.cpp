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
// Description: complex (state-space) modal analysis on the real undamped
// modes of the last `eigen` run, and the `complexEigen` command. See
// LadrunoComplexEigen.h for the linearization, sign conventions and
// references.

#include "LadrunoComplexEigen.h"
#include "LadrunoDampingAssembler.h"

#include <Matrix.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Vector.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <NodeIter.h>
#include <elementAPI.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#ifdef _WIN32

extern "C" int DGGEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                     double *BETA, double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);

#else

extern "C" int dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                      double *BETA, double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);

#endif

namespace {

// || (lambda^2 Mt + lambda Ct + Kt) z ||_2 for a complex (lambda, z).
double quadResidual(const Matrix &Mt, const Matrix &Ct, const Matrix &Kt,
                    double lr, double li,
                    const std::vector<double> &zRe,
                    const std::vector<double> &zIm)
{
    const int p = Mt.noRows();
    // lambda^2 = (lr^2 - li^2) + i (2 lr li)
    const double l2r = lr * lr - li * li;
    const double l2i = 2.0 * lr * li;
    double norm2 = 0.0;
    for (int i = 0; i < p; i++) {
        double rr = 0.0, ri = 0.0;
        for (int j = 0; j < p; j++) {
            const double m = Mt(i, j), c = Ct(i, j), k = Kt(i, j);
            const double ar = l2r * m + lr * c + k;   // Re of the (i,j) operator
            const double ai = l2i * m + li * c;       // Im of the (i,j) operator
            rr += ar * zRe[j] - ai * zIm[j];
            ri += ar * zIm[j] + ai * zRe[j];
        }
        norm2 += rr * rr + ri * ri;
    }
    return std::sqrt(norm2);
}

} // namespace

int LadrunoComplexEigen::solveReducedPencil(const Matrix &Mt, const Matrix &Ct,
                                            const Matrix &Kt,
                                            std::vector<ComplexMode> &modes,
                                            double tol)
{
    modes.clear();

    const int p = Mt.noRows();
    if (p <= 0 || Mt.noCols() != p ||
        Ct.noRows() != p || Ct.noCols() != p ||
        Kt.noRows() != p || Kt.noCols() != p) {
        opserr << "LadrunoComplexEigen::solveReducedPencil - Mt/Ct/Kt must be "
               << "square with a common size; got " << Mt.noRows() << "x"
               << Mt.noCols() << ", " << Ct.noRows() << "x" << Ct.noCols()
               << ", " << Kt.noRows() << "x" << Kt.noCols() << "\n";
        return -1;
    }
    if (tol <= 0.0)
        tol = 1.0e-8;

    // Build the 2p x 2p pair (-B, A) column-major for LAPACK:
    //   (-B) = [-Kt 0; 0 Mt],   A = [Ct Mt; Mt 0],   (-B) w = lambda A w.
    const int n = 2 * p;
    std::vector<double> negB(static_cast<size_t>(n) * n, 0.0);
    std::vector<double> Amat(static_cast<size_t>(n) * n, 0.0);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            const size_t colTop = static_cast<size_t>(j) * n;       // columns 0..p-1
            const size_t colBot = static_cast<size_t>(j + p) * n;   // columns p..2p-1
            negB[colTop + i]     = -Kt(i, j);
            negB[colBot + (i+p)] =  Mt(i, j);
            Amat[colTop + i]     =  Ct(i, j);
            Amat[colBot + i]     =  Mt(i, j);
            Amat[colTop + (i+p)] =  Mt(i, j);
        }
    }

    std::vector<double> alphaR(n), alphaI(n), beta(n);
    std::vector<double> vr(static_cast<size_t>(n) * n);
    double vlDummy = 0.0;

    char jobvl = 'N';
    char jobvr = 'V';
    int  N = n, ldA = n, ldB = n, ldvl = 1, ldvr = n, info = 0;

    // workspace query, then solve
    int lwork = -1;
    double workQuery = 0.0;
#ifdef _WIN32
    DGGEV(&jobvl, &jobvr, &N, negB.data(), &ldA, Amat.data(), &ldB,
          alphaR.data(), alphaI.data(), beta.data(),
          &vlDummy, &ldvl, vr.data(), &ldvr, &workQuery, &lwork, &info);
#else
    dggev_(&jobvl, &jobvr, &N, negB.data(), &ldA, Amat.data(), &ldB,
           alphaR.data(), alphaI.data(), beta.data(),
           &vlDummy, &ldvl, vr.data(), &ldvr, &workQuery, &lwork, &info);
#endif
    if (info != 0) {
        opserr << "LadrunoComplexEigen::solveReducedPencil - dggev workspace "
               << "query failed, INFO = " << info << "\n";
        return -2;
    }
    lwork = static_cast<int>(workQuery);
    if (lwork < 8 * n)
        lwork = 8 * n;
    std::vector<double> work(static_cast<size_t>(lwork));

#ifdef _WIN32
    DGGEV(&jobvl, &jobvr, &N, negB.data(), &ldA, Amat.data(), &ldB,
          alphaR.data(), alphaI.data(), beta.data(),
          &vlDummy, &ldvl, vr.data(), &ldvr, work.data(), &lwork, &info);
#else
    dggev_(&jobvl, &jobvr, &N, negB.data(), &ldA, Amat.data(), &ldB,
           alphaR.data(), alphaI.data(), beta.data(),
           &vlDummy, &ldvl, vr.data(), &ldvr, work.data(), &lwork, &info);
#endif
    if (info < 0) {
        opserr << "LadrunoComplexEigen::solveReducedPencil - invalid argument "
               << -info << " passed to LAPACK dggev\n";
        return -2;
    }
    if (info > 0) {
        opserr << "LadrunoComplexEigen::solveReducedPencil - LAPACK dggev "
               << "failed to converge, INFO = " << info << " (of " << n
               << " roots)\n";
        return -3;
    }

    // scale for the beta ~ 0 (infinite eigenvalue) test
    double alphaScale = 0.0, betaScale = 0.0;
    for (int k = 0; k < n; k++) {
        alphaScale = std::max(alphaScale,
                              std::hypot(alphaR[k], alphaI[k]));
        betaScale  = std::max(betaScale, std::fabs(beta[k]));
    }
    const double betaTol = tol * std::max(betaScale, alphaScale);

    // First pass: recover lambda_k, find the frequency scale for the
    // rigid-mode (omega0 ~ 0) test.
    std::vector<double> lamRe(n, 0.0), lamIm(n, 0.0);
    std::vector<bool>   infinite(n, false);
    double omegaScale = 0.0;
    for (int k = 0; k < n; k++) {
        if (std::fabs(beta[k]) <= betaTol) {
            infinite[k] = true;
            continue;
        }
        lamRe[k] = alphaR[k] / beta[k];
        lamIm[k] = alphaI[k] / beta[k];
        omegaScale = std::max(omegaScale, std::hypot(lamRe[k], lamIm[k]));
    }
    const double omegaTol = tol * std::max(omegaScale, 1.0);

    // Second pass: classify, collapse conjugate pairs, extract z.
    // LAPACK's own pair convention drives the storage layout: alphaI(k) != 0
    // marks a conjugate pair stored in consecutive columns k (real part of v,
    // eigenvalue with POSITIVE imaginary part) and k+1 (imag part of v, the
    // conjugate). The engineering tolerance omegaTol only picks the reported
    // kind, never the column layout.
    for (int k = 0; k < n; k++) {
        ComplexMode mode;
        mode.zRe.assign(static_cast<size_t>(p), 0.0);
        mode.zIm.assign(static_cast<size_t>(p), 0.0);

        if (infinite[k]) {
            // pathological (Mt singular): tag, don't divide
            mode.omega0 = 0.0;  mode.omegaD = 0.0;  mode.zeta = 0.0;
            mode.lambdaRe = 0.0; mode.lambdaIm = 0.0;
            mode.kind = RIGID;   mode.resid = 0.0;
            for (int i = 0; i < p; i++)
                mode.zRe[static_cast<size_t>(i)] =
                    vr[static_cast<size_t>(k) * n + i];
            modes.push_back(mode);
            continue;
        }

        const double lr = lamRe[k], li = lamIm[k];
        const double mag = std::hypot(lr, li);

        if (alphaI[k] != 0.0 && k + 1 < n) {
            // conjugate pair: report once (the +Im representative),
            // consume both columns
            mode.lambdaRe = lr;
            mode.lambdaIm = std::fabs(li);
            mode.omega0 = mag;
            mode.omegaD = mode.lambdaIm;
            mode.zeta   = (mag > 0.0) ? -lr / mag : 0.0;
            // near-critically-damped pairs (|Im| below the engineering
            // tolerance) are tagged OVERDAMPED but keep the pair's vector
            mode.kind = (mode.lambdaIm > omegaTol) ? UNDERDAMPED : OVERDAMPED;
            // v = VR(:,k) + i VR(:,k+1) pairs with lambda = alpha(k)/beta(k).
            // LAPACK fixes the pair layout by alphaI(k) > 0, NOT by the sign
            // of lambda: when beta(k) < 0 the column-k eigenvalue has
            // NEGATIVE imaginary part, so pairing its vector with the
            // reported +Im lambda requires the conjugate (negate zIm).
            const double imSign = (li < 0.0) ? -1.0 : 1.0;
            for (int i = 0; i < p; i++) {
                mode.zRe[static_cast<size_t>(i)] =
                    vr[static_cast<size_t>(k) * n + i];
                mode.zIm[static_cast<size_t>(i)] =
                    imSign * vr[static_cast<size_t>(k + 1) * n + i];
            }
            k++;   // skip the conjugate column
        } else {
            // exactly-real root: rigid (lambda ~ 0) or overdamped,
            // one entry per root
            mode.lambdaRe = lr;  mode.lambdaIm = 0.0;
            mode.omegaD = 0.0;
            if (mag <= omegaTol) {
                mode.kind = RIGID;
                mode.omega0 = 0.0;
                mode.zeta = 0.0;
            } else {
                mode.kind = OVERDAMPED;
                mode.omega0 = mag;
                mode.zeta = -lr / mag;   // +1 decaying, -1 growing
            }
            for (int i = 0; i < p; i++)
                mode.zRe[static_cast<size_t>(i)] =
                    vr[static_cast<size_t>(k) * n + i];
        }

        // phase normalization: largest-|z_j| entry made real-positive
        int    jmax = 0;
        double amax = -1.0;
        for (int j = 0; j < p; j++) {
            const double a = std::hypot(mode.zRe[static_cast<size_t>(j)],
                                        mode.zIm[static_cast<size_t>(j)]);
            if (a > amax) { amax = a; jmax = j; }
        }
        if (amax > 0.0) {
            const double cr = mode.zRe[static_cast<size_t>(jmax)] / amax;
            const double ci = mode.zIm[static_cast<size_t>(jmax)] / amax;
            for (int j = 0; j < p; j++) {          // z <- z * conj(c)
                const double zr = mode.zRe[static_cast<size_t>(j)];
                const double zi = mode.zIm[static_cast<size_t>(j)];
                mode.zRe[static_cast<size_t>(j)] =  zr * cr + zi * ci;
                mode.zIm[static_cast<size_t>(j)] = -zr * ci + zi * cr;
            }
        }

        mode.resid = quadResidual(Mt, Ct, Kt, mode.lambdaRe, mode.lambdaIm,
                                  mode.zRe, mode.zIm);
        modes.push_back(mode);
    }

    std::sort(modes.begin(), modes.end(),
              [](const ComplexMode &a, const ComplexMode &b) {
                  if (a.omega0 != b.omega0)
                      return a.omega0 < b.omega0;
                  return a.lambdaRe < b.lambdaRe;
              });

    return 0;
}

int
LadrunoComplexEigen::solveFromDomain(Domain &theDomain, int numModes,
                                     std::vector<ComplexMode> &modes,
                                     double tol, bool closedForm)
{
    modes.clear();

    // presence probe FIRST: getEigenvalues() exit(-1)s when never set
    const int pAvail = theDomain.getNumEigenvalues();
    if (pAvail < 1) {
        opserr << "LadrunoComplexEigen::solveFromDomain - no eigenvalues in "
               << "the domain; run 'eigen <N>' before 'complexEigen'\n";
        return -10;
    }
    const Vector &eigenvalues = theDomain.getEigenvalues();
    const int p = (numModes <= 0) ? pAvail : numModes;
    if (p > pAvail) {
        opserr << "LadrunoComplexEigen::solveFromDomain - requested "
               << p << " modes but the last eigen run produced only "
               << pAvail << "\n";
        return -11;
    }

    int rc;
    if (!closedForm) {
        // element-by-element projection of the M and C the model carries.
        // Damping that does not flow through getDamp() is warned about, not
        // absorbed: modalDamping lives in the integrator, so the assembled C
        // cannot see it. Element 'damping' objects (folded into
        // getResistingForceIncInertia, e.g. ElasticBeam2d theDamping) and
        // HHT-alpha numerical damping are likewise invisible; no generic
        // probe exists for those.
        if (theDomain.getModalDampingFactors() != 0)
            opserr << "WARNING LadrunoComplexEigen: modalDamping is active; "
                   << "it is applied by the integrator, NOT by element "
                   << "getDamp(), so the assembled projection cannot include "
                   << "it - the reported zeta EXCLUDES the modal damping\n";
        Matrix Mt(p, p), Ct(p, p), Kt(p, p);
        if (LadrunoDampingAssembler::projectReduced(theDomain, p,
                                                    Mt, Ct, Kt) < 0)
            return -14;
        rc = solveReducedPencil(Mt, Ct, Kt, modes, tol);
        if (rc < 0)
            return rc;
        return pushComplexShapes(theDomain, p, modes);
    }

    double aM, bK, bK0, bKc;
    theDomain.getRayleighDampingFactors(aM, bK, bK0, bKc);

    // region or element Rayleigh factors ('region ... -rayleigh', per-element
    // '-rayleigh') are written straight onto the elements without touching
    // the domain copy (MeshRegion::setRayleighDampingFactors), so the closed
    // form below would use the global factors only. Detect and warn; the
    // default assembled path is the one that handles scoped damping.
    {
        bool warned = false;
        Element *elePtr;
        ElementIter &theElemIter = theDomain.getElements();
        while (!warned && (elePtr = theElemIter()) != 0) {
            const Vector f = elePtr->getRayleighDampingFactors();
            if (f.Size() >= 4 &&
                (f(0) != aM || f(1) != bK || f(2) != bK0 || f(3) != bKc)) {
                opserr << "WARNING LadrunoComplexEigen: element "
                       << elePtr->getTag() << " carries Rayleigh factors that "
                       << "differ from the last global 'rayleigh' call "
                       << "(region-/element-scoped damping). The closed-form "
                       << "projection uses the GLOBAL factors only - the "
                       << "reported zeta ignores the scoped damping. Drop "
                       << "'-closedForm': the default assembled path "
                       << "captures scoped damping exactly\n";
                warned = true;
            }
        }
    }

    if (bK0 != 0.0 || bKc != 0.0) {
        opserr << "LadrunoComplexEigen::solveFromDomain - betaKinit/betaKcomm "
               << "Rayleigh terms are set; the closed form covers only "
               << "alphaM/betaK (K_init/K_commit are not proportional to the "
               << "eigen K). Drop '-closedForm': the default assembled path "
               << "handles them - refusing rather than misreporting zeta\n";
        return -12;
    }

    // closed form: since K phi_a = omega_a^2 M phi_a, in ANY
    // normalization Phi^T M Phi = diag(g_a) and Phi^T K Phi = diag(g_a w_a^2),
    // so after the g_a^{-1/2} re-normalization the reduced pencil is exactly
    //   Mt = I,  Kt = diag(w_a^2),  Ct = diag(aM + bK w_a^2)
    // with no matrix ever assembled. Negative computed eigenvalues (possible
    // from a shifted / indefinite prior state) are refused - the closed form
    // presumes a genuine undamped spectrum.
    Matrix Mt(p, p), Ct(p, p), Kt(p, p);
    for (int a = 0; a < p; a++) {
        const double w2 = eigenvalues(a);
        if (w2 < 0.0) {
            opserr << "LadrunoComplexEigen::solveFromDomain - eigenvalue "
                   << a + 1 << " is negative (" << w2 << "); the Rayleigh "
                   << "closed form needs a genuine undamped spectrum\n";
            return -13;
        }
        Mt(a, a) = 1.0;
        Kt(a, a) = w2;
        Ct(a, a) = aM + bK * w2;
    }

    rc = solveReducedPencil(Mt, Ct, Kt, modes, tol);
    if (rc < 0)
        return rc;
    return pushComplexShapes(theDomain, p, modes);
}

int
LadrunoComplexEigen::pushComplexShapes(Domain &theDomain, int p,
                                       const std::vector<ComplexMode> &modes)
{
    // psi_k = Phi z_k, distributed node-wise: this node's psi rows are
    // nodeEig(:, 0..p-1) * z_k. Nodes without stored real eigenvectors
    // (fully fixed) hold no complex shapes either: their psi is zero.
    //
    // Phase convention: the largest-|psi| physical entry of each mode is
    // made real-positive. The kernel normalizes z (reduced space) for -qz
    // reproducibility, but psi = Phi z does not inherit that, so pass 1 finds each mode's largest-magnitude physical entry and
    // pass 2 writes the rotated psi * conj(psi_max)/|psi_max|.
    const int nModes = static_cast<int>(modes.size());
    if (nModes < 1)
        return 0;

    std::vector<double> bestMag(static_cast<size_t>(nModes), 0.0);
    std::vector<double> rotRe(static_cast<size_t>(nModes), 1.0);
    std::vector<double> rotIm(static_cast<size_t>(nModes), 0.0);

    for (int pass = 0; pass < 2; pass++) {
        Node *nodePtr;
        NodeIter &theNodeIter = theDomain.getNodes();
        while ((nodePtr = theNodeIter()) != 0) {
            const int stored = nodePtr->getNumEigenvectors();
            if (stored < p)
                continue;   // guarded earlier for analysis nodes; skip others
            const Matrix &vecs = nodePtr->getEigenvectors();  // ndf x stored
            const int ndf = nodePtr->getNumberDOF();
            if (pass == 1 && nodePtr->setNumComplexEigenvectors(nModes) < 0)
                return -20;
            Vector re(ndf), im(ndf);
            for (int k = 0; k < nModes; k++) {
                const ComplexMode &mk = modes[static_cast<size_t>(k)];
                const size_t ks = static_cast<size_t>(k);
                for (int i = 0; i < ndf; i++) {
                    double sr = 0.0, si = 0.0;
                    for (int a = 0; a < p; a++) {
                        const double phi = vecs(i, a);
                        sr += phi * mk.zRe[static_cast<size_t>(a)];
                        si += phi * mk.zIm[static_cast<size_t>(a)];
                    }
                    if (pass == 0) {
                        const double mag = std::hypot(sr, si);
                        if (mag > bestMag[ks]) {
                            bestMag[ks] = mag;
                            rotRe[ks] = sr / mag;    // e^{i arg psi_max}
                            rotIm[ks] = si / mag;
                        }
                    } else {
                        // psi * conj(rot): rotates psi_max onto the real axis
                        re(i) = sr * rotRe[ks] + si * rotIm[ks];
                        im(i) = -sr * rotIm[ks] + si * rotRe[ks];
                    }
                }
                if (pass == 1 &&
                    nodePtr->setComplexEigenvector(k + 1, re, im) < 0)
                    return -21;
            }
        }
    }
    return 0;
}

// complexEigen command. Two entries, one solver:
//   complexEigen <numModes> <-numModes $p> <-tol $eps> <-closedForm>
//     domain path: project M and C onto the modes of the last `eigen` run
//     (all of them by default) and solve the reduced pencil.
//   complexEigen -qz $p <p*p Mt> <p*p Ct> <p*p Kt> <-tol $eps>   (row-major)
//     solve an explicitly given reduced pencil (verification entry).
// Both return 7 doubles per reported mode:
//   {omega0, omegaD, zeta, Re(lambda), Im(lambda), kind, resid}
// with kind 0 = underdamped (one entry per conjugate pair), 1 = overdamped,
// 2 = rigid, and resid = ||(lambda^2 Mt + lambda Ct + Kt) z||, which ties
// each reported lambda to its own eigenvector.
int OPS_complexEigen(void)
{
    std::vector<LadrunoComplexEigen::ComplexMode> cmodes;
    int numdata = 1;

    // peek the first token (if any) to pick the entry. OPS_GetStringFromAll,
    // not OPS_GetString: under Python a bare integer argument
    // (complexEigen(2)) is not returned by OPS_GetString.
    char optBuf[64];
    const char* opt = nullptr;
    if (OPS_GetNumRemainingInputArgs() > 0)
        opt = OPS_GetStringFromAll(optBuf, 64);

    if (opt != nullptr && strcmp(opt, "-qz") == 0) {

        int p = 0;
        if (OPS_GetIntInput(&numdata, &p) < 0 || p < 1 || p > 500) {
            opserr << "WARNING complexEigen -qz: cannot read a valid pencil size p "
                   << "(1..500)\n";
            return -1;
        }

        const int pp = p * p;
        if (OPS_GetNumRemainingInputArgs() < 3 * pp) {
            opserr << "WARNING complexEigen -qz " << p << ": expected " << 3 * pp
                   << " matrix entries (row-major Mt, Ct, Kt), got "
                   << OPS_GetNumRemainingInputArgs() << "\n";
            return -1;
        }

        Matrix Mt(p, p), Ct(p, p), Kt(p, p);
        Matrix* mats[3] = {&Mt, &Ct, &Kt};
        for (int m = 0; m < 3; m++) {
            std::vector<double> flat(static_cast<size_t>(pp));
            numdata = pp;
            if (OPS_GetDoubleInput(&numdata, flat.data()) < 0) {
                opserr << "WARNING complexEigen -qz: cannot read matrix " << m
                       << " entries\n";
                return -1;
            }
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                    (*mats[m])(i, j) = flat[static_cast<size_t>(i) * p + j];
        }

        double tol = 1.0e-8;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char* flag = OPS_GetString();
            if (strcmp(flag, "-tol") == 0 &&
                OPS_GetNumRemainingInputArgs() > 0) {
                numdata = 1;
                if (OPS_GetDoubleInput(&numdata, &tol) < 0) {
                    opserr << "WARNING complexEigen -qz: cannot read -tol value\n";
                    return -1;
                }
            } else {
                opserr << "WARNING complexEigen -qz: unknown option '" << flag
                       << "'\n";
                return -1;
            }
        }

        if (LadrunoComplexEigen::solveReducedPencil(Mt, Ct, Kt, cmodes, tol) < 0) {
            opserr << "WARNING complexEigen -qz: reduced-pencil solve failed\n";
            return -1;
        }

    } else {

        // domain path. Default: assembled element-by-element projection
        // (scoped Rayleigh, betaKinit/betaKcomm, damper materials and nodal
        // alphaM all captured); -closedForm: global Rayleigh closed form.
        int numModes = 0;   // 0 = all modes from the last eigen run
        double tol = 1.0e-8;
        bool closedForm = false;
        while (true) {
            if (opt == nullptr) {
                if (OPS_GetNumRemainingInputArgs() <= 0) break;
                opt = OPS_GetStringFromAll(optBuf, 64);
            }
            if (strcmp(opt, "-closedForm") == 0) {
                closedForm = true;
            } else if (strcmp(opt, "-numModes") == 0 &&
                       OPS_GetNumRemainingInputArgs() > 0) {
                numdata = 1;
                if (OPS_GetIntInput(&numdata, &numModes) < 0 || numModes < 1) {
                    opserr << "WARNING complexEigen: cannot read a valid "
                           << "-numModes value\n";
                    return -1;
                }
            } else if (strcmp(opt, "-tol") == 0 &&
                       OPS_GetNumRemainingInputArgs() > 0) {
                numdata = 1;
                if (OPS_GetDoubleInput(&numdata, &tol) < 0) {
                    opserr << "WARNING complexEigen: cannot read -tol value\n";
                    return -1;
                }
            } else {
                // bare positive integer = numModes (as in `eigen`)
                char* end = nullptr;
                const long v = strtol(opt, &end, 10);
                if (end != nullptr && *end == '\0' && v > 0) {
                    numModes = static_cast<int>(v);
                } else {
                    opserr << "WARNING complexEigen: unknown option '" << opt
                           << "'\n";
                    return -1;
                }
            }
            opt = nullptr;
        }

        Domain* theDomain = OPS_GetDomain();
        if (theDomain == 0) {
            opserr << "WARNING complexEigen: no domain\n";
            return -1;
        }
        if (LadrunoComplexEigen::solveFromDomain(*theDomain, numModes,
                                                 cmodes, tol, closedForm) < 0)
            return -1;
    }

    std::vector<double> out;
    out.reserve(cmodes.size() * 7);
    for (const LadrunoComplexEigen::ComplexMode &mk : cmodes) {
        out.push_back(mk.omega0);
        out.push_back(mk.omegaD);
        out.push_back(mk.zeta);
        out.push_back(mk.lambdaRe);
        out.push_back(mk.lambdaIm);
        out.push_back(static_cast<double>(mk.kind));
        out.push_back(mk.resid);
    }

    int numOut = static_cast<int>(out.size());
    if (OPS_SetDoubleOutput(&numOut, out.data(), false) < 0) {
        opserr << "WARNING complexEigen: failed to set output\n";
        return -1;
    }
    return 0;
}
