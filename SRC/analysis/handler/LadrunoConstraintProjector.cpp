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
// References:
//   J.O. Hallquist, LS-DYNA Theory Manual, Livermore Software Technology
//   Corporation, 2006, Sections 28.1 and 28.2.
//
// Description: implementation of LadrunoConstraintProjector (see the header).

#include <LadrunoConstraintProjector.h>
#include <LinearSOE.h>
#include <DiagonalSOE.h>
#include <OPS_Globals.h>
#include <cmath>
#include <vector>
#include <limits>

// Condition-number limits for L^T M L of a group. The matrix is symmetric
// positive definite for positive masses and a full-rank L; above the refuse
// limit the projection would amplify round-off into the acceleration.
#define LADRUNO_LTML_COND_REFUSE 1.0e12
#define LADRUNO_LTML_COND_WARN   1.0e8

// Eigenvalues of a small symmetric matrix A by the cyclic Jacobi method.
// Fills eig (unsorted); only used for the condition number of L^T M L,
// whose size is the number of retained DOFs of a group (typically <= 6).
// Returns 0 on success.
static int
ladrunoSymEigJacobi(const Matrix &A, Vector &eig)
{
    int n = A.noRows();
    if (n <= 0 || A.noCols() != n) return -1;
    eig = Vector(n);
    if (n == 1) { eig(0) = A(0, 0); return 0; }

    std::vector<double> a((size_t)n * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[(size_t)i * n + j] = A(i, j);

    const int maxSweeps = 100;
    for (int sweep = 0; sweep < maxSweeps; sweep++) {
        double off = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                off += a[(size_t)i * n + j] * a[(size_t)i * n + j];
        if (off <= 1.0e-30) break;

        for (int p = 0; p < n; p++) {
            for (int q = p + 1; q < n; q++) {
                double apq = a[(size_t)p * n + q];
                if (apq == 0.0) continue;
                double app = a[(size_t)p * n + p];
                double aqq = a[(size_t)q * n + q];
                // t is the smaller root of t^2 + 2 theta t - 1 = 0; the plane
                // rotation G (G_pq = s, G_qp = -s, diagonal c) zeroes the
                // (p,q) entry of G^T A G
                double theta = (aqq - app) / (2.0 * apq);
                double t = (theta >= 0.0 ? 1.0 : -1.0) /
                           (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                double c = 1.0 / std::sqrt(t * t + 1.0);
                double s = t * c;
                // A G (columns p, q)
                for (int k = 0; k < n; k++) {
                    double g = a[(size_t)k * n + p];
                    double h = a[(size_t)k * n + q];
                    a[(size_t)k * n + p] = c * g - s * h;
                    a[(size_t)k * n + q] = s * g + c * h;
                }
                // G^T (A G) (rows p, q)
                for (int k = 0; k < n; k++) {
                    double g = a[(size_t)p * n + k];
                    double h = a[(size_t)q * n + k];
                    a[(size_t)p * n + k] = c * g - s * h;
                    a[(size_t)q * n + k] = s * g + c * h;
                }
            }
        }
    }
    for (int i = 0; i < n; i++) eig(i) = a[(size_t)i * n + i];
    return 0;
}

LadrunoConstraintProjector::LadrunoConstraintProjector()
  : groups(), massReady(false), projectICs(false), icTol(1.0e-8)
{
}

LadrunoConstraintProjector::~LadrunoConstraintProjector()
{
}

int
LadrunoConstraintProjector::enforceIC(Vector &u, Vector &v) const
{
    if (projectICs) {
        this->snapICs(u, v);
        return 0;
    }
    return this->checkIC(u, v, icTol);
}

void
LadrunoConstraintProjector::addGroup(const ID &allEqn, int nRet, const Matrix &L,
                                     const Vector &delta, const ID &fixedEqn)
{
    Group g;
    g.allEqn = allEqn;
    g.nRet   = nRet;
    g.L      = L;
    g.delta  = delta;
    g.fixedEqn = fixedEqn;
    g.m      = Vector(allEqn.Size());
    g.LtML   = Matrix(nRet, nRet);
    groups.push_back(g);
    massReady = false;
}

int
LadrunoConstraintProjector::buildMass(LinearSOE *theSOE)
{
    DiagonalSOE *theDiag = dynamic_cast<DiagonalSOE *>(theSOE);
    if (theDiag == 0) {
        opserr << "LadrunoConstraintProjector::buildMass() - constraints LadrunoProjection "
                  "requires system Diagonal (lumped mass); the current system is not a "
                  "DiagonalSOE.\n";
        return -1;
    }

    // The assembled diagonal is the lumped mass the integrator inverts; it is
    // read before the solver replaces A[i] by 1/A[i].
    const double *A = theDiag->getDiagonalA();
    int neq = theDiag->getNumEqn();
    if (A == 0) {
        opserr << "LadrunoConstraintProjector::buildMass() - the system has not been sized.\n";
        return -1;
    }

    for (int gi = 0; gi < (int)groups.size(); gi++) {
        Group &g = groups[gi];
        int nrows = g.allEqn.Size();
        for (int r = 0; r < nrows; r++) {
            int e = g.allEqn(r);
            if (e < 0 || e >= neq) {
                opserr << "LadrunoConstraintProjector::buildMass() - equation " << e
                       << " out of range [0," << neq << ") in group " << gi << ".\n";
                return -1;
            }
            g.m(r) = A[e];
            if (g.m(r) <= 0.0) {
                opserr << "LadrunoConstraintProjector::buildMass() - massless DOF "
                          "(equation " << e << ", lumped mass " << g.m(r) << ") in a "
                          "constraint group. The projection handler keeps every tied DOF "
                          "in the equation set, so each needs a nonzero lumped mass. Add "
                          "translational mass, and for rotational ties (rigidLink -beam, "
                          "rigidDiaphragm) a small rotational mass (about 0.01-0.1% of the "
                          "translational mass of the node), or use constraints "
                          "Transformation or Penalty.\n";
                return -1;
            }
        }

        // L^T diag(m) L
        int nRet = g.nRet;
        g.LtML.Zero();
        for (int p = 0; p < nRet; p++) {
            for (int q = 0; q < nRet; q++) {
                double s = 0.0;
                for (int r = 0; r < nrows; r++)
                    s += g.L(r, p) * g.m(r) * g.L(r, q);
                g.LtML(p, q) = s;
            }
        }

        // Condition number lambda_max / lambda_min of L^T M L. A nearly
        // singular matrix (a retained direction with almost no mass, nearly
        // redundant constraints) is refused; an ill-conditioned one is noted.
        if (nRet > 0) {
            Vector eig;
            if (ladrunoSymEigJacobi(g.LtML, eig) < 0) {
                opserr << "LadrunoConstraintProjector::buildMass() - eigenvalue estimate "
                          "failed for L^T M L in group " << gi << ".\n";
                return -1;
            }
            double lmax = eig(0), lmin = eig(0);
            for (int k = 1; k < eig.Size(); k++) {
                if (eig(k) > lmax) lmax = eig(k);
                if (eig(k) < lmin) lmin = eig(k);
            }
            double cond = (lmin > 0.0) ? (lmax / lmin)
                                       : std::numeric_limits<double>::infinity();
            if (lmin <= 0.0 || cond > LADRUNO_LTML_COND_REFUSE) {
                opserr << "LadrunoConstraintProjector::buildMass() - near-singular L^T M L in "
                          "group " << gi << " (condition number " << cond << ", minimum "
                          "eigenvalue " << lmin << "). A retained direction has almost no "
                          "mass, or the constraints are redundant. Correct the constraint "
                          "set, add mass, or use constraints Transformation or Penalty.\n";
                return -1;
            }
            if (cond > LADRUNO_LTML_COND_WARN)
                opserr << "LadrunoConstraintProjector::buildMass() - NOTE L^T M L in group " << gi
                       << " is ill-conditioned (condition number " << cond << "); the projected "
                          "acceleration may lose precision. Check for nearly redundant "
                          "constraints or widely different tied masses.\n";
        }
    }

    massReady = true;
    return 0;
}

void
LadrunoConstraintProjector::project(Vector &a)
{
    if (!massReady) {
        opserr << "LadrunoConstraintProjector::project() - mass not built; buildMass() "
                  "must run after the mass-only formTangent and before project().\n";
        return;
    }

    for (int gi = 0; gi < (int)groups.size(); gi++) {
        Group &g = groups[gi];

        // constrained DOFs whose retained DOFs are all fixed: a_c = 0
        for (int i = 0; i < g.fixedEqn.Size(); i++)
            a(g.fixedEqn(i)) = 0.0;

        int nrows = g.allEqn.Size();
        int nRet  = g.nRet;
        if (nRet == 0 || nrows == 0)
            continue;

        // rhs = L^T M a_raw
        Vector rhs(nRet);
        rhs.Zero();
        for (int p = 0; p < nRet; p++) {
            double s = 0.0;
            for (int r = 0; r < nrows; r++)
                s += g.L(r, p) * g.m(r) * a(g.allEqn(r));
            rhs(p) = s;
        }

        // a_r = (L^T M L)^-1 rhs, solved on a copy
        Vector ar(nRet);
        Matrix tmp(g.LtML);
        if (tmp.Solve(rhs, ar) < 0) {
            opserr << "LadrunoConstraintProjector::project() - solve failed in group "
                   << gi << "\n";
            continue;
        }

        // a_group = L a_r, written to the retained and the constrained rows
        for (int r = 0; r < nrows; r++) {
            double aproj = 0.0;
            for (int p = 0; p < nRet; p++)
                aproj += g.L(r, p) * ar(p);
            a(g.allEqn(r)) = aproj;
        }
    }
}

int
LadrunoConstraintProjector::checkIC(const Vector &u, const Vector &v, double icTol) const
{
    int nViol = 0;

    for (int gi = 0; gi < (int)groups.size(); gi++) {
        const Group &g = groups[gi];
        int nrows = g.allEqn.Size();
        int nRet  = g.nRet;

        // tolerance scaled by the displacement and velocity magnitude of the group
        double umax = 0.0, vmax = 0.0;
        for (int r = 0; r < nrows; r++) {
            double au = std::fabs(u(g.allEqn(r)));  if (au > umax) umax = au;
            double av = std::fabs(v(g.allEqn(r)));  if (av > vmax) vmax = av;
        }
        double tolU = icTol * (1.0 + umax);
        double tolV = icTol * (1.0 + vmax);

        // constrained rows: u_c - C u_r = delta, v_c - C v_r = 0
        for (int j = 0; j < nrows - nRet; j++) {
            int row = nRet + j;
            int ec  = g.allEqn(row);
            double Cur = 0.0, Cvr = 0.0;
            for (int p = 0; p < nRet; p++) {
                Cur += g.L(row, p) * u(g.allEqn(p));
                Cvr += g.L(row, p) * v(g.allEqn(p));
            }
            double ru = u(ec) - Cur - g.delta(j);
            double rv = v(ec) - Cvr;
            if (std::fabs(ru) > tolU || std::fabs(rv) > tolV) {
                opserr << "LadrunoConstraintProjector::checkIC() - the initial state "
                          "violates a projected constraint at equation " << ec
                       << " (displacement residual " << ru << ", velocity residual " << rv
                       << "). The constrained DOF and its retained DOFs must start on "
                          "the constraint manifold; correct the initial conditions or use "
                          "-projectICs.\n";
                nViol++;
            }
        }

        // constrained DOFs whose retained DOFs are all fixed start at zero
        for (int i = 0; i < g.fixedEqn.Size(); i++) {
            int ec = g.fixedEqn(i);
            if (std::fabs(u(ec)) > tolU || std::fabs(v(ec)) > tolV) {
                opserr << "LadrunoConstraintProjector::checkIC() - constrained DOF "
                          "(equation " << ec << ") tied to fixed DOFs must start at zero "
                          "(displacement " << u(ec) << ", velocity " << v(ec) << ").\n";
                nViol++;
            }
        }
    }

    return nViol;
}

void
LadrunoConstraintProjector::snapICs(Vector &u, Vector &v) const
{
    for (int gi = 0; gi < (int)groups.size(); gi++) {
        const Group &g = groups[gi];
        int nrows = g.allEqn.Size();
        int nRet  = g.nRet;

        for (int j = 0; j < nrows - nRet; j++) {
            int row = nRet + j;
            int ec  = g.allEqn(row);
            double Cur = 0.0, Cvr = 0.0;
            for (int p = 0; p < nRet; p++) {
                Cur += g.L(row, p) * u(g.allEqn(p));
                Cvr += g.L(row, p) * v(g.allEqn(p));
            }
            u(ec) = Cur + g.delta(j);
            v(ec) = Cvr;
        }
        for (int i = 0; i < g.fixedEqn.Size(); i++) {
            u(g.fixedEqn(i)) = 0.0;
            v(g.fixedEqn(i)) = 0.0;
        }
    }
}
