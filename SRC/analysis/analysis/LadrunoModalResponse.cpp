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
// Description: modal-superposition response commands (modalResponseHistory,
// frequencyResponse, steadyStateDynamics, randomResponse).
// See LadrunoModalResponse.h for the formulation and references.

#include <LadrunoModalResponse.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <DomainModalProperties.h>
#include <TimeSeries.h>
#include <Node.h>
#include <NodeIter.h>
#include <LoadPattern.h>
#include <NodalLoad.h>
#include <NodalLoadIter.h>
#include <ElementalLoadIter.h>
#include <SP_ConstraintIter.h>
#include <classTags.h>
#include <Vector.h>
#include <Matrix.h>
#include <elementAPI.h>
#include <cmath>
#include <cstring>
#include <complex>
#include <fstream>
#include <iomanip>
#include <limits>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// shared modal-damping coefficient d_a = c_a/m_a (see header).
double
LadrunoModalDamping::coeff(int mode0, double w) const
{
    switch (kind) {
    case RAYLEIGH:
        return a0 + a1 * w * w;              // == 2 xi w for w>0; a0 for a rigid mode
    case MODAL_LIST: {
        double x = 0.0;
        if (mode0 >= 0 && mode0 < static_cast<int>(xis.size()))
            x = xis[mode0];
        return 2.0 * x * w;                  // rigid mode -> 0
    }
    case UNIFORM:
    default:
        return 2.0 * xi * w;                 // rigid mode -> 0
    }
}

// ---------------------------------------------------------------------------
// Exact PWL SDOF recurrence blocks in (w, d) form.
//   q'' + d q' + w^2 q = f(t),  f linear on the step.
//   [q1;qd1] = A [q0;qd0] + B [f0;f1],   exact for piecewise-linear f.
// Unified over all damping branches by Delta = w^2 - (d/2)^2, and covering the
// rigid-body Rayleigh case (w=0, d>0).  Agrees to ~1e-14 with the
// first-order-hold matrix exponential of Van Loan (1978) over 40 (w,d) cases
// and the rigid-Rayleigh steady state.
// ---------------------------------------------------------------------------
void
LadrunoModalResponse::recurrenceBlocks(double w, double d, double dt,
                                       double A[2][2], double B[2][2])
{
    const double CRIT_TOL = 1.0e-8; // relative closeness to the critical root
    const double half_d = 0.5 * d;

    // --- homogeneous propagator A -----------------------------------------
    double p11, p12, p21, p22;
    if (w == 0.0) {
        if (d == 0.0) {
            // double integrator
            p11 = 1.0; p12 = dt; p21 = 0.0; p22 = 1.0;
        } else {
            // roots {0, -d}: rigid mode with mass-proportional damping
            const double e = std::exp(-d * dt);
            p11 = 1.0; p12 = (1.0 - e) / d; p21 = 0.0; p22 = e;
        }
    } else {
        const double e = std::exp(-half_d * dt);
        const double Delta = w * w - half_d * half_d; // w^2 - (d/2)^2
        double c, sd;                                 // sd = sin(wd dt)/wd
        if (std::abs(Delta) <= CRIT_TOL * (w * w)) {
            c = 1.0; sd = dt;                          // critically damped limit
        } else if (Delta > 0.0) {
            const double wd = std::sqrt(Delta);        // underdamped
            c = std::cos(wd * dt); sd = std::sin(wd * dt) / wd;
        } else {
            const double wh = std::sqrt(-Delta);       // overdamped
            c = std::cosh(wh * dt); sd = std::sinh(wh * dt) / wh;
        }
        p11 = e * (c + half_d * sd);
        p12 = e * sd;
        p21 = e * (-w * w * sd);
        p22 = e * (c - half_d * sd);
    }
    A[0][0] = p11; A[0][1] = p12;
    A[1][0] = p21; A[1][1] = p22;

    // --- load block B -----------------------------------------------------
    // particular ramp solution q_p, cancelled to REST at t=0 by the
    // homogeneous part, then propagated by A to t=dt.
    if (w == 0.0 && d == 0.0) {
        B[0][0] = dt * dt / 3.0; B[0][1] = dt * dt / 6.0;
        B[1][0] = dt / 2.0;      B[1][1] = dt / 2.0;
        return;
    }
    for (int j = 0; j < 2; ++j) {
        const double f0 = (j == 0) ? 1.0 : 0.0;
        const double f1 = (j == 0) ? 0.0 : 1.0;
        const double r  = (f1 - f0) / dt;
        double qp0, qpd0, qpDt, qpdDt;
        if (w == 0.0) {
            // w=0, d>0:  v'=q', v' + d v = f0 + r t
            //   v_p(t) = (f0 + r t)/d - r/d^2
            //   q_p(t) = (f0 t + r t^2/2)/d - r t/d^2      (q_p(0)=0)
            const double invd = 1.0 / d, invd2 = invd * invd;
            qp0   = 0.0;
            qpd0  = f0 * invd - r * invd2;
            qpDt  = (f0 * dt + r * dt * dt * 0.5) * invd - r * dt * invd2;
            qpdDt = (f0 + r * dt) * invd - r * invd2;
        } else {
            // w>0, any d:  q_p(t) = (f0 + r t)/w^2 - d r / w^4
            const double w2 = w * w, invw2 = 1.0 / w2, invw4 = invw2 * invw2;
            qp0   = f0 * invw2 - d * r * invw4;
            qpd0  = r * invw2;
            qpDt  = (f0 + r * dt) * invw2 - d * r * invw4;
            qpdDt = r * invw2;
        }
        const double qh0  = -qp0;
        const double qhd0 = -qpd0;
        const double qh_dt  = p11 * qh0 + p12 * qhd0;
        const double qhd_dt = p21 * qh0 + p22 * qhd0;
        B[0][j] = qpDt  + qh_dt;
        B[1][j] = qpdDt + qhd_dt;
    }
}

// ---------------------------------------------------------------------------
// Shared -load helper: raw-eigenvector modal load sums
//   psiTP_raw[a] = sum_nodes sum_i evec(i, modes[a]) * P_node(i)
// from a load pattern holding ONLY plain NodalLoads (reference values; the
// pattern's own TimeSeries is deliberately not read: the caller supplies the
// time/frequency dependence).  The Vscale factor is applied by the callers
// together with the recovery psi so the pair matches the basis in which
// DomainModalProperties computed generalizedMasses(); the result is then
// independent of the eigenvector normalization.  Anything that would drop
// part of the excitation is refused.  `modes` = active 0-based indices.
// ---------------------------------------------------------------------------
static int
LadrunoAssembleModalLoadShape(Domain* domain, int patternTag, int ndf,
                              const std::vector<int>& modes, const char* cmd,
                              std::vector<double>& psiTP)
{
    const int nm = static_cast<int>(modes.size());
    LoadPattern* pat = domain->getLoadPattern(patternTag);
    if (pat == 0) {
        opserr << cmd << " - -load: no loadPattern with tag "
               << patternTag << ".\n";
        return -1;
    }
    // only plain nodal loads define the spatial shape P: refuse patterns
    // carrying element loads or sp constraints rather than silently dropping
    // part of the excitation.
    if (pat->getElementalLoads()() != 0 || pat->getSPs()() != 0) {
        opserr << cmd << " - -load: pattern " << patternTag
               << " carries elemental loads and/or sp constraints; only "
                  "plain nodal loads define the load shape P (build a "
                  "separate pattern holding just the nodal loads).\n";
        return -1;
    }
    psiTP.assign(nm, 0.0);
    bool anyLoad = false;
    NodalLoadIter& nlIter = pat->getNodalLoads();
    NodalLoad* nl;
    while ((nl = nlIter()) != 0) {
        // getData() returns the reference load vector only for a PLAIN
        // NodalLoad; derived types (e.g. thermal actions) overload it with
        // other semantics; refuse them.
        if (nl->getClassTag() != LOAD_TAG_NodalLoad) {
            opserr << cmd << " - -load: pattern " << patternTag
                   << " holds a non-plain nodal load (classTag "
                   << nl->getClassTag() << "); only plain 'load' nodal "
                      "loads are supported.\n";
            return -1;
        }
        Node* lnode = domain->getNode(nl->getNodeTag());
        if (lnode == 0 || lnode->getNumEigenvectors() < 1) {
            // a load on a node with NO eigenvector (orphan / outside the
            // eigen analysis) would silently vanish from the modal
            // projection: refuse it. (A fully-fixed node HAS a zero
            // eigenvector and legitimately contributes nothing.)
            opserr << cmd << " - -load: nodal load on node "
                   << nl->getNodeTag() << " which has no eigenvector "
                      "(unknown node, or not part of the last eigen "
                      "analysis).\n";
            return -1;
        }
        int ltype = 0;
        const Vector& Pn = nl->getData(ltype);
        const Matrix& lev = lnode->getEigenvectors();
        const int lndf = lev.noRows();
        const int nd = std::min(std::min(lndf, Pn.Size()), ndf);
        anyLoad = true;
        for (int a = 0; a < nm; ++a) {
            double s = 0.0;
            for (int i = 0; i < nd; ++i) {
                // pressure-DOF skip, same heuristic as the recovery
                if (ndf == 6 && lndf == 4 && i == 3)
                    continue;
                s += lev(i, modes[a]) * Pn(i);
            }
            psiTP[a] += s;
        }
    }
    if (!anyLoad) {
        opserr << cmd << " - -load: pattern " << patternTag
               << " has no nodal loads (the load shape P is empty).\n";
        return -1;
    }
    return 0;
}

// ---------------------------------------------------------------------------
LadrunoModalResponse::LadrunoModalResponse(AnalysisModel* theModel,
                                           TimeSeries*    baseAccel,
                                           int            direction,
                                           double         dt,
                                           int            nsteps,
                                           double         t0)
    : m_model(theModel)
    , m_baseAccel(baseAccel)
    , m_direction(direction)
    , m_dt(dt)
    , m_nsteps(nsteps)
    , m_t0(t0)
    , m_damp_kind(DAMP_UNIFORM)
    , m_xi(0.0)
    , m_a0(0.0), m_a1(0.0)
{
}

LadrunoModalResponse::~LadrunoModalResponse() {}

void LadrunoModalResponse::setDampingUniform(double xi)
{ m_damp_kind = DAMP_UNIFORM; m_xi = xi; }

void LadrunoModalResponse::setDampingRayleigh(double a0, double a1)
{ m_damp_kind = DAMP_RAYLEIGH; m_a0 = a0; m_a1 = a1; }

void LadrunoModalResponse::setDampingModalList(const std::vector<double>& xis)
{ m_damp_kind = DAMP_MODAL_LIST; m_xi_list = xis; }

void LadrunoModalResponse::setModeSubset(const std::vector<int>& modes_1based)
{
    m_modes.clear();
    for (int m : modes_1based) m_modes.push_back(m - 1); // to 0-based
}

double LadrunoModalResponse::dampingCoeff(int mode0, double w) const
{
    switch (m_damp_kind) {
    case DAMP_RAYLEIGH:
        // c/m = a0 + a1 w^2   (== 2 xi w for w>0; == a0 for a rigid mode)
        return m_a0 + m_a1 * w * w;
    case DAMP_MODAL_LIST: {
        double xi = 0.0;
        if (mode0 >= 0 && mode0 < static_cast<int>(m_xi_list.size()))
            xi = m_xi_list[mode0];
        return 2.0 * xi * w; // rigid mode -> 0
    }
    case DAMP_UNIFORM:
    default:
        return 2.0 * m_xi * w; // rigid mode -> 0
    }
}

int LadrunoModalResponse::check(DomainModalProperties& mp)
{
    Domain* domain = m_model->getDomainPtr();
    const Vector& ev = domain->getEigenvalues();
    if (ev.Size() < 1) {
        opserr << "LadrunoModalResponse - No eigenvalues; call 'eigen' first.\n";
        return -1;
    }
    if (ev.Size() != mp.eigenvalues().Size()) {
        opserr << "LadrunoModalResponse - modalProperties is stale vs the last "
                  "eigen run; call 'modalProperties' right after 'eigen'.\n";
        return -1;
    }
    double tol = std::max(1.0e-15, 1.0e-12 * ev.Norm());
    for (int i = 0; i < ev.Size(); ++i)
        if (std::abs(ev(i) - mp.eigenvalues()(i)) > tol) {
            opserr << "LadrunoModalResponse - eigenvalue mismatch between the "
                      "domain and DomainModalProperties (re-run modalProperties "
                      "after eigen).\n";
            return -1;
        }
    return 0;
}

int LadrunoModalResponse::analyze()
{
    if (m_model == 0 || m_model->getDomainPtr() == 0) {
        opserr << "LadrunoModalResponse - no AnalysisModel/Domain.\n";
        return -1;
    }
    Domain* domain = m_model->getDomainPtr();

    if (m_baseAccel == 0) {
        opserr << "LadrunoModalResponse - an excitation time series is "
                  "required (-baseAccel <tsTag>, or -series <tsTag> with "
                  "-load).\n";
        return -1;
    }
    if (m_dt <= 0.0 || m_nsteps < 1) {
        opserr << "LadrunoModalResponse - need dt>0 and nsteps>=1.\n";
        return -1;
    }

    DomainModalProperties mp;
    if (domain->getModalProperties(mp) < 0) {
        opserr << "LadrunoModalResponse - eigen and modalProperties have not "
                  "been called.\n";
        return -1;
    }
    if (check(mp) < 0)
        return -1;

    const int num_eigen = domain->getEigenvalues().Size();
    const int ndf = mp.totalMass().Size();
    int exdof = 0;
    if (!m_useLoad) {
        if (m_direction < 1 || m_direction > ndf) {
            opserr << "LadrunoModalResponse - direction (" << m_direction
                   << ") must be in 1.." << ndf << ".\n";
            return -1;
        }
        exdof = m_direction - 1;
    }

    // active mode list (default: all)
    std::vector<int> modes = m_modes;
    if (modes.empty()) {
        modes.resize(num_eigen);
        for (int i = 0; i < num_eigen; ++i) modes[i] = i;
    }
    for (int m : modes)
        if (m < 0 || m >= num_eigen) {
            opserr << "LadrunoModalResponse - requested mode " << (m + 1)
                   << " is out of range 1.." << num_eigen << ".\n";
            return -1;
        }
    // reject duplicate modes: the recovery sums per-mode contributions, so a
    // repeated index would silently double-count that mode's response.
    {
        std::vector<int> seen = modes;
        std::sort(seen.begin(), seen.end());
        if (std::adjacent_find(seen.begin(), seen.end()) != seen.end()) {
            opserr << "LadrunoModalResponse - -modes contains a duplicate index "
                      "(each mode may appear at most once).\n";
            return -1;
        }
    }
    if (m_damp_kind == DAMP_MODAL_LIST) {
        int maxmode = 0;
        for (int m : modes) maxmode = std::max(maxmode, m);
        if (static_cast<int>(m_xi_list.size()) <= maxmode) {
            opserr << "LadrunoModalResponse - -modalDamp list has "
                   << static_cast<int>(m_xi_list.size()) << " ratios but mode "
                   << (maxmode + 1) << " is used.\n";
            return -1;
        }
    }
    const int nm = static_cast<int>(modes.size());

    // negative-eigenvalue guard: a small negative lambda is numerical noise on a
    // rigid-body mode (treat as w=0); a genuinely negative one means the eigen
    // run is ill-posed and the modal transient would be meaningless: refuse
    // rather than silently reinterpret it as rigid.
    double maxAbsLambda = 0.0;
    for (int m : modes)
        maxAbsLambda = std::max(maxAbsLambda, std::abs(mp.eigenvalues()(m)));
    const double negTol = std::max(1.0e-12, 1.0e-8 * maxAbsLambda);

    // -load channel: raw-eigenvector modal load sums via the shared helper
    std::vector<double> psiTP;
    if (m_useLoad) {
        if (LadrunoAssembleModalLoadShape(domain, m_patternTag, ndf, modes,
                                          "modalResponseHistory", psiTP) < 0)
            return -1;
    }

    // per-mode precompute.  fcoef is the modal-load coefficient:
    //   base accel:  f_a(t) = fcoef_a * ug''(t),  fcoef_a = -Gamma_a
    //   nodal load:  f_a(t) = fcoef_a * s(t),     fcoef_a = (psi_a^T P)/m~_a
    std::vector<double> w(nm), fcoef(nm), Vscale(nm), dcoef(nm);
    std::vector<double> Amat(nm * 4), Bmat(nm * 4);
    std::vector<double> q(nm, 0.0), qd(nm, 0.0), qdd(nm, 0.0);
    for (int a = 0; a < nm; ++a) {
        const int mode0 = modes[a];
        const double lambda = mp.eigenvalues()(mode0);
        if (lambda < -negTol) {
            opserr << "LadrunoModalResponse - mode " << (mode0 + 1)
                   << " has a negative eigenvalue (" << lambda << "); the eigen "
                      "result is ill-posed for a modal transient.\n";
            return -1;
        }
        const double wa = (lambda > 0.0) ? std::sqrt(lambda) : 0.0; // clamp rigid
        w[a] = wa;
        Vscale[a] = mp.eigenVectorScaleFactors()(mode0);
        if (m_useLoad) {
            const double ma = mp.generalizedMasses()(mode0);
            if (!(ma > 0.0)) {
                opserr << "LadrunoModalResponse - mode " << (mode0 + 1)
                       << " has a non-positive generalized mass (" << ma
                       << "); the modal load f_a = psi^T P / m_a is "
                          "undefined.\n";
                return -1;
            }
            fcoef[a] = (psiTP[a] * Vscale[a]) / ma;
        } else {
            fcoef[a] = -mp.modalParticipationFactors()(mode0, exdof);
        }
        dcoef[a] = dampingCoeff(mode0, wa);
        double A[2][2], B[2][2];
        recurrenceBlocks(wa, dcoef[a], m_dt, A, B);
        Amat[a * 4 + 0] = A[0][0]; Amat[a * 4 + 1] = A[0][1];
        Amat[a * 4 + 2] = A[1][0]; Amat[a * 4 + 3] = A[1][1];
        Bmat[a * 4 + 0] = B[0][0]; Bmat[a * 4 + 1] = B[0][1];
        Bmat[a * 4 + 2] = B[1][0]; Bmat[a * 4 + 3] = B[1][1];
    }

    // --- station loop: reconstruct u,udot,uddot and COMMIT one step each ---
    // station 0 is the initial rest state at t0; then nsteps advances.
    for (int n = 0; n <= m_nsteps; ++n) {
        const double t = m_t0 + n * m_dt;
        const double ug = m_baseAccel->getFactor(t);

        // modal accelerations q'' from the ODE (for the accel recovery)
        for (int a = 0; a < nm; ++a) {
            const double f_a = fcoef[a] * ug;
            qdd[a] = f_a - dcoef[a] * qd[a] - w[a] * w[a] * q[a];
        }

        // physical reconstruction  u = sum_a (phi_a Vscale) q_a  (+ vel/accel)
        Node* node;
        NodeIter& theNodes = domain->getNodes();
        while ((node = theNodes()) != 0) {
            const Matrix& evec = node->getEigenvectors();
            const int node_ndf = evec.noRows();
            if (node_ndf < 1) continue;

            Vector ud(node_ndf), uv(node_ndf), ua(node_ndf);
            ud.Zero(); uv.Zero(); ua.Zero();

            const int nd = std::min(node_ndf, ndf);
            for (int i = 0; i < nd; ++i) {
                // skip a pressure DOF (3D U-P node: node_ndf==4, i==3), same
                // heuristic as ResponseSpectrumAnalysis. NOTE: the 2D U-P case
                // (node_ndf==3, pressure at i==2) is ambiguous with a 2D beam
                // node (Ux,Uy,Rz) by DOF count alone, so, as in RSA, it is not
                // special-cased. Skipped/unmodeled DOFs are left at 0 by the
                // full-vector overwrite below (correct for a from-rest modal run).
                if (ndf == 6 && node_ndf == 4 && i == 3)
                    continue;
                double sd = 0.0, sv = 0.0, sa = 0.0;
                for (int a = 0; a < nm; ++a) {
                    const double psi = evec(i, modes[a]) * Vscale[a];
                    sd += psi * q[a];
                    sv += psi * qd[a];
                    sa += psi * qdd[a];
                }
                ud(i) = sd; uv(i) = sv; ua(i) = sa;
            }
            node->setTrialDisp(ud);
            node->setTrialVel(uv);
            node->setTrialAccel(ua);
        }

        // advance domain clock, propagate to elements, commit + record
        m_model->setCurrentDomainTime(t);
        if (m_model->updateDomain() < 0) {
            opserr << "LadrunoModalResponse - updateDomain failed at step "
                   << n << ".\n";
            return -1;
        }
        if (m_model->commitDomain() < 0) {
            opserr << "LadrunoModalResponse - commitDomain failed at step "
                   << n << ".\n";
            return -1;
        }

        // advance the modal SDOFs from station n to n+1
        if (n < m_nsteps) {
            const double t_next = m_t0 + (n + 1) * m_dt;
            const double ug_next = m_baseAccel->getFactor(t_next);
            for (int a = 0; a < nm; ++a) {
                const double f0 = fcoef[a] * ug;
                const double f1 = fcoef[a] * ug_next;
                const double q0 = q[a], qd0 = qd[a];
                q[a]  = Amat[a * 4 + 0] * q0 + Amat[a * 4 + 1] * qd0
                      + Bmat[a * 4 + 0] * f0 + Bmat[a * 4 + 1] * f1;
                qd[a] = Amat[a * 4 + 2] * q0 + Amat[a * 4 + 3] * qd0
                      + Bmat[a * 4 + 2] * f0 + Bmat[a * 4 + 3] * f1;
            }
        }
    }

    return 0;
}

// ---------------------------------------------------------------------------
// Command entry point (openseespy + Tcl share this via the interpreter
// wrappers, mirroring OPS_ResponseSpectrumAnalysis).
//
//   modalResponseHistory -dt $dt -nsteps $n
//       (-baseAccel $tsTag -dir $dir | -load $patternTag -series $tsTag)
//       (-damp $xi | -rayleigh $a0 $a1 | -modalDamp $xi1 $xi2 ...)
//       [-modes $m1 $m2 ...] [-t0 $t0]
//
//   Excitation is exactly ONE of:
//     -baseAccel $tsTag -dir k   uniform base accel ug''(t) (RELATIVE response)
//     -load $patternTag -series $tsTag   nodal forces P(t) = s(t)*P, P from
//         the pattern's plain NodalLoads (its own TimeSeries is IGNORED), s(t)
//         the -series time series sampled at the stations (ABSOLUTE response)
//
// Requires a prior `eigen N` and `modalProperties`.
// ---------------------------------------------------------------------------
int
OPS_LadrunoModalResponseHistory(void)
{
    AnalysisModel* theModel = *OPS_GetAnalysisModel();
    if (theModel == 0 || theModel->getDomainPtr() == 0) {
        opserr << "modalResponseHistory - no AnalysisModel/Domain available.\n";
        return -1;
    }

    double dt = 0.0, t0 = 0.0;
    int nsteps = 0, dir = 0, loadPatternTag = 0;
    bool haveLoad = false;
    TimeSeries* ts = 0;          // -baseAccel series
    TimeSeries* seriesTs = 0;    // -series scalar s(t) for -load
    int dampKind = -1;              // 0 uniform, 1 rayleigh, 2 modal-list
    double xi = 0.0, a0 = 0.0, a1 = 0.0;
    std::vector<double> xilist;
    std::vector<int> modes;
    int numData = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();

        if (strcmp(opt, "-dt") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&numData, &dt) < 0) {
                opserr << "modalResponseHistory - -dt requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-nsteps") == 0 || strcmp(opt, "-numSteps") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&numData, &nsteps) < 0) {
                opserr << "modalResponseHistory - -nsteps requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-t0") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&numData, &t0) < 0) {
                opserr << "modalResponseHistory - -t0 requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-dir") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&numData, &dir) < 0) {
                opserr << "modalResponseHistory - -dir requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-baseAccel") == 0) {
            int tag;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&numData, &tag) < 0) {
                opserr << "modalResponseHistory - -baseAccel requires a timeSeries tag.\n"; return -1;
            }
            ts = OPS_getTimeSeries(tag);
            if (ts == 0) {
                opserr << "modalResponseHistory - no timeSeries with tag " << tag << ".\n"; return -1;
            }
        }
        else if (strcmp(opt, "-load") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&numData, &loadPatternTag) < 0) {
                opserr << "modalResponseHistory - -load requires a loadPattern tag.\n"; return -1;
            }
            haveLoad = true;
        }
        else if (strcmp(opt, "-series") == 0) {
            int tag;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&numData, &tag) < 0) {
                opserr << "modalResponseHistory - -series requires a timeSeries tag.\n"; return -1;
            }
            seriesTs = OPS_getTimeSeries(tag);
            if (seriesTs == 0) {
                opserr << "modalResponseHistory - no timeSeries with tag " << tag << ".\n"; return -1;
            }
        }
        else if (strcmp(opt, "-damp") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&numData, &xi) < 0) {
                opserr << "modalResponseHistory - -damp requires a value.\n"; return -1;
            }
            // refuse negative damping: a negative ratio is unphysical (it
            // grows the modal response) and is almost always a typo. xi==0
            // stays legal (undamped modal history).
            if (xi < 0.0) {
                opserr << "modalResponseHistory - -damp xi must be >= 0 (got "
                       << xi << ").\n"; return -1;
            }
            dampKind = 0;
        }
        else if (strcmp(opt, "-rayleigh") == 0) {
            double v[2]; int nd = 2;
            if (OPS_GetNumRemainingInputArgs() < 2 || OPS_GetDoubleInput(&nd, v) < 0) {
                opserr << "modalResponseHistory - -rayleigh requires a0 a1.\n"; return -1;
            }
            a0 = v[0]; a1 = v[1]; dampKind = 1;
        }
        else if (strcmp(opt, "-modalDamp") == 0 || strcmp(opt, "-modalDamping") == 0) {
            xilist.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double item; int one = 1;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetDoubleInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1); // unconsume the non-numeric token
                    break;
                }
                xilist.push_back(item);
            }
            if (xilist.empty()) {
                opserr << "modalResponseHistory - -modalDamp requires >=1 ratio.\n"; return -1;
            }
            // every per-mode ratio must be >= 0.
            for (size_t j = 0; j < xilist.size(); ++j) {
                if (xilist[j] < 0.0) {
                    opserr << "modalResponseHistory - -modalDamp entry " << (int)j
                           << " = " << xilist[j] << " must be >= 0.\n"; return -1;
                }
            }
            dampKind = 2;
        }
        else if (strcmp(opt, "-modes") == 0) {
            modes.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int item; int one = 1;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                modes.push_back(item);
            }
            if (modes.empty()) {
                opserr << "modalResponseHistory - -modes requires >=1 mode index.\n"; return -1;
            }
        }
        else {
            opserr << "modalResponseHistory - unknown option '" << opt << "'.\n";
            return -1;
        }
    }

    // required arguments: exactly one excitation channel
    if (haveLoad && ts != 0) {
        opserr << "modalResponseHistory - -baseAccel and -load are mutually "
                  "exclusive (one excitation channel per run).\n"; return -1;
    }
    if (haveLoad) {
        if (seriesTs == 0) {
            opserr << "modalResponseHistory - -load requires -series <tsTag> "
                      "(the scalar s(t) multiplying the pattern's load shape).\n";
            return -1;
        }
        if (dir != 0)
            opserr << "WARNING modalResponseHistory - -dir is ignored with "
                      "-load (the pattern's nodal loads define the shape).\n";
    } else {
        if (seriesTs != 0) {
            opserr << "modalResponseHistory - -series is only meaningful with "
                      "-load (use -baseAccel <tsTag> for base motion).\n";
            return -1;
        }
        if (ts == 0) {
            opserr << "modalResponseHistory - an excitation channel is required: "
                      "-baseAccel <tsTag> -dir <1..ndf>, or -load <patternTag> "
                      "-series <tsTag>.\n"; return -1;
        }
        if (dir < 1) {
            opserr << "modalResponseHistory - -dir <1..ndf> is required.\n"; return -1;
        }
    }
    if (dt <= 0.0 || nsteps < 1) {
        opserr << "modalResponseHistory - -dt >0 and -nsteps >=1 are required.\n"; return -1;
    }
    if (dampKind < 0) {
        opserr << "modalResponseHistory - a damping channel is required "
                  "(-damp | -rayleigh | -modalDamp).\n"; return -1;
    }

    LadrunoModalResponse mrh(theModel, haveLoad ? seriesTs : ts,
                             haveLoad ? 0 : dir, dt, nsteps, t0);
    if (haveLoad)           mrh.setLoadPattern(loadPatternTag);
    if (dampKind == 0)      mrh.setDampingUniform(xi);
    else if (dampKind == 1) mrh.setDampingRayleigh(a0, a1);
    else                    mrh.setDampingModalList(xilist);
    if (!modes.empty())     mrh.setModeSubset(modes);

    return mrh.analyze();
}

// ===========================================================================
// modal frequency response / steady-state dynamics
// ===========================================================================
LadrunoFrequencyResponse::LadrunoFrequencyResponse(AnalysisModel* theModel)
    : m_model(theModel)
    , m_ex(EX_BASE_ACCEL)
    , m_direction(0)
    , m_patternTag(0)
    , m_amp(1.0)
    , m_fmin(0.0), m_fmax(0.0)
    , m_nf(0)
    , m_sweep(SWEEP_LIN)
    , m_nodeTag(0)
    , m_dof(0)
    , m_resp(RESP_DISP)
{
}

LadrunoFrequencyResponse::~LadrunoFrequencyResponse() {}

void LadrunoFrequencyResponse::setModeSubset(const std::vector<int>& modes_1based)
{
    m_modes.clear();
    for (int m : modes_1based) m_modes.push_back(m - 1); // to 0-based
}

// Build the Hz sample grid.  LIN: uniform in f.  LOG: geometric.  BIASED: a
// LIN base grid PLUS a fine cluster around each in-band modal frequency f_a =
// w_a/2pi so sharp (low-damping) resonant peaks are resolved.
int
LadrunoFrequencyResponse::buildGrid(std::vector<double>& freqs) const
{
    freqs.clear();
    if (m_nf < 1 || m_fmin < 0.0)
        return -1;

    if (m_nf == 1) { freqs.push_back(m_fmin); return 0; }

    // a multi-point sweep needs a real interval; fmax==fmin with nf>1 would emit
    // nf identical rows (LIN/LOG don't de-dup). Use nf==1 for a single frequency.
    if (m_fmax <= m_fmin)
        return -1;

    if (m_sweep == SWEEP_LOG) {
        if (m_fmin <= 0.0) return -1; // log needs a positive lower bound
        const double r = std::log(m_fmax / m_fmin) / (m_nf - 1);
        for (int i = 0; i < m_nf; ++i)
            freqs.push_back(m_fmin * std::exp(r * i));
        return 0;
    }

    // LIN base grid (also the base for BIASED)
    const double df = (m_fmax - m_fmin) / (m_nf - 1);
    for (int i = 0; i < m_nf; ++i)
        freqs.push_back(m_fmin + df * i);

    if (m_sweep == SWEEP_BIASED && m_model && m_model->getDomainPtr()) {
        const Vector& ev = m_model->getDomainPtr()->getEigenvalues();
        // active modes for the cluster (default: all extracted)
        std::vector<int> modes = m_modes;
        if (modes.empty())
            for (int i = 0; i < ev.Size(); ++i) modes.push_back(i);
        const int NCLUST = 15; // points per side within the window
        for (int m : modes) {
            if (m < 0 || m >= ev.Size()) continue;
            const double lam = ev(m);
            if (lam <= 0.0) continue; // rigid mode: no resonant peak to resolve
            const double fa = std::sqrt(lam) / (2.0 * M_PI);
            if (fa < m_fmin || fa > m_fmax) continue;
            const double halfw = 0.05 * fa; // +/-5% window around the peak
            for (int k = -NCLUST; k <= NCLUST; ++k) {
                const double f = fa + halfw * (double(k) / NCLUST);
                if (f >= m_fmin && f <= m_fmax) freqs.push_back(f);
            }
        }
        std::sort(freqs.begin(), freqs.end());
        freqs.erase(std::unique(freqs.begin(), freqs.end()), freqs.end());
    }
    return 0;
}

int
LadrunoFrequencyResponse::run(std::vector<std::vector<double> >& rows, bool ampOnly)
{
    rows.clear();
    if (m_model == 0 || m_model->getDomainPtr() == 0) {
        opserr << "frequencyResponse - no AnalysisModel/Domain.\n";
        return -1;
    }
    Domain* domain = m_model->getDomainPtr();

    const Vector& ev = domain->getEigenvalues();
    if (ev.Size() < 1) {
        opserr << "frequencyResponse - no eigenvalues; call 'eigen' first.\n";
        return -1;
    }
    DomainModalProperties mp;
    if (domain->getModalProperties(mp) < 0) {
        opserr << "frequencyResponse - eigen and modalProperties have not been "
                  "called.\n";
        return -1;
    }
    if (ev.Size() != mp.eigenvalues().Size()) {
        opserr << "frequencyResponse - modalProperties is stale vs the last eigen "
                  "run; call 'modalProperties' right after 'eigen'.\n";
        return -1;
    }
    // element-wise staleness check (same as modalResponseHistory): a re-run of
    // `eigen` with the SAME mode count but a changed model leaves a same-size
    // DomainModalProperties that the size test alone accepts; the FRF would then
    // be built from
    // stale w/Gamma/Vscale while the -biased grid reads the NEW domain eigenvalues
    // (inconsistent, silent). Refuse unless the snapshot matches the domain.
    {
        const double tol = std::max(1.0e-15, 1.0e-12 * ev.Norm());
        for (int i = 0; i < ev.Size(); ++i)
            if (std::abs(ev(i) - mp.eigenvalues()(i)) > tol) {
                opserr << "frequencyResponse - eigenvalue mismatch between the domain "
                          "and DomainModalProperties (re-run modalProperties after "
                          "eigen).\n";
                return -1;
            }
    }
    const int num_eigen = ev.Size();
    const int ndf = mp.totalMass().Size();
    int exdof = 0;
    if (m_ex == EX_BASE_ACCEL) {
        if (m_direction < 1 || m_direction > ndf) {
            opserr << "frequencyResponse - -dir (" << m_direction
                   << ") must be in 1.." << ndf << ".\n";
            return -1;
        }
        exdof = m_direction - 1;
    }

    // response node/dof
    Node* rnode = domain->getNode(m_nodeTag);
    if (rnode == 0) {
        opserr << "frequencyResponse - response node " << m_nodeTag
               << " not found.\n";
        return -1;
    }
    // probe with the non-exiting accessor first: Node::getEigenvectors() exit(0)s
    // the whole process when unset (a fully-constrained node, or one outside the
    // eigen analysis), which for a user-chosen response node would be a silent kill.
    if (rnode->getNumEigenvectors() < 1) {
        opserr << "frequencyResponse - response node " << m_nodeTag
               << " has no eigenvector (fully constrained, or not in the last "
                  "eigen analysis); pick a free node.\n";
        return -1;
    }
    const Matrix& revec = rnode->getEigenvectors();
    const int rndf = revec.noRows();
    if (m_dof < 1 || m_dof > rndf) {
        opserr << "frequencyResponse - -dof (" << m_dof << ") must be in 1.." << rndf
               << " for node " << m_nodeTag << ".\n";
        return -1;
    }
    const int rdof = m_dof - 1;

    // active mode list (default: all), 0-based, no duplicates
    std::vector<int> modes = m_modes;
    if (modes.empty()) {
        modes.resize(num_eigen);
        for (int i = 0; i < num_eigen; ++i) modes[i] = i;
    }
    for (int m : modes)
        if (m < 0 || m >= num_eigen) {
            opserr << "frequencyResponse - requested mode " << (m + 1)
                   << " is out of range 1.." << num_eigen << ".\n";
            return -1;
        }
    {
        std::vector<int> seen = modes;
        std::sort(seen.begin(), seen.end());
        if (std::adjacent_find(seen.begin(), seen.end()) != seen.end()) {
            opserr << "frequencyResponse - -modes contains a duplicate index.\n";
            return -1;
        }
    }
    if (m_damp.kind == LadrunoModalDamping::MODAL_LIST) {
        int maxmode = 0;
        for (int m : modes) maxmode = std::max(maxmode, m);
        if (static_cast<int>(m_damp.xis.size()) <= maxmode) {
            opserr << "frequencyResponse - -modalDamp list has "
                   << static_cast<int>(m_damp.xis.size()) << " ratios but mode "
                   << (maxmode + 1) << " is used.\n";
            return -1;
        }
    }
    const int nm = static_cast<int>(modes.size());

    // negative-eigenvalue guard (same policy as modalResponseHistory)
    double maxAbsLambda = 0.0;
    for (int m : modes)
        maxAbsLambda = std::max(maxAbsLambda, std::abs(mp.eigenvalues()(m)));
    const double negTol = std::max(1.0e-12, 1.0e-8 * maxAbsLambda);

    // -load channel: assemble the RAW-eigenvector modal load sums
    //   psiTP_raw[a] = sum_nodes sum_i evec(i, mode_a) * P_node(i)
    // via the shared helper (Vscale applied with the recovery psi below so the
    // pair matches the basis generalizedMasses() was computed in).
    std::vector<double> psiTP;
    if (m_ex == EX_LOAD) {
        if (LadrunoAssembleModalLoadShape(domain, m_patternTag, ndf, modes,
                                          "frequencyResponse", psiTP) < 0)
            return -1;
    }

    // per-mode precompute: w, d, and the scalar recovery weight
    //   base accel:  c_a = psi_a(node,dof) * (-Gamma_a)
    //   nodal load:  c_a = psi_a(node,dof) * (psi_a^T P) / m~_a
    // with psi_a = evec*Vscale and m~_a = generalizedMasses()(a) (same basis).
    std::vector<double> w(nm), d(nm), cw(nm);
    for (int a = 0; a < nm; ++a) {
        const int mode0 = modes[a];
        const double lambda = mp.eigenvalues()(mode0);
        if (lambda < -negTol) {
            opserr << "frequencyResponse - mode " << (mode0 + 1)
                   << " has a negative eigenvalue (" << lambda << ").\n";
            return -1;
        }
        const double wa = (lambda > 0.0) ? std::sqrt(lambda) : 0.0;
        w[a] = wa;
        d[a] = m_damp.coeff(mode0, wa);
        const double Vscale = mp.eigenVectorScaleFactors()(mode0);
        const double psi    = revec(rdof, mode0) * Vscale;
        if (m_ex == EX_LOAD) {
            const double ma = mp.generalizedMasses()(mode0);
            if (!(ma > 0.0)) {
                opserr << "frequencyResponse - mode " << (mode0 + 1)
                       << " has a non-positive generalized mass (" << ma
                       << "); the modal load f_a = psi^T P / m_a is "
                          "undefined.\n";
                return -1;
            }
            cw[a] = psi * (psiTP[a] * Vscale) / ma;
        } else {
            const double Gamma = mp.modalParticipationFactors()(mode0, exdof);
            cw[a] = psi * (-Gamma);
        }
    }

    std::vector<double> freqs;
    if (buildGrid(freqs) < 0) {
        opserr << "frequencyResponse - bad frequency sweep "
                  "(need nf>=1, 0<=fmin<=fmax; log needs fmin>0).\n";
        return -1;
    }

    // sweep
    const std::complex<double> I(0.0, 1.0);
    for (double f : freqs) {
        const double Om = 2.0 * M_PI * f;
        std::complex<double> uc(0.0, 0.0);
        for (int a = 0; a < nm; ++a) {
            const std::complex<double> denom(w[a] * w[a] - Om * Om, Om * d[a]);
            uc += cw[a] / denom;               // psi*(-Gamma)*H_a(Om)
        }
        // response-quantity multiplier (relative), then excitation amplitude
        if (m_resp == RESP_VEL)        uc *= I * Om;
        else if (m_resp == RESP_ACCEL) uc *= -(Om * Om);
        uc *= m_amp;

        std::vector<double> row;
        row.push_back(f);
        if (ampOnly) {
            row.push_back(std::abs(uc));
        } else {
            row.push_back(uc.real());
            row.push_back(uc.imag());
        }
        rows.push_back(row);
    }

    // optional file output
    if (!m_outfile.empty()) {
        std::ofstream os(m_outfile.c_str());
        if (!os) {
            opserr << "frequencyResponse - cannot open output file '"
                   << m_outfile.c_str() << "'.\n";
            return -1;
        }
        os << std::scientific << std::setprecision(15);
        for (const std::vector<double>& r : rows) {
            for (size_t j = 0; j < r.size(); ++j)
                os << (j ? " " : "") << r[j];
            os << "\n";
        }
    }
    return 0;
}

// ---------------------------------------------------------------------------
// Shared parser for frequencyResponse / steadyStateDynamics.
//
//   frequencyResponse -freq fmin fmax nf [-lin|-log|-biased]
//       (-baseAccel -dir $dir | -load $patternTag) [-amp $a]
//       (-damp $xi | -rayleigh $a0 $a1 | -modalDamp $xi1 ...)
//       -node $tag -dof $dof [-resp disp|vel|accel] [-modes ...] [-out $file]
//
//   steadyStateDynamics ...same... (reports |response| instead of Re/Im)
//
//   Excitation is exactly ONE of:
//     -baseAccel -dir k   uniform base accel amp*e^{+iOm t} (relative response)
//     -load $patternTag   nodal forces amp*P*e^{+iOm t}, P = the pattern's
//                         plain NodalLoad reference values (the pattern's own
//                         TimeSeries is IGNORED; absolute response)
//
// Requires a prior `eigen N` and `modalProperties`.  Returns the table to the
// interpreter via OPS_SetDoubleListsOutput (and writes -out if given).
// ---------------------------------------------------------------------------
static int
OPS_FreqResponseImpl(bool ampOnly, const char* cmd)
{
    AnalysisModel* theModel = *OPS_GetAnalysisModel();
    if (theModel == 0 || theModel->getDomainPtr() == 0) {
        opserr << cmd << " - no AnalysisModel/Domain available.\n";
        return -1;
    }

    double fmin = 0.0, fmax = 0.0, amp = 1.0;
    int nf = 0, dir = 0, nodeTag = 0, dof = 0, loadPatternTag = 0;
    LadrunoFrequencyResponse::SweepKind sweep = LadrunoFrequencyResponse::SWEEP_LIN;
    LadrunoFrequencyResponse::RespKind  resp  = LadrunoFrequencyResponse::RESP_DISP;
    bool haveFreq = false, haveBase = false, haveLoad = false;
    bool haveDamp = false, haveNode = false;
    LadrunoModalDamping damp;
    std::vector<int> modes;
    std::string outfile;
    int one = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();

        if (strcmp(opt, "-freq") == 0) {
            double fv[2]; int nd = 2;
            if (OPS_GetNumRemainingInputArgs() < 3 || OPS_GetDoubleInput(&nd, fv) < 0
                    || OPS_GetIntInput(&one, &nf) < 0) {
                opserr << cmd << " - -freq requires fmin fmax nf.\n"; return -1;
            }
            fmin = fv[0]; fmax = fv[1]; haveFreq = true;
        }
        else if (strcmp(opt, "-lin") == 0)    sweep = LadrunoFrequencyResponse::SWEEP_LIN;
        else if (strcmp(opt, "-log") == 0)    sweep = LadrunoFrequencyResponse::SWEEP_LOG;
        else if (strcmp(opt, "-biased") == 0) sweep = LadrunoFrequencyResponse::SWEEP_BIASED;
        else if (strcmp(opt, "-baseAccel") == 0) haveBase = true;
        else if (strcmp(opt, "-load") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &loadPatternTag) < 0) {
                opserr << cmd << " - -load requires a loadPattern tag.\n"; return -1;
            }
            haveLoad = true;
        }
        else if (strcmp(opt, "-dir") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &dir) < 0) {
                opserr << cmd << " - -dir requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-amp") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&one, &amp) < 0) {
                opserr << cmd << " - -amp requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-node") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &nodeTag) < 0) {
                opserr << cmd << " - -node requires a tag.\n"; return -1;
            }
            haveNode = true;
        }
        else if (strcmp(opt, "-dof") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &dof) < 0) {
                opserr << cmd << " - -dof requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-resp") == 0 || strcmp(opt, "-response") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << cmd << " - -resp requires disp|vel|accel.\n"; return -1;
            }
            const char* rname = OPS_GetString();
            if (strcmp(rname, "disp") == 0)       resp = LadrunoFrequencyResponse::RESP_DISP;
            else if (strcmp(rname, "vel") == 0)   resp = LadrunoFrequencyResponse::RESP_VEL;
            else if (strcmp(rname, "accel") == 0) resp = LadrunoFrequencyResponse::RESP_ACCEL;
            else {
                opserr << cmd << " - unknown -resp '" << rname
                       << "' (use disp|vel|accel).\n"; return -1;
            }
        }
        else if (strcmp(opt, "-damp") == 0) {
            double xi;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&one, &xi) < 0) {
                opserr << cmd << " - -damp requires a value.\n"; return -1;
            }
            // refuse negative damping: d_a = 2 xi w_a < 0 silently conjugates
            // the FRF (H -> H*), identical |H| with the phase sign flipped, so a
            // mistyped '-0.05' passes every magnitude check. xi==0 is legal (the
            // undamped resonance is a genuine singularity).
            if (xi < 0.0) {
                opserr << cmd << " - -damp xi must be >= 0 (got " << xi
                       << "); negative damping conjugates the FRF (same magnitude, "
                          "flipped phase). Use xi=0 for the undamped case.\n";
                return -1;
            }
            damp.kind = LadrunoModalDamping::UNIFORM; damp.xi = xi; haveDamp = true;
        }
        else if (strcmp(opt, "-rayleigh") == 0) {
            double v[2]; int nd = 2;
            if (OPS_GetNumRemainingInputArgs() < 2 || OPS_GetDoubleInput(&nd, v) < 0) {
                opserr << cmd << " - -rayleigh requires a0 a1.\n"; return -1;
            }
            damp.kind = LadrunoModalDamping::RAYLEIGH; damp.a0 = v[0]; damp.a1 = v[1];
            haveDamp = true;
        }
        else if (strcmp(opt, "-modalDamp") == 0 || strcmp(opt, "-modalDamping") == 0) {
            damp.xis.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double item;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetDoubleInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                damp.xis.push_back(item);
            }
            if (damp.xis.empty()) {
                opserr << cmd << " - -modalDamp requires >=1 ratio.\n"; return -1;
            }
            // every per-mode ratio must be >= 0 (same FRF conjugation by
            // negative damping as -damp).
            for (size_t j = 0; j < damp.xis.size(); ++j) {
                if (damp.xis[j] < 0.0) {
                    opserr << cmd << " - -modalDamp entry " << (int)j << " = "
                           << damp.xis[j] << " must be >= 0 (negative damping "
                              "conjugates the FRF).\n";
                    return -1;
                }
            }
            damp.kind = LadrunoModalDamping::MODAL_LIST; haveDamp = true;
        }
        else if (strcmp(opt, "-modes") == 0) {
            modes.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int item;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                modes.push_back(item);
            }
            if (modes.empty()) {
                opserr << cmd << " - -modes requires >=1 mode index.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-out") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << cmd << " - -out requires a filename.\n"; return -1;
            }
            outfile = OPS_GetString();
        }
        else {
            opserr << cmd << " - unknown option '" << opt << "'.\n"; return -1;
        }
    }

    if (!haveFreq) { opserr << cmd << " - -freq fmin fmax nf is required.\n"; return -1; }
    if (haveBase && haveLoad) {
        opserr << cmd << " - -baseAccel and -load are mutually exclusive "
                  "(one excitation channel per sweep).\n"; return -1;
    }
    if (haveLoad) {
        if (dir != 0)
            opserr << "WARNING " << cmd << " - -dir is ignored with -load "
                      "(the pattern's nodal loads define the shape).\n";
    }
    else if (!haveBase || dir < 1) {
        opserr << cmd << " - an excitation channel is required: "
                  "-baseAccel -dir <1..ndf>, or -load <patternTag>.\n"; return -1;
    }
    if (!haveNode || dof < 1) {
        opserr << cmd << " - -node <tag> -dof <1..ndf> are required.\n"; return -1;
    }
    if (!haveDamp) {
        opserr << cmd << " - a damping channel is required "
                  "(-damp | -rayleigh | -modalDamp).\n"; return -1;
    }

    // a -biased grid places k=0 exactly on each in-band
    // resonance f_a=w_a/2pi (f = fa + halfw*(k/NCLUST), k=0 => f=fa). With an
    // exactly-zero d_a the denominator w_a^2-Om^2+i Om d_a is real and hits 0
    // there => an inf/NaN FRF row. That is an honest singularity (not guarded),
    // but -biased produces it by construction, so warn once: a zero-damping
    // resonance-biased sweep is almost never intended.
    if (sweep == LadrunoFrequencyResponse::SWEEP_BIASED) {
        bool zeroDamp = false;
        if (damp.kind == LadrunoModalDamping::UNIFORM)       zeroDamp = (damp.xi == 0.0);
        else if (damp.kind == LadrunoModalDamping::RAYLEIGH) zeroDamp = (damp.a0 == 0.0 && damp.a1 == 0.0);
        else if (damp.kind == LadrunoModalDamping::MODAL_LIST) {
            for (size_t j = 0; j < damp.xis.size(); ++j)
                if (damp.xis[j] == 0.0) { zeroDamp = true; break; }
        }
        if (zeroDamp)
            opserr << "WARNING " << cmd << " - -biased with zero damping samples "
                      "EXACTLY on an undamped resonance (f_a=w_a/2pi) => inf/NaN "
                      "FRF row at that mode. Use a positive damping ratio or a "
                      "-lin/-log grid.\n";
    }

    LadrunoFrequencyResponse fr(theModel);
    if (haveLoad) fr.setLoadPattern(loadPatternTag, amp);
    else          fr.setBaseAccel(dir, amp);
    fr.setSweep(fmin, fmax, nf, sweep);
    fr.setDamping(damp);
    fr.setResponse(nodeTag, dof, resp);
    if (!modes.empty())  fr.setModeSubset(modes);
    if (!outfile.empty()) fr.setOutputFile(outfile.c_str());

    std::vector<std::vector<double> > rows;
    if (fr.run(rows, ampOnly) < 0)
        return -1;

    OPS_SetDoubleListsOutput(rows);
    return 0;
}

int OPS_LadrunoFrequencyResponse(void)  { return OPS_FreqResponseImpl(false, "frequencyResponse"); }
int OPS_LadrunoSteadyStateDynamics(void){ return OPS_FreqResponseImpl(true,  "steadyStateDynamics"); }

// ===========================================================================
// stationary random response (PSD -> RMS)
// ===========================================================================
LadrunoRandomResponse::LadrunoRandomResponse(AnalysisModel* theModel)
    : m_model(theModel)
    , m_psd(0)
    , m_ex(LadrunoFrequencyResponse::EX_BASE_ACCEL)
    , m_direction(0)
    , m_patternTag(0)
    , m_fmin(0.0), m_fmax(0.0)
    , m_nf(0)
    , m_sweep(LadrunoFrequencyResponse::SWEEP_LIN)
    , m_nodeTag(0)
    , m_dof(0)
    , m_resp(LadrunoFrequencyResponse::RESP_DISP)
    , m_duration(0.0)
{
}

LadrunoRandomResponse::~LadrunoRandomResponse() {}

int
LadrunoRandomResponse::run(std::vector<std::vector<double> >& rows, Stats& out)
{
    rows.clear();
    out = Stats();

    if (m_psd == 0) {
        opserr << "randomResponse - an input-PSD time series is required "
                  "(-inputPSD).\n";
        return -1;
    }
    // RMS is a band INTEGRAL: a single-frequency "grid" has zero measure, so
    // unlike frequencyResponse (where nf==1 is a legitimate point query)
    // refuse it.
    if (m_nf < 2) {
        opserr << "randomResponse - the sweep needs nf >= 2 points to carry a "
                  "band integral (got nf=" << m_nf << ").\n";
        return -1;
    }

    // The complex FRF sweep does ALL the shared validation (eigen/modalProperties
    // freshness, node/dof, mode subset, damping-list length, negative eigenvalue)
    // and emits H_x(f) rows {f, Re, Im}, the exact operator the PSD propagates
    // through.  amp stays 1: the input PSD carries the excitation scale.
    LadrunoFrequencyResponse fr(m_model);
    if (m_ex == LadrunoFrequencyResponse::EX_LOAD)
        fr.setLoadPattern(m_patternTag, 1.0);
    else
        fr.setBaseAccel(m_direction, 1.0);
    fr.setSweep(m_fmin, m_fmax, m_nf, m_sweep);
    fr.setDamping(m_damp);
    fr.setResponse(m_nodeTag, m_dof, m_resp);
    if (!m_modes_1based.empty())
        fr.setModeSubset(m_modes_1based);

    std::vector<std::vector<double> > frf;
    if (fr.run(frf, false) < 0)
        return -1; // (the FRF layer already reported why, prefixed
                   //  "frequencyResponse -"; randomResponse rides it)

    // zero-damped in-band mode -> the variance integral diverges; refuse rather
    // than integrate a finite-but-meaningless (or inf-poisoned, under -biased)
    // sample of a divergent integrand.  Out-of-band zero-damped modes keep a
    // finite in-band integrand and pass.
    {
        Domain* domain = m_model->getDomainPtr();
        const Vector& ev = domain->getEigenvalues();
        std::vector<int> modes0;
        if (m_modes_1based.empty())
            for (int i = 0; i < ev.Size(); ++i) modes0.push_back(i);
        else
            for (int m : m_modes_1based) modes0.push_back(m - 1);
        for (int mode0 : modes0) {
            const double lam = ev(mode0);
            const double wa = (lam > 0.0) ? std::sqrt(lam) : 0.0;
            const double da = m_damp.coeff(mode0, wa);
            const double fa = wa / (2.0 * M_PI);
            const bool inBand = (fa >= m_fmin && fa <= m_fmax);
            // a RIGID mode's FRF diverges at its own f=0 for ANY damping
            // (denominator -Om^2 + i Om d -> 0 as Om -> 0), so d>0 does not
            // rescue it: with f=0 in the band the variance integral diverges
            // regardless of the Rayleigh a0 the transient legitimately carries
            // there.  (A da<=0 test alone would let a w=0, -rayleigh a0>0 mode
            // through to an inf/NaN RMS.)
            if (inBand && wa == 0.0) {
                opserr << "randomResponse - mode " << (mode0 + 1) << " is a "
                          "rigid-body mode (w=0) and f=0 lies inside the sweep "
                          "band: its response-variance integral diverges for "
                          "any damping. Use fmin > 0 or exclude the mode "
                          "(-modes).\n";
                return -1;
            }
            if (inBand && da <= 0.0) {
                opserr << "randomResponse - mode " << (mode0 + 1) << " (f="
                       << fa << " Hz) lies inside the sweep band with ZERO "
                          "damping: its response-variance integral diverges. "
                          "Give the mode a positive damping ratio or exclude "
                          "its resonance from the band.\n";
                return -1;
            }
        }
    }

    // propagate: G_xx(f) = |H_x(f)|^2 * G_in(f), G_in sampled at f [Hz]
    rows.reserve(frf.size());
    for (const std::vector<double>& r : frf) {
        const double f   = r[0];
        const double gin = m_psd->getFactor(f);
        if (gin < 0.0) {
            opserr << "randomResponse - the input PSD is negative at f=" << f
                   << " Hz (G=" << gin << "); a PSD must be >= 0 everywhere "
                      "on the sweep.\n";
            rows.clear();
            return -1;
        }
        const double h2 = r[1] * r[1] + r[2] * r[2];
        std::vector<double> row;
        row.push_back(f);
        row.push_back(gin);
        row.push_back(h2 * gin);
        rows.push_back(row);
    }

    // band moments (trapezoid on the sorted sweep grid)
    double m0 = 0.0, m2 = 0.0;
    for (size_t i = 1; i < rows.size(); ++i) {
        const double f0 = rows[i - 1][0], f1 = rows[i][0];
        const double g0 = rows[i - 1][2], g1 = rows[i][2];
        const double df = f1 - f0;
        m0 += 0.5 * df * (g0 + g1);
        m2 += 0.5 * df * (f0 * f0 * g0 + f1 * f1 * g1);
    }
    out.m0  = m0;
    out.m2  = m2;
    out.rms = std::sqrt(m0);
    if (m0 > 0.0)
        out.nu0 = std::sqrt(m2 / m0);

    // Davenport expected peak over a duration (optional)
    if (m_duration > 0.0) {
        const double nT = out.nu0 * m_duration;
        if (nT <= 1.0) {
            opserr << "WARNING randomResponse - nu0*T = " << nT << " <= 1: the "
                      "Davenport peak factor needs many upcrossings (nu0*T >> 1); "
                      "the peak entry is reported as NaN.\n";
        } else {
            if (nT < 2.0)
                opserr << "WARNING randomResponse - nu0*T = " << nT << " is "
                          "barely above 1: the Davenport asymptotic factor is "
                          "unreliable this close to its validity edge (needs "
                          "nu0*T >> 1); treat the peak estimate as indicative "
                          "only.\n";
            const double lnt = std::sqrt(2.0 * std::log(nT));
            out.peak = (lnt + 0.5772 / lnt) * out.rms;
            out.hasPeak = true;
        }
    }

    // optional file output: {f, G_in, G_xx} rows
    if (!m_outfile.empty()) {
        std::ofstream os(m_outfile.c_str());
        if (!os) {
            opserr << "randomResponse - cannot open output file '"
                   << m_outfile.c_str() << "'.\n";
            return -1;
        }
        os << std::scientific << std::setprecision(15);
        for (const std::vector<double>& r : rows) {
            for (size_t j = 0; j < r.size(); ++j)
                os << (j ? " " : "") << r[j];
            os << "\n";
        }
    }
    return 0;
}

// ---------------------------------------------------------------------------
// Command entry point.
//
//   randomResponse -freq fmin fmax nf [-lin|-log|-biased]
//       (-baseAccel -dir $dir | -load $patternTag) -inputPSD $tsTag
//       (-damp $xi | -rayleigh $a0 $a1 | -modalDamp $xi1 ...)
//       -node $tag -dof $dof [-resp disp|vel|accel] [-modes ...]
//       [-out $file] [-stats] [-duration $T]
//
//   -inputPSD: ONE-SIDED PSD G(f) in Hz supplied as a timeSeries sampled at f
//       (e.g. Path with f->G points, or Constant).  With -baseAccel it is the
//       PSD of the base acceleration; with -load it is the PSD of the scalar
//       s(t) multiplying the pattern's nodal-load shape P (P(t)=s(t)*P, all
//       loads fully correlated).
//   Returns the RMS (scalar).  With -stats: the list {rms, nu0, m0, m2}.
//   -duration T IMPLIES the stats list form (even without -stats) and appends
//   a 5th entry E[peak], NaN when the estimate is invalid (nu0*T <= 1), so
//   the shape is always stable at 5.  -out writes {f, G_in, G_xx}.
//
// Requires a prior `eigen N` and `modalProperties`.
// ---------------------------------------------------------------------------
int
OPS_LadrunoRandomResponse(void)
{
    const char* cmd = "randomResponse";
    AnalysisModel* theModel = *OPS_GetAnalysisModel();
    if (theModel == 0 || theModel->getDomainPtr() == 0) {
        opserr << cmd << " - no AnalysisModel/Domain available.\n";
        return -1;
    }

    double fmin = 0.0, fmax = 0.0, duration = 0.0;
    int nf = 0, dir = 0, nodeTag = 0, dof = 0, loadPatternTag = 0;
    LadrunoFrequencyResponse::SweepKind sweep = LadrunoFrequencyResponse::SWEEP_LIN;
    LadrunoFrequencyResponse::RespKind  resp  = LadrunoFrequencyResponse::RESP_DISP;
    bool haveFreq = false, haveBase = false, haveLoad = false;
    bool haveDamp = false, haveNode = false;
    bool wantStats = false;
    TimeSeries* psd = 0;
    LadrunoModalDamping damp;
    std::vector<int> modes;
    std::string outfile;
    int one = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();

        if (strcmp(opt, "-freq") == 0) {
            double fv[2]; int nd = 2;
            if (OPS_GetNumRemainingInputArgs() < 3 || OPS_GetDoubleInput(&nd, fv) < 0
                    || OPS_GetIntInput(&one, &nf) < 0) {
                opserr << cmd << " - -freq requires fmin fmax nf.\n"; return -1;
            }
            fmin = fv[0]; fmax = fv[1]; haveFreq = true;
        }
        else if (strcmp(opt, "-lin") == 0)    sweep = LadrunoFrequencyResponse::SWEEP_LIN;
        else if (strcmp(opt, "-log") == 0)    sweep = LadrunoFrequencyResponse::SWEEP_LOG;
        else if (strcmp(opt, "-biased") == 0) sweep = LadrunoFrequencyResponse::SWEEP_BIASED;
        else if (strcmp(opt, "-baseAccel") == 0) haveBase = true;
        else if (strcmp(opt, "-load") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &loadPatternTag) < 0) {
                opserr << cmd << " - -load requires a loadPattern tag.\n"; return -1;
            }
            haveLoad = true;
        }
        else if (strcmp(opt, "-dir") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &dir) < 0) {
                opserr << cmd << " - -dir requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-inputPSD") == 0) {
            int tag;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &tag) < 0) {
                opserr << cmd << " - -inputPSD requires a timeSeries tag.\n"; return -1;
            }
            psd = OPS_getTimeSeries(tag);
            if (psd == 0) {
                opserr << cmd << " - no timeSeries with tag " << tag << ".\n"; return -1;
            }
        }
        else if (strcmp(opt, "-node") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &nodeTag) < 0) {
                opserr << cmd << " - -node requires a tag.\n"; return -1;
            }
            haveNode = true;
        }
        else if (strcmp(opt, "-dof") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetInt(&one, &dof) < 0) {
                opserr << cmd << " - -dof requires a value.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-resp") == 0 || strcmp(opt, "-response") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << cmd << " - -resp requires disp|vel|accel.\n"; return -1;
            }
            const char* rname = OPS_GetString();
            if (strcmp(rname, "disp") == 0)       resp = LadrunoFrequencyResponse::RESP_DISP;
            else if (strcmp(rname, "vel") == 0)   resp = LadrunoFrequencyResponse::RESP_VEL;
            else if (strcmp(rname, "accel") == 0) resp = LadrunoFrequencyResponse::RESP_ACCEL;
            else {
                opserr << cmd << " - unknown -resp '" << rname
                       << "' (use disp|vel|accel).\n"; return -1;
            }
        }
        else if (strcmp(opt, "-damp") == 0) {
            double xi;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&one, &xi) < 0) {
                opserr << cmd << " - -damp requires a value.\n"; return -1;
            }
            // negative damping is refused: |H|^2 hides the conjugation
            // entirely, a negative xi gives an identical RMS, so it could not
            // be caught downstream.
            if (xi < 0.0) {
                opserr << cmd << " - -damp xi must be >= 0 (got " << xi << ").\n";
                return -1;
            }
            damp.kind = LadrunoModalDamping::UNIFORM; damp.xi = xi; haveDamp = true;
        }
        else if (strcmp(opt, "-rayleigh") == 0) {
            double v[2]; int nd = 2;
            if (OPS_GetNumRemainingInputArgs() < 2 || OPS_GetDoubleInput(&nd, v) < 0) {
                opserr << cmd << " - -rayleigh requires a0 a1.\n"; return -1;
            }
            damp.kind = LadrunoModalDamping::RAYLEIGH; damp.a0 = v[0]; damp.a1 = v[1];
            haveDamp = true;
        }
        else if (strcmp(opt, "-modalDamp") == 0 || strcmp(opt, "-modalDamping") == 0) {
            damp.xis.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double item;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetDoubleInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                damp.xis.push_back(item);
            }
            if (damp.xis.empty()) {
                opserr << cmd << " - -modalDamp requires >=1 ratio.\n"; return -1;
            }
            for (size_t j = 0; j < damp.xis.size(); ++j) {
                if (damp.xis[j] < 0.0) {
                    opserr << cmd << " - -modalDamp entry " << (int)j << " = "
                           << damp.xis[j] << " must be >= 0.\n";
                    return -1;
                }
            }
            damp.kind = LadrunoModalDamping::MODAL_LIST; haveDamp = true;
        }
        else if (strcmp(opt, "-modes") == 0) {
            modes.clear();
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int item;
                int rem_before = OPS_GetNumRemainingInputArgs();
                if (OPS_GetIntInput(&one, &item) < 0) {
                    if (OPS_GetNumRemainingInputArgs() < rem_before)
                        OPS_ResetCurrentInputArg(-1);
                    break;
                }
                modes.push_back(item);
            }
            if (modes.empty()) {
                opserr << cmd << " - -modes requires >=1 mode index.\n"; return -1;
            }
        }
        else if (strcmp(opt, "-out") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << cmd << " - -out requires a filename.\n"; return -1;
            }
            outfile = OPS_GetString();
        }
        else if (strcmp(opt, "-stats") == 0) wantStats = true;
        else if (strcmp(opt, "-duration") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDouble(&one, &duration) < 0) {
                opserr << cmd << " - -duration requires a value.\n"; return -1;
            }
            if (duration <= 0.0) {
                opserr << cmd << " - -duration must be > 0 (got " << duration
                       << ").\n"; return -1;
            }
            wantStats = true; // the peak estimate rides the stats output
        }
        else {
            opserr << cmd << " - unknown option '" << opt << "'.\n"; return -1;
        }
    }

    if (!haveFreq) { opserr << cmd << " - -freq fmin fmax nf is required.\n"; return -1; }
    if (haveBase && haveLoad) {
        opserr << cmd << " - -baseAccel and -load are mutually exclusive "
                  "(one excitation channel per run).\n"; return -1;
    }
    if (haveLoad) {
        if (dir != 0)
            opserr << "WARNING " << cmd << " - -dir is ignored with -load "
                      "(the pattern's nodal loads define the shape).\n";
    }
    else if (!haveBase || dir < 1) {
        opserr << cmd << " - an excitation channel is required: "
                  "-baseAccel -dir <1..ndf>, or -load <patternTag>.\n"; return -1;
    }
    if (psd == 0) {
        opserr << cmd << " - -inputPSD <tsTag> is required.\n"; return -1;
    }
    if (!haveNode || dof < 1) {
        opserr << cmd << " - -node <tag> -dof <1..ndf> are required.\n"; return -1;
    }
    if (!haveDamp) {
        opserr << cmd << " - a damping channel is required "
                  "(-damp | -rayleigh | -modalDamp).\n"; return -1;
    }

    LadrunoRandomResponse rr(theModel);
    if (haveLoad) rr.setLoadPattern(loadPatternTag);
    else          rr.setBaseAccel(dir);
    rr.setInputPSD(psd);
    rr.setSweep(fmin, fmax, nf, sweep);
    rr.setDamping(damp);
    rr.setResponse(nodeTag, dof, resp);
    if (!modes.empty())   rr.setModeSubset(modes);
    if (!outfile.empty()) rr.setOutputFile(outfile.c_str());
    if (duration > 0.0)   rr.setDuration(duration);

    std::vector<std::vector<double> > rows;
    LadrunoRandomResponse::Stats st;
    if (rr.run(rows, st) < 0)
        return -1;

    if (wantStats) {
        std::vector<double> vals;
        vals.push_back(st.rms);
        vals.push_back(st.nu0);
        vals.push_back(st.m0);
        vals.push_back(st.m2);
        // with -duration the return shape is always 5 entries: a peak that
        // could not be estimated (nu0*T <= 1, e.g. a zero response) is NaN,
        // never dropped (a 4-or-5 shape depending on the response would break
        // the caller's unpacking at run time).
        if (duration > 0.0)
            vals.push_back(st.hasPeak
                               ? st.peak
                               : std::numeric_limits<double>::quiet_NaN());
        int n = static_cast<int>(vals.size());
        if (OPS_SetDoubleOutput(&n, vals.data(), false) < 0) {
            opserr << cmd << " - failed to set the stats output.\n"; return -1;
        }
    } else {
        int n = 1;
        double rms = st.rms;
        if (OPS_SetDoubleOutput(&n, &rms, true) < 0) {
            opserr << cmd << " - failed to set the RMS output.\n"; return -1;
        }
    }
    return 0;
}
