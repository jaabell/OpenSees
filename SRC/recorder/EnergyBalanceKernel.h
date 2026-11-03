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
//   T. Belytschko, W.K. Liu, B. Moran (2000). Nonlinear Finite Elements for
//     Continua and Structures, Wiley, Section 6.2.3 (energy balance check).
//   C.-M. Uang, V.V. Bertero (1990). Evaluation of seismic energy in
//     structures. Earthquake Engineering and Structural Dynamics 19(1),
//     77-90. doi:10.1002/eqe.4290190108
//
// Description: the per-entity energy terms and the per-scope work integrals
// of the structural-dynamics energy balance used by EnergyBalanceRecorder.
// Header-only, namespace ebkernel.
//
// Six components, in this order: KE, IE, DW, ULW, RES, ERR where
//   KE  = 1/2 v^T M v              (instantaneous; element + nodal mass)
//   IE  = int F_resisting^T v dt   (trapezoidal)
//   DW  = int v^T C v dt           (trapezoidal; element + nodal damping)
//   ULW = int v^T P_ext dt         (trapezoidal; Node::getUnbalancedLoad)
//   RES = ULW - (KE + IE + DW)     (closure residual)
//   ERR = 100 * RES / E_ref        (E_ref = running max of
//                                   |KE|+|IE|+|DW|+|ULW|, per scope)
//
// Element terms:
//   * the element mass is copied before getDamp() is called, because many
//     elements return the same static scratch Matrix from getMass() and
//     getDamp() (and the base Element::getDamp() calls getMass());
//   * velocities are gathered node by node with zero fill;
//   * empty (0x0) mass or damping matrices are skipped;
//   * IE uses the element resisting force F dotted with v.

#ifndef EnergyBalanceKernel_h
#define EnergyBalanceKernel_h

#include <Domain.h>
#include <Node.h>
#include <Element.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>

#include <LadrunoEnergyChannels.h>

#include <cmath>

namespace ebkernel {

static const int NUM_ENERGY_COMPONENTS = 6;

// Decide whether an element's damping matrix C describes the damping forces
// it adds outside getResistingForce(). Returns true, and the power v.g of the
// actual extra force g = F_incInertia - F - M a in gRate, when g differs from
// C v by more than round-off, i.e. when part of C is the damping tangent of a
// material whose viscous force is already inside F (or the element damping
// force is otherwise not C v). Returns false, keeping v^T C v, when C is zero,
// v is zero, or g equals C v to round-off.
// Must be called after C and F have been read: getResistingForceIncInertia()
// overwrites the scratch storage that many elements return from getDamp() and
// getResistingForce().
inline bool dampingOutsideResistingForce(Element *ele, const Matrix &C,
                                         const Vector &F, const Matrix &M,
                                         bool haveMass, const Vector &vel,
                                         int numDOF, double &gRate)
{
    gRate = 0.0;
    if (F.Size() != numDOF)
        return false;

    bool vZero = true, cZero = true;
    for (int i = 0; i < numDOF && vZero; ++i)
        if (vel(i) != 0.0) vZero = false;
    if (vZero)
        return false;
    for (int i = 0; i < numDOF && cZero; ++i)
        for (int j = 0; j < numDOF; ++j)
            if (C(i, j) != 0.0) { cZero = false; break; }
    if (cZero)
        return false;

    static Vector Fs, Cv, Ma, acc;
    if (Fs.Size() != numDOF) {
        Fs.resize(numDOF); Cv.resize(numDOF); Ma.resize(numDOF); acc.resize(numDOF);
    }
    // Copy F and form C v while their storage is still valid.
    for (int i = 0; i < numDOF; ++i) {
        Fs(i) = F(i);
        double s = 0.0;
        for (int j = 0; j < numDOF; ++j)
            s += C(i, j) * vel(j);
        Cv(i) = s;
    }
    Ma.Zero();
    if (haveMass) {
        Node **nodes = ele->getNodePtrs();
        const int nn = ele->getNumExternalNodes();
        int cnt = 0;
        for (int i = 0; i < nn && cnt < numDOF; ++i) {
            const Vector &a = nodes[i]->getAccel();
            for (int j = 0; j < a.Size() && cnt < numDOF; ++j)
                acc(cnt++) = a(j);
        }
        for (; cnt < numDOF; ++cnt)
            acc(cnt) = 0.0;
        Ma.addMatrixVector(0.0, M, acc, 1.0);
    }

    const Vector &Finc = ele->getResistingForceIncInertia();
    if (Finc.Size() != numDOF)
        return false;

    double scale = 0.0, diff = 0.0, rate = 0.0;
    for (int i = 0; i < numDOF; ++i) {
        const double g = Finc(i) - Fs(i) - Ma(i);
        scale = std::fmax(scale, std::fmax(std::fmax(std::fabs(Finc(i)), std::fabs(Fs(i))),
                                           std::fmax(std::fabs(Ma(i)), std::fabs(Cv(i)))));
        diff = std::fmax(diff, std::fabs(g - Cv(i)));
        rate += g * vel(i);
    }
    if (diff <= 1.0e-12 * scale)
        return false;
    gRate = rate;
    return true;
}

// One element's instantaneous energy terms:
//   ke = 1/2 v^T M v   (kinetic energy)
//   ie = F^T v         (internal power; the caller integrates over dt)
//   dw = damping power (the caller integrates over dt)
struct ElementEnergyTerms {
    double ke, ie, dw;
};

inline void elementEnergyTerms(Element *ele, ElementEnergyTerms &t)
{
    t.ke = 0.0; t.ie = 0.0; t.dw = 0.0;
    // A Subdomain (the ShadowSubdomains on P0 of a PartitionedDomain) has DOFs
    // and external nodes but no node pointers; its energy is accounted on its
    // own process.
    if (ele->isSubdomain() || ele->getNodePtrs() == 0)
        return;

    static Matrix M;      // copy of the mass (see the file comment)
    static Vector vel;
    M = ele->getMass();
    const Matrix &C = ele->getDamp();
    const Vector &F = ele->getResistingForce();
    const int numExternalNodes = ele->getNumExternalNodes();
    const int numDOF = ele->getNumDOF();
    Node **elenodes = ele->getNodePtrs();
    if (vel.Size() < numDOF)
        vel.resize(numDOF);

    int cnt = 0;
    for (int i = 0; i < numExternalNodes && cnt < numDOF; ++i) {
        const Vector &node_vel = elenodes[i]->getVel();
        for (int j = 0; j < node_vel.Size() && cnt < numDOF; ++j)
            vel(cnt++) = node_vel(j);
    }
    for (; cnt < numDOF; ++cnt)
        vel(cnt) = 0.0;

    const bool haveMass = (M.noRows() == numDOF && M.noCols() == numDOF);
    const bool haveDamp = (C.noRows() == numDOF && C.noCols() == numDOF);
    double kinetic = 0.0;
    if (haveMass)
        for (int i = 0; i < numDOF; ++i)
            for (int j = 0; j < numDOF; ++j)
                kinetic += 0.5 * M(i, j) * vel(i) * vel(j);
    t.ke = kinetic;

    // IE is formed before the damping term because the check below calls
    // getResistingForceIncInertia(), which many elements write into the same
    // scratch vector getResistingForce() returned.
    if (F.Size() == numDOF) {
        double dot = 0.0;
        for (int i = 0; i < numDOF; ++i)
            dot += F(i) * vel(i);
        t.ie = dot;
    }

    if (haveDamp) {
        double dr = 0.0;
        for (int i = 0; i < numDOF; ++i)
            for (int j = 0; j < numDOF; ++j)
                dr += C(i, j) * vel(i) * vel(j);

        // DW is the power of the element damping forces that are not inside
        // getResistingForce() (IE already integrates F.v). For Rayleigh
        // damping those forces are C v and v^T C v is right. An element whose
        // getDamp() also returns the damping tangent of a rate-dependent
        // material (ZeroLength, Truss, TwoNodeLink or bearings with Viscous,
        // ViscousDamper, Parallel(..., Viscous), py/qz/tz dashpots) already
        // has that material's viscous force inside F, so v^T C v would count
        // the dashpot a second time (RES = -dissipation for a linear
        // dashpot). The force the element adds beyond F and inertia is
        //     g = F_incInertia - F - M a ;
        // when g equals C v to round-off (every Rayleigh-only element) v^T C v
        // is kept, otherwise DW uses v.g.
        double gRate = 0.0;
        if (dampingOutsideResistingForce(ele, C, F, M, haveMass, vel, numDOF, gRate))
            t.dw = gRate;
        else
            t.dw = dr;
    }
}

// Nodal contributions: lumped-mass kinetic energy, nodal damping
// (Node::getDamp, alphaM times the nodal mass) and the work rate of the
// external nodal load. M is fully consumed before getDamp() is called because
// Node::getMass() and Node::getDamp() can share the same scratch matrix.
inline void addNodeEnergy(Node *node,
                          double &kinetic,
                          double &damping_rate,
                          double &unbalanced_rate)
{
    const Vector &vel = node->getVel();
    const int n = vel.Size();

    const Matrix &M = node->getMass();
    if (M.noRows() == n && M.noCols() == n)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                kinetic += 0.5 * M(i, j) * vel(i) * vel(j);

    const Matrix &C = node->getDamp();
    if (C.noRows() == n && C.noCols() == n)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                damping_rate += C(i, j) * vel(i) * vel(j);

    const Vector &Pext = node->getUnbalancedLoad();
    if (Pext.Size() == n)
        unbalanced_rate += vel ^ Pext;
}

// Per-scope accumulator: the trapezoidal work integrals and the closure and
// ERR arithmetic of the six-column layout.
struct EnergyAccumulator {
    double internal_energy, damping_work, unbalanced_load_work;
    double prev_internal_rate, prev_damping_rate, prev_unbalanced_rate;
    double eref;          // running max of |KE|+|IE|+|DW|+|ULW| for this scope

    EnergyAccumulator() { reset(); }

    void reset() {
        internal_energy = 0.0;
        damping_work = 0.0;
        unbalanced_load_work = 0.0;
        prev_internal_rate = 0.0;
        prev_damping_rate = 0.0;
        prev_unbalanced_rate = 0.0;
        eref = 0.0;
    }

    // Integrate one step. The first recorded step uses the rectangle rule
    // (there is no previous rate yet; skipping the increment would leave a
    // constant offset in RES), later steps the trapezoidal rule.
    // Fills out[6] = {KE, IE, DW, ULW, RES, ERR}.
    void step(double dT, bool firstRecord,
              double ke, double internal_rate, double damping_rate, double unbalanced_rate,
              double out[NUM_ENERGY_COMPONENTS])
    {
        if (firstRecord) {
            internal_energy      += internal_rate   * dT;
            damping_work         += damping_rate    * dT;
            unbalanced_load_work += unbalanced_rate * dT;
        } else {
            internal_energy      += 0.5 * (prev_internal_rate   + internal_rate)   * dT;
            damping_work         += 0.5 * (prev_damping_rate    + damping_rate)    * dT;
            unbalanced_load_work += 0.5 * (prev_unbalanced_rate + unbalanced_rate) * dT;
        }

        prev_internal_rate   = internal_rate;
        prev_damping_rate    = damping_rate;
        prev_unbalanced_rate = unbalanced_rate;

        const double sum = ke + internal_energy + damping_work;
        const double res = unbalanced_load_work - sum;
        // The energy scale is the running max of the summed magnitudes; it
        // does not collapse to zero when KE and IE cancel (free vibration).
        const double mag = std::fabs(ke) + std::fabs(internal_energy)
                         + std::fabs(damping_work) + std::fabs(unbalanced_load_work);
        if (mag > eref) eref = mag;
        const double err = (eref > 1.0e-16) ? 100.0 * res / eref : 0.0;

        out[0] = ke;
        out[1] = internal_energy;
        out[2] = damping_work;
        out[3] = unbalanced_load_work;
        out[4] = res;
        out[5] = err;
    }
};

// ---------------------------------------------------------------------------
// Channel-aware accumulator with an element/nodal split (-v2 layout).
//
// Differences from EnergyAccumulator:
//   * KE and DW are carried as element and nodal parts. Under MPI the element
//     parts reduce correctly by summation across ranks (elements are uniquely
//     partitioned); the nodal parts are counted on every rank that holds a
//     shared boundary node and are kept in their own columns so an offline
//     reduction can treat them separately.
//   * Reads the energy channels of LadrunoEnergyChannels.h: ABSORB_LEAK L
//     (the resisting-sign injection work contained in the raw IE integral),
//     LNVD_WORK and MODAL_WORK (dissipation applied by the integrator,
//     invisible to element and node getDamp()). The hourglass total Hg is
//     the sum of the element "hourglassEnergy" responses. All four are
//     cumulative; a baseline is taken at the first record and deltas are
//     used, so a new recorder starts at zero.
//   * Balance:  IE_shown = IE_raw - L - Hg
//               E_inject = -L                (source, positive = energy in)
//               RES = (ULW + E_inject)
//                     - (KE + IE_shown + DW + E_lnvd + E_modal + E_hg)
//     L and Hg cancel in RES (moving them between terms cannot fake
//     closure); E_lnvd and E_modal close what the six-column layout leaves
//     in RES. For a viscous hourglass control Hg is cumulative dissipation;
//     for a stiffness control it is the stored stabilization energy, which
//     can decrease on unloading.
//
// Output slot order (fixed; the recorder skips undeclared channel slots when
// writing): KE_ele KE_nod IE DW_ele DW_nod ULW E_inject E_lnvd E_modal E_hg
// RES ERR.
// ---------------------------------------------------------------------------
static const int NUM_V2_SLOTS = 12;
enum V2Slot {
    V2_KE_ELE = 0, V2_KE_NOD, V2_IE, V2_DW_ELE, V2_DW_NOD, V2_ULW,
    V2_E_INJECT, V2_E_LNVD, V2_E_MODAL, V2_E_HG, V2_RES, V2_ERR
};

struct EnergyAccumulatorV2 {
    double internal_energy;                       // raw IE (may contain L)
    double damping_work_ele, damping_work_nod;
    double unbalanced_load_work;
    double prev_internal_rate;
    double prev_damping_rate_ele, prev_damping_rate_nod;
    double prev_unbalanced_rate;
    double eref;                                  // running max, this scope
    bool   baselined;
    double base_leak, base_lnvd;                  // channel totals at first record
    double base_modal, base_hg;

    EnergyAccumulatorV2() { reset(); }

    void reset() {
        internal_energy = 0.0;
        damping_work_ele = 0.0;
        damping_work_nod = 0.0;
        unbalanced_load_work = 0.0;
        prev_internal_rate = 0.0;
        prev_damping_rate_ele = 0.0;
        prev_damping_rate_nod = 0.0;
        prev_unbalanced_rate = 0.0;
        eref = 0.0;
        baselined = false;
        base_leak = 0.0;
        base_lnvd = 0.0;
        base_modal = 0.0;
        base_hg = 0.0;
    }

    // Same integration rule as EnergyAccumulator::step(). Channel totals are
    // read from the registry by the caller and passed in, so the arithmetic
    // here does not depend on the registry.
    void step(double dT, bool firstRecord,
              double ke_ele, double ke_nod,
              double internal_rate,
              double damping_rate_ele, double damping_rate_nod,
              double unbalanced_rate,
              double leak_total, double lnvd_total,
              double modal_total, double hg_total,
              double out[NUM_V2_SLOTS])
    {
        if (firstRecord) {
            internal_energy      += internal_rate    * dT;
            damping_work_ele     += damping_rate_ele * dT;
            damping_work_nod     += damping_rate_nod * dT;
            unbalanced_load_work += unbalanced_rate  * dT;
        } else {
            internal_energy      += 0.5 * (prev_internal_rate    + internal_rate)    * dT;
            damping_work_ele     += 0.5 * (prev_damping_rate_ele + damping_rate_ele) * dT;
            damping_work_nod     += 0.5 * (prev_damping_rate_nod + damping_rate_nod) * dT;
            unbalanced_load_work += 0.5 * (prev_unbalanced_rate  + unbalanced_rate)  * dT;
        }

        prev_internal_rate    = internal_rate;
        prev_damping_rate_ele = damping_rate_ele;
        prev_damping_rate_nod = damping_rate_nod;
        prev_unbalanced_rate  = unbalanced_rate;

        if (!baselined) {
            base_leak = leak_total;
            base_lnvd = lnvd_total;
            base_modal = modal_total;
            base_hg = hg_total;
            baselined = true;
        }
        const double L  = leak_total - base_leak;    // resisting-sign leak
        const double Dl = lnvd_total - base_lnvd;    // LNVD dissipation (>= 0)
        const double Dm = modal_total - base_modal;  // modal dissipation (>= 0)
        const double Hg = hg_total - base_hg;        // hourglass part of IE

        const double ke      = ke_ele + ke_nod;
        const double dw      = damping_work_ele + damping_work_nod;
        const double ie_disp = internal_energy - L - Hg;
        const double e_inj   = -L;

        const double res = (unbalanced_load_work + e_inj)
                         - (ke + ie_disp + dw + Dl + Dm + Hg);

        const double mag = std::fabs(ke) + std::fabs(ie_disp) + std::fabs(dw)
                         + std::fabs(unbalanced_load_work)
                         + std::fabs(e_inj) + std::fabs(Dl)
                         + std::fabs(Dm) + std::fabs(Hg);
        if (mag > eref) eref = mag;
        const double err = (eref > 1.0e-16) ? 100.0 * res / eref : 0.0;

        out[V2_KE_ELE]   = ke_ele;
        out[V2_KE_NOD]   = ke_nod;
        out[V2_IE]       = ie_disp;
        out[V2_DW_ELE]   = damping_work_ele;
        out[V2_DW_NOD]   = damping_work_nod;
        out[V2_ULW]      = unbalanced_load_work;
        out[V2_E_INJECT] = e_inj;
        out[V2_E_LNVD]   = Dl;
        out[V2_E_MODAL]  = Dm;
        out[V2_E_HG]     = Hg;
        out[V2_RES]      = res;
        out[V2_ERR]      = err;
    }
};

} // namespace ebkernel

#endif // EnergyBalanceKernel_h
