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
// Description: named energy channels read by EnergyBalanceRecorder.
//
// The recorder computes kinetic energy, internal work, damping work and
// external nodal-load work from element and node quantities (getMass,
// getDamp, getResistingForce, getUnbalancedLoad). Two classes of energy
// escape those quantities:
//   * closure gaps: energy applied inside the integrator solve (modal
//     damping, local non-viscous damping) reaches no element or node
//     quantity, so it ends up in the closure residual RES;
//   * attribution gaps: energy that is read but booked in the wrong term.
//     The incident-wave injection of an absorbing boundary is the main case:
//     the injection force is returned by the element's getResistingForce(),
//     so the internal-work integral IE absorbs the injection work with
//     resisting-force sign. RES can then close while IE is wrong.
//
// The recorder holds only a Domain and has no handle on the integrator or on
// element internals, so producers publish into this registry and the
// recorder reads it at record time.
//
// Semantics:
//   * Channels store cumulative energy (units of work), not rates and not
//     per-step increments. Producers integrate at their own cadence (per
//     commit for elements, per sub-step for integrators) and add increments.
//     A recorder that samples every N steps then reads a running total,
//     which does not alias, whereas trapezoidal integration of a published
//     rate at recorder cadence would.
//   * Totals are monotone over the lifetime of a model and are reset by
//     Domain::clearAll() (the wipe path of every interpreter), see
//     resetOnWipe(). There is no per-step clear. The recorder snapshots each
//     channel at its first record() and reports the change since then, so a
//     recorder created during an analysis starts from zero.
//   * declare(c) marks a channel as relevant without adding energy (an
//     absorbing element declares in setDomain even if no incident motion is
//     ever applied). The recorder writes a column for every declared
//     channel; undeclared channels cost nothing and print nothing.
//
// Channel conventions:
//   * ABSORB_LEAK stores L, the sum over absorbing elements of the
//     cumulative resisting-sign injection work int R_inj^T v dt, which is
//     exactly the quantity contained in the recorder's raw IE integral. The
//     recorder derives both corrections from it: IE_shown = IE_raw - L and
//     E_inject = -L. Moving L between terms cannot change RES (it cancels);
//     it makes IE and E_inject correct.
//   * LNVD_WORK stores the cumulative local non-viscous damping dissipation
//     sum alpha*|r_i|*|v_i|*dt >= 0 (a true sink).
//   * MODAL_WORK stores the cumulative modal-damping dissipation
//     v^T C_modal v integrated per commit (a true sink).
//
// Single-threaded, like the rest of the analysis core. Header-only.

#ifndef LadrunoEnergyChannels_h
#define LadrunoEnergyChannels_h

#include <limits>
#include <OPS_Globals.h>

namespace Ladruno {

class EnergyChannelRegistry {
public:
    enum Channel {
        ABSORB_LEAK = 0,   // cumulative resisting-sign injection work L
        LNVD_WORK   = 1,   // cumulative local non-viscous damping dissipation
        MODAL_WORK  = 2,   // cumulative modal-damping dissipation
        NUM_CHANNELS = 3
    };

    static EnergyChannelRegistry &instance() {
        static EnergyChannelRegistry theRegistry;
        return theRegistry;
    }

    // Mark a channel relevant (column written) without adding energy.
    void declare(Channel c) {
        if (c >= 0 && c < NUM_CHANNELS)
            isDeclared[c] = true;
    }

    // Accumulate dE (units of work) into a channel; declares it as a side
    // effect so a producer cannot publish into an invisible channel.
    // A non-finite dE (a diverging analysis) is discarded with a warning
    // printed once per process: a single NaN would otherwise make the channel
    // total, and every RES computed from it, NaN until the next wipe. The
    // diverged run's own balance is already non-finite through the element
    // and node terms, so no information is lost.
    void addEnergy(Channel c, double dE) {
        if (c >= 0 && c < NUM_CHANNELS) {
            isDeclared[c] = true;
            if (!(dE == dE) || dE == std::numeric_limits<double>::infinity()
                            || dE == -std::numeric_limits<double>::infinity()) {
                if (!warnedNonFinite) {
                    warnedNonFinite = true;
                    opserr << "WARNING EnergyChannelRegistry: discarded a "
                              "non-finite energy increment (diverged analysis?); "
                              "printed once\n";
                }
                return;
            }
            total[c] += dE;
        }
    }

    bool declared(Channel c) const {
        return (c >= 0 && c < NUM_CHANNELS) ? isDeclared[c] : false;
    }

    bool anyDeclared() const {
        for (int i = 0; i < NUM_CHANNELS; ++i)
            if (isDeclared[i])
                return true;
        return false;
    }

    // Cumulative total since the last wipe; consumers subtract their own
    // first-record baseline.
    double energy(Channel c) const {
        return (c >= 0 && c < NUM_CHANNELS) ? total[c] : 0.0;
    }

    // Called from Domain::clearAll(). wipe() destroys every producer and
    // consumer, so it is the natural zero point. Resetting there also keeps
    // the totals at full precision for the next model: a large total left by
    // a previous model (e.g. 1e300 before an overflow) would absorb the next
    // model's small increments (total + dE == total). Declared flags reset as
    // well, so the channel columns describe the current model only.
    void resetOnWipe() {
        for (int i = 0; i < NUM_CHANNELS; ++i) {
            total[i] = 0.0;
            isDeclared[i] = false;
        }
    }

private:
    EnergyChannelRegistry() {
        for (int i = 0; i < NUM_CHANNELS; ++i) {
            total[i] = 0.0;
            isDeclared[i] = false;
        }
        warnedNonFinite = false;
    }
    EnergyChannelRegistry(const EnergyChannelRegistry &);
    EnergyChannelRegistry &operator=(const EnergyChannelRegistry &);

    double total[NUM_CHANNELS];
    bool isDeclared[NUM_CHANNELS];
    bool warnedNonFinite;
};

} // namespace Ladruno

#endif // LadrunoEnergyChannels_h
