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
//   Wriggers, P. (2006). Computational Contact Mechanics, 2nd ed. Springer,
//     Berlin. doi:10.1007/978-3-540-32609-0
//   Laursen, T.A. (2002). Computational Contact and Impact Mechanics.
//     Springer, Berlin. doi:10.1007/978-3-662-04864-1
//   Hallquist, J.O. (2006). LS-DYNA Theory Manual, Livermore Software
//     Technology Corporation, Section 26 (Contact-Impact Algorithm).
//
// Description: LadrunoContactDomain, the contact engine attached to a Domain.
// It owns the contact surfaces and the contact definitions, and the
// path-dependent pair state, so that the assembly-side adapters
// (LadrunoContactFE, rebuilt in every ConstraintHandler::handle()) are
// stateless views.
//
// Lifetime: a Domain owns at most one engine (none means no cost and stock
// behaviour). The engine survives domainChanged() (which clears the
// AnalysisModel, not the Domain) and is deleted by Domain::clearAll() (wipe)
// and by ~Domain. Domain::commit(), revertToLastCommit() and revertToStart()
// call the matching engine methods, so path-dependent contact state commits
// and reverts with the rest of the model whatever the integrator.

#ifndef LadrunoContactDomain_h
#define LadrunoContactDomain_h

#include <vector>
#include <map>
#include <set>

class LadrunoContactSurface;
class Domain;

class LadrunoContactDomain
{
  public:
    LadrunoContactDomain();
    ~LadrunoContactDomain();

    // --- model definition (contactSurface / contact commands) ---
    int addSurface(LadrunoContactSurface *surf);          // takes ownership
    LadrunoContactSurface *getSurface(int tag) const;

    // Reconcile the declared contact with a domain the analysis has been
    // shrinking (`remove node`, element removal). Called at the start of every
    // handle(): removed nodes are pruned from every surface, and an interaction
    // whose master or slave surface is left empty is retired; both are reported
    // on opserr. Returns the number of interactions retired by this call.
    int pruneRemovedNodes(Domain *dom);

    // A node-to-segment contact between a MASTER_SEGMENTS surface and a
    // SLAVE_NODES surface; the handler builds one LadrunoContactFE per
    // (slave node, candidate master segment) pair.
    //   outward:  optional direction toward the slave's allowed half-space,
    //             used to orient the derived segment normal; null = derived
    //             per pair from the slave and segment reference positions.
    //   knAuto:   kn is sized per pair at handle() time from the stiffness of
    //             the solid element that owns the master segment (kn is then
    //             a placeholder).
    //   cellFrac: broad-phase cell size as a multiple of the median segment
    //             diagonal (1.0 default; a huge value gives one bucket, i.e.
    //             brute force).
    //   consistentNormal: add the geometric part kn gN d^2gN/du^2 of the
    //             consistent normal tangent (curved or large-sliding masters).
    //   kt, mu, cohesion, tauMax: penalty-regularized friction, tangential
    //             penalty kt and cone cap = min(mu N + cohesion, tauMax)
    //             (tauMax <= 0: no cap); all of mu, cohesion, tauMax <= 0 is
    //             frictionless.
    //   consistentTan: assemble the non-symmetric consistent friction tangent
    //             (false: the symmetric tangent, safe on any linear solver).
    //   muc:      viscous normal stabilization, a dashpot muc on the normal gap
    //             rate while in contact (0: none).
    int addContact(int tag, int masterSurfTag, int slaveSurfTag,
                   double kn, double kt, double mu, const double *outward = 0,
                   bool knAuto = false, double cellFrac = 1.0,
                   bool consistentTan = false, bool consistentNormal = false,
                   double cohesion = 0.0, double tauMax = 0.0, double muc = 0.0);
    int getNumContacts(void) const { return (int)theContacts.size(); }

    struct Contact {
        // Set when the analysis removed every node of one of the surfaces of
        // this interaction. A retired interaction is skipped by handle() with a
        // named message; it cannot hide a typo, because a node tag that does
        // not exist is refused when the surface is declared.
        bool retired = false;
        int tag, masterSurfTag, slaveSurfTag;
        double kn, kt, mu;
        double cohesion;        // cone intercept c
        double tauMax;          // shear cap of the cone (<= 0: none)
        bool knAuto;            // kn sized from the master solid element
        bool hasOutward;        // an explicit orientation direction was given
        double outward[3];      // orientation direction toward the allowed half-space
        double cellFrac;        // bucket-sort cell = cellFrac * median segment diagonal
        bool consistentNormal;  // assemble the geometric normal tangent
        bool consistentTan;     // non-symmetric consistent friction tangent
        double muc;             // viscous normal stabilization coefficient (0: none)
        bool hasFriction(void) const { return mu > 0.0 || cohesion > 0.0 || tauMax > 0.0; }
    };
    const Contact &getContact(int i) const { return theContacts[i]; }

    // --- rigid analytical plane (point p0, outward unit normal n) against a
    //     slave node set; one adapter per slave node. ---
    struct RigidPlane {
        bool retired = false;   // set when the analysis removed every slave node
        int tag, slaveSurfTag;
        double p0[3], n[3];     // n stored normalized
        double kn;
        double muc;             // viscous normal stabilization coefficient (0: none)
    };
    int addRigidPlane(int tag, int slaveSurfTag,
                      const double p0[3], const double n[3], double kn, double muc = 0.0);
    int getNumRigidPlanes(void) const { return (int)theRigidPlanes.size(); }
    const RigidPlane &getRigidPlane(int i) const { return theRigidPlanes[i]; }

    // --- node-to-segment normal force readout (ladrunoContactForce) ---
    // Each segment adapter writes its normal penalty force tn = kn <-gap>+ into a
    // (contactTag, slaveTag, segIndex) slot at every residual evaluation, and
    // writes 0 when the pair is not active, so after convergence the slots hold
    // the converged forces of the pairs that are in contact. getNtsForce sums the
    // slots of one slave node. The slots are cleared at every handle(). No effect
    // on the assembled residual or tangent.
    void   clearNtsForce(void);
    void   setNtsForce(int contactTag, int slaveTag, int segIndex, double tn);
    double getNtsForce(int slaveTag) const;   // sum over (contactTag, slaveTag, *)

    // --- friction pair state, keyed (contactTag, slaveNodeTag, global master
    //     segment ordinal), which is stable across adapter rebuilds. The adapter
    //     reads the committed plastic slip gpT and the engagement origin gT0,
    //     runs the return map and writes gpTtrial; commit() promotes the trial
    //     and revertToLastCommit() drops it. ---
    struct FrictionState {
        double gpT[3];        // committed plastic slip from the engagement configuration
        double gpTtrial[3];   // trial plastic slip (written at every residual evaluation)
        double gT0[3];        // tangential origin at engagement (captured once)
        bool   engaged;       // gT0 has been captured
        // gT0 and engaged are set during a residual evaluation, so a rejected
        // step must not keep an origin captured at its trial configuration:
        // commit() copies them here and revertToLastCommit() restores them.
        double gT0committed[3];
        bool   engagedCommitted;
        FrictionState() : engaged(false), engagedCommitted(false) {
            for (int d = 0; d < 3; d++)
                gpT[d] = gpTtrial[d] = gT0[d] = gT0committed[d] = 0.0;
        }
    };
    // the slot of a pair, created zeroed (engaged = false) on first use
    FrictionState &getOrCreateFrictionState(int contactTag, int slaveTag, int segIndex);

    // The engine survives domainChanged() and later analyses, so slots of pairs
    // that no longer exist (re-meshing, element removal) are collected: the
    // handler calls frictionGCBegin(), frictionGCMark() for every frictional
    // pair it builds, and frictionGCEnd(), which erases the unmarked slots.
    void frictionGCBegin(void);
    void frictionGCMark(int contactTag, int slaveTag, int segIndex);
    void frictionGCEnd(void);
    int  getNumFrictionStates(void) const { return (int)theFrictionStates.size(); }

    // --- lifecycle (driven by Domain::commit / revertToLastCommit / revertToStart) ---
    // augmenting: the Domain's within-step augmentation flag, passed for
    // contact formulations that update multipliers on commit.
    int commit(bool augmenting = false);
    int revertToLastCommit(void);
    // Drop all path-dependent state so that Domain::revertToStart() returns
    // contact to the state of a fresh run; the definitions are kept.
    int revertToStart(void);
    int getNumCommits(void) const { return numCommits; }
    int getNumReverts(void) const { return numReverts; }

    // One-time diagnostics, latched per (contact tag, topic) on the engine, so
    // they reset with the model (wipe deletes the engine) and a second contact
    // with the same condition is still reported. True exactly once per key.
    enum WarnTopic {
        WARN_VISC_STATIC = 0        // -visc under a static integrator: dashpot disabled
    };
    bool warnOnce(int contactTag, int topic);

  private:
    std::vector<LadrunoContactSurface *> theSurfaces;   // owned
    std::vector<Contact> theContacts;
    std::vector<RigidPlane> theRigidPlanes;
    // The contact tag is the leading key of every path-state store, so tags
    // must be unique across all contact definitions (checked by every add*).
    bool contactTagInUse(int tag) const;
    int numCommits;
    int numReverts;

    // (contactTag, slaveTag, segIndex) key of the pair state stores.
    struct PairKey {
        int c, s, g;
        bool operator<(const PairKey &o) const {
            if (c != o.c) return c < o.c;
            if (s != o.s) return s < o.s;
            return g < o.g;
        }
    };
    std::map<PairKey, double> theNtsForce;              // per-pair normal force readout
    std::map<PairKey, FrictionState> theFrictionStates;
    std::set<PairKey> liveKeys;                         // collection scratch (per handle())
    std::set<std::pair<int, int> > theWarnOnce;         // (contactTag, WarnTopic) latches
};

#endif
