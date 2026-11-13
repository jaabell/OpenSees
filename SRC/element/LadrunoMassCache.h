/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno (ADR-77 T2/G2 extension): shared per-instance element mass-matrix
// cache -- the LadrunoBrick pattern (which itself mirrors the Ki idiom)
// factored into a reusable helper so the other fork solids do not each
// re-implement the guard logic.
//
// CONTRACT. The element calls lookup() at the top of getMass() with
//   sig    : every MUTABLE scalar input of its mass formation OTHER than
//            nodal coordinates -- per-GP material getRho(), any element-level
//            rho override, thickness, ... . Anything fixed at construction
//            (massType/cMass flags, quadrature) may be omitted.
//   nodes  : the element's node-pointer array; coordinates are read via
//            getCrds() and guarded here, so setNodeCoord invalidates.
// On a hit it returns the cached Matrix*. On a miss (first call, or any
// guard value changed) it returns 0; the element forms its mass exactly as
// before and calls fill(), which stores a byte-copy plus the guards.
//
// SCOPE OF THE GUARANTEE (ADR-77 review wave): the guards guarantee the
// cached matrix is bit-identical to what the element's OWN formation would
// produce right now -- no more. Elements whose formation reads a setDomain
// snapshot (BezierT*'s controlPts, SolidShell's X, the SSP J-tables) re-form
// from that snapshot after setNodeCoord: the guard trips, the re-formed mass
// still reflects old geometry, cached and uncached stay byte-identical. A
// guard match is NOT proof the matrix matches current node coordinates.
// Also: a hit's Matrix* is owned here and freed on the next invalidation --
// callers must consume it immediately, never hold it across another
// getMass()/setEnabled/recvSelf (the vanilla class-static never dangled).
// And the -noMassCache escape is per-instance and NOT serialized: broker-
// built remote copies (SP/MP/DDM) default to enabled, so the A/B escape is
// effective on the building rank only.
//
// WHY GUARDS, NOT ASSUMPTIONS (ADR-76 App. A.4, made concrete by the brick's
// acceptance test): rho IS mutable mid-analysis on materials that register
// setParameter "rho" (e.g. LadrunoJ2), and coordinates via setNodeCoord.
// Invariance is therefore CHECKED on every call -- a handful of double
// compares (~ns) against a formation that costs 100s of ns to us. A guard
// mismatch invalidates and re-fills, so a stale answer is unreachable rather
// than merely unlikely. This is what separates this cache from the ADR-76
// R2 design that was withdrawn.
//
// BIT-IDENTITY (G-BYTE): the cache returns a byte-copy of the element's own
// formation with every mutable input guarded, so cached and uncached runs
// are bit-identical by construction. Proven for the brick by an 11-check
// acceptance battery; the shared-helper integrations are gated by the same
// on/off byte-equality tests (g2_mass_cache_ext_verify.py).
//
// Also retires one ADR-75b 5.4 hazard per element on the cached path: these
// getMass() implementations all return CLASS-STATIC scratch matrices; a hit
// returns the per-instance copy instead.
//
// NOT serialized (the T7/brick policy): the cache re-forms on first use
// after recvSelf. Not copyable. ~(matrix + guards) heap bytes per element.

#ifndef LadrunoMassCache_h
#define LadrunoMassCache_h

#include <Matrix.h>
#include <Node.h>
#include <Vector.h>

class LadrunoMassCache {
public:
  LadrunoMassCache() : M(0), guard(0), nGuard(0), enabled(true) {}
  ~LadrunoMassCache() { invalidate(); delete [] guard; }

  void setEnabled(bool e) { enabled = e; if (!e) invalidate(); }
  bool isEnabled() const { return enabled; }
  void invalidate() { if (M != 0) { delete M; M = 0; } }
  // the matrix fill() just stored (0 when disabled or uncacheable) -- lets an
  // element return its per-instance copy on the MISS path too (WP-124: LadrunoBrick)
  const Matrix *cached() const { return M; }

  // 0 = miss (form + fill); non-null = hit (return *it from getMass)
  const Matrix *lookup(const double *sig, int nsig,
                       Node **nodes, int nnode, int ndm) {
    if (!enabled || M == 0)
      return 0;
    if (nsig + nnode * ndm != nGuard) {   // shape changed: treat as miss
      invalidate();
      return 0;
    }
    int k = 0;
    for (int i = 0; i < nsig; i++)
      if (sig[i] != guard[k++]) { invalidate(); return 0; }
    for (int a = 0; a < nnode; a++) {
      // a warm cache can outlive the nodes (setDomain(0) on element removal /
      // collapse-recorder flows) -- the uncached formation never touches the
      // node array on those paths, so the cache must not either (review wave)
      if (nodes[a] == 0) { invalidate(); return 0; }
      const Vector &crd = nodes[a]->getCrds();
      for (int j = 0; j < ndm; j++)
        if (crd(j) != guard[k++]) { invalidate(); return 0; }
    }
    return M;
  }

  void fill(const Matrix &m, const double *sig, int nsig,
            Node **nodes, int nnode, int ndm) {
    if (!enabled)
      return;
    invalidate();
    for (int a = 0; a < nnode; a++)
      if (nodes[a] == 0)                  // pre-setDomain getMass (broker /
        return;                           // direct-C++): can't guard geometry,
                                          // stay uncached (review wave)
    M = new Matrix(m);                    // byte-copy of the fresh formation
    const int n = nsig + nnode * ndm;
    if (n != nGuard) {
      delete [] guard;
      guard = new double[n];
      nGuard = n;
    }
    int k = 0;
    for (int i = 0; i < nsig; i++)
      guard[k++] = sig[i];
    for (int a = 0; a < nnode; a++) {
      const Vector &crd = nodes[a]->getCrds();
      for (int j = 0; j < ndm; j++)
        guard[k++] = crd(j);
    }
  }

private:
  LadrunoMassCache(const LadrunoMassCache &);              // not copyable
  LadrunoMassCache &operator=(const LadrunoMassCache &);

  Matrix *M;
  double *guard;
  int nGuard;
  bool enabled;
};

#endif
