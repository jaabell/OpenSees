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
//   Dassault Systemes Simulia. Abaqus Theory Guide, Section 5.1.1
//     (contact formulations): unit normals computed at every master node by
//     averaging the adjacent segment normals give a smoothly varying normal
//     field N(X).
//   Wriggers, P. (2006). Computational Contact Mechanics, 2nd ed. Springer,
//     Berlin. doi:10.1007/978-3-540-32609-0 (smoothing of faceted surfaces).
//
// Description: LadrunoContactNormalField, header-only averaged nodal-normal
// smoothing of a faceted master surface for node-to-segment contact. It uses
// no OpenSees types (raw int/double arrays and the STL only).
//
// The per-segment faceted normal of LadrunoContactProjection::normalOriented
// takes its sign per (slave, segment) pair from the slave's reference side, so
// a slave sliding over a convex ridge can see the sign flip, and the normal
// jumps across segment junctions. The averaged field shares each node's normal
// between its incident segments (C0-continuous across junctions, not C1), and
// the outward sense is one global surface datum captured once.
//
// Three stages:
//   propagateOrientation(): flood-fill manifold orientation, a per-segment
//     winding sign sigma_s that makes the windings coherent from a seed.
//     Topological (connectivity only), so the caller can cache it. Refuses a
//     non-manifold edge (more than 2 segments), a non-orientable surface (sign
//     conflict on a cycle), a disconnected multi-shell surface and a closed
//     surface.
//   nodalNormals(): per node, the area-weighted (Newell) sum of the coherently
//     wound adjacent facet normals with the global outward sign, scattered to a
//     per-segment-per-node layout.
//   smoothNormal(): interpolation of a segment's nodal normals with the shape
//     functions, renormalized; false on a degenerate blend.
//
// No tangent terms are provided here: a consumer assembles the symmetric
// frozen-field tangent kn B^T B.

#ifndef LadrunoContactNormalField_h
#define LadrunoContactNormalField_h

#include <cmath>
#include <vector>
#include <map>

namespace LadrunoContactNormalField {

// ------------------------------------------------------------------ small vec3
inline double nf_dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline double nf_norm(const double a[3]) { return std::sqrt(nf_dot(a, a)); }

// ----------------------------------------------------------------- status codes
// Orientation outcome. Anything but OK: the surface is refused for smoothing
// and the caller keeps the faceted normal.
enum Status {
    OK             = 0,
    NON_MANIFOLD   = 1,   // an edge is shared by more than 2 segments
    DISCONNECTED   = 2,   // not all segments reachable from the seed (multi-shell)
    NON_ORIENTABLE = 3,   // a sign conflict on a cycle (Moebius strip, folded surface)
    CLOSED         = 4    // no boundary edges (a closed shell): one global sign cannot
                          // serve both sides; the supported target is an open patch
};

// Area-weighted facet normal of one segment by Newell's method: the magnitude
// is the facet area and the direction is the face normal in node-winding order
// (exact for a planar facet, well defined for a warped quad). P = nps x 3,
// row-major. Winding-dependent by design; propagateOrientation() makes the
// windings coherent.
inline void newellAreaNormal(int nps, const double *P, double out[3]) {
    out[0] = out[1] = out[2] = 0.0;
    for (int i = 0; i < nps; i++) {
        const double *a = P + (size_t)i * 3;
        const double *b = P + (size_t)((i + 1) % nps) * 3;
        out[0] += (a[1] - b[1]) * (a[2] + b[2]);
        out[1] += (a[2] - b[2]) * (a[0] + b[0]);
        out[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    for (int d = 0; d < 3; d++) out[d] *= 0.5;   // the Newell sum is 2 area n; halve it
}

// ------------------------------------------------ flood-fill manifold orientation
// Fill sigma[nSeg] (each +1/-1) with a coherent winding from segment 0
// (sigma = +1), propagating across shared edges: in a coherently oriented
// manifold two adjacent faces traverse their shared edge in opposite node
// order. mTags = flat nSeg*nps connectivity (segment s spans
// mTags[s*nps .. +nps)). nps = 3 (tri) or 4 (quad). sigma is meaningful only
// on OK. The global outward sign is applied later by nodalNormals().
inline int propagateOrientation(const int *mTags, int nSeg, int nps, int *sigma) {
    if (mTags == 0 || sigma == 0 || nSeg <= 0 || (nps != 3 && nps != 4))
        return NON_ORIENTABLE;   // degenerate input

    // Undirected edge {lo,hi} -> list of (segment, dir), dir = +1 if the segment
    // traverses lo->hi in winding order, -1 if hi->lo.
    struct EdgeUse { int seg; int dir; };
    std::map<std::pair<int,int>, std::vector<EdgeUse> > edges;
    for (int s = 0; s < nSeg; s++)
        for (int k = 0; k < nps; k++) {
            int a = mTags[s*nps + k];
            int b = mTags[s*nps + (k + 1) % nps];
            if (a == b) return NON_ORIENTABLE;            // collapsed edge
            int lo = a < b ? a : b, hi = a < b ? b : a;
            EdgeUse u; u.seg = s; u.dir = (a < b) ? +1 : -1;
            edges[std::make_pair(lo, hi)].push_back(u);
        }
    // Adjacency: for each interior (2-segment) edge, the relative sign that makes
    // the pair coherent. 1-segment edges are boundary edges of an open patch.
    struct Adj { int nbr; int relSign; };
    std::vector<std::vector<Adj> > adj((size_t)nSeg);
    int nBoundary = 0;
    for (std::map<std::pair<int,int>, std::vector<EdgeUse> >::const_iterator it = edges.begin();
         it != edges.end(); ++it) {
        const std::vector<EdgeUse> &use = it->second;
        if (use.size() > 2) return NON_MANIFOLD;
        if (use.size() == 1) nBoundary++;
        if (use.size() == 2) {
            // coherent means opposite traversal; relSign = -(dir1 dir2):
            //   same direction     => -1 (the neighbour must flip),
            //   opposite direction => +1 (already coherent).
            int rel = -(use[0].dir * use[1].dir);
            Adj a0; a0.nbr = use[1].seg; a0.relSign = rel; adj[use[0].seg].push_back(a0);
            Adj a1; a1.nbr = use[0].seg; a1.relSign = rel; adj[use[1].seg].push_back(a1);
        }
    }
    // Depth-first sweep from segment 0; sigma[nbr] = sigma[cur] relSign. A conflict
    // on a visited segment means non-orientable; unvisited segments afterwards
    // mean a disconnected surface.
    for (int s = 0; s < nSeg; s++) sigma[s] = 0;   // 0 = unvisited
    std::vector<int> stack;
    sigma[0] = +1; stack.push_back(0);
    int visited = 1;
    while (!stack.empty()) {
        int cur = stack.back(); stack.pop_back();
        for (size_t e = 0; e < adj[cur].size(); e++) {
            int nbr = adj[cur][e].nbr;
            int want = sigma[cur] * adj[cur][e].relSign;
            if (sigma[nbr] == 0) { sigma[nbr] = want; visited++; stack.push_back(nbr); }
            else if (sigma[nbr] != want) return NON_ORIENTABLE;   // cycle sign conflict
        }
    }
    if (visited != nSeg) return DISCONNECTED;
    // A coherently wound, connected manifold with no boundary edge is closed
    // (a box or a tube); a single facet is open.
    if (nBoundary == 0 && nSeg > 1) return CLOSED;
    return OK;
}

// ------------------------------------------------ per-segment shared-edge mask
// Fill sharedEdge[nSeg*nps] (row-major [seg][localEdge]) with 1 where local edge
// k of segment s (from local node k to (k+1)%nps, the traversal used by
// propagateOrientation) is shared with another segment, else 0 (a free edge of
// an open patch). Topological, so the caller can cache it with sigma.
//
// Consumer: the facet-ownership guard of
// LadrunoContactProjection::evalSegmentSmooth. At a sharp convex ridge a slave
// pressed into one facet can have its closest point on the adjacent facet land
// on the shared edge with a large spurious penetration; a projection on a
// shared edge with a large gap is rejected (the neighbour owns it), one on a
// free edge is kept.
inline void segmentSharedEdges(const int *mTags, int nSeg, int nps, int *sharedEdge) {
    if (mTags == 0 || sharedEdge == 0 || nSeg <= 0 || (nps != 3 && nps != 4)) return;
    std::map<std::pair<int,int>, int> useCount;   // undirected edge {lo,hi} -> number of uses
    for (int s = 0; s < nSeg; s++)
        for (int k = 0; k < nps; k++) {
            int a = mTags[s*nps + k], b = mTags[s*nps + (k + 1) % nps];
            int lo = a < b ? a : b, hi = a < b ? b : a;
            useCount[std::make_pair(lo, hi)]++;
        }
    for (int s = 0; s < nSeg; s++)
        for (int k = 0; k < nps; k++) {
            int a = mTags[s*nps + k], b = mTags[s*nps + (k + 1) % nps];
            int lo = a < b ? a : b, hi = a < b ? b : a;
            sharedEdge[s*nps + k] = (useCount[std::make_pair(lo, hi)] >= 2) ? 1 : 0;
        }
}

// True if the projection (xi,eta) lies within edgeTol (parametric) of a shared
// segment edge. sharedEdge[k] = 1 for the local edge from shape node k to
// (k+1)%nps, the node ordering of LadrunoContactProjection::shape(). Quad edges:
// k0 eta=-1, k1 xi=+1, k2 eta=+1, k3 xi=-1. Tri edges: k0 eta=0, k1 xi+eta=1,
// k2 xi=0. sharedEdge == 0: no guard (false). edgeTol is calibrated on the quad
// parametric span of 2; the tri span is 1, so the tri band uses edgeTol/2 to
// keep the same fraction of the facet.
inline bool onSharedInteriorEdge(int nps, double xi, double eta,
                                 const int *sharedEdge, double edgeTol) {
    if (sharedEdge == 0) return false;
    if (nps == 4) {
        if (sharedEdge[0] && eta <= -1.0 + edgeTol) return true;   // k0: eta = -1
        if (sharedEdge[1] && xi  >=  1.0 - edgeTol) return true;   // k1: xi = +1
        if (sharedEdge[2] && eta >=  1.0 - edgeTol) return true;   // k2: eta = +1
        if (sharedEdge[3] && xi  <= -1.0 + edgeTol) return true;   // k3: xi = -1
        return false;
    }
    double triTol = 0.5 * edgeTol;                                 // tri span 1 vs quad span 2
    if (sharedEdge[0] && eta <= triTol)            return true;    // k0: eta = 0
    if (sharedEdge[1] && xi + eta >= 1.0 - triTol) return true;    // k1: xi + eta = 1
    if (sharedEdge[2] && xi <= triTol)             return true;    // k2: xi = 0
    return false;
}

// ------------------------------------------------------- area-weighted nodal normals
// Build the nodal-normal field from current coordinates and the cached sigma_s
// with the global outward sign G, and scatter it to segNodalNorm[nSeg*nps*3]
// (row-major [seg][node][xyz]):
//
//   m_a (unnormalized) = sum_{s containing a} sigma_s newellAreaNormal(s)
//   n_a                = G normalize(m_a)
//
// G (+1/-1) is decided once by the caller (voteSign) and passed in frozen, so
// the field cannot flip sign from one configuration to the next; only the
// magnitudes and directions track the deformation. A node whose |m_a| is
// negligible relative to its largest incident facet contribution (a sharp fold
// where coherent facets cancel) is left zero, so a query that touches it falls
// back through smoothNormal()'s degeneracy test. Returns the number of
// degenerate (zero) unique nodal normals.
inline int nodalNormals(const int *mTags, int nSeg, int nps, const double *segCoords,
                        const int *sigma, double G, double *segNodalNorm) {
    // accumulate per unique node tag, tracking the largest single-facet contribution
    std::map<int, std::vector<double> > acc;   // tag -> [x,y,z]
    std::map<int, double> accMax;              // tag -> max |sigma_s newell_s| over incident segments
    for (int s = 0; s < nSeg; s++) {
        double fn[3];
        newellAreaNormal(nps, segCoords + (size_t)s*nps*3, fn);
        double sgn = (double)sigma[s];
        double c[3] = { sgn*fn[0], sgn*fn[1], sgn*fn[2] };   // coherently oriented area normal
        double cmag = nf_norm(c);
        for (int k = 0; k < nps; k++) {
            int tag = mTags[s*nps + k];
            std::vector<double> &v = acc[tag];
            if (v.empty()) v.assign(3, 0.0);
            for (int d = 0; d < 3; d++) v[d] += c[d];
            if (cmag > accMax[tag]) accMax[tag] = cmag;
        }
    }
    // normalize each unique node normal (apply G), count degeneracies
    int nDeg = 0;
    std::map<int, std::vector<double> > unit;   // tag -> unit n_a (G applied); zero if degenerate
    for (std::map<int, std::vector<double> >::iterator it = acc.begin(); it != acc.end(); ++it) {
        double m[3] = { it->second[0], it->second[1], it->second[2] };
        double mn = nf_norm(m);
        std::vector<double> u(3, 0.0);
        // relative degeneracy test: a coherent nodal sum is of the order of its
        // incident facet areas; far below that, the facets cancel (a cusp)
        double ref = accMax[it->first];
        if (mn > 1e-9 * (ref + 1e-300)) {
            for (int d = 0; d < 3; d++) u[d] = G * m[d] / mn;
        } else {
            nDeg++;   // left zero: queries touching this node fall back to the faceted normal
        }
        unit[it->first] = u;
    }
    // scatter to the per-segment-per-node layout
    for (int s = 0; s < nSeg; s++)
        for (int k = 0; k < nps; k++) {
            const std::vector<double> &u = unit[mTags[s*nps + k]];
            double *dst = segNodalNorm + ((size_t)s*nps + k) * 3;
            for (int d = 0; d < 3; d++) dst[d] = u[d];
        }
    return nDeg;
}

// The coherent area-weighted vote vector sum_s sigma_s newellAreaNormal(s), the
// aggregate oriented normal of the surface. The global sign is sign(vote.seed);
// its confidence is |cos(vote, seed)| in [0,1] (near 0: the seed is nearly
// perpendicular to the field and the sign is ill-conditioned, as for a two-sided
// or saddle-shaped master). Returns the sign (+1/-1); fills *conf if non-null.
inline double voteSign(int nSeg, int nps, const double *segCoords, const int *sigma,
                       const double seed[3], double *conf) {
    double vote[3] = {0.0, 0.0, 0.0};
    for (int s = 0; s < nSeg; s++) {
        double fn[3];
        newellAreaNormal(nps, segCoords + (size_t)s*nps*3, fn);
        for (int d = 0; d < 3; d++) vote[d] += (double)sigma[s] * fn[d];
    }
    double vn = nf_norm(vote), sn = nf_norm(seed), vs = nf_dot(vote, seed);
    if (conf != 0) *conf = (vn > 1e-300 && sn > 1e-300) ? (vs >= 0 ? vs : -vs) / (vn * sn) : 0.0;
    return (vs < 0.0) ? -1.0 : +1.0;
}

// ------------------------------------------------------------- smooth normal query
// n = normalize(sum_i N[i] nodalNorm[i]) over a segment's nps nodes. N[] = the
// master shape functions at the projection (LadrunoContactProjection::shape);
// nodalNorm = the segment's nodal normals (one row of nodalNormals()). Returns
// false on a degenerate blend (antiparallel corner normals across a folded
// element, or a node left zero by nodalNormals()); the caller then uses the
// faceted normal for this query.
inline bool smoothNormal(int nps, const double N[4], const double nodalNorm[4][3], double n[3]) {
    double s[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < nps; i++)
        for (int d = 0; d < 3; d++) s[d] += N[i] * nodalNorm[i][d];
    double sn = nf_norm(s);
    if (sn < 1e-12) return false;            // degenerate blend
    for (int d = 0; d < 3; d++) n[d] = s[d] / sn;
    return true;
}

} // namespace LadrunoContactNormalField

#endif
