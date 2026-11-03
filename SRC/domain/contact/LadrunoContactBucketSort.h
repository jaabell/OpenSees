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
//   Hallquist, J.O. (2006). LS-DYNA Theory Manual, Livermore Software
//     Technology Corporation, Section 26 (Contact-Impact Algorithm):
//     bucket sorting, Eqs. (26.37)-(26.44).
//
// Description: LadrunoContactBucketSort, a header-only spatial bucket-sort
// broad phase for node-to-segment contact. It uses no OpenSees types (raw
// double arrays and the STL only).
//
// A grid is built over the master segments once; for each slave node only the
// segments registered in its 27 neighbour buckets are returned, instead of
// all nSeg segments (brute force, O(nSlave nSeg)). The narrow phase
// (projection, gap and penalty in LadrunoContactFE) is unchanged, so the
// result is identical for every kept pair. Each segment is registered in
// every bucket its corner nodes span, so the kept set is a superset of every
// geometrically near pair and no active pair is lost.
//
// Scope: the grid is built in ConstraintHandler::handle() (a domainChanged
// event) from reference coordinates. For static and small-motion contact the
// candidate set covers every pair that activates, and the result is
// bitwise-identical to brute force. Large sliding (a slave that migrates
// farther than the search band from its cell) needs the broad phase rebuilt;
// a new analysis after large motion rebuilds it.

#ifndef LadrunoContactBucketSort_h
#define LadrunoContactBucketSort_h

#include <cmath>
#include <vector>
#include <algorithm>

namespace LadrunoContactBucketSort {

class Grid {
  public:
    // segCoords: flat nSeg*nps*3 (master segment node REFERENCE coords, row-major
    //   [seg][node][xyz]). lmaxFrac: cell = lmaxFrac * median segment diagonal
    //   (1.0 default; 1.5 for distorted meshes). clipPct: runaway-node guard; 0
    //   disables it, else clamp the bbox to [clipPct, 100-clipPct] percentiles of
    //   the segment centroids (e.g. 1.0 -> [1%,99%]) so one fly-off node can't
    //   collapse all real segments into a single cell.
    Grid(int nSeg, int nps, const double *segCoords,
         double lmaxFrac = 1.0, double clipPct = 1.0)
      : nSeg_(nSeg), nps_(nps), NX_(1), NY_(1), NZ_(1),
        xmin_(0), xmax_(0), ymin_(0), ymax_(0), zmin_(0), zmax_(0),
        lmax_(1.0), guardFired_(false), stampTick_(0)
    {
        if (nSeg_ <= 0 || nps_ <= 0 || segCoords == 0) return;

        // --- segment centroids + diagonals ---
        std::vector<double> cx(nSeg_), cy(nSeg_), cz(nSeg_);
        std::vector<double> diags(nSeg_);
        for (int s = 0; s < nSeg_; s++) {
            const double *P = segCoords + (size_t)s * nps_ * 3;
            double sx = 0, sy = 0, sz = 0;
            for (int i = 0; i < nps_; i++) { sx += P[i*3+0]; sy += P[i*3+1]; sz += P[i*3+2]; }
            cx[s] = sx / nps_; cy[s] = sy / nps_; cz[s] = sz / nps_;
            double dmax = 0.0;
            for (int i = 0; i < nps_; i++)
                for (int j = i+1; j < nps_; j++) {
                    double dx = P[i*3+0]-P[j*3+0], dy = P[i*3+1]-P[j*3+1], dz = P[i*3+2]-P[j*3+2];
                    double d = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (d > dmax) dmax = d;
                }
            diags[s] = dmax;
        }
        std::vector<double> sd(diags);
        std::sort(sd.begin(), sd.end());
        double medDiag = sd[nSeg_ / 2];
        lmax_ = std::max(1e-12, lmaxFrac * medDiag);   // degenerate-segment floor

        // --- bbox of centroids, with the optional runaway percentile clip ---
        clip(cx, clipPct, xmin_, xmax_);
        clip(cy, clipPct, ymin_, ymax_);
        clip(cz, clipPct, zmin_, zmax_);

        // --- grid resolution, capped NX*NY*NZ <= min(nSeg, 5000) (Eq. 26.37) ---
        // Per-axis pre-clamp to 5000.0 in double before the int cast: the total is capped
        // at min(nSeg,5000) cells anyway, and with clipPct=0 a diverging node can make
        // (range/lmax_) exceed the int range; a double to int overflow is undefined.
        // On LLP64 Windows `long` is 32-bit, so the product test uses long long to stay
        // exact (a wrapped product would request a huge allocation in grid_.assign()).
        NX_ = std::max(1, (int)std::min(5000.0, (xmax_ - xmin_) / lmax_) + 1);
        NY_ = std::max(1, (int)std::min(5000.0, (ymax_ - ymin_) / lmax_) + 1);
        NZ_ = std::max(1, (int)std::min(5000.0, (zmax_ - zmin_) / lmax_) + 1);
        long long cap = std::min((long long)nSeg_, (long long)5000);
        while ((long long)NX_ * NY_ * NZ_ > cap) {
            int m = std::max(NX_, std::max(NY_, NZ_));
            if (NX_ == m && NX_ > 1) NX_--;
            else if (NY_ == m && NY_ > 1) NY_--;
            else if (NZ_ > 1) NZ_--;
            else break;
        }

        // --- segment-span registration (LS-DYNA type-13, Eq. 26.44): each segment in
        //     every bucket its corner nodes span (not centroid-only; type-4 is
        //     deficient, misses true-near pairs for large/skewed segments). ---
        grid_.assign((size_t)NX_ * NY_ * NZ_, std::vector<int>());
        for (int s = 0; s < nSeg_; s++) {
            const double *P = segCoords + (size_t)s * nps_ * 3;
            int imin = NX_, imax = 1, jmin = NY_, jmax = 1, kmin = NZ_, kmax = 1;
            for (int i = 0; i < nps_; i++) {
                int px = clampP(P[i*3+0], xmin_, xmax_, NX_);
                int py = clampP(P[i*3+1], ymin_, ymax_, NY_);
                int pz = clampP(P[i*3+2], zmin_, zmax_, NZ_);
                imin = std::min(imin, px); imax = std::max(imax, px);
                jmin = std::min(jmin, py); jmax = std::max(jmax, py);
                kmin = std::min(kmin, pz); kmax = std::max(kmax, pz);
            }
            for (int ix = imin; ix <= imax; ix++)
                for (int iy = jmin; iy <= jmax; iy++)
                    for (int iz = kmin; iz <= kmax; iz++)
                        grid_[linId(ix, iy, iz)].push_back(s);
        }
        stamp_.assign(nSeg_, 0);
        stampTick_ = 0;
    }

    // Fill out[] (caller-sized >= nSeg) with the de-duplicated candidate segment
    // indices for slavePt (the union over its 27 neighbour buckets, Eq. 26.42); returns
    // the count. De-dup via a sparse-set stamp (no per-call alloc/clear).
    int candidates(const double slavePt[3], int *out) const
    {
        if (nSeg_ <= 0 || out == 0) return 0;
        int PX = clampP(slavePt[0], xmin_, xmax_, NX_);
        int PY = clampP(slavePt[1], ymin_, ymax_, NY_);
        int PZ = clampP(slavePt[2], zmin_, zmax_, NZ_);
        ++stampTick_;
        int n = 0;
        for (int ix = std::max(1, PX-1); ix <= std::min(NX_, PX+1); ix++)
            for (int iy = std::max(1, PY-1); iy <= std::min(NY_, PY+1); iy++)
                for (int iz = std::max(1, PZ-1); iz <= std::min(NZ_, PZ+1); iz++) {
                    const std::vector<int> &b = grid_[linId(ix, iy, iz)];
                    for (size_t t = 0; t < b.size(); t++) {
                        int sid = b[t];
                        if (stamp_[sid] != stampTick_) { stamp_[sid] = stampTick_; out[n++] = sid; }
                    }
                }
        return n;
    }

    int numSegments(void) const { return nSeg_; }
    bool runawayGuardFired(void) const { return guardFired_; }

  private:
    int nSeg_, nps_, NX_, NY_, NZ_;
    double xmin_, xmax_, ymin_, ymax_, zmin_, zmax_, lmax_;
    bool guardFired_;
    std::vector<std::vector<int>> grid_;
    mutable std::vector<long long> stamp_;  // sparse-set de-dup (cache, hence mutable)
    mutable long long stampTick_;           // wide so the per-query ++ never wraps

    int linId(int px, int py, int pz) const { return (px-1) + (py-1)*NX_ + (pz-1)*NX_*NY_; }

    // bucket index in [1,N] for coordinate v in [lo,hi]; clamps out-of-range to the
    // edge cell (a slave beyond the master bbox lands in the nearest boundary cell,
    // and the +-1 neighbour search still sees the boundary segments).
    static int clampP(double v, double lo, double hi, int N)
    {
        if (hi <= lo) return 1;
        int p = (int)(N * (v - lo) / (hi - lo)) + 1;
        return std::min(std::max(p, 1), N);
    }

    // bbox of vals; if pct>0 clamp to [pct, 100-pct] percentiles (runaway guard) and
    // set guardFired_ when the clamp actually moved a bound.
    void clip(std::vector<double> &vals, double pct, double &mn, double &mx)
    {
        double rawMin = *std::min_element(vals.begin(), vals.end());
        double rawMax = *std::max_element(vals.begin(), vals.end());
        if (pct <= 0.0) { mn = rawMin; mx = rawMax; return; }
        std::vector<double> s(vals);
        std::sort(s.begin(), s.end());
        int n = (int)s.size();
        int loI = std::min(n-1, std::max(0, (int)(n * pct / 100.0)));
        int hiI = std::min(n-1, std::max(0, (int)(n * (100.0 - pct) / 100.0)));
        mn = s[loI]; mx = s[hiI];
        if (mn > rawMin || mx < rawMax) guardFired_ = true;
    }
};

} // namespace LadrunoContactBucketSort

#endif
