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
//   Hinton, E., Rock, T. and Zienkiewicz, O.C. (1976). A note on mass
//     lumping and related processes in the finite element method.
//     Earthquake Engineering and Structural Dynamics, 4(3), 245-249.
//     doi:10.1002/eqe.4290040305
//
// Description: mass-conserving diagonal lumping of a consistent element mass
// matrix (Hinton, Rock and Zienkiewicz 1976, "HRZ"). Per translational
// direction d the consistent diagonal entries are scaled so that they sum to
// the element mass in that direction:
//
//   S_d = sum_{a in d} Mc_aa ,   m_d = sum_{a in d} sum_j Mc_aj ,
//   alpha_d = m_d / S_d ,        M_aa = alpha_d Mc_aa .
//
// m_d (the sum of the full row sums over the direction-d DOFs) equals the
// element mass in direction d when the consistent mass has no net coupling
// into direction d; this holds for global-frame translational blocks and for
// symmetric beam masses, where the transverse-rotational coupling terms cancel
// in the sum. The lump is positive (Mc_aa > 0 for a properly integrated
// consistent mass), conserves the mass per direction, and gives rotational
// DOFs a nonzero inertia: a rotational DOF is scaled by the mean of the
// element's translational factors alpha_d, which stays defined when the
// factors differ by direction (consistent beam mass: 140 axial, 156
// transverse).
//
// On geometrically regular elements (square quadrilateral, parallelepiped
// hexahedron, bar) all Mc_aa of a direction are equal and the lump equals the
// row sum. On distorted and higher-order elements the two differ per node,
// and the row sum can produce zero or negative translational entries and zero
// rotational inertia, while the HRZ lump stays positive and conserves mass.

#ifndef LadrunoMassLumping_h
#define LadrunoMassLumping_h

#include <Matrix.h>
#include <ID.h>

namespace Ladruno {

// Return codes of the HRZ lumper.
enum HRZStatus {
  HRZ_OK            = 0,  // HRZ lump applied
  HRZ_FELLBACK_DIAG = 1,  // a guard failed: mdiag holds the consistent diagonal
                          //   (not mass-conserving; the caller should warn)
  HRZ_BAD_INPUT     = -1  // n <= 0 or null buffers
};

// Core kernel. Mc is the n x n consistent mass packed row-major
// (Mc[i*n + j]); dir[i] is the spatial direction 0..ndm-1 of a translational
// DOF, or -1 for a rotational (or other non-translational) DOF. Fills
// mdiag[i] with the HRZ diagonal.
//
// Guards (any failure falls back to the consistent diagonal):
//   - every translational Mc_aa > 0,
//   - every per-direction S_d > 0 and m_d > 0,
//   - every translational M_aa > 0.
inline int hrzLumpRaw(const double *Mc, const int *dir, int n, double *mdiag)
{
  if (n <= 0 || Mc == 0 || dir == 0 || mdiag == 0)
    return HRZ_BAD_INPUT;

  // ndm = largest translational direction index + 1 (at most 3)
  int ndm = 0;
  for (int i = 0; i < n; ++i) {
    int d = dir[i];
    if (d >= 0 && d + 1 > ndm) ndm = d + 1;
  }
  if (ndm > 3) ndm = 3;

  // per-direction S_d (sum of the consistent diagonal) and m_d (sum of the
  // full row sums, the element mass in that direction)
  double S[3]  = {0.0, 0.0, 0.0};
  double Mt[3] = {0.0, 0.0, 0.0};
  bool ok = true;

  for (int i = 0; i < n && ok; ++i) {
    int d = dir[i];
    if (d < 0 || d >= ndm) continue;          // rotational DOFs handled below
    double mii = Mc[i * n + i];
    if (!(mii > 0.0)) { ok = false; break; }
    S[d] += mii;
    double rs = 0.0;
    for (int j = 0; j < n; ++j) rs += Mc[i * n + j];
    Mt[d] += rs;
  }

  if (ok)
    for (int d = 0; d < ndm; ++d)
      if (!(S[d] > 0.0) || !(Mt[d] > 0.0)) { ok = false; break; }

  if (ok) {
    double alpha[3] = {0.0, 0.0, 0.0};
    double alphaSum = 0.0;
    for (int d = 0; d < ndm; ++d) { alpha[d] = Mt[d] / S[d]; alphaSum += alpha[d]; }
    double alphaMean = (ndm > 0) ? alphaSum / ndm : 0.0;

    for (int i = 0; i < n; ++i) {
      int d = dir[i];
      double mii = Mc[i * n + i];
      if (d >= 0 && d < ndm) {
        mdiag[i] = alpha[d] * mii;
        if (!(mdiag[i] > 0.0)) { ok = false; break; }
      } else {
        // Non-translational DOF. A positive inertia is scaled by the mean
        // translational factor. A zero or negative entry (shell drilling DOF
        // without rotational mass, a pressure DOF with a compressibility
        // term) is passed through unchanged instead of failing the element,
        // as the consistent-diagonal lump would do.
        mdiag[i] = (mii > 0.0) ? alphaMean * mii : mii;
      }
    }
  }

  if (!ok) {
    for (int i = 0; i < n; ++i) mdiag[i] = Mc[i * n + i];
    return HRZ_FELLBACK_DIAG;
  }
  return HRZ_OK;
}

// OpenSees wrapper: copies the consistent mass M into a row-major buffer and
// calls the kernel. dofDir(i) is the direction 0..ndm-1 of a translational
// DOF or -1 otherwise.
inline int hrzLump(const Matrix &M, const ID &dofDir, double *mdiag, int n)
{
  if (n <= 0 || M.noRows() != n || M.noCols() != n || dofDir.Size() != n)
    return HRZ_BAD_INPUT;
  double *Mc  = new double[n * n];
  int    *dir = new int[n];
  for (int i = 0; i < n; ++i) {
    dir[i] = dofDir(i);
    for (int j = 0; j < n; ++j) Mc[i * n + j] = M(i, j);
  }
  int rc = hrzLumpRaw(Mc, dir, n, mdiag);
  delete[] Mc;
  delete[] dir;
  return rc;
}

} // namespace Ladruno

#endif
