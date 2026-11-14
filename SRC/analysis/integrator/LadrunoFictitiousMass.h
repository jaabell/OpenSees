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
//   Papadrakakis, M. (1981). A method for the automatic evaluation of the
//     dynamic relaxation parameters. Computer Methods in Applied Mechanics
//     and Engineering, 25(1), 35-48. doi:10.1016/0045-7825(81)90066-9
//   Underwood, P. (1983). Dynamic relaxation. In T. Belytschko and
//     T.J.R. Hughes (eds.), Computational Methods for Transient Analysis,
//     North-Holland, Amsterdam, pp. 245-265.
//
// Description: integrator-owned diagonal fictitious mass shared by the
// static integrators that regularize or replace the tangent with a diagonal
// matrix (LadrunoArcLength -stabilize and LadrunoDynamicRelaxation).
//
// The diagonal M* is built from an element probe, not from
// Element::getMass(): getMass() is zero on models without density, which are
// the usual case for static analysis. Three modes:
//   0 gershgorin  m_i = f * sum_j |K_ij|, the row sums of the assembled
//                 element tangent stiffness. By the Gershgorin circle
//                 theorem the largest eigenvalue of M*^-1 K is at most
//                 max_i (sum_j |K_ij|) / m_i = 1/f, so f = dt^2/4 bounds
//                 omega_max * dt <= 2 (Papadrakakis 1981; Underwood 1983).
//                 f <= 0 means no prefactor (f = 1).
//   1 lumped      m_i = massScale * sum_j |M_ij| (row sums of the element
//                 mass; needs a nonzero density).
//   2 unity       m_i = 1.
// Any zero or negative diagonal (a free DOF with no stiffness or mass) is
// raised to 1e-8 times the largest diagonal so that M* stays invertible.

#ifndef LadrunoFictitiousMass_h
#define LadrunoFictitiousMass_h

#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Element.h>
#include <math.h>

namespace Ladruno {

// Fill Mstar (already sized to the number of equations) with the diagonal
// fictitious mass of the given mode. Returns 0 on success, -1 without a model.
inline int
buildGershgorinDiagonal(AnalysisModel *theModel, Vector &Mstar, int mode,
                        double massScale, double dt2Quarter)
{
  if (theModel == 0) return -1;
  int size = Mstar.Size();
  Mstar.Zero();

  if (mode == 2) {                       // unity
    for (int i = 0; i < size; i++) Mstar(i) = 1.0;
    return 0;
  }

  // gershgorin or lumped: accumulate the absolute row sums of the element K
  // (or M) at the element equation numbers
  FE_EleIter &theEles = theModel->getFEs();
  FE_Element *elePtr;
  while ((elePtr = theEles()) != 0) {
    Element *e = elePtr->getElement();
    if (e == 0) continue;                // constraint FE without an element
    const ID &id = elePtr->getID();
    const Matrix &Ke = (mode == 1) ? e->getMass() : e->getTangentStiff();
    int n = id.Size();
    if (Ke.noRows() != n) continue;
    for (int i = 0; i < n; i++) {
      int loc = id(i);
      if (loc < 0) continue;
      double rs = 0.0;
      for (int j = 0; j < n; j++) rs += fabs(Ke(i, j));
      Mstar(loc) += rs;
    }
  }

  if (mode == 0) {
    double f = (dt2Quarter > 0.0) ? dt2Quarter : 1.0;
    if (f != 1.0)
      for (int i = 0; i < size; i++) Mstar(i) *= f;
  } else {
    for (int i = 0; i < size; i++) Mstar(i) *= massScale;
  }

  // floor any zero or negative diagonal
  double mmax = 0.0;
  for (int i = 0; i < size; i++) if (Mstar(i) > mmax) mmax = Mstar(i);
  double floor = (mmax > 0.0) ? mmax * 1.0e-8 : 1.0;
  for (int i = 0; i < size; i++)
    if (!(Mstar(i) > 0.0)) Mstar(i) = floor;

  return 0;
}

// Add fac*M* to the diagonal of the SOE matrix. With zeroFirst the matrix is
// zeroed before, so A = fac*M* and a solve is a diagonal inverse (dynamic
// relaxation); without it fac*M* regularizes the tangent already assembled
// (stabilized static step).
inline void
addDiagonalToSOE(LinearSOE *theSOE, const Vector &Mstar, int size, double fac,
                 bool zeroFirst)
{
  if (theSOE == 0) return;
  if (zeroFirst) theSOE->zeroA();
  static Matrix m1(1, 1);
  static ID id1(1);
  for (int i = 0; i < size; i++) {
    m1(0, 0) = fac * Mstar(i);
    id1(0) = i;
    theSOE->addA(m1, id1);
  }
}

} // namespace Ladruno

#endif
