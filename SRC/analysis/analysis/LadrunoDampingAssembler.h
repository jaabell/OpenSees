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
//   Veletsos, A.S. and Ventura, C.E. (1986). Modal analysis of
//     non-classically damped linear systems. Earthquake Engineering and
//     Structural Dynamics, 14(2), 217-243. doi:10.1002/eqe.4290140205
//
// Description: reduced-operator projection for complex modal analysis
// (LadrunoComplexEigen).
//
//   Projects the model's damping and mass operators onto the real undamped
//   modal basis Phi of the last `eigen` run without assembling a global
//   matrix: accumulate element by element
//
//       Ct += Phi_e^T C_e Phi_e     (C_e = Element::getDamp(): Rayleigh
//                                    element factors and damper materials,
//                                    the C of the transient analysis)
//       Mt += Phi_e^T M_e Phi_e     (M_e = Element::getMass())
//
//   plus the nodal contributions (lumped nodal mass and its alphaM Rayleigh
//   term via Node::getDamp()/getMass()). Phi_e is gathered from the nodal
//   eigenvector storage (Node::getEigenvectors) in the element's node/DOF
//   ordering: no DOF numbering, no system of equations, no integrator.
//
//   The stiffness projection is not assembled: K phi_a = w_a^2 M phi_a gives
//   Kt = Mt diag(w^2) exactly (symmetrized to remove round-off), for any
//   basis normalization and for repeated eigenvalues (vectors of distinct
//   eigenvalues are M-orthogonal; within a repeated eigenspace Lambda is a
//   multiple of I, so Mt and Lambda commute). Passing the full projected Mt
//   (not an assumed identity) to the pencil solver is what keeps the
//   repeated-eigenvalue case exact without re-orthonormalization.
//
//   For a model with global Rayleigh damping only, the element and node
//   getDamp() sums reproduce the closed form Ct = diag(alphaM + betaK w^2);
//   region or element Rayleigh factors and betaKinit/betaKcomm terms are
//   captured exactly because getDamp() already carries them per element.
//
//   Contract:
//   - C is what getDamp() returns. Damping that does not flow through
//     getDamp() is invisible: modalDamping (applied by the integrator; the
//     driver warns), element 'damping' objects (folded into
//     getResistingForceIncInertia), HHT-alpha numerical damping, and any
//     element whose -doRayleigh flag is left at its default 0 (Truss and
//     the zeroLength family default to off).
//   - Mt and Ct are read from the current model while Lambda is the
//     spectrum stored by the last eigen run. Editing the model after eigen
//     (mass, stiffness, state) makes them inconsistent, and because
//     Kt = Mt Lambda reuses the stored spectrum, the reported omega0 still
//     equals the old frequencies. Re-run eigen after any model change (the
//     same convention as modalProperties).
//   - Eigen must have run under a constraint handler that distributes
//     eigenvectors to constrained DOFs (Transformation does; 'constraints
//     Plain' with MP constraints writes zero rows at the retained-node
//     dependents and would misproject).

#ifndef LadrunoDampingAssembler_h
#define LadrunoDampingAssembler_h

class Domain;
class Matrix;

class LadrunoDampingAssembler
{
  public:
    // Project the model's M and C onto the first p modes of the last eigen
    // run and synthesize Kt from the eigenvalues (see the note above). Mt,
    // Ct, Kt must be p x p (they are zeroed here). Returns 0 on success,
    // < 0 on error (no or short spectrum, a node with fewer stored modes
    // than p, or an element whose damping/mass dimension does not match its
    // nodes' DOFs).
    static int projectReduced(Domain &theDomain, int p,
                              Matrix &Mt, Matrix &Ct, Matrix &Kt);
};

#endif
