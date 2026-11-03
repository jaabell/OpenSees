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
//   Puso, M.A. (2004). A 3D mortar method for solid mechanics. International
//     Journal for Numerical Methods in Engineering, 59(3), 315-336.
//     doi:10.1002/nme.865
//   Wohlmuth, B.I. (2000). A mortar finite element method using dual spaces
//     for the Lagrange multiplier. SIAM Journal on Numerical Analysis, 38(3),
//     989-1012. doi:10.1137/S0036142999350929
//   Popp, A., Wohlmuth, B.I., Gee, M.W. and Wall, W.A. (2012). Dual quadratic
//     mortar finite element methods for 3D finite deformation contact. SIAM
//     Journal on Scientific Computing, 34(4), B421-B446.
//     doi:10.1137/110848190
//   Dassault Systemes Simulia. Abaqus Analysis User's Guide, "Mesh tie
//     constraints" (*TIE, node-to-surface) and "Shell-to-solid coupling".
//   Hinton, E., Rock, T. and Zienkiewicz, O.C. (1976). A note on mass lumping
//     and related processes in the finite element method. Earthquake
//     Engineering and Structural Dynamics, 4(3), 245-249.
//     doi:10.1002/eqe.4290040305
//
// Description: LadrunoTie, a setup-time generator of surface ties between
// non-conforming meshes. A tie is a permanent, non-sliding bond, so the
// pairing is computed once, from the reference coordinates, and written as
// ordinary EQ_Constraints, one per slave node per tied DOF:
//
//        u_s = sum_k P_sk u_k            (1*u_s - sum_k P_sk u_k = 0)
//
// The generator writes no enforcement code: constraints Lagrange or Penalty
// enforce the rows (the Transformation and Auto handlers do not apply
// EQ_Constraints). Four transfers P are available:
//
//   collocation (default)  each slave node is projected onto its nearest
//                          master facet (bucket-sort broad phase, closest-
//                          point projection); P_sk are the master shape
//                          functions at the projection.
//   -hermite               ndf-6 shell edge ties: the shell-normal
//                          translation and the slope along the interface are
//                          reconstructed with cubic Hermite functions of the
//                          master edge (w, theta) pairs.
//   -mortar                integral (weak) tie of the slave surface:
//                          P = D^-1 M with the mortar matrices of Puso (2004),
//                          optionally with the dual basis of Wohlmuth (2000)
//                          (-dual), which makes D diagonal and P local.
//   -shellSolid            an ndf-6 shell edge (master) tied to an ndf-3 solid
//                          face (slaves) by the linearized rigid plane-section
//                          map u_s = sum_j N_j (u_j + theta_j x d).
//
// Every row ties a slave DOF to master DOFs only (no chained constraints).

#ifndef LadrunoTie_h
#define LadrunoTie_h

class Domain;
class ID;

namespace LadrunoTie {

// Collocation tie. Pair each slave node with its nearest master facet,
// evaluate the master shape functions N_i at the projection, and emit one
// EQ_Constraint per slave node per tied DOF (u_s = sum_i N_i u_{m,i}).
//
//   slaveNodes        slave node tags (each tied to ONE facet).
//   masterFacetNodes  flat node tags, nps per facet, row-major [facet][node].
//   nps               nodes per master facet (3 tri-3, 4 quad-4).
//   dofs              1-based DOFs to tie. Empty: 1..ndm, plus the rotations
//                     4,5,6 for a 3D ndf-6 node.
//   tolFrac           conforming tolerance: a slave that projects farther than
//                     tolFrac * facet size off the master surface is refused
//                     (slave nodes are not moved onto the master surface).
//   hermite           rotation-consistent cubic Hermite w-theta transfer for
//                     ndf-6 shell EDGE (butt-joint) ties. The linear transfer
//                     leaves an O(h^2) error on the transverse displacement w
//                     when the curvature varies ALONG the interface (w
//                     quadratic along the tie line). With hermite the
//                     shell-normal translation and the slope along the
//                     interface of each slave node are a cubic Hermite
//                     combination of the master edge nodes' (w, theta_t)
//                     pairs: exact through cubic w. In-plane translations and
//                     the other rotations keep the linear transfer. Requires
//                     ndf-6 nodes on both sides, the full 1..6 DOF set, and
//                     each slave on a master facet EDGE (or corner). Assumes
//                     Kirchhoff kinematics on the tie line (slope = rotation):
//                     a shear-deformable state carries an O(gamma h) tie
//                     error, gamma being the transverse shear angle.
//
// Returns the number of EQ_Constraints emitted (>= 0), or -1 on a refusal:
// slave and master surfaces not node-disjoint, a slave listed twice, a slave
// off the master surface, a Hermite precondition not met, a missing node or
// DOF, or (only under a handler that requires it) a tied DOF without positive
// lumped mass.
int generate(Domain *theDomain,
             const ID &slaveNodes,
             const ID &masterFacetNodes, int nps,
             const ID &dofs, double tolFrac,
             bool hermite = false);

// Integral mortar tie. Ties the slave SURFACE to the master surface in the
// weak sense: the mortar matrices over the clipped overlap
//        D_IJ = INT N_I^s N_J^s dGamma   (slave-slave)
//        M_IK = INT N_I^s phi_K^m dGamma (slave-master)
// are assembled with LadrunoMortarKernel::integratePair, condensed once to
//        u_s = P u_m ,  P = D^-1 M
// and emitted as one EQ_Constraint per slave node per tied DOF. Inverting D
// globally makes every row retain MASTER nodes only.
//
//   slaveFacetNodes   flat tags, npsS per slave facet (3, 4 or 8).
//   masterFacetNodes  flat tags, npsM per master facet (3, 4, 6 or 8).
//                     Quadratic facets: corners first, flat, midside nodes at
//                     the edge midpoints. tri-6 slave facets are refused: the
//                     integrals of their corner shape functions are zero
//                     (Popp et al. 2012), which makes the dual basis and the
//                     coverage measures degenerate.
//   outward           direction that orients the facet normals; if 0, the
//                     average master facet normal is used.
//   dual              dual (biorthogonal) slave basis (Wohlmuth 2000): D is
//                     diagonal, so P is local (each slave retains only the
//                     masters under its own facets). Built per slave facet as
//                     A^e = diag(row sums of D^e) (D^e)^-1. Linear
//                     completeness is preserved (row-sum lumping of D would
//                     lose it).
//
// A (slave facet, master facet) pair whose mean separation exceeds the
// conforming tolerance is excluded from the integrals: on a closed or folded
// surface the far side projects onto the auxiliary plane of a near facet,
// and including it would tie the slave to the opposite side of the surface.
//
// Returns the number of EQ_Constraints emitted (>= 0), or -1 on a refusal:
// surfaces not node-disjoint, a slave facet not covered by the master
// surface, a slave facet off the master surface, an overlapping master
// surface, an invalid quadratic facet, an ill-conditioned D, or (only under a
// handler that requires it) a tied DOF without positive lumped mass.
int generateMortar(Domain *theDomain,
                   const ID &slaveFacetNodes, int npsS,
                   const ID &masterFacetNodes, int npsM,
                   const ID &dofs, double tolFrac, const double *outward,
                   bool dual = false);

// Shell-to-solid plane-section tie. Ties an ndf-6 shell EDGE (the master, a
// polyline of segment node pairs) to an ndf-3 solid FACE (the slaves). Each
// slave is projected onto the closest point of the edge polyline, with a
// FROZEN arm d = x_s - xbar(xi), and three translation rows are emitted:
//
//        u_{s,i} = sum_{j=A,B} N_j(xi) [ u_{j,i} + (theta_j x d)_i ]   i = 1..3
//
// the linearized rigid plane-section map (the kinematics of a shell-to-solid
// coupling). The arm runs along the shell normal, so the drilling (theta.n)
// coefficient of theta x d is ~zero and is dropped by the 1e-12 coefficient
// filter: the shell drilling DOF stays free. The tied DOFs are the solid
// translations {1,2,3}. Small rotations only (frozen arm); the
// through-thickness Poisson stretch is suppressed at the seam (a boundary
// layer proportional to nu*t*strain, exact at nu = 0); a shear-deformable
// shell leaves an O(gamma t) local misfit while the resultant transfer is
// exact. Collocation only.
//
//   slaveNodes       solid FACE node tags.
//   masterEdgeNodes  flat shell-edge segment node PAIRS, 2 per segment
//                    (a1 b1 a2 b2 ..); kinks allowed (per-node frozen arms).
//   thickness        shell thickness t for the arm check
//                    |d| <= (0.5 + slack) t; <= 0 skips the check (with a
//                    note).
//   tolFrac          slack for the arm check (at least 0.02).
//
// Returns the number of EQ_Constraints emitted (>= 0), or -1 on a refusal:
// duplicate slave, slave on the master edge, master edge node with ndf < 6,
// zero-length segment, slave outside the polyline, the arm check, or (only
// under a handler that requires it) a tied DOF without positive lumped mass.
int generateShellSolid(Domain *theDomain,
                       const ID &slaveNodes,
                       const ID &masterEdgeNodes,
                       double thickness, double tolFrac);

} // namespace LadrunoTie

#endif
