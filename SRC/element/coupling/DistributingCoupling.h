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

// Written: N. Mora Bowen, P. Palacios, J.A. Abell
// Created: 2026
//
// Description: DistributingCoupling, a distributing (interpolation)
// coupling (RBE3-type) enforced by penalty or augmented Lagrangian.
//
// The motion of a reference node R (translations and rotations) is the
// weighted least-squares rigid-body fit of the translations of N
// independent nodes I_i with weights w_i. Equivalently, a force or moment
// applied at R is distributed to the independent nodes as a statically
// equivalent force set, and the coupling adds no stiffness to the
// independent nodes against deformation of the set. A typical use is the
// connection of a 6-DOF beam or shell node to a face of 3-DOF solid nodes.
//
// With W = sum w_i, the weighted centroid x_c = sum w_i x_i / W,
// r_i = x_i - x_c and the weighted position inertia
//     I_c = sum w_i [ (r_i . r_i) 1 - r_i (x) r_i ],
// the constraint is
//
//     g_t = u_R + theta_R x (x_c - x_R) - (1/W) sum w_i u_i      (ndm rows)
//     g_r = P theta_R - I_c^+ sum w_i (r_i x u_i)                 (nrot rows)
//
// I_c^+ is the spectral pseudo-inverse of I_c and P the projector on its
// range: rotation axes about which the independent set has no spread
// (collinear or coincident nodes) are left unconstrained instead of
// producing a singular operator. The transport term theta_R x (x_c - x_R)
// makes the constraint exact for a reference node offset from the
// centroid. g = B u with constant B, and
//
//     t = D g (+ lambda),   P = B^T t,   K = B^T D B,
//     D = diag(K_t I, K_r I).
//
// The force on independent node i is (w_i/W) t_t + w_i (alpha x r_i) with
// alpha = I_c^+ t_r, the transpose of the kinematic map, so the net force
// and moment of the distributed set equal the reference load for any
// penalty value. With -enforce al the multipliers are updated once per
// committed step, lambda <- lambda + D g. The default rotational penalty is
// K_r = K_t * sum w_i |r_i|^2 / W.
//
// The same kinematics are the Nastran RBE3 element and the Abaqus
// *COUPLING, *DISTRIBUTING constraint:
//   MSC Nastran Quick Reference Guide, bulk data entry RBE3.
//   Abaqus Analysis User's Guide, section "Coupling constraints"
//     (distributing coupling).
// Penalty and augmented Lagrangian enforcement of linear constraints:
//   P. Wriggers, Computational Contact Mechanics, 2nd ed., Springer, 2006,
//     chapter 6.
//   J.C. Simo, T.A. Laursen, An augmented Lagrangian treatment of contact
//     problems involving friction, Computers & Structures 42 (1992) 97-116.
//
// Usage:
//   element DistributingCoupling tag refNode N i1 ... iN
//       <-w w1 ... wN> <-k Kt | -k auto -host eleTag <-kAlpha a>>
//       <-kr Kr> <-enforce penalty|al> <-absolute>
//
// Constraints on use: small rotations (linear kinematics, B is constant);
// the reference node needs ndm+nrot DOFs (6 in 3D, 3 in 2D); the coupling
// carries no mass and no Rayleigh damping.

#ifndef DistributingCoupling_h
#define DistributingCoupling_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

class Node;
class Channel;
class FEM_ObjectBroker;
class Response;

class DistributingCoupling : public Element
{
 public:
  DistributingCoupling(int tag, int ndm, int refNode, const ID &indepNodes,
                       const Vector &weights, double kt, double kr, bool krUser,
                       int enforce, double kAlpha, int hostEleTag, bool ktAuto,
                       bool initGapCapture = true);
  DistributingCoupling();
  ~DistributingCoupling();

  const char *getClassType(void) const { return "DistributingCoupling"; }

  // domain
  int getNumExternalNodes(void) const;
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  int getNumDOF(void);
  void setDomain(Domain *theDomain);

  // state
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  // matrices and forces
  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getMass(void);
  const Matrix &getDamp(void);
  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);

  // Rayleigh factors are not applied to a constraint element: they are
  // stored as zero, so the element contributes no damping.
  int setRayleighDampingFactors(double alphaM, double betaK,
                                double betaK0, double betaKc);

  // parallel and database
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInfo);

 private:
  int ndm;                  // 2 or 3
  int nrot;                 // rotation components: 3 (3D) or 1 (2D)
  int nIndep;               // number of independent nodes N
  ID connectedNodes;        // [ refNode, indep_1 .. indep_N ]
  Vector weights;           // w_i

  double Kt;                // translational penalty
  double Kr;                // rotational penalty
  bool krUser;              // -kr given (else K_r = K_t * l^2)
  bool ktAuto;              // -k auto: K_t from the -host element stiffness
  double kAlpha;            // multiplier for -k auto
  int hostEleTag;           // representative host element, -1 if none
  bool ktResolved;          // K_t / K_r resolved (transient)
  bool hostMissWarned;      // missing -host reported (transient)

  int enforce;              // 0 = penalty, 1 = augmented Lagrangian
  Vector lambda;            // translational multipliers (ndm)
  Vector lambda_r;          // rotational multipliers (nrot)

  // geometry, resolved at setDomain
  bool valid;
  double W;                 // sum w_i
  double ell2;              // sum w_i |r_i|^2 / W
  Vector dRef;              // x_c - x_R
  Matrix rvec;              // r_i (nIndep x ndm)
  Matrix Icplus;            // pseudo-inverse of I_c (nrot x nrot)
  Matrix Pproj;             // projector on the range of I_c (nrot x nrot)
  int nKept;                // number of constrained reference rotation axes

  int nDOF;                 // sum of node ndf
  ID nodeNdf;               // ndf of each node (1+N)
  ID dofOffset;             // element DOF offset of each node (1+N)

  Matrix *B;                // gap operator (ndm+nrot) x nDOF

  bool initGapCapture;      // false with -absolute
  bool g0Computed;
  Vector g0;                // initial translational gap (ndm)
  Vector gr0;               // initial rotational gap (nrot)

  Node **theNodes;
  Matrix *K;
  Vector *P;
  Matrix *M0;               // zero mass and damping matrix (nDOF x nDOF)

  void allocate(void);
  void resolveGeometry(void);
  void buildB(void);
  void resolvePenalty(void);
  void computeGap(Vector &g);
  void captureInitialGap(void);
  double crossOp(int i, int r, int j) const;
  double transOp(int k, int r) const;
};

#endif
