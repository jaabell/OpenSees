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
// Description: KinematicCoupling, a rigid kinematic coupling (RBE2-type)
// enforced by penalty or augmented Lagrangian.
//
// A reference node R (translations and, optionally, rotations) rigidly
// drives a set of N slave nodes S_i. Each slave follows the small-rotation
// rigid-body motion of R:
//
//     u_i     = u_R + theta_R x d_i,      d_i = x_i - x_R
//     theta_i = theta_R                   (where a slave rotation is tied)
//
// The constraint is written as one scalar gap per tied slave component:
//
//     g_i^t = u_i - u_R - theta_R x d_i        (translation rows)
//     g_i^r = theta_i - theta_R                (rotation rows)
//
// g = B u with a constant operator B (built once from the initial
// coordinates). The tie force and tangent are
//
//     t = D g (+ lambda),   P = B^T t,   K = B^T D B,
//     D = diag(K_t on translation rows, K_r on rotation rows).
//
// P is self-equilibrated for any penalty value. With -enforce al the
// multipliers are updated once per committed step (Uzawa iteration),
// lambda <- lambda + D g, which drives the gap toward zero over successive
// steps at a moderate penalty.
//
// The default rotational penalty is K_r = K_t * l^2, with l^2 the larger of
// the mean and the maximum |d_i|^2 (unit length when every slave is
// coincident with R), so that translation and rotation rows have the same
// order of stiffness.
//
// The same kinematics are the Nastran RBE2 element and the Abaqus
// *COUPLING, *KINEMATIC constraint:
//   MSC Nastran Quick Reference Guide, bulk data entry RBE2.
//   Abaqus Analysis User's Guide, section "Coupling constraints"
//     (kinematic coupling).
// Penalty and augmented Lagrangian enforcement of linear constraints:
//   P. Wriggers, Computational Contact Mechanics, 2nd ed., Springer, 2006,
//     chapter 6.
//   J.C. Simo, T.A. Laursen, An augmented Lagrangian treatment of contact
//     problems involving friction, Computers & Structures 42 (1992) 97-116.
//
// Usage:
//   element KinematicCoupling tag refNode N s1 ... sN
//       <-dof c1 ... cK> <-k Kt | -k auto -host eleTag <-kAlpha a>>
//       <-kr Kr> <-enforce penalty|al> <-absolute>
//
// Constraints on use: small rotations (linear kinematics, B is constant);
// the coupling carries no mass and no Rayleigh damping.

#ifndef KinematicCoupling_h
#define KinematicCoupling_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

class Node;
class Channel;
class FEM_ObjectBroker;
class Response;

class KinematicCoupling : public Element
{
 public:
  KinematicCoupling(int tag, int ndm, int refNode, const ID &slaveNodes,
                    const ID &dofSel, double kt, double kr, bool krUser,
                    int enforce, double kAlpha, int hostEleTag, bool ktAuto,
                    bool initGapCapture = true);
  KinematicCoupling();
  ~KinematicCoupling();

  const char *getClassType(void) const { return "KinematicCoupling"; }

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
  int nSlave;               // number of slave nodes N
  ID connectedNodes;        // [ refNode, slave_1 .. slave_N ]
  ID dofSel;                // requested 1-based slave components (empty = all)

  double Kt;                // translational penalty
  double Kr;                // rotational penalty
  bool krUser;              // -kr given (else K_r = K_t * l^2)
  bool ktAuto;              // -k auto: K_t from the -host element stiffness
  double kAlpha;            // multiplier for -k auto
  int hostEleTag;           // representative host element, -1 if none
  bool ktResolved;          // K_t / K_r resolved (transient)
  bool hostMissWarned;      // missing -host reported (transient)
  double ell2;              // rotation length scale l^2

  int enforce;              // 0 = penalty, 1 = augmented Lagrangian
  Vector lambdaAL;          // per-row multipliers (size nGap)

  bool hasRefRot;           // R carries rotation DOFs

  // geometry and gap layout, resolved at setDomain
  bool valid;
  Matrix dvec;              // d_i = x_i - x_R (nSlave x ndm)
  int nGap;                 // number of tied scalar components
  ID gapNode;               // node slot (1..N) of each row
  ID gapDof;                // node-local DOF index of each row
  ID gapIsRot;              // 1 for a rotation row

  int nDOF;                 // sum of node ndf
  ID nodeNdf;               // ndf of each node (1+N)
  ID dofOffset;             // element DOF offset of each node (1+N)

  Matrix *B;                // gap operator nGap x nDOF

  bool initGapCapture;      // false with -absolute
  bool g0Computed;
  Vector g0;                // initial gap (size nGap)

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
  double rowPenalty(int row) const { return gapIsRot(row) ? Kr : Kt; }
  double transOp(int i, int k, int r) const;
};

#endif
