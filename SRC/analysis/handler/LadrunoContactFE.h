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
//   Puso, M.A. and Laursen, T.A. (2004). A mortar segment-to-segment contact
//     method for large deformation solid mechanics. Computer Methods in Applied
//     Mechanics and Engineering, 193(6-8), 601-629. doi:10.1016/j.cma.2003.10.010
//   Simo, J.C. and Laursen, T.A. (1992). An augmented Lagrangian treatment of
//     contact problems involving friction. Computers & Structures, 42(1),
//     97-116. doi:10.1016/0045-7949(92)90540-G
//   Ericson, C. (2005). Real-Time Collision Detection. Morgan Kaufmann,
//     San Francisco, Section 5.1.9. doi:10.1201/b14581
//
// Description: LadrunoContactFE, an FE_Element adapter that adds contact
// contributions to the assembly without a Domain Element behind it
// (myEle == 0). getResidual() and getTangent() read the trial displacements of
// the nodes, as an Element would, and assemble the contact force and
// stiffness. getTangent() goes through Integrator::formEleTangent(), so the
// integrator decides what is assembled (an explicit integrator that forms
// only M and C never sees the contact stiffness). The adapter is a stateless
// view: path-dependent pair state lives on the Domain-owned
// LadrunoContactDomain, so adapters can be rebuilt in every handle().
//
// Modes:
//   RIGID_PLANE  one slave node against a rigid analytical plane, penalty
//            normal, connectivity {slave}, ndof = ndm.
//   SEGMENT  node-to-segment penalty contact of one slave node against one
//            master segment (tri3 or quad4), connectivity {slave, segment
//            nodes}, ndof = 3 (1 + nps).
//   MORTAR   mortar segment-to-segment contact of one slave facet against one
//            master facet (Puso and Laursen 2004): the overlap of the two
//            facets on the auxiliary plane of the master facet is clipped and
//            integrated with a triangle rule, giving the mortar matrices D
//            (slave x slave), M (slave x master) and the weighted gap gt.
//            Connectivity {slave facet nodes, master facet nodes},
//            ndof = 3 (npsS + npsM).
//   EDGE_EDGE  edge-to-edge penalty contact of one slave edge against one
//            master edge, for facets that meet edge-on: the closest points of
//            the two segments (Ericson 2005, Section 5.1.9), the common
//            perpendicular as normal and a body-fixed sign
//            (LadrunoEdgeKernel). Connectivity {slave edge nodes, master edge
//            nodes}, ndof = 12.
//
// 2D (ndm == 2) is a parameterization of the SEGMENT and MORTAR modes, with
// the geometry of LadrunoContact2DKernel and two DOFs per node:
//   SEGMENT  a slave node against a 2-node master line (nps = 2), or against
//            one master vertex (nps = 1, a degenerate segment) at a concave
//            corner or at an open end of the master chain;
//   MORTAR   a 2-node slave segment against a 2-node master segment.
//
// With myEle == 0 the FE_Element getResidual()/getTangent() and the force
// helpers have no element to call, so this class owns its buffers and
// overrides them.

#ifndef LadrunoContactFE_h
#define LadrunoContactFE_h

#include <FE_Element.h>
#include <Vector.h>
#include <Matrix.h>

class Integrator;
class Node;
class Domain;

class LadrunoContactFE : public FE_Element
{
  public:
    // Rigid analytical plane against ONE slave node. Connectivity = {slave}
    // (its first ndm translational DOFs). p0 = a point on the plane, n = the
    // outward unit normal (toward the slave's allowed half-space), kn = penalty.
    // contactTag = the contactPlane tag. muc = viscous normal stabilization
    // coefficient (0: none).
    LadrunoContactFE(int tag, Node *slaveNode, int ndm,
                     const double p0[3], const double n[3], double kn,
                     Domain *theDomain = 0, int contactTag = 0, double muc = 0.0);
    // Node-to-segment penalty contact against ONE master segment (tri3 or
    // quad4, nps nodes). Connectivity = {slave} + {segment nodes}. The normal is
    // derived at each evaluation and oriented toward orientDir (the slave's
    // allowed half-space), so it does not depend on the node winding.
    //   theDomain, contactTag, segIndex: the engine lookup (fetched again at
    //     every evaluation, because wipe deletes the engine) and the pair key
    //     (contact tag, slave tag, global master-segment ordinal).
    //   consistentNormal: add the geometric normal tangent kn gN d^2gN/du^2
    //     (false: the main term kn B^T B only, exact for a flat fixed master).
    //   kt, mu, cohesion, tauMax: penalty-regularized Coulomb/Tresca friction,
    //     cone cap = min(mu N + cohesion, tauMax). The plastic slip lives on the
    //     engine, keyed by the pair. mu, cohesion and tauMax all <= 0: the
    //     frictionless path, no engine state is touched.
    //   consistentTan: assemble the non-symmetric consistent friction tangent
    //     (false: the symmetric one, safe on any solver). The friction tangent
    //     enters only through addKtToTang/addKiToTang, i.e. implicit integrators.
    //   muc: viscous normal stabilization coefficient (0: none).
    LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes, int nps, double kn,
                     const double orientDir[3], double kt = 0.0, double mu = 0.0,
                     Domain *theDomain = 0, int contactTag = 0, int segIndex = 0,
                     bool consistentTan = false, bool consistentNormal = false,
                     double cohesion = 0.0, double tauMax = 0.0, double muc = 0.0);
    // Mortar contact of ONE slave facet against ONE master facet (tri3 or quad4
    // each). Connectivity = {slave facet nodes} + {master facet nodes}.
    //   epsN: normal penalty of the augmented Lagrangian.
    //   orientDir: direction toward the slave's allowed half-space (orients the
    //     master facet normal).
    //   contactTag, slaveFacetIndex: the engine keys; theDomain is the engine
    //     lookup (null: penalty only, no multiplier).
    //   ngp: order of the triangle rule on the clipped overlap.
    //   isTie: a mesh tie (permanent bond of the full relative displacement,
    //     epsN is then the tie penalty).
    //   mu, epsT, cohesion, tauMax: Coulomb/Tresca friction on the cone
    //     min(mu N + cohesion, tauMax) with tangential penalty epsT; all of mu,
    //     cohesion, tauMax <= 0: frictionless, no engine state is touched.
    //   consistentTan: the non-symmetric consistent friction tangent.
    //   masterFacetIndex: ordinal of the master facet (the friction state key).
    LadrunoContactFE(int tag, Node **slaveNodes, int npsS, Node **masterNodes, int npsM,
                     double epsN, const double orientDir[3], int contactTag,
                     int slaveFacetIndex, Domain *theDomain, int ngp = 2,
                     bool isTie = false, double mu = 0.0, double epsT = 0.0,
                     double cohesion = 0.0, double tauMax = 0.0, bool consistentTan = false,
                     int masterFacetIndex = 0);
    // Edge-to-edge penalty contact of ONE slave edge (sNodeA, sNodeB) against
    // ONE master edge (mNodeA, mNodeB). Connectivity = {sa, sb, ma, mb}.
    //   epsN: normal penalty; orientDir: reference direction of the first sign
    //     capture (toward the slave's allowed side).
    //   contactTag, theDomain: the engine keys and lookup (the sign, friction
    //     and multiplier state is per edge pair on the engine).
    //   mu, kt, cohesion, tauMax, consistentTan: friction, as for the other
    //     modes; all of mu, cohesion, tauMax <= 0: frictionless.
    //   edgeAlm: use the committed multiplier of the pair, p = min(0,
    //     lambdaN + epsN gN) (augmented Lagrangian, updated by commit()).
    LadrunoContactFE(int tag, Node *sNodeA, Node *sNodeB, Node *mNodeA, Node *mNodeB,
                     double epsN, const double orientDir[3], int contactTag,
                     Domain *theDomain, double mu = 0.0, double kt = 0.0,
                     double cohesion = 0.0, double tauMax = 0.0, bool consistentTan = false,
                     bool edgeAlm = false);
    // 2D node-to-segment contact (SEGMENT mode, ndm2 == 2). nps == 2: the slave
    // against the master line segNodes[0] -> segNodes[1]; nps == 1: the slave
    // against the vertex segNodes[0] (connectivity {slave, vertex}).
    //   sigma: the orientation of the master surface (+1/-1), decided once per
    //     interface by the handler; Lref: the minimum initial segment length of
    //     the master surface (the relative gauge of the kernel).
    //   prevFar, nextFar: the far nodes of the neighbouring segments in surface
    //     order (null at an open end). On a segment pair a previous neighbour
    //     makes the previous segment own a slave that both project in bounds;
    //     on a vertex pair both give the corner (or, with one of them null,
    //     the end cap). They steer the active set only and are not part of the
    //     connectivity; they are stored by tag and looked up at every
    //     evaluation, so a removed node reads as an open end.
    //   kt, mu, cohesion, tauMax, consistentTan: scalar friction along the
    //     tangent perp(n); muc: viscous normal stabilization.
    LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes, int nps, int ndm2,
                     double kn, double sigma, double Lref, Node *prevFar, Node *nextFar,
                     double kt, double mu, bool consistentTan, double muc,
                     Domain *theDomain, int contactTag, int segIndex,
                     double cohesion = 0.0, double tauMax = 0.0);
    // 2D mortar contact (MORTAR mode, ndm2 == 2) of a 2-node slave segment
    // against a 2-node master segment. sigma and Lref as above (Lref over both
    // surfaces); epsN, epsT, cohesion and tauMax already multiplied by the
    // thickness. ntsCoDeclared: a node-to-segment contact covers both slave
    // nodes on the same master, so a pair too nearly perpendicular to be
    // integrated is left to it silently (else it is reported once).
    LadrunoContactFE(int tag, Node **slaveNodes, int npsS, Node **masterNodes, int npsM,
                     int ndm2, double epsN, double sigma, double Lref, int contactTag,
                     int slaveFacetIndex, Domain *theDomain, bool isTie, double mu,
                     double epsT, double cohesion, double tauMax, bool consistentTan,
                     int masterFacetIndex, bool ntsCoDeclared);
    ~LadrunoContactFE();

    // Mortar normal-gap shift: gbar_I -> gbar_I + gapOffset, and with adjust
    // the weighted gap of this pair at the reference configuration is
    // subtracted (only where |gbar_ref| <= adjustTol when adjustTol > 0). The
    // reference gap is integrated at the first evaluation.
    void setMortarGapShift(double gapOffset, bool adjust, double adjustTol);

    // Smoothed mortar contact law: smoothN > 0 (a length g0) replaces the kink
    // of p = min(0, epsN gbar) by a C1 quadratic onset over |gbar| < g0 and
    // fades the friction traction in over the same band; 0 < smoothT < 1 rounds
    // the stick/slip corner (LadrunoFrictionKernel::frictionReturnMapSmooth).
    void setMortarSmoothing(double smoothN, double smoothT);

    // Fill myID with the first ndm equations of each node (the translations).
    // FE_Element::setID() copies each node's whole DOF_Group ID, which on a
    // node with ndf > ndm (a beam or shell node, a u-p node) overflows the
    // ndm-per-node layout of the adapter. The extra DOFs are not coupled.
    int setID(void);

    const Vector &getResidual(Integrator *theIntegrator);
    const Matrix &getTangent(Integrator *theIntegrator);

    void zeroTangent(void);
    void addKtToTang(double fact = 1.0);  // += fact K_c when the pair is active
    void addKiToTang(double fact = 1.0);  // initial stiffness: the penalty K_c, so
                                          // Newton -initial keeps the contact stiffness
    void addCtoTang(double fact = 1.0);   // viscous stabilization muc B^T B (else none)
    void addMtoTang(double fact = 1.0);   // contact pairs carry no mass

    // The FE_Element versions need an element; reachable off the assembly
    // path (e.g. modal damping), so they return a zero vector here.
    const Vector &getTangForce(const Vector &x, double fact = 1.0);
    const Vector &getK_Force(const Vector &x, double fact = 1.0);
    const Vector &getKi_Force(const Vector &x, double fact = 1.0);
    const Vector &getC_Force(const Vector &x, double fact = 1.0);
    const Vector &getM_Force(const Vector &x, double fact = 1.0);

    // Soft-constraint penalty (-soft / -edgeSoft) of this pair, s > 0: under
    // an integrator that assembles no contact stiffness the penalty is
    // s 4 m_eff / dt^2 (see LadrunoContactDomain::setSoftScale); 0: off.
    void setSoftScale(double s) { softScale = (s > 0.0) ? s : 0.0; }

  private:
    // gap n.(x_s - p0) at the slave's current trial position; < 0 = penetration
    double rigidPlaneGap(void) const;

    // Project the slave onto the master segment at the current trial
    // configuration and derive the oriented normal, gap and shape weights.
    // Returns true when active (penetrating and in-bounds); then fills gap (< 0),
    // n[3], N[nps] and the gap operator B[ndof] = [ n^T | -N_i n^T ] over
    // [slave, segment nodes].
    // If gTvec != 0 it also fills the tangential relative displacement at the
    // contact point (the slip measure of the friction return map).
    bool segmentActive(double &gap, double n[3], double N[4], double *B,
                       double *gTvec = 0) const;

    // The -visc dashpot is active when muc > 0 and the integrator is not a
    // StaticIntegrator. Under a static integrator the nodal trial velocities are
    // state left by an earlier transient stage, not rates, and the damping
    // tangent is never assembled, so the dashpot is disabled there (reported
    // once per contact).
    bool viscousActive(Integrator *theIntegrator) const;

    // true when the pair carries friction (mu, cohesion or tauMax > 0)
    bool hasFriction(void) const { return mu > 0.0 || cohesion > 0.0 || tauMax > 0.0; }

    // Scatter the 3 x 3 friction block K_ss (LadrunoFrictionKernel::
    // frictionTangentBlock) into `tang` with G = [I | -N_i I]: block (a,b) is
    // w_a w_b K_ss, w = [1, -N_0, ..., -N_nps].
    void addFrictionTang(double fact, const double n[3], const double N[4], double tn,
                         const double gTeff[3], const double gpT[3], bool consistent);

    // Geometric part of the consistent normal tangent, kn gN d^2gN/du^2
    // (LadrunoContactProjection::normalGeomTangent), added to `tang`.
    void addNormalGeomTang(double fact);

    // Integrate the mortar facet pair at the current trial configuration
    // (LadrunoMortarKernel::integratePair). Returns false when the overlap is
    // empty or degenerate; otherwise fills D, M, the weighted gap g (one per
    // slave node) and the master facet normal n oriented by orientDir.
    bool mortarActive(double D[4][4], double M[4][4], double g[4], double n[3]) const;

    // Mortar tangent K = epsN Bt^T diag(act_I / a_I) Bt (x) (n (x) n),
    // Bt = [D, -M] over [slave | master] nodes, on the active set of the
    // residual. The variation of D, M and n with the displacements is not
    // included (exact for flat facets that do not slide).
    void addMortarTang(double fact, bool initialStiff = false);

    // Mesh-tie force: for every slave node of the pair the weighted relative
    // displacement r_I = sum_J D_IJ u_s,J - sum_K M_IK u_m,K (from the trial
    // displacements: the bond holds the as-built configuration), the traction
    // t_I = lambdaTie_I + epsTie r_I / a_I (all three components, no clamp),
    // scattered as the normal force: f^s_K = -sum_I D_KI t_I,
    // f^m_L = +sum_I M_IL t_I.
    void addMortarTieForce(const double D[4][4], const double M[4][4],
                           class LadrunoContactDomain *cd);

    // Mortar friction force (the tangent is in addMortarTang). For every slave
    // node in normal contact (p_normal[I] < 0): the weighted relative
    // displacement r_I, its tangential part over a_I as the slip (positions
    // carry no tangential information after the projection), measured from the
    // engagement origin and offset by lambdaT / epsT (augmented trial traction
    // lambdaT + epsT (gTeff - gpT)), the return map of LadrunoFrictionKernel
    // with N = -p_normal[I], and the scatter +D tF on the slave facet and
    // -M^T tF on the master facet.
    // chi (null: 1): the onset weights of the smoothed law, scaling the traction.
    void addMortarFriction(const double D[4][4], const double M[4][4], const double n[3],
                           const double p_normal[4], class LadrunoContactDomain *cd,
                           const double *chi = 0);

    // Smoothed normal law at the augmented argument pr = lambda + epsN gbar,
    // x = -pr, S = epsN g0: P = 0 (x <= -S), (x + S)^2 / (4 S) (|x| < S),
    // x (x >= S); dP = dP/dx; chi = C1 smoothstep of the friction onset over
    // the same band and dchi = dchi/dx.
    void mortarSmoothLaw(double pr, double &P, double &dP, double &chi, double &dchi) const;

    // Edge-to-edge geometry at the current trial configuration
    // (LadrunoEdgeKernel::closestPtSegSeg). Returns false unless the edges are
    // non-parallel, of non-zero length and cross strictly inside both (a pair
    // near a vertex belongs to the other lanes). Fills the signed gap gN, the
    // oriented common-perpendicular normal n, the closest-point parameters s
    // and t and the gap operator B[4][3] = [(1-s) n, s n, -(1-t) n, -t n].
    // committedSign (+1/-1) is applied as given; 0 orients n by orientDir and
    // returns the chosen sign in *outSign for the caller to keep.
    bool edgeGeom(double &gN, double n[3], double &s, double &t, double B[4][3],
                  int committedSign, int *outSign) const;

    // Tangential slip of the edge pair at its closest points: the tangential
    // part of the relative displacement (1-s) u_a0 + s u_a1 - (1-t) u_b0 - t u_b1
    // (the closest-point construction makes the relative position purely
    // normal, so positions carry no tangential information).
    void edgeSlip(double s, double t, const double n[3], double gT[3]) const;

    // Step size of the soft penalty: < 0 when the integrator assembles a
    // stiffness for this adapter (the configured penalty applies), else the
    // step size (0 when not known yet).
    double explicitTimeStep(Integrator *theIntegrator);
    void softWarning(void);   // once per contact: the soft penalty could not be sized
    // B M^-1 B^T of the gap operator of a segment pair along dir (N == 0: the
    // slave only, a rigid plane) and of an edge pair
    double gapModeInvMass(const double dir[3], const double N[4]) const;
    double edgeGapModeInvMass(const double dir[3], double s, double t) const;
    // the normal and stick penalties of this evaluation: the soft values under
    // an explicit integrator, the configured kn and kt otherwise
    double softKn(Integrator *theIntegrator, const double n[3], const double N[4]);
    double softKt(Integrator *theIntegrator, const double n[3], const double N[4]);
    double softKnEdge(Integrator *theIntegrator, const double n[3], double s, double t);
    double softKtEdge(Integrator *theIntegrator, const double n[3], double s, double t);
    // the segment-based soft penalty of a 3D mortar pair (residual only)
    void addSoft2Penalty(double dt);

    // Edge-to-edge tangent: epsN B^T B on the active set of the residual plus
    // the friction block; initialStiff uses the stick block.
    void addEdgeTang(double fact, bool initialStiff);

    // 2D node-to-segment narrow phase (SEGMENT, ndm == 2). nps == 2: the
    // projection on the segment, refused when the previous segment in surface
    // order also projects in bounds (it owns the slave). nps == 1 with both
    // neighbours: the concave vertex, active when both adjacent segments
    // refuse and the slave is in the wedge (aPrev > 0, aNext < 0); with one
    // neighbour: the end cap, active past the open end within 5 percent of
    // Lref. The side sign of a vertex pair is committedSide when non-zero,
    // else the live one, returned in *liveSideOut for the caller to commit.
    // Fills gap (< 0), n (n[2] = 0), the weights N (segment: shape functions,
    // vertex: {1, 0}), B over [slave xy | node xy ...] and, if gTout != 0, the
    // scalar tangential relative displacement along perp(n). Returns true only
    // when this adapter owns the slave and it penetrates.
    bool segment2DActive(double &gap, double n[3], double N[2], double B[6],
                         double committedSide, double *liveSideOut,
                         double *gTout = 0) const;
    // the committed side sign of a 2D vertex pair (0: none or not captured)
    double vtx2DCommittedSide(void) const;
    Node *nts2dPrevFarNode(void) const;   // null: no neighbour or removed
    Node *nts2dNextFarNode(void) const;

    // 2D friction tangent G^T K G, G = [I2 | -N_i I2], with the 2 x 2 block
    // Kss t (x) t + dTN t (x) n from the scalar kernel (t = perp(n)).
    void addFrictionTang2D(double fact, const double n[3], const double N2[2], double tn,
                           double gTeff, double gpT, bool consistent);

    // 2D mortar pair at the current trial configuration
    // (LadrunoContact2DKernel::mortarSegmentPair2D): D, M, g and the master
    // normal; refusalOut and cosOut are always written.
    bool mortarActive2D(double D[2][2], double M[2][2], double g[2], double n[2],
                        int &refusalOut, double &cosOut) const;
    // 2D mortar friction force (scalar per slave node along perp(n)), tie
    // force and tangent, the 2D forms of addMortarFriction, addMortarTieForce
    // and addMortarTang
    void addMortarFriction2D(const double D[2][2], const double M[2][2], const double n[2],
                             const double p_normal[2], class LadrunoContactDomain *cd);
    void addMortarTie2D(const double D[2][2], const double M[2][2],
                        class LadrunoContactDomain *cd);
    void addMortarTang2D(double fact, bool initialStiff = false);

    Vector resid;
    Matrix tang;

    enum Mode { RIGID_PLANE = 1, SEGMENT = 2, MORTAR = 3, EDGE_EDGE = 4 };
    Mode mode;
    Node *theSlave;
    int ndm;
    double kn;

    // RIGID_PLANE binding
    double planeP0[3];
    double planeN[3];

    // SEGMENT binding
    Node *segNode[4];   // master segment nodes (tri3: 3, quad4: 4)
    int nps;            // nodes per segment
    double orientDir[3];// fixed direction toward the slave's allowed half-space

    Domain *theDomain;  // engine lookup (the engine pointer is not cached: wipe deletes it)
    int contactTag;     // pair key: contact definition tag ...
    int segIndex;       // ... and the global master-segment ordinal (stable across rebuilds)
    bool consistentNormal;

    // friction
    double kt;          // tangential penalty
    double mu;          // Coulomb coefficient
    double cohesion;    // cone intercept
    double tauMax;      // shear cap of the cone (<= 0: none)
    bool consistentTan; // include the non-symmetric d_TN (x) n coupling

    // viscous normal stabilization: force -muc gdot B, gdot = B v, while in
    // contact; damping tangent muc B^T B through addCtoTang
    double muc;

    // MORTAR binding (kn holds epsN)
    Node *mortarSlave[4];   // slave facet nodes
    Node *mortarMaster[4];  // master facet nodes
    int npsS, npsM;         // nodes per slave / master facet
    int slaveFacetIndex;    // ordinal of the slave facet in its surface (stable across rebuilds)
    int mortarOrder;        // order of the triangle rule on the clipped overlap
    bool isTie;             // mesh tie
    int masterFacetIndex;   // ordinal of the master facet in its surface (stable across rebuilds)

    // mortar gap shift (setMortarGapShift); off unless an option was given
    bool   hasGapShift;
    double gapShiftOffset;
    bool   gapShiftAdjust;
    double gapShiftAdjTol;          // 0: adjust every node of the pair
    mutable bool   gapShiftRefReady;
    mutable double gapShiftRef[4];  // -gbar_I at the reference configuration

    // smoothed law (setMortarSmoothing); 0: off
    double smoothBand;              // g0
    double smoothSlip;              // r

    // EDGE_EDGE binding (kn holds epsN)
    Node *edgeNode[4];      // slave edge nodes A, B, master edge nodes A, B
    bool  edgeAlm;          // use the committed multiplier of the pair

    // 2D bindings (set only by the 2D constructors)
    double nts2dSigma      = 0.0;   // orientation of the master surface
    double nts2dLref       = 0.0;   // minimum initial master segment length
    bool   nts2dHasPrevFar = false; // a previous neighbour was given ...
    int    nts2dPrevFarTag = 0;     // ... with this node tag
    bool   nts2dHasNextFar = false;
    int    nts2dNextFarTag = 0;
    double mortar2dSigma   = 0.0;
    double mortar2dLref    = 0.0;   // minimum initial segment length of both surfaces
    bool   mortar2dNtsCoDeclared = false;
    bool   mortar2dPerpWarned    = false;   // the perpendicular-pair report was made

    // soft penalty (setSoftScale) and the stiffness probe of explicitTimeStep()
    double softScale         = 0.0;
    bool   stiffnessProbe    = false;   // add*ToTang record instead of assembling
    bool   probeSawStiffness = false;
};

#endif
