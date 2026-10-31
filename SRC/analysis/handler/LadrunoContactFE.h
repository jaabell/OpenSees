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
    ~LadrunoContactFE();

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

    // Current time step of the explicit integrator driving the analysis, for
    // a penalty sized from the Courant limit (the -soft option). Returns 0 when
    // the integrator does not report one; the -soft option is refused at the
    // command level while no integrator does.
    static double explicitTimeStep(Integrator *theIntegrator);

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

    Vector resid;
    Matrix tang;

    enum Mode { RIGID_PLANE = 1, SEGMENT = 2 };
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
};

#endif
