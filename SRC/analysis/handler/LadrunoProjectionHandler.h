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
//   J.O. Hallquist, LS-DYNA Theory Manual, Livermore Software Technology
//   Corporation, 2006, Sections 28.1 and 28.2.
//
// Description: constraint handler for explicit dynamics that enforces the
// multi-point constraints (MP_Constraint, EQ_Constraint) by projecting the
// solved acceleration onto the constraint manifold instead of eliminating
// the constrained DOFs.
//
// Assembly is that of PlainHandler: every node DOF keeps its own equation
// and its own diagonal mass, and homogeneous single-point constraints are
// removed from the equation set. The projection is carried out by a
// LadrunoConstraintProjector, which is passed to an integrator implementing
// LadrunoProjectionConsumer (CentralDifferenceLadruno, ExplicitBathe).
//
// handle() builds the connected groups of constrained DOFs (union-find over
// (node, dof) pairs), refuses constraint chains, DOFs constrained twice,
// a tied DOF that is also fixed, consistent mass and massless tied DOFs, and
// records the prescribed-motion DOFs. doneNumberingDOF() stores the groups
// with their equation numbers in the projector and passes it on.
//
// Limits: the constraint matrix is assembled once, so the lever arm of
// rigidLink -beam and rigidDiaphragm is frozen at the initial geometry (a
// note is printed when a retained node has rotated more than 0.1 rad); every
// tied DOF, rotations included, needs a nonzero lumped mass; constraints
// must not cross partition boundaries.
//
// Command: constraints LadrunoProjection <-verbose> <-projectICs> <-icTol $tol>

#ifndef LadrunoProjectionHandler_h
#define LadrunoProjectionHandler_h

#include <ConstraintHandler.h>
#include <vector>
#include <set>
#include <utility>

class LadrunoConstraintProjector;
class FE_Element;
class DOF_Group;
class Node;
class SP_Constraint;

class LadrunoProjectionHandler : public ConstraintHandler
{
  public:
    LadrunoProjectionHandler();
    LadrunoProjectionHandler(bool verbose, bool projectICs, double icTol);
    ~LadrunoProjectionHandler();

    int handle(const ID *nodesNumberedLast = 0);
    int doneNumberingDOF(void);
    int applyLoad(void);          // imposes the prescribed-motion displacements
    void clearAll(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  protected:

  private:
    // one constrained DOF and the retained DOFs it is tied to
    struct SlaveRec {
        int vtx;                                 // vertex of the constrained DOF
        std::vector<int> masterVtx;              // vertices of its retained DOFs
        std::vector<double> coeff;               // entries of the C row
        double delta;                            // Uc0(i) - sum_j C(i,j) Ur0(j)
    };
    // a connected group of constrained DOFs (filled at handle())
    struct RawGroup {
        std::vector<int> masterVtx;              // retained DOF vertices (columns of L)
        std::vector<SlaveRec> slaves;            // constrained rows of L
    };

    // A DOF with a non-homogeneous SP or imposedMotion. It is removed from
    // the equation set like a fixed DOF; applyLoad() imposes its displacement
    // each step, as TransformationDOF_Group::enforceSPs() does. ImposedMotionSP
    // sets the velocity and acceleration itself; a plain SP supplies only the
    // displacement.
    struct PrescribedDOF {
        int nodeTag;
        int dof;
        SP_Constraint *sp;
    };
    std::vector<PrescribedDOF> prescribedDOFs;
    std::set<std::pair<int,int> > prescribedKey;   // (nodeTag, dof) of every prescribed DOF

    // A constrained DOF tied only to prescribed retained DOFs is fully
    // determined by the prescribed motion: u_c = sum_k C_k u_k + delta,
    // v_c = sum_k C_k v_k, a_c = sum_k C_k a_k. It is removed from the
    // equation set and its state is imposed in applyLoad() after that of its
    // retained DOFs. A DOF tied to both free and prescribed retained DOFs is
    // refused at handle().
    struct DerivedSlave {
        int nodeTag, dof;
        std::vector<int> masterNode, masterDof;
        std::vector<double> coeff;
        double delta;                            // Uc0 - sum_k C_k Ur0_k
    };
    std::vector<DerivedSlave> derivedSlaves;
    std::set<std::pair<int,int> > derivedKey;

    // Retained node whose rotation enters a constrained translation through
    // the frozen lever arm of C. applyLoad() prints one note per node when
    // the rotation since handle() exceeds the small-rotation limit.
    struct RotMonitor {
        int nodeTag;
        std::vector<int> rotDof;     // rotational DOFs (2D: {2}; 3D: {3,4,5})
        std::vector<double> theta0;  // rotation at handle()
        bool warned;
    };
    std::vector<RotMonitor> rotMonitors;

    // a vertex is one (node, dof) pair
    std::vector<int> vtxNode;
    std::vector<int> vtxDof;
    std::vector<RawGroup> theGroups;

    LadrunoConstraintProjector *theProjector;    // owned
    bool verbose;
    bool projectICs;
    double icTol;

    int buildGroups(void);
    int classifyDerivedSlaves(int &countDOF);
    int massGuard(void);
    void flagRotMonitor(int masterNode, int masterDof, int slaveNode, int slaveDof);
};

#endif
