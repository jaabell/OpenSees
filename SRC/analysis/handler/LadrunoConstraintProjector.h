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
//   Corporation, 2006, Sections 28.1 (nodal constraint sets) and 28.2
//   (linear constraint equations).
//
// Description: mass-orthogonal projection of the explicit acceleration onto
// the manifold of linear multi-point constraints.
//
// For each connected group of constrained DOFs the relation between all the
// DOFs of the group and its retained (free) DOFs is stored as
//     a_group = L a_r ,   L = [ I ; C ]
// (rows ordered retained first, then constrained). The enforced acceleration
// is the projection of the raw acceleration in the metric of the lumped mass
// matrix M:
//     a_proj = L (L^T M L)^-1 L^T M a_raw
// which minimises |a - a_raw|_M on the constraint subspace. L^T M (a_proj -
// a_raw) = 0, so the generalised momentum of every retained direction is
// conserved. For equalDOF (C = I) the result is the mass-weighted mean
// acceleration of Section 28.1 of Hallquist (2006); for a general linear
// constraint it is the form [L^T M L]^-1 L^T F of Section 28.2.
//
// The lumped mass of each equation is read from the assembled DiagonalSOE,
// which is the mass the integrator inverts and includes element mass as well
// as nodal mass. buildMass() must run between the integrator's mass-only
// formTangent() and its solve(), because DiagonalDirectSolver overwrites the
// diagonal with its reciprocal when it factors.
//
// C is assembled once from MP_Constraint::getConstraint(); for rigidLink
// -beam and rigidDiaphragm the lever arm is frozen at the initial geometry
// (small rotation). Every DOF in a group keeps its own equation and needs a
// nonzero lumped mass.

#ifndef LadrunoConstraintProjector_h
#define LadrunoConstraintProjector_h

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class LinearSOE;

class LadrunoConstraintProjector
{
  public:
    LadrunoConstraintProjector();
    ~LadrunoConstraintProjector();

    // Called by the handler once per group, with equation numbers. The rows
    // of L are the nRet retained rows (identity) followed by the constrained
    // rows (C). allEqn holds the equation number of each row (all >= 0).
    // delta holds the offset Uc0 - C Ur0 of each constrained row, used by the
    // initial-condition check. fixedEqn lists constrained equations whose
    // retained DOFs are all fixed; their acceleration is set to zero.
    void addGroup(const ID &allEqn, int nRet, const Matrix &L,
                  const Vector &delta, const ID &fixedEqn);

    int numGroups(void) const { return (int)groups.size(); }

    // Reads diag(M) from the assembled DiagonalSOE and forms L^T M L for
    // each group. Returns -1 with a message for a SOE that is not a
    // DiagonalSOE, a massless DOF in a group or a numerically singular
    // L^T M L.
    int buildMass(LinearSOE *theSOE);

    // In-place projection of the acceleration vector a (size numEqn).
    void project(Vector &a);

    // Compliance of (u, v) with every group in the multi-point constraint
    // increment convention u_c - C u_r = delta, v_c - C v_r = 0. Returns the
    // number of violated constrained rows and prints one message per row.
    int checkIC(const Vector &u, const Vector &v, double icTol) const;

    // Sets the constrained components of u and v to the values implied by
    // the retained DOFs (u_c = C u_r + delta, v_c = C v_r; zero for the
    // constrained DOFs whose retained DOFs are all fixed).
    void snapICs(Vector &u, Vector &v) const;

    // Initial-condition policy. enforceIC() is the call the integrator makes
    // in domainChanged(): with projectICs it snaps (u, v) and returns 0,
    // otherwise it returns the number of violations of checkIC().
    void setICPolicy(bool projectICs_, double icTol_) { projectICs = projectICs_; icTol = icTol_; }
    int enforceIC(Vector &u, Vector &v) const;

    bool isMassReady(void) const { return massReady; }
    void invalidateMass(void) { massReady = false; }

  private:
    struct Group {
        ID allEqn;      // equation of each row of L (retained first), all >= 0
        int nRet;       // number of retained rows and columns (top block of L = I)
        Matrix L;       // (nRet + nCon) x nRet
        Vector delta;   // nCon offsets of the constrained rows
        Vector m;       // lumped mass of each row (filled by buildMass)
        Matrix LtML;    // nRet x nRet
        ID fixedEqn;    // constrained equations set to zero (all retained DOFs fixed)
        Group() : nRet(0) {}
    };
    std::vector<Group> groups;
    bool massReady;
    bool projectICs;
    double icTol;
};

#endif
