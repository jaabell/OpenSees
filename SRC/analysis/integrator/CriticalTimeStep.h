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
//   Belytschko, T., Liu, W.K., Moran, B. and Elkhodary, K.I. (2014).
//     Nonlinear Finite Elements for Continua and Structures, 2nd ed.,
//     Wiley, Ch. 6 (stability of the central difference method and the
//     element eigenvalue inequality).
//   Hughes, T.J.R. (1987). The Finite Element Method: Linear Static and
//     Dynamic Finite Element Analysis, Prentice-Hall, Ch. 9.
//
// Description: per-element critical time step for explicit integrators.
// The governing step is the minimum over the elements of the local domain of
// 2/omega_max, with omega_max^2 the largest eigenvalue of the element pencil
// K_e v = lambda M_e v and M_e the lumped (diagonal) element mass. By the
// element eigenvalue inequality the largest element eigenvalue bounds the
// largest eigenvalue of the assembled lumped-mass system from above, so the
// estimate is conservative for an unconstrained mesh. With Rayleigh damping
// the central-difference limit reduces to
//
//   dt_cr = (2/omega_max) (sqrt(1 + xi^2) - xi),
//   xi = (alphaM/omega_max + betaK omega_max) / 2,
//
// with betaK the sum of the three stiffness-proportional coefficients.
//
// The estimate ignores constraints (equalDOF, rigid diaphragms, MP
// constraints) and pure nodal mass, so on constrained or nodal-mass-only
// models it can be non-binding or absent (no element contributes).

#ifndef CriticalTimeStep_h
#define CriticalTimeStep_h

class AnalysisModel;
class Element;
class Matrix;
class Vector;

// How the element mass is lumped to the diagonal of the element pencil.
enum class CTSLumping {
    RowSum,   // (0) row sum onto the diagonal. Conserves mass, but can give
              //   zero or negative entries on higher-order elements and zero
              //   rotational inertia on beams and shells.
    Diagonal, // (1) diagonal of the consistent mass, M_ii. Positive for a
              //   consistent mass, so it avoids the rotational defect of the
              //   row sum, but it does not conserve mass.
    HRZ       // (2) Hinton-Rock-Zienkiewicz lump (LadrunoMassLumping.h): the
              //   consistent diagonal scaled per translational direction to
              //   the element mass; rotational DOFs take the mean factor.
              //   Positive, mass-conserving, nonzero rotational inertia.
};

// Result of one scan over the elements of the (local) domain.
struct CTSResult {
    double damped_dt;      // central-difference limit including element
                           //   Rayleigh damping; +inf if nothing contributed
    double undamped_dt;    // 2/omega_max; +inf if nothing contributed
    int    damped_tag;     // governing element tag (damped); -1 if none
    int    undamped_tag;   // governing element tag (undamped); -1 if none
    int    n_contributing; // elements that produced a finite positive estimate
    int    n_scanned;      // elements iterated
};

// Scan all elements of the domain of theModel and return the governing step.
//   useTangent : the pencil uses getTangentStiff() (current state) instead of
//                getInitialStiff(), so a refresh can track stiffening or
//                softening (contact, plasticity).
//   lumping    : how the element mass is lumped to the diagonal.
//
// An element that reports a positive Element::getExplicitCriticalTimeStep()
// contributes that value and skips the eigensolve.
//
// Numerics: the lumped mass is diagonal and usually positive, so the
// symmetric-definite driver DSYGVX is tried first, asking for the largest
// eigenvalue only; if the mass is not positive definite (massless DOFs) the
// general driver DGGEV is used, with a relative threshold on beta that
// discards the near-infinite eigenvalues of massless modes. In a parallel
// build the result is reduced (MPI_MIN) across processes.
//
// The mass is copied to the diagonal before the stiffness is fetched: many
// elements return a reference to one shared static matrix from getMass(),
// getInitialStiff() and getTangentStiff(), so holding two of those references
// at once would alias.
CTSResult computeCriticalTimeStep(AnalysisModel *theModel,
                                  bool useTangent,
                                  CTSLumping lumping);

// Lump the consistent mass M of one element into mdiag (length M.noRows())
// according to lumping. M is fully copied before returning. For HRZ the
// direction of each DOF is taken from the node layout of the element: the
// first ndm DOFs of each node are translational and the element matrix is
// assumed node by node.
void lumpElementMass(Element *ele, const Matrix &M, CTSLumping lumping,
                     double *mdiag);

// Largest generalized eigenvalue lambda_max = omega_max^2 of K v = lambda M v
// for one element, with M the diagonal mdiag of length n. useTangent selects
// getTangentStiff() over getInitialStiff(). Returns < 0 on a size mismatch or
// an eigensolver failure. The self-reported step of the element is not
// considered here.
double elementLambdaMax(Element *ele, bool useTangent,
                        const double *mdiag, int n);

// Total stiffness-proportional Rayleigh coefficient of an element. The
// element damping is C = alphaM M + betaK K + betaK0 K_init + betaKc K_commit;
// the three beta terms produce the same high-frequency damping (equal at the
// initial state, where K = K_init = K_commit), so all three shrink the
// explicit step and all enter the damped estimate. Each coefficient is
// clamped at 0 before the sum.
double stiffnessRayleighBeta(Element *ele);

// true if every entry of v is finite. Divergence checks on solved vectors use
// this rather than Vector::pNorm(0): the max-compare in pNorm skips NaN
// entries, so only +/-Inf would be detected.
bool vectorIsFinite(const Vector &v);

#endif
