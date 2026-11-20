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
//   Foss, K.A. (1958). Co-ordinates which uncouple the equations of motion
//     of damped linear dynamic systems. Journal of Applied Mechanics,
//     25(3), 361-364. doi:10.1115/1.4011828
//   Veletsos, A.S. and Ventura, C.E. (1986). Modal analysis of
//     non-classically damped linear systems. Earthquake Engineering and
//     Structural Dynamics, 14(2), 217-243. doi:10.1002/eqe.4290140205
//   Tisseur, F. and Meerbergen, K. (2001). The quadratic eigenvalue
//     problem. SIAM Review, 43(2), 235-286.
//     doi:10.1137/S0036144500381988
//
// Description: complex (state-space) modal analysis of a non-classically
// damped model, reduced onto the real undamped modes of the last `eigen`
// run.
//
//   Reduced pencil. With Phi the p retained real modes, the damped free
//   vibration u = Phi z e^{lambda t} gives the quadratic eigenproblem
//
//       (lambda^2 Mt + lambda Ct + Kt) z = 0,   Mt,Ct,Kt  p x p,
//
//   Mt = Phi^T M Phi, Ct = Phi^T C Phi, Kt = Phi^T K Phi. Ct is diagonal
//   only for classical (e.g. Rayleigh) damping; a localized damper couples
//   the modes and the eigenvalues become complex with mode shapes whose
//   entries have different phases.
//
//   Linearization (Foss 1958; Tisseur and Meerbergen 2001, the
//   symmetric form): with w = {z ; lambda z} the quadratic problem is
//
//       B w = -lambda A w,   A = [Ct Mt; Mt 0],  B = [Kt 0; 0 -Mt].
//
//   LAPACK dggev receives the pair (-B, A), (-B) w = lambda A w, so its
//   eigenvalue is lambda itself. A is nonsingular whenever Mt is, so an
//   infinite eigenvalue (beta ~ 0) only arises from a singular Mt; it is
//   tagged, never divided through.
//
//   Projection of C. The default path projects element by element,
//   Ct = sum_e Phi_e^T C_e Phi_e with C_e = Element::getDamp() plus the
//   nodal alphaM term (LadrunoDampingAssembler), so region or element
//   Rayleigh factors, betaKinit/betaKcomm and damper materials are seen
//   exactly as the transient analysis sees them. The -closedForm path uses
//   the global Rayleigh factors only: since K phi_a = omega_a^2 M phi_a,
//   Ct = diag(alphaM + betaK omega_a^2) exactly in the mass-normalized
//   basis, with no matrix assembled.
//
//   Per mode: omega0 = |lambda|, omegaD = Im(lambda), zeta =
//   -Re(lambda)/|lambda| (Veletsos and Ventura 1986). The physical complex
//   shape psi = Phi z is stored on the nodes, the largest-|psi| entry made
//   real and positive, and recorded by the Node recorder as
//   complexEigenRe<k> / complexEigenIm<k>.

#ifndef LadrunoComplexEigen_h
#define LadrunoComplexEigen_h

#include <vector>

class Matrix;
class Domain;

class LadrunoComplexEigen
{
  public:
    enum ModeKind {
        UNDERDAMPED = 0,  // conjugate pair, reported once (+Im representative)
        OVERDAMPED  = 1,  // real non-oscillatory root, reported per root
        RIGID       = 2   // lambda ~ 0, or infinite eigenvalue (beta ~ 0)
    };

    struct ComplexMode {
        double omega0;    // undamped natural circular frequency  |lambda|
        double omegaD;    // damped circular frequency  Im(lambda) (>= 0)
        double zeta;      // modal damping ratio  -Re(lambda)/|lambda|
        double lambdaRe;  // Re(lambda)
        double lambdaIm;  // Im(lambda) of the reported (+Im) representative
        ModeKind kind;
        double resid;     // ||(lambda^2 Mt + lambda Ct + Kt) z||_2, quality metric
        std::vector<double> zRe;  // reduced eigenvector z (top block), real part
        std::vector<double> zIm;  //                                    imag part
    };

    // Solve the reduced quadratic pencil. Mt, Ct, Kt must be p x p with the
    // same p. Underdamped conjugate pairs are collapsed to one physical mode
    // (the +Im root); overdamped/rigid real roots are each reported. Modes
    // are returned sorted by ascending omega0 (rigid first). z is
    // phase-normalized: the largest-|z_j| entry is made real-positive.
    // tol scales the beta ~ 0 (infinite-eigenvalue) and Im ~ 0
    // (non-oscillatory) classification tests.
    // Returns 0 on success; < 0 on dimension mismatch or LAPACK failure
    // (dggev INFO is reported via opserr).
    static int solveReducedPencil(const Matrix &Mt, const Matrix &Ct,
                                  const Matrix &Kt,
                                  std::vector<ComplexMode> &modes,
                                  double tol = 1.0e-8);

    // Domain-coupled solve on the last `eigen` result. numModes <= 0 means
    // "all available real modes"; a positive numModes uses only the first
    // numModes of them.
    //
    // Default (closedForm = false): project M and C element by element via
    // LadrunoDampingAssembler. This captures per-element and region
    // Rayleigh factors, betaKinit/betaKcomm, damper materials (Viscous in a
    // zeroLength element, etc.) and nodal alphaM damping: the C the
    // transient analysis assembles.
    //
    // closedForm = true: the diagonal Rayleigh closed form
    // Ct = diag(alphaM + betaK w^2) from the global factors of the last
    // `rayleigh` command. Refuses betaKinit/betaKcomm (not proportional to
    // the eigen K) and warns when element or region factors differ from
    // the global ones.
    //
    // Both refuse when no prior eigen result is available or numModes
    // exceeds it (returns < 0 with an opserr message).
    static int solveFromDomain(Domain &theDomain, int numModes,
                               std::vector<ComplexMode> &modes,
                               double tol = 1.0e-8,
                               bool closedForm = false);

    // Distribute the complex (phased) mode shapes psi_k = Phi z_k onto the
    // nodes (Node::setComplexEigenvector), one column per reported mode
    // (underdamped pairs collapsed, so the count can differ from p). Called
    // by solveFromDomain on success.
    static int pushComplexShapes(Domain &theDomain, int p,
                                 const std::vector<ComplexMode> &modes);
};

#endif
