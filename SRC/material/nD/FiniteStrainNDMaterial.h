/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// Description: Abstract base class for FINITE-STRAIN material adaptors — the
// "seam 3" contract between an F-delivering finite-strain element and an
// UNCHANGED small-strain NDMaterial.
//
// This class adds ONE thing to the NDMaterial interface: a deformation-gradient
// entry point, setTrialF(F). A concrete subclass (e.g. LogStrainNDMaterial) wraps
// an inner small-strain NDMaterial and lifts it to finite strain by the spatial
// multiplicative / logarithmic-strain technique of de Souza Neto, Perić & Owen,
// *Computational Methods for Plasticity* (2008), Box 14.3 (MATISU):
//
//     bᵉ_trial = F_Δ · bᵉ_n · F_Δᵀ                       (14.94)
//     εᵉ_trial = ½ ln[bᵉ_trial]   (Eulerian Hencky strain, 14.30/14.104)
//     → drive the UNCHANGED small-strain return map with εᵉ_trial as "strain"
//       (14.89–90 share the format of the infinitesimal scheme 7.25)
//     → inner material returns Kirchhoff stress τ and the small-strain tangent D
//     → Cauchy stress  σ = J⁻¹ τ        (J = det F)        (Box 14.3, iv)
//     → spatial consistent tangent  a   (§14.5, eq. 14.99)
//
// CONTRACT (locked jointly with the solid-transformation-wrapper team — the
// element side, LadrunoBrick + SolidTransformation::finite, classTag 33002):
//
//   * The ELEMENT calls setTrialF(F) with the TOTAL deformation gradient F (3×3)
//     at the Gauss point. The adaptor computes F_Δ = F·F_n⁻¹ internally — the
//     element does NOT pass increments (use setTrialFincr only if explicitly
//     supported). The element NEVER calls setTrialStrain(...) on this material.
//
//   * getStress()  returns the CAUCHY stress σ as a 6-vector in the canonical
//     3D Voigt order {00, 11, 22, 01, 12, 20} with ENGINEERING shear
//     (γ_ij = 2 ε_ij). The element assembles the SPATIAL internal force
//     ∫ bᵀ σ dv on the current configuration.
//
//   * The spatial CONSTITUTIVE modulus is c_ijkl = (1/2J)[D:L:B] — the MATERIAL
//     part only of the eq.(14.99) spatial tangent (NOT a PK2/reference tangent;
//     the formulation is spatial multiplicative, NOT total-Lagrangian).
//
//     SEAM-3 TANGENT CONTRACT (revised 2026-06-01 after the element FD-tangent
//     gate; supersedes the earlier "∫BᵀcB + ∫GᵀΣG" wording):
//       - c is NOT minor-symmetric in (k,l) (the factor B, eq.14.102, is not
//         k↔l symmetric), so it does NOT fit a 6×6 Voigt matrix. getTangent()
//         (6×6) is therefore the LOSSY minor-symmetric projection — usable for
//         output, NOT for a consistent element tangent.
//       - The element gets the FULL tensor via getSpatialTangentTensor(c), forms
//         the consistent spatial tangent  a_ijkl = c_ijkl − σ_il δ_jk  (the
//         −σ_il δ_jk geometric/initial-stress term is element-owned, built from
//         the Cauchy stress σ), and assembles
//             K_{(a,i)(b,k)} = ∫ (∂Nₐ/∂x_j) a_ijkl (∂N_b/∂x_l) dv  (full gradients).
//         a_ijkl yields a SYMMETRIC K only when the FULL non-symmetric c and the
//         FULL −σδ are both kept (they cancel) — collapsing c to 6×6 breaks it.
//       - The element still owns the geometric term; the material owns the
//         constitutive tensor. The arbiter is the element's FD-tangent test.
//
//   * getType() must report "ThreeDimensional"; getOrder() must return 6.
//
//   * STATE OWNERSHIP: the adaptor owns the committed elastic left Cauchy–Green
//     tensor bᵉ_n (equivalently εᵉ_n) per material instance (one per Gauss
//     point), and WRAPS the inner material's commitState / revertToLastCommit /
//     revertToStart / sendSelf / recvSelf. bᵉ_n is NOT a plastic multiplier and
//     NOT a total strain; commitState sets bᵉ_n ← exp[2 εᵉ_{n+1}].
//
//   * DEGENERACY GUARANTEE: update() must return finite, non-NaN σ and tangent
//     under repeated principal stretches (uniaxial → 2-fold, hydrostatic /
//     equibiaxial → 3-fold, any axisymmetric strain). The 0/0 term
//     (ln bᵢ − ln bⱼ)/(bᵢ − bⱼ) in the log-map tangent must use the analytic
//     limit (Miehe 1998; dSNPO App. A.5), not ε-perturbation.
//
// See Ladruno_implementation/09_finite_strain_material_wrapper.md for the full
// plan, the GREEN/YELLOW/RED inner-material matrix, and the theory.

#ifndef FiniteStrainNDMaterial_h
#define FiniteStrainNDMaterial_h

#include <NDMaterial.h>
#include <OPS_Globals.h>

class Matrix;
class Vector;

class FiniteStrainNDMaterial : public NDMaterial
{
 public:
  FiniteStrainNDMaterial(int tag, int classTag) : NDMaterial(tag, classTag) {}
  FiniteStrainNDMaterial() : NDMaterial() {}
  virtual ~FiniteStrainNDMaterial() {}

  // ===================================================================== //
  //  THE FINITE-STRAIN SEAM (seam 3)                                      //
  // ===================================================================== //

  // Set the trial state from the TOTAL deformation gradient F (3×3) at the GP.
  // Returns 0 on success, < 0 on failure (e.g. det F ≤ 0, non-convergence).
  // This is THE entry point a finite-strain element must call.
  virtual int setTrialF(const Matrix &F) = 0;

  // Optional increment form: F_Δ maps config n onto config n+1 (F_{n+1}=F_Δ F_n).
  // Default: unsupported — a finite-strain element should pass total F.
  virtual int setTrialFincr(const Matrix &Fincr) {
    opserr << "FiniteStrainNDMaterial::setTrialFincr() not supported by this "
              "material; pass the total F via setTrialF().\n";
    return -1;
  }

  // ===================================================================== //
  //  SPATIAL RESPONSE — what the element consumes                         //
  // ===================================================================== //

  // CAUCHY stress σ (6-vector, Voigt {00,11,22,01,12,20}, engineering shear).
  virtual const Vector &getStress(void) = 0;
  // SPATIAL CONSTITUTIVE modulus c = (1/2J)[D:L:B] as a 6×6. WARNING: c is NOT
  // minor-symmetric in its second index pair (k,l), so this 6×6 is the LOSSY
  // minor-symmetric projection — it CANNOT carry the full modulus. Kept for
  // small-strain-style consumers/output only. A finite-strain ELEMENT building a
  // consistent tangent MUST use getSpatialTangentTensor() below.  // Ladruno
  virtual const Matrix &getTangent(void) = 0;

  // FULL 4th-order spatial constitutive modulus c_ijkl = (1/2J)[D:L:B], all 81
  // components (NOT minor-symmetric in (k,l) — see above). This is the element's
  // tangent channel: the element forms the consistent spatial tangent
  //     a_ijkl = c_ijkl − σ_il δ_jk
  // (the −σ_il δ_jk is the element-owned geometric/initial-stress term) and
  // assembles K_{(a,i)(b,k)} = ∫ (∂Nₐ/∂x_j) a_ijkl (∂N_b/∂x_l) dv. Returns 0 on
  // success; the default −1 means "unsupported" so element misuse fails loudly.
  // (The 6×6 getTangent() above is insufficient because it loses the (k,l)
  // antisymmetric content on which the K-symmetric cancellation with −σδ relies.)
  virtual int getSpatialTangentTensor(double c[3][3][3][3]) { return -1; }

  // Current Jacobian J = det F (> 0). τ = J·σ recovers Kirchhoff if needed.
  virtual double getJ(void) const = 0;

  // ===================================================================== //
  //  The small-strain entry point is DISABLED for finite-strain materials //
  // ===================================================================== //

  // A finite-strain material is driven by F, not by a strain vector. These
  // default to an error so element misuse is caught loudly rather than silently
  // feeding a linear strain into the multiplicative machinery.
  virtual int setTrialStrain(const Vector & /*v*/) override {
    opserr << "FiniteStrainNDMaterial: setTrialStrain() is invalid — this is a "
              "finite-strain material; the element must call setTrialF(F).\n";
    return -1;
  }
  virtual int setTrialStrain(const Vector & /*v*/, const Vector & /*r*/) override {
    opserr << "FiniteStrainNDMaterial: setTrialStrain() is invalid — this is a "
              "finite-strain material; the element must call setTrialF(F).\n";
    return -1;
  }

  // NOTE: getType() / getOrder() / commitState() / revertToLastCommit() /
  // revertToStart() / getCopy() / sendSelf() / recvSelf() remain pure-virtual
  // (from NDMaterial) — the concrete adaptor implements them, owning bᵉ_n and
  // wrapping the inner material's state. getStrain() may return the Hencky strain.

 protected:

 private:
};

#endif
