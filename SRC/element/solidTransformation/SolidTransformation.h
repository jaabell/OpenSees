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
// Description: SolidTransformation — the geometry-method layer for solid
// (continuum) elements. It is the solid-element analogue of OpenSees'
// beam-column CrdTransf / CorotCrdTransf (which exist for FRAMES only): a small
// interface that decouples an element's KINEMATICS (linear / corotational /
// finite strain) from its constitutive + integration CORE, so the core stays
// oblivious to which geometry method wraps it.
//
// AXIS ORTHOGONALITY. This is a *geometry-method* axis — linear ⊂ corot ⊂
// finite — orthogonal to the element *formulation* axis (std / bbar / uri /
// eas). The formulation (which OWNS B-matrix construction) runs unchanged inside
// whichever geometry method is selected.
//
// DIVISION OF LABOUR (the load-bearing design decision).
//   * The FORMULATION owns B-matrix construction and assembles the local
//     (core-frame) internal force f and stiffness K. The transformation NEVER
//     builds B.
//   * The TRANSFORMATION owns exactly four things:
//       - getStrainMeasure()                  — what strain the element feeds
//                                               the material core;
//       - localizeDisp()      (seam 2)        — global nodal disp → core frame;
//       - globalizeForce/Stiff (seam 3, +K_geo) — core f/K → global DOFs, plus
//                                               the geometric stiffness the
//                                               method implies;
//       - the recorder frame  (frameTimeVarying / getCurrentFrame).
//
// THE THREE METHODS.
//   linear  (v1) — identity. uCore = uGlobal; small engineering strain ε = B·u;
//                  no K_geo; frame fixed. Reuses the full small-strain material
//                  library and every formulation with ZERO behaviour change —
//                  routing a small-strain element through linear must reproduce
//                  the direct kernel to ~1e-12 (the seam proof).
//   corot   (v2) — large displacement + large ROTATION, small material strain.
//                  Strips the rigid rotation R (iterative Higham polar) so the
//                  core sees small strain in the de-rotated frame; adds the
//                  corotational K_geo (EICR). Still reuses the full small-strain
//                  material library.
//   finite  (v3) — large STRAIN. getStrainMeasure()==DeformationGradient: the
//                  element computes the full F and the material is the seam-3
//                  FiniteStrainNDMaterial adaptor (setTrialF(F) → Cauchy σ +
//                  spatial tangent). See SRC/material/nD/FiniteStrainNDMaterial.h.
//   hypo    (v4) — large STRAIN, rate-form (ADR 79). A MARKER like finite: the
//                  element integrates the Hughes–Winget midpoint objective
//                  strain increment in the unrotated (Green–Naghdi) material
//                  frame and drives an UNCHANGED small-strain material by
//                  accumulate-and-setTrialStrain — so getStrainMeasure() is
//                  honestly SmallStrain and all seam methods are identity.
//                  Per-GP machinery lives in the element + LadrunoHypoKernel.h.
//
// SERIALIZATION. The wrapper is NOT a MovableObject in v1. A consuming element
// serializes only the method id (getMethodID(); Linear=0) and rebuilds the
// transformation in recvSelf via SolidTransformation::create(methodID).
//
// The wrapper is ELEMENT-AGNOSTIC: it touches no classTags.h entry and knows
// nothing of any particular element. See Ladruno_implementation/
// solid_transformation_wrapper.md and 09_ladruno_brick.md.

#ifndef SolidTransformation_h
#define SolidTransformation_h

class Matrix;
class Vector;

class SolidTransformation
{
 public:
  // Which strain quantity the element computes and hands the material core.
  enum class StrainMeasure {
    SmallStrain,          // linear, corot — element feeds ε to setTrialStrain
    DeformationGradient   // finite       — element feeds F to setTrialF
  };

  // Method ids — stable across versions; serialized by the consuming element.
  static const int METHOD_LINEAR = 0;
  static const int METHOD_COROT  = 1;
  static const int METHOD_FINITE = 2;
  static const int METHOD_HYPO   = 3;   // Ladruno (ADR 79): rate-form UL marker

  SolidTransformation() {}
  virtual ~SolidTransformation() {}

  // ----- identity / bookkeeping ----------------------------------------- //
  virtual int getMethodID() const = 0;            // Linear=0, Corot=1, Finite=2
  virtual const char *getName() const = 0;        // "linear" / "corot" / "finite"
  virtual SolidTransformation *getCopy() const = 0;

  // ----- strain measure -------------------------------------------------- //
  virtual StrainMeasure getStrainMeasure() const = 0;

  // ----- (0) refresh geometry state, ONCE per evaluation, before the Gauss
  //           loop. refCrds / curCrds are (numNodes x 3): reference and current
  //           nodal coordinates. linear ignores them; corot/finite build the
  //           frame / deformation gradient from them.
  virtual int update(int numNodes, const Matrix &refCrds,
                     const Matrix &curCrds) = 0;

  // ----- (seam 2) localize global nodal displacements into the core frame.
  //           uGlobal / uCore are (numNodes*ndf). linear: identity. uCore MAY
  //           alias uGlobal (identity writes nothing).
  virtual int localizeDisp(const Vector &uGlobal, Vector &uCore) const = 0;

  // ----- (seam 3) globalize the core internal force / stiffness back to global
  //           DOFs, adding the geometric stiffness the method implies. linear:
  //           identity (+0). fCore is supplied to globalizeStiff because the
  //           corotational K_geo depends on the current internal force (EICR),
  //           not just on stress. Outputs MAY alias their core inputs.
  virtual int globalizeForce(const Vector &fCore, Vector &fGlobal) const = 0;
  virtual int globalizeStiff(const Matrix &kCore, const Vector &fCore,
                             Matrix &kGlobal) const = 0;

  // ----- recorder frame (element contract §B.2). Is the local frame
  //           time-varying, and what is the current 3x3 basis (columns = the
  //           local axes)? linear: false / identity.
  virtual bool frameTimeVarying() const = 0;
  virtual void getCurrentFrame(Matrix &R) const = 0;

  // ----- optional per-step state. corot/finite may carry a committed frame or
  //           birth gradient F0; linear is stateless, so these default to
  //           no-ops and a consuming element need not route them in v1.
  virtual int commitState() { return 0; }
  virtual int revertToLastCommit() { return 0; }
  virtual int revertToStart() { return 0; }

  // Factory used by a consuming element's recvSelf to rebuild a transformation
  // from a serialized method id. Defined in SolidTransformationLinear.cpp;
  // future methods register their case there. Returns 0 on an unknown id.
  static SolidTransformation *create(int methodID);
};

#endif
