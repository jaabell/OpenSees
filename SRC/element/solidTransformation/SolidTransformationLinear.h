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
// Description: SolidTransformationLinear — the IDENTITY geometry method (v1).
//
// Every seam is a pass-through:
//   getStrainMeasure() == SmallStrain   (element computes ε = B·u directly)
//   update()           — does nothing (reference coords are ignored)
//   localizeDisp()      — uCore = uGlobal   (no rigid rotation stripped)
//   globalizeForce()    — fGlobal = fCore   (no transform)
//   globalizeStiff()    — kGlobal = kCore   (no geometric stiffness, K_geo = 0)
//   frameTimeVarying() == false; getCurrentFrame() == 3x3 identity
//
// Routing a small-strain element through this transformation must reproduce the
// direct-kernel result to ~1e-12 — that bit-for-bit equality IS the proof that
// the seam is wired correctly (so corot/finite can later drop in on the same
// interface with no element rewrite).

#ifndef SolidTransformationLinear_h
#define SolidTransformationLinear_h

#include <SolidTransformation.h>

class SolidTransformationLinear : public SolidTransformation
{
 public:
  SolidTransformationLinear();
  virtual ~SolidTransformationLinear();

  int getMethodID() const { return METHOD_LINEAR; }
  const char *getName() const { return "linear"; }
  SolidTransformation *getCopy() const;

  StrainMeasure getStrainMeasure() const { return StrainMeasure::SmallStrain; }

  int update(int numNodes, const Matrix &refCrds, const Matrix &curCrds);

  int localizeDisp(const Vector &uGlobal, Vector &uCore) const;

  int globalizeForce(const Vector &fCore, Vector &fGlobal) const;
  int globalizeStiff(const Matrix &kCore, const Vector &fCore,
                     Matrix &kGlobal) const;

  bool frameTimeVarying() const { return false; }
  void getCurrentFrame(Matrix &R) const;
};

#endif
