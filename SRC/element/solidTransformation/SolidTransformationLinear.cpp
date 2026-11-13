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
// SolidTransformationLinear — identity geometry method (v1). See the header.

#include <SolidTransformationLinear.h>
#include <Matrix.h>
#include <Vector.h>

SolidTransformationLinear::SolidTransformationLinear()
  : SolidTransformation()
{
}

SolidTransformationLinear::~SolidTransformationLinear()
{
}

SolidTransformation *
SolidTransformationLinear::getCopy() const
{
  return new SolidTransformationLinear();
}

// Identity: there is no geometry state to refresh. The reference / current
// coordinates are unused (corot/finite consume them to build the frame / F).
int
SolidTransformationLinear::update(int /*numNodes*/, const Matrix & /*refCrds*/,
                                  const Matrix & /*curCrds*/)
{
  return 0;
}

// Identity: the core frame IS the global frame — no rigid rotation is stripped.
// Guard against aliasing so the element may legitimately pass uCore == uGlobal.
int
SolidTransformationLinear::localizeDisp(const Vector &uGlobal, Vector &uCore) const
{
  if (&uCore != &uGlobal)
    uCore = uGlobal;
  return 0;
}

// Identity: fGlobal = fCore (no transform, no geometric force contribution).
int
SolidTransformationLinear::globalizeForce(const Vector &fCore, Vector &fGlobal) const
{
  if (&fGlobal != &fCore)
    fGlobal = fCore;
  return 0;
}

// Identity: kGlobal = kCore (no transform, K_geo = 0). fCore is ignored; corot
// uses it for the EICR geometric stiffness.
int
SolidTransformationLinear::globalizeStiff(const Matrix &kCore, const Vector & /*fCore*/,
                                          Matrix &kGlobal) const
{
  if (&kGlobal != &kCore)
    kGlobal = kCore;
  return 0;
}

// Identity: the local frame is the global axes, fixed for all time.
void
SolidTransformationLinear::getCurrentFrame(Matrix &R) const
{
  R.Zero();
  R(0, 0) = 1.0;
  R(1, 1) = 1.0;
  R(2, 2) = 1.0;
}

// ---------------------------------------------------------------------------
// Factory (declared on the base) — rebuilds a transformation from a serialized
// method id. v1 ships only the identity; corot/finite register their cases here
// as they land.
// ---------------------------------------------------------------------------
SolidTransformation *
SolidTransformation::create(int methodID)
{
  switch (methodID) {
  case METHOD_LINEAR:
    return new SolidTransformationLinear();
  // jaabell/OpenSees port: only the small-strain linear method is brought
  // over; corot, finite and hypo are not, and create() returns 0 for them.
  default:
    return 0;
  }
}
