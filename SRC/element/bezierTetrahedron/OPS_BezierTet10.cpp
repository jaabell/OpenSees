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

// Authors: Nicolas Mora Bowen, Patricio Palacios, Jose Abell, Guppi (Ladruño)
// Created: 05/2026
//
// Factory function for the BezierTet10 element.
// Parses user input from both Tcl and Python interfaces.
//
// Usage:
//   element BezierTet10 $tag $nd1 $nd2 $nd3 $nd4 $nd5 $nd6 $nd7 $nd8 $nd9 $nd10
//                       $matTag
//                       <-bbar> <-cMass> <-rho $r>
//                       <-bodyForce $b1 $b2 $b3> <-pressure $p>
//                       <-geom linear|corot|finite>
//
// Required arguments:
//   $tag        - unique element tag
//   $nd1-$nd10  - node tags (4 corners + 6 mid-edge), TenNodeTetrahedron order:
//                 v1 v2 v3 v4, e12 e23 e13 e14 e34 e24
//   $matTag     - 3D NDMaterial tag
//
// Optional arguments:
//   -bbar               - use B-bar formulation (near-incompressibility)
//   -cMass              - consistent mass (default is lumped ρVe/10, Eq. 57)
//   -rho $r             - mass density (else taken from the material)
//   -bodyForce $b1..$b3 - body force per unit volume (rampable via SelfWeight)
//   -pressure $p        - volume "pressure" hack acting in +z (as BezierTri6)
//   -geom linear|corot|finite - geometry method (default linear = small strain;
//                         corot = large rotation / small strain, EICR via the
//                         SolidTransformation layer; finite = large-strain
//                         updated-Lagrangian, requires a finite-strain material
//                         driven by setTrialF(F), e.g. nDMaterial LogStrain).
//                         -geom finite -bbar = F-bar (near-incompressible cure;
//                         GENERALLY UNSYMMETRIC tangent → use system FullGeneral).
//                         pressure is unsupported under finite (rejected at parse).
//   -fbar centroid|mean_dilatation - F-bar variant (default centroid). Only with
//                         -bbar -geom finite. centroid = single centroid J₀ (dSNPO
//                         eq 15.5); mean_dilatation = volume-averaged J̄ (reduces to
//                         the small-strain mean-dilatation bbar).
//
// Example (Python):
//   ops.nDMaterial('ElasticIsotropic', 1, 1000.0, 0.3)
//   ops.element('BezierTet10', 1, 1,2,3,4,5,6,7,8,9,10, 1, '-bbar')
//
// LOADS AND BCs — CONTROL-VALUE SEMANTICS (read before porting a tet10 deck):
//   Element DOFs are Bernstein CONTROL values. At the four VERTICES the
//   quadratic Bernstein basis is interpolatory, so vertex DOFs, vertex point
//   loads and vertex BCs mean exactly what they mean on TenNodeTetrahedron.
//   At the six MID-EDGE nodes they do NOT: a nodal force there is conjugate to
//   the (non-interpolatory) mid-edge Bernstein function, and a nonzero sp() is
//   a control value, not the displacement at that point (homogeneous fixes are
//   equivalent — zero control values = zero field).
//   Consequence for surface tractions: LAGRANGE-consistent nodal forces (the
//   tet10 rule — vertices ~0, midsides ~q·A/3) applied to this element
//   represent an OSCILLATORY traction, exact in resultant but locally spiking;
//   elastically benign, it can silently drive surface Gauss points into
//   yield/apex regimes the intended traction never reaches (TIMs T2, 2026-08).
//   The BERNSTEIN-consistent rule is simpler: every quadratic Bernstein face
//   function integrates to A/6, so a uniform traction q on a face is q·A/6 on
//   EACH of its six nodes (general t: f_a = ∫ t·B_a dΓ, Bernstein B_a).

#include "BezierTet10.h"

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>
#include <FiniteStrainNDMaterial.h>  // Ladruno — -geom finite material guards
#include <SolidTransformation.h>   // Ladruno — geometry-method ids (linear/corot/finite)

#include <string.h>

// Built into OpenSees statically (registered in OpenSeesElementCommands.cpp via
// `void* OPS_BezierTet10();`), so use plain C++ linkage — NOT the extern "C"
// OPS_Export DLL style, which would not resolve against that declaration.
void *OPS_BezierTet10()
{
    // Required: tag nd1..nd10 matTag = 12 remaining args (eleType consumed)
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();
    if (numRemainingArgs < 12) {
        opserr << "WARNING insufficient arguments\n"
               << "Want: element BezierTet10 tag nd1 nd2 nd3 nd4 nd5 nd6 nd7 "
                  "nd8 nd9 nd10 matTag <-bbar> <-cMass> <-rho r> "
                  "<-bodyForce b1 b2 b3> <-pressure p>\n";
        return 0;
    }

    // ─── Parse integer data: tag + 10 node tags ───────────────
    int iData[11];
    int numData = 11;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data for BezierTet10\n";
        return 0;
    }

    int tag = iData[0];
    bool massCacheOn = true;   // Ladruno (ADR-77 G2 ext): default on, guard-checked

    // ─── Parse material tag ───────────────────────────────────
    int matTag;
    numData = 1;
    if (OPS_GetIntInput(&numData, &matTag) != 0) {
        opserr << "WARNING invalid matTag for BezierTet10 " << tag << "\n";
        return 0;
    }

    NDMaterial *theMat = OPS_getNDMaterial(matTag);
    if (theMat == 0) {
        opserr << "WARNING NDMaterial " << matTag
               << " not found for BezierTet10 " << tag << "\n";
        return 0;
    }

    // ─── Parse optional arguments ─────────────────────────────
    double rho = 0.0;
    double pressure = 0.0;
    double b1 = 0.0, b2 = 0.0, b3 = 0.0;
    bool useBbar = false;
    bool cMass = false;
    int geomMethodID = SolidTransformation::METHOD_LINEAR;   // -geom (Ladruno)
    int fbarMode = BezierTet10::FBAR_CENTROID;               // -fbar (Ladruno)

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *option = OPS_GetString();

        if (strcmp(option, "-bbar") == 0 || strcmp(option, "-Bbar") == 0) {
            useBbar = true;
        }
        else if (strcmp(option, "-fbar") == 0 || strcmp(option, "-fBar") == 0) {
            // Ladruno: F-bar variant (only meaningful with -bbar -geom finite).
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING -fbar needs a value for BezierTet10 " << tag
                       << " (use centroid|mean_dilatation)\n";
                return 0;
            }
            const char *fb = OPS_GetString();
            if (strcmp(fb, "centroid") == 0)
                fbarMode = BezierTet10::FBAR_CENTROID;
            else if (strcmp(fb, "mean_dilatation") == 0 ||
                     strcmp(fb, "mean") == 0 || strcmp(fb, "meandilatation") == 0)
                fbarMode = BezierTet10::FBAR_MEAN;
            else {
                opserr << "WARNING unknown -fbar '" << fb << "' for BezierTet10 "
                       << tag << " (use centroid|mean_dilatation)\n";
                return 0;
            }
        }
        else if (strcmp(option, "-cMass") == 0 || strcmp(option, "-cmass") == 0) {
            cMass = true;
        }
        else if (strcmp(option, "-rho") == 0 || strcmp(option, "-mass") == 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &rho) != 0) {
                opserr << "WARNING invalid rho for BezierTet10 " << tag << "\n";
                return 0;
            }
        }
        else if (strcmp(option, "-pressure") == 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &pressure) != 0) {
                opserr << "WARNING invalid pressure for BezierTet10 " << tag << "\n";
                return 0;
            }
        }
        else if (strcmp(option, "-bodyForce") == 0 || strcmp(option, "-body") == 0) {
            double bData[3];
            numData = 3;
            if (OPS_GetDoubleInput(&numData, bData) != 0) {
                opserr << "WARNING invalid body force for BezierTet10 " << tag << "\n";
                return 0;
            }
            b1 = bData[0];
            b2 = bData[1];
            b3 = bData[2];
        }
        else if (strcmp(option, "-noMassCache") == 0) {
            massCacheOn = false;   // Ladruno (ADR-77 G2 ext): A/B escape
        }
        else if (strcmp(option, "-geom") == 0 || strcmp(option, "-geometry") == 0) {
            // Ladruno: geometry method. linear (default), corot, or finite.
            // finite = large-strain updated-Lagrangian; the material guards
            // below enforce that it is paired with a FiniteStrainNDMaterial.
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING -geom needs a value for BezierTet10 " << tag
                       << " (use linear|corot|finite)\n";
                return 0;
            }
            const char *g = OPS_GetString();
            if (strcmp(g, "linear") == 0)
                geomMethodID = SolidTransformation::METHOD_LINEAR;
            else if (strcmp(g, "corot") == 0 || strcmp(g, "corotational") == 0
                     || strcmp(g, "finite") == 0) {
                // jaabell/OpenSees port: only the small-strain linear geometry
                // was brought over from Ladruno (no corot/finite transformations)
                opserr << "WARNING -geom " << g << " is not available in this build of "
                       << "BezierTet10 " << tag << " (only -geom linear)\n";
                return 0;
            }
            else {
                opserr << "WARNING unknown -geom '" << g << "' for BezierTet10 "
                       << tag << " (use linear|corot|finite)\n";
                return 0;
            }
        }
        else {
            opserr << "WARNING unknown option '" << option
                   << "' for BezierTet10 " << tag << "\n";
        }
    }

    // Ladruno: v1 corot does not yet support the +z pressure hack (its
    // fixed-direction-vs-co-rotating-load behaviour under large rotation is
    // unvalidated). Reject the combination at parse time.
    if (geomMethodID == SolidTransformation::METHOD_COROT && pressure != 0.0) {
        opserr << "WARNING BezierTet10 " << tag
               << ": -geom corot does not support -pressure in v1\n";
        return 0;
    }

    // Ladruno: -geom finite (updated-Lagrangian) guards. std (plain F) or bbar
    // (F-bar, dSNPO §15.1) are both valid; the +z pressure hack is unvalidated
    // under finite. It requires a finite-strain material driven by setTrialF(F)
    // (e.g. nDMaterial LogStrain). Mirrors OPS_LadrunoBrick.cpp.
    if (geomMethodID == SolidTransformation::METHOD_FINITE) {
        if (pressure != 0.0) {
            opserr << "WARNING BezierTet10 " << tag
                   << ": -geom finite does not support -pressure in v1\n";
            return 0;
        }
        if (dynamic_cast<FiniteStrainNDMaterial *>(theMat) == 0) {
            opserr << "WARNING BezierTet10 " << tag
                   << ": -geom finite requires a finite-strain NDMaterial driven "
                      "by setTrialF(F) (e.g. nDMaterial LogStrain); material "
                   << matTag << " is not one\n";
            return 0;
        }
        // F-bar (bbar + finite) has a GENERALLY UNSYMMETRIC tangent (dSNPO eq
        // 15.10); a symmetric solver silently drops the coupling and breaks Newton.
        // Advise once per process (not per element) to keep large meshes quiet.
        if (useBbar) {
            static bool fbarSolverAdvised = false;
            if (!fbarSolverAdvised) {
                fbarSolverAdvised = true;
                opserr << "BezierTet10: -geom finite -bbar (F-bar) has a GENERALLY "
                          "UNSYMMETRIC tangent (dSNPO eq 15.10) — use an unsymmetric "
                          "solver (e.g. 'system FullGeneral', 'UmfPack', or "
                          "'SparseGEN'); a symmetric system (BandSPD/ProfileSPD) drops "
                          "the coupling term and may not converge.\n";
            }
        }
        else if (fbarMode != BezierTet10::FBAR_CENTROID) {
            // -fbar without -bbar is a no-op (F-bar is off); warn so it is not
            // silently ignored.
            opserr << "WARNING BezierTet10 " << tag
                   << ": -fbar has no effect without -bbar (F-bar is off)\n";
        }
    }

    // Ladruno: converse guard — a finite-strain NDMaterial (e.g. LogStrain) is
    // driven ONLY by setTrialF(F). Under -geom linear|corot the element uses the
    // small-strain setTrialStrain() path, which a FiniteStrainNDMaterial disables
    // (returns -1, sets no state) — leaving a zero-stress phantom. Reject it.
    if (geomMethodID != SolidTransformation::METHOD_FINITE &&
        dynamic_cast<FiniteStrainNDMaterial *>(theMat) != 0) {
        opserr << "WARNING BezierTet10 " << tag
               << ": a finite-strain NDMaterial (e.g. nDMaterial LogStrain) "
                  "requires -geom finite; it cannot be driven by -geom "
                  "linear|corot\n";
        return 0;
    }

    // ─── Create element ───────────────────────────────────────
    Element *theEle = new BezierTet10(tag,
                                      iData[1], iData[2], iData[3], iData[4],
                                      iData[5], iData[6], iData[7], iData[8],
                                      iData[9], iData[10],
                                      *theMat, rho, b1, b2, b3,
                                      useBbar, cMass, pressure, geomMethodID, fbarMode);

    if (theEle == 0) {
        opserr << "WARNING ran out of memory creating BezierTet10 " << tag << "\n";
        return 0;
    }

    ((BezierTet10 *)theEle)->setMassCache(massCacheOn);   // Ladruno (ADR-77 G2 ext): transient, unserialized
    return theEle;
}
