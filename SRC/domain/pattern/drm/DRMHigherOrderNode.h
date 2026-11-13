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

// Ladruno (ADR-88): optional interface a solid Element MAY implement to tell
// H5DRMLoadPattern how to reconstruct the free-field motion at one of its own
// SECONDARY (non-corner) local nodes -- e.g. the mid-edge node of a quadratic
// serendipity/Bezier element such as BezierTet10 or LadrunoBrick20 -- from a
// linear combination of its own PRIMARY (corner) local nodes, instead of
// requiring a direct match against a real H5DRM station (which the ShakerMaker
// grid never provides for a mid-edge point).
//
// This is exact, not an approximation, whenever the true free-field motion
// varies linearly across the element -- the same assumption a straight-sided
// quadratic element (BezierTet10 v1) already makes about its own geometry, and
// consistent with sizing the mesh to resolve the wavelength of interest. See
// the design rationale in the CERN ATLAS DRM project,
// notebooks/01_calibration_drm_soil/drm_load_pattern/
// 86_ladruno_h5drm_higher_order_elements_adr_pxp.md.
//
// Elements that do NOT implement this interface behave exactly as before this
// ADR: every one of their nodes must match a real H5DRM station directly, or
// the element is excluded from the DRM boundary set.
#ifndef DRMHigherOrderNode_h
#define DRMHigherOrderNode_h

#include <vector>

class DRMHigherOrderNode
{
public:
    virtual ~DRMHigherOrderNode() {}

    // localNode: 0-based local node index, in the SAME order as
    // Element::getExternalNodes(). Returns false if localNode is a PRIMARY
    // node (no interpolation -- H5DRMLoadPattern must match it to a real
    // station, exactly as before this interface existed). Returns true if
    // localNode is a SECONDARY node, filling primaryLocalNodes (0-based local
    // indices of this same element's PRIMARY nodes) and weights (same length,
    // summing to 1.0) with the linear combination used to reconstruct
    // displacement/acceleration at that secondary node from its parents' real
    // station data.
    virtual bool getDRMInterpolation(int localNode,
                                      std::vector<int>& primaryLocalNodes,
                                      std::vector<double>& weights) const = 0;
};

#endif // DRMHigherOrderNode_h
