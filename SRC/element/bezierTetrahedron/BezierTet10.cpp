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
// Description: 10-node quadratic Bézier tetrahedral element for 3D analysis.
//
// Reference:
//   Kadapa, C. "Novel quadratic Bézier triangular and tetrahedral elements
//   using existing mesh generators." Int. J. Numer. Methods Eng.,
//   2019; 117(5):543-573. doi:10.1002/nme.5967 (§5)
//
// ═══════════════════════════════════════════════════════════════════
//  MATHEMATICAL FOUNDATION
// ═══════════════════════════════════════════════════════════════════
//
//  Shape functions: Quadratic Bernstein polynomials on the tetrahedron
//  with barycentric coordinates (L1,L2,L3,L4), L4 = 1-L1-L2-L3.
//
//    N1=L1²  N2=L2²  N3=L3²  N4=L4²              (4 vertices)
//    N5=2L1L2  N6=2L2L3  N7=2L1L3                (mid-edges 1-2, 2-3, 1-3)
//    N8=2L1L4  N9=2L3L4  N10=2L2L4               (mid-edges 1-4, 3-4, 2-4)
//
//  Node order matches OpenSees TenNodeTetrahedron (jaabell/Larenas), so
//  Lagrange T10 / Gmsh meshes feed directly with only the basis swapped.
//
//  Key properties:
//    - Partition of unity: Σ Nₐ = 1
//    - Nonnegativity: Nₐ ≥ 0  ⇒  all-positive lumped mass (explicit dynamics)
//    - Equal integration: ∫ Nₐ dΩ = Vₑ/10 for all a (vertices and mid-edges)
//
//  Voigt order {xx, yy, zz, xy, yz, zx} == Kadapa Eq. 26-27 == TenNodeTet.
//  Isoparametric map uses CONTROL POINTS: x(ξ) = Σ Nₐ(ξ)·Pₐ (P = X if straight).
// ═══════════════════════════════════════════════════════════════════

#include "BezierTet10.h"

#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <LadrunoElementShell.h>   // Ladruno (WP-124): shared Element-contract helpers
#include <Vector.h>
#include <ID.h>
#include <Domain.h>
#include <OPS_Globals.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Renderer.h>
#include <SolidTransformation.h>         // Ladruno — geometry-method layer
#include <SolidTransformationLinear.h>   // Ladruno — default/fallback (identity)
#include <FiniteStrainNDMaterial.h>      // Ladruno — -geom finite: setTrialF(F) seam

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// ═══════════════════════════════════════════════════════════════════
//  STATIC DATA
// ═══════════════════════════════════════════════════════════════════

int BezierTet10::numInstances = 0;

Matrix BezierTet10::K_return(NELD, NELD);
Matrix BezierTet10::M_return(NELD, NELD);
Vector BezierTet10::P_return(NELD);

// ─── 4-Point Tetrahedral Gauss Quadrature (degree 2) ──────────
// Symmetric rule; alpha/beta == TenNodeTetrahedron. Each weight = 1/24,
// so the weights sum to 1/6 (the reference tetrahedron volume).
//   a = (5 + 3√5)/20 ≈ 0.585410,  b = (5 - √5)/20 ≈ 0.138197
// Stored as the 3 free barycentric coords (L1,L2,L3); L4 = 1-L1-L2-L3.
const double BezierTet10::GP4_L[][3] = {
    {0.585410196624968, 0.138196601125011, 0.138196601125011},
    {0.138196601125011, 0.585410196624968, 0.138196601125011},
    {0.138196601125011, 0.138196601125011, 0.585410196624968},
    {0.138196601125011, 0.138196601125011, 0.138196601125011}
};
const double BezierTet10::GP4_w[] = {
    1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0
};

// ─── 1D Gauss-Legendre (4-pt) on [0,1] ────────────────────────
// Collapsed (Duffy) into a 4×4×4 rule for the consistent mass: exact to
// degree 7, which covers the degree-6 integrand the (1-a)² Duffy Jacobian
// produces from the degree-4 product NₐNᵦ (ADR O12).
const double BezierTet10::GL4_t[] = {
    0.069431844202974, 0.330009478207572,
    0.669990521792428, 0.930568155797026
};
const double BezierTet10::GL4_w[] = {
    0.173927422568727, 0.326072577431273,
    0.326072577431273, 0.173927422568727
};

// ─── Edge → (vertexA, vertexB) for mid-edge nodes 5..10 (idx 4..9) ──
// Matches TenNodeTetrahedron edge order (jaabell/Larenas N9↔N10 swap).
const int BezierTet10::edgeV[6][2] = {
    {0, 1},   // node 5  : edge (1-2)
    {1, 2},   // node 6  : edge (2-3)
    {0, 2},   // node 7  : edge (1-3)
    {0, 3},   // node 8  : edge (1-4)
    {2, 3},   // node 9  : edge (3-4)
    {1, 3}    // node 10 : edge (2-4)
};

// Ladruno (ADR-88): H5DRM free-field interpolation -- localNode 0..3 are the 4
// corners (PRIMARY, no interpolation); localNode 4..9 are the 6 mid-edge nodes
// (SECONDARY), reusing the SAME edgeV table above verbatim -- no separate
// topology data to keep in sync. Weight 0.5/0.5 is exact under the element's
// own v1 straight-sided assumption (see the ADR-88 design doc, 10ter, on why
// this does not conflict with the Bernstein/blossom non-interpolatory basis).
bool BezierTet10::getDRMInterpolation(int localNode,
                                       std::vector<int>& primaryLocalNodes,
                                       std::vector<double>& weights) const
{
    if (localNode < 4 || localNode > 9)
        return false; // corner node -- must match a real H5DRM station directly

    const int edgeIdx = localNode - 4;
    primaryLocalNodes.assign({edgeV[edgeIdx][0], edgeV[edgeIdx][1]});
    weights.assign({0.5, 0.5});
    return true;
}


// ═══════════════════════════════════════════════════════════════════
//  CONSTRUCTORS AND DESTRUCTOR
// ═══════════════════════════════════════════════════════════════════

BezierTet10::BezierTet10(int tag,
                         int nd1, int nd2, int nd3, int nd4, int nd5,
                         int nd6, int nd7, int nd8, int nd9, int nd10,
                         NDMaterial &m, double r,
                         double bx, double by, double bz,
                         bool bbar, bool cmass, double press,
                         int geomMethodID, int fbar)
  : Element(tag, ELE_TAG_BezierTet10),
    theMaterial(0),
    connectedExternalNodes(NEN),
    rho(r), pressure(press),
    useBbar(bbar), cMass(cmass), fbarMode(fbar),   // Ladruno — F-bar variant
    Q(NELD), applyLoad(0), Ki(0), theGeom(0)   // Ladruno — geometry method
{
    // Ladruno: geometry-method layer (linear default / corot). Linear is the
    // identity wrapper so std/bbar run bit-for-bit unchanged.
    theGeom = SolidTransformation::create(geomMethodID);
    if (theGeom == 0)
        theGeom = new SolidTransformationLinear();   // safe fallback (unknown id)

    // One-time attribution print (Abell pattern) — Ladruño banner + authors
    if (numInstances == 0) {
        numInstances++;
        opserr << "\n"
" ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄\n"
"███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███\n"
"███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███\n"
"███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███\n"
"███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███\n"
"███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███\n"
"███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███\n"
"█████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀\n"
"▀                                   ███    ███\n";
        opserr << "BezierTet10 - Quadratic Bezier Tetrahedron\n"
               << "  Authors: Nicolas Mora Bowen, Patricio Palacios, "
                  "Jose Abell, Guppi (Ladruño)\n"
               << "  Ref: Kadapa, IJNME 2019; 117(5):543-573. "
                  "doi:10.1002/nme.5967\n";
    }

    // Store body forces
    b[0] = bx; b[1] = by; b[2] = bz;
    appliedB[0] = appliedB[1] = appliedB[2] = 0.0;

    // Store connectivity
    int nd[NEN] = {nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8, nd9, nd10};
    for (int i = 0; i < NEN; i++)
        connectedExternalNodes(i) = nd[i];

    for (int i = 0; i < NEN; i++) {
        theNodes[i] = 0;
        controlPts[i][0] = controlPts[i][1] = controlPts[i][2] = 0.0;
    }

    // Allocate 3D materials at the NGAUSS Gauss points
    theMaterial = new NDMaterial*[NGAUSS];
    for (int i = 0; i < NGAUSS; i++) {
        theMaterial[i] = m.getCopy("ThreeDimensional");
        if (theMaterial[i] == 0) {
            opserr << "BezierTet10::BezierTet10 - failed to get a 3D copy "
                      "of the material\n";
            exit(-1);
        }
    }
}


BezierTet10::BezierTet10()
  : Element(0, ELE_TAG_BezierTet10),
    theMaterial(0),
    connectedExternalNodes(NEN),
    rho(0.0), pressure(0.0),
    useBbar(false), cMass(false), fbarMode(FBAR_CENTROID),   // Ladruno — F-bar variant
    Q(NELD), applyLoad(0), Ki(0), theGeom(0)   // Ladruno — geometry method
{
    // Ladruno: broker/database reconstruction path — recvSelf rebuilds the real
    // method from the serialized id; default to the identity wrapper meanwhile.
    theGeom = new SolidTransformationLinear();

    b[0] = b[1] = b[2] = 0.0;
    appliedB[0] = appliedB[1] = appliedB[2] = 0.0;

    for (int i = 0; i < NEN; i++) {
        theNodes[i] = 0;
        controlPts[i][0] = controlPts[i][1] = controlPts[i][2] = 0.0;
    }

}


BezierTet10::~BezierTet10()
{
    if (theMaterial != 0) {
        for (int i = 0; i < NGAUSS; i++)
            if (theMaterial[i])
                delete theMaterial[i];
        delete[] theMaterial;
    }
    if (Ki != 0)
        delete Ki;
    if (theGeom != 0)        // Ladruno — geometry-method layer
        delete theGeom;
}


// ═══════════════════════════════════════════════════════════════════
//  INFORMATION METHODS
// ═══════════════════════════════════════════════════════════════════

int BezierTet10::getNumExternalNodes(void) const { return NEN; }
const ID &BezierTet10::getExternalNodes(void) { return connectedExternalNodes; }
Node **BezierTet10::getNodePtrs(void) { return theNodes; }
int BezierTet10::getNumDOF(void) { return NELD; }


// ═══════════════════════════════════════════════════════════════════
//  LIFECYCLE METHODS
// ═══════════════════════════════════════════════════════════════════

void BezierTet10::setDomain(Domain *theDomain)
{
    if (theDomain == 0) {
        for (int i = 0; i < NEN; i++)
            theNodes[i] = 0;
        return;
    }

    for (int i = 0; i < NEN; i++) {
        theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
        if (theNodes[i] == 0) {
            opserr << "BezierTet10::setDomain -- node "
                   << connectedExternalNodes(i) << " does not exist\n";
            return;
        }
        if (theNodes[i]->getNumberDOF() != NDOF) {
            opserr << "BezierTet10::setDomain -- node "
                   << connectedExternalNodes(i)
                   << " has " << theNodes[i]->getNumberDOF()
                   << " DOFs, expected " << NDOF << "\n";
            return;
        }
    }

    computeControlPoints();

    // ─── v1 straight-sided guard (ADR 06_bezier_tet10 D9′) ────────
    //  v1 is validated for STRAIGHT-SIDED elements only: all 6 mid-edge
    //  nodes at the edge midpoints ⇒ control points = nodes (P = X).
    //  A curved mesh breaks the interpolatory Dirichlet BCs and under-
    //  integrates BᵀDB. Warn; curved support + Eq. 14 is a follow-up.
    for (int e = 0; e < 6; e++) {
        const Vector &Xc1 = theNodes[edgeV[e][0]]->getCrds();
        const Vector &Xc2 = theNodes[edgeV[e][1]]->getCrds();
        const Vector &Xm  = theNodes[4 + e]->getCrds();
        double L2 = 0.0, d2 = 0.0;
        for (int k = 0; k < 3; k++) {
            double mid = 0.5 * (Xc1(k) + Xc2(k));
            d2 += (Xm(k) - mid) * (Xm(k) - mid);
            L2 += (Xc2(k) - Xc1(k)) * (Xc2(k) - Xc1(k));
        }
        if (L2 > 0.0 && d2 > 1.0e-12 * L2) {   // > ~1e-6 relative deviation
            opserr << "WARNING BezierTet10 " << this->getTag()
                   << " - mid-edge node " << connectedExternalNodes(4 + e)
                   << " is off the edge midpoint (curved element). v1 is "
                      "validated for straight-sided meshes only.\n";
        }
    }

    this->DomainComponent::setDomain(theDomain);
}


int BezierTet10::commitState()
{
    int retVal = 0;
    if ((retVal = this->Element::commitState()) != 0)
        opserr << "BezierTet10::commitState() - failed in base class\n";
    for (int i = 0; i < NGAUSS; i++)
        retVal += theMaterial[i]->commitState();
    retVal += theGeom->commitState();   // Ladruno (corot is stateless: no-op)
    return retVal;
}


int BezierTet10::revertToLastCommit()
{
    int retVal = 0;
    for (int i = 0; i < NGAUSS; i++)
        retVal += theMaterial[i]->revertToLastCommit();
    retVal += theGeom->revertToLastCommit();   // Ladruno (corot stateless)
    return retVal;
}


int BezierTet10::revertToStart()
{
    int retVal = 0;
    for (int i = 0; i < NGAUSS; i++)
        retVal += theMaterial[i]->revertToStart();
    retVal += theGeom->revertToStart();   // Ladruno (corot stateless)
    if (theNodes[0] != 0)
        computeControlPoints();
    return retVal;
}


int BezierTet10::update()
{
    // Ladruno: -geom finite (updated-Lagrangian) drives the material via
    // setTrialF(F) instead of the small-strain setTrialStrain path below.
    if (this->isFinite())
        return this->updateFinite();

    // Ladruno (seam 0+2): refresh the geometry method and localize the trial
    // displacement into the core frame. For -geom linear this is the identity
    // (uCore == uGlobal), so the strain below is bit-for-bit the direct kernel;
    // for -geom corot it is the de-rotated (Rᵀ) displacement. Strains ε = B·uCore
    // (or B̄·uCore) are then pushed to the materials.
    const Vector &u = this->computeLocalDisp();

    int ret = 0;

    double dN_avg[3][NEN];
    double volume;
    if (useBbar)
        computeVolumeAveragedDerivatives(dN_avg, volume);

    for (int gp = 0; gp < NGAUSS; gp++) {
        double B[NSTRESS][NELD];
        double factor = formBAtGauss(gp, dN_avg, B);
        if (factor < 0.0) {
            opserr << "BezierTet10::update() - degenerate Jacobian at GP " << gp
                   << " (element " << this->getTag() << ")\n";
            return -1;
        }

        Vector strain(NSTRESS);
        for (int i = 0; i < NSTRESS; i++) {
            double sum = 0.0;
            for (int j = 0; j < NELD; j++)
                sum += B[i][j] * u(j);
            strain(i) = sum;
        }

        ret += theMaterial[gp]->setTrialStrain(strain);
    }

    return ret;
}


// ═══════════════════════════════════════════════════════════════════
//  STIFFNESS AND FORCE
// ═══════════════════════════════════════════════════════════════════

const Matrix &BezierTet10::getTangentStiff()
{
    // Ladruno: -geom finite — full updated-Lagrangian consistent tangent
    // K = ∫ (∂Nₐ/∂x_j) a_ijkl (∂N_b/∂x_l) dv (geometric term already folded into
    // a_ijkl = c − σδ). No globalize seam — identity for finite.
    if (this->isFinite()) {
        static Vector fScratch(NELD);
        this->formResidAndTangentFinite(1, fScratch, &K_return);
        return K_return;
    }

    // K = Σᵢ wᵢ |Jᵢ| · Bᵢᵀ Dᵢ Bᵢ   (B → B̄ for B-bar), assembled in the CORE
    // frame in ONE pass that also yields the core internal force, then globalized
    // (seam 3). Ladruno: refresh the geometry method so K and the fCore used for
    // K_geo share one fresh R.
    this->computeLocalDisp();

    // one Gauss pass → core K (into K_return) AND the core internal force fCore.
    static Vector fCore(NELD);
    this->formCore(1, fCore, &K_return);

    // seam 3: globalize the core-frame K and add the corotational geometric
    // stiffness K_geo, which depends on fCore — the SAME force globalizeForce
    // rotates in getResistingForce (both via formCore), so they match by
    // construction. Identity (kGlobal=kCore, K_geo=0) for -geom linear.  // Ladruno
    theGeom->globalizeStiff(K_return, fCore, K_return);

    return K_return;
}


const Matrix &BezierTet10::getInitialStiff()
{
    // Ladruno: NO isFinite() branch — and that is intentional (mirrors
    // LadrunoBrick::getInitialStiff). At the undeformed reference state F = I,
    // σ = 0, and the spatial modulus equals the material's initial tangent, so
    // the finite tangent reduces EXACTLY to this small-strain ∫BᵀD₀B. The
    // converged solution comes from getTangentStiff → formResidAndTangentFinite;
    // Ki here only seeds initial-guess operators (algorithm('Initial'), Krylov,
    // Rayleigh βK0, eigen-on-Ki), so the small-strain reference tangent is the
    // correct "initial" stiffness and never enters the equilibrium path.
    if (Ki != 0)
        return *Ki;

    K_return.Zero();

    double dN_avg[3][NEN];
    double volume;
    if (useBbar)
        computeVolumeAveragedDerivatives(dN_avg, volume);

    for (int gp = 0; gp < NGAUSS; gp++) {
        double B[NSTRESS][NELD];
        double factor = formBAtGauss(gp, dN_avg, B);
        // Abandon (return), not skip-GP: a partial K is never correct, and a
        // degenerate reference geometry is permanent — return the abandoned
        // K_return WITHOUT caching it as Ki (mirrors the formCore guard).
        if (factor < 0.0) {
            opserr << "BezierTet10::getInitialStiff - degenerate Jacobian at GP "
                   << gp << " (element " << this->getTag() << ")\n";
            return K_return;
        }

        // Full BᵀD₀B — getInitialTangent() may be unsymmetric for the same
        // materials whose getTangent() is (see formCore); no mirror shortcut.
        const Matrix &D = theMaterial[gp]->getInitialTangent();
        double DB[NSTRESS][NELD];
        for (int k = 0; k < NSTRESS; k++)
            for (int Jc = 0; Jc < NELD; Jc++) {
                double sum = 0.0;
                for (int l = 0; l < NSTRESS; l++)
                    sum += D(k, l) * B[l][Jc];
                DB[k][Jc] = sum;
            }
        for (int I = 0; I < NELD; I++)
            for (int Jc = 0; Jc < NELD; Jc++) {
                double sum = 0.0;
                for (int k = 0; k < NSTRESS; k++)
                    sum += B[k][I] * DB[k][Jc];
                K_return(I, Jc) += factor * sum;
            }
    }

    // seam 3 (reference config): pin the geometry frame to R = I by refreshing
    // theGeom with cur == ref, THEN globalize with a zero core force. Without the
    // update-to-reference, globalizeStiff would reuse a STALE current-config R if
    // getInitialStiff is queried mid-analysis. R = I is deterministic, so caching
    // Ki below stays valid. Identity for -geom linear.  // Ladruno
    if (theNodes[0] != 0) {
        static Matrix refC(NEN, 3);
        for (int i = 0; i < NEN; i++) {
            const Vector &X = theNodes[i]->getCrds();
            for (int d = 0; d < 3; d++)
                refC(i, d) = X(d);
        }
        theGeom->update(NEN, refC, refC);          // Rmat = I
        static Vector zeroF(NELD);
        zeroF.Zero();
        theGeom->globalizeStiff(K_return, zeroF, K_return);
    }

    return LadrunoShell::cacheKi(Ki, K_return);
}


const Matrix &BezierTet10::getMass()
{
    // Ladruno (ADR-77 G2 ext): per-instance mass cache -- see LadrunoMassCache.h
    // for the contract and the brick's measured motivation. Signature = every
    // mutable scalar input: the element rho override + the NGAUSS material
    // rhos (mutable via setParameter "rho"). Coordinates guarded inside.
    // cMass/bbar are fixed at construction and deliberately omitted. The
    // 64-GP collapsed-Duffy consistent mass is the single most expensive
    // formation among the fork solids, which is why this element is in the
    // extension's scope at all.
    double mcSig[1 + NGAUSS];
    mcSig[0] = rho;
    for (int i = 0; i < NGAUSS; i++)
        mcSig[1 + i] = theMaterial[i]->getRho();
    if (const Matrix *Mc = massCache.lookup(mcSig, 1 + NGAUSS, theNodes, NEN, 3))
        return *Mc;

    // Lumped (default) or consistent (-cMass). Density: element rho
    // overrides; otherwise the material's own density.
    M_return.Zero();

    double rhoEff = rho;
    if (rhoEff == 0.0) {
        double sum = 0.0;
        for (int i = 0; i < NGAUSS; i++)
            sum += theMaterial[i]->getRho();
        rhoEff = sum / NGAUSS;
    }
    if (rhoEff == 0.0) {
        massCache.fill(M_return, mcSig, 1 + NGAUSS, theNodes, NEN, 3);   // Ladruno (ADR-77 G2 ext)
        return M_return;
    }

    if (!cMass) {
        // ─── Lumped mass (Kadapa Eq. 57) ──────────────────────
        //   Mᵉ = (ρ Vₑ / 10) · diag[1₁₀, 1₁₀, 1₁₀]
        //  Each of the 10 nodes gets an EQUAL share ρVₑ/10 — exact for the
        //  quadratic Bernstein tet (∫Nₐ = Vₑ/10 ∀a) and ALL POSITIVE, the
        //  reason this element exists for explicit dynamics.
        double Ve = this->computeVolume();
        double m  = rhoEff * Ve / 10.0;
        for (int a = 0; a < NEN; a++) {
            M_return(3*a,     3*a)     = m;
            M_return(3*a + 1, 3*a + 1) = m;
            M_return(3*a + 2, 3*a + 2) = m;
        }
        massCache.fill(M_return, mcSig, 1 + NGAUSS, theNodes, NEN, 3);   // Ladruno (ADR-77 G2 ext)
        return M_return;
    }

    // ─── Consistent mass:  M = ∫_Ω ρ Nᵤᵀ Nᵤ dΩ ────────────────
    //  4×4×4 collapsed-Duffy rule (degree 7) integrates the degree-6
    //  integrand exactly for straight edges. All entries ≥ 0 (Nₐ ≥ 0).
    for (int i = 0; i < 4; i++) {
        double a = GL4_t[i];
        for (int j = 0; j < 4; j++) {
            double bb = GL4_t[j];
            for (int k = 0; k < 4; k++) {
                double cc = GL4_t[k];
                double L1 = a;
                double L2 = (1.0 - a) * bb;
                double L3 = (1.0 - a) * (1.0 - bb) * cc;
                double jac = (1.0 - a) * (1.0 - a) * (1.0 - bb);
                double w = GL4_w[i] * GL4_w[j] * GL4_w[k] * jac;

                double N[NEN];
                shapeFunctions(L1, L2, L3, N);

                double dN[3][NEN];
                shapeDerivatives(L1, L2, L3, dN);

                double J[3][3];
                double dN_dx[3][NEN];
                double detJ = computeJacobian(dN, J, dN_dx);

                double factor = w * fabs(detJ) * rhoEff;

                for (int p = 0; p < NEN; p++) {
                    for (int q = 0; q < NEN; q++) {
                        double val = factor * N[p] * N[q];
                        M_return(3*p,     3*q)     += val;
                        M_return(3*p + 1, 3*q + 1) += val;
                        M_return(3*p + 2, 3*q + 2) += val;
                    }
                }
            }
        }
    }

    massCache.fill(M_return, mcSig, 1 + NGAUSS, theNodes, NEN, 3);   // Ladruno (ADR-77 G2 ext)
    return M_return;
}


void BezierTet10::zeroLoad()
{
    Q.Zero();
    applyLoad = 0;
    appliedB[0] = appliedB[1] = appliedB[2] = 0.0;
}


int BezierTet10::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    // Supports the standard self-weight (gravity) load so the body force can
    // be ramped by a load pattern / time series (same idiom as FourNodeQuad).
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_SelfWeight) {
        applyLoad = 1;
        appliedB[0] += loadFactor * data(0) * b[0];
        appliedB[1] += loadFactor * data(1) * b[1];
        appliedB[2] += loadFactor * data(2) * b[2];
        return 0;
    }

    opserr << "BezierTet10::addLoad - load type " << type
           << " unknown for element " << this->getTag()
           << " (only SelfWeight is supported)\n";
    return -1;
}


int BezierTet10::addInertiaLoadToUnbalance(const Vector &accel)
{
    const Matrix &M = this->getMass();
    bool hasMass = false;
    for (int i = 0; i < NELD && !hasMass; i++)
        if (M(i, i) != 0.0) hasMass = true;
    if (!hasMass)
        return 0;

    static Vector a(NELD);
    // Q -= M × R·a_g  (the OpenSees convention: getResistingForce() subtracts Q,
    // so the unbalance gains -M·R·a_g). Was `+1.0` until WP-117, which drove
    // the element mass with -a_g under UniformExcitation. Matches
    // TenNodeTetrahedron / LadrunoBrick / every other fork element.  // Ladruno
    return LadrunoShell::addGroundInertia(Q, M, theNodes, NEN, NDOF, false, accel, a,
                                          "BezierTet10");
}


const Vector &BezierTet10::getResistingForce()
{
    // F = globalize(∫_Ω Bᵀσ dΩ)  −  body force  −  pressure  −  Q.
    //
    // Ladruno (corot load-frame contract): the INTERNAL force ∫Bᵀσ is assembled
    // in the core frame and rotated to global by globalizeForce (seam 3). The
    // fixed-direction EXTERNAL loads — body force (gravity), the +z pressure
    // hack, and Q (applied / inertia) — are applied in the GLOBAL frame AFTER
    // globalizeForce so the rotation R never co-rotates them. This keeps the
    // fCore fed to globalizeStiff equal to pure ∫Bᵀσ (so K_geo carries no
    // spurious external-load term) and makes the f/K paths trivially consistent.
    // For -geom linear globalizeForce is the identity and this reproduces the
    // direct kernel.
    //
    // Ladruno: -geom finite assembles ∫Bᵀσ dv directly on the current config
    // (Cauchy σ; no localize/globalize — identity for finite). Both paths then
    // share the fixed-direction external-load tail (body force / pressure / Q)
    // applied in the GLOBAL frame below.
    static Vector fInt(NELD);
    if (this->isFinite()) {
        this->formResidAndTangentFinite(0, fInt, 0);   // current-config ∫Bᵀσ dv
        P_return = fInt;
    } else {
        this->computeLocalDisp();                     // seam 0+2: refresh frame
        this->formCore(0, fInt, 0);                   // core-frame ∫Bᵀσ (no tangent)
        theGeom->globalizeForce(fInt, fInt);          // seam 3: → global
        P_return = fInt;
    }

    // fixed-direction external loads, GLOBAL frame
    double bx = (applyLoad == 1) ? appliedB[0] : b[0];
    double by = (applyLoad == 1) ? appliedB[1] : b[1];
    double bz = (applyLoad == 1) ? appliedB[2] : b[2];
    if (bx != 0.0 || by != 0.0 || bz != 0.0 || pressure != 0.0) {
        for (int gp = 0; gp < NGAUSS; gp++) {
            double w = GP4_w[gp];
            double N[NEN];
            shapeFunctions(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], N);

            double dN[3][NEN];
            shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);

            double J[3][3];
            double dN_dx[3][NEN];
            double detJ = computeJacobian(dN, J, dN_dx);
            double factor = w * fabs(detJ);

            for (int a = 0; a < NEN; a++) {
                P_return(3*a)     -= factor * N[a] * bx;
                P_return(3*a + 1) -= factor * N[a] * by;
                P_return(3*a + 2) -= factor * N[a] * bz;
                if (pressure != 0.0)   // +z volume hack (mirrors BezierTri6)
                    P_return(3*a + 2) -= factor * N[a] * pressure;
            }
        }
    }

    P_return.addVector(1.0, Q, -1.0);
    return P_return;
}


// ─── Geometry-method seam helpers (Ladruno) ───────────────────────────
// computeLocalDisp: seam 0+2 — refresh theGeom from the current geometry and
// return the localized (core-frame) trial displacement. Identity for -geom
// linear (uCore == uGlobal). Mirrors LadrunoBrick::computeLocalDisp.
const Vector &BezierTet10::computeLocalDisp(void)
{
    static Matrix refCrds(NEN, 3), curCrds(NEN, 3);
    static Vector uGlobal(NELD), uCore(NELD);

    for (int i = 0; i < NEN; i++) {
        const Vector &X = theNodes[i]->getCrds();
        const Vector &u = theNodes[i]->getTrialDisp();
        for (int d = 0; d < 3; d++) {
            refCrds(i, d)     = X(d);
            curCrds(i, d)     = X(d) + u(d);
            uGlobal(3 * i + d) = u(d);
        }
    }

    theGeom->update(NEN, refCrds, curCrds);
    theGeom->localizeDisp(uGlobal, uCore);
    return uCore;
}

// formBAtGauss: the single guarded strain-displacement assembly shared by every
// Gauss loop (update / formCore / getInitialStiff). Returns w·|detJ|, or −1 if
// the element is degenerate (detJ == 0, where computeJacobian leaves dN_dx unset
// — so the degenerate check MUST live here, not be duplicated/forgotten per loop).
double BezierTet10::formBAtGauss(int gp, const double dN_avg[3][NEN],
                                 double B[NSTRESS][NELD]) const
{
    double dN[3][NEN];
    shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);

    double J[3][3];
    double dN_dx[3][NEN];
    double detJ = computeJacobian(dN, J, dN_dx);
    if (fabs(detJ) <= 0.0)
        return -1.0;

    if (useBbar)
        computeBBarMatrix(dN_dx, dN_avg, B);
    else
        computeBMatrix(dN_dx, B);

    return GP4_w[gp] * fabs(detJ);
}

// formCore: ONE Gauss pass → core-frame internal force fInt = ∫Bᵀσ dΩ (always)
// and, when tangFlag, the core tangent K = ∫BᵀDB dΩ. NO body force / pressure /
// Q. getResistingForce calls formCore(0,…) and getTangentStiff formCore(1,…), so
// the fCore globalizeForce rotates equals the fCore globalizeStiff uses for K_geo
// BY CONSTRUCTION — and the tangent path no longer runs a second Bᵀσ loop.
void BezierTet10::formCore(int tangFlag, Vector &fInt, Matrix *K)
{
    fInt.Zero();
    if (tangFlag && K != 0)
        K->Zero();

    double dN_avg[3][NEN];
    double volume;
    if (useBbar)
        computeVolumeAveragedDerivatives(dN_avg, volume);

    for (int gp = 0; gp < NGAUSS; gp++) {
        double B[NSTRESS][NELD];
        double factor = formBAtGauss(gp, dN_avg, B);
        // Abandon the pass (return), not skip-GP (continue): a partially-
        // assembled fInt/K is never correct. For a straight-sided tet detJ
        // is constant so one bad GP means the element is degenerate
        // everywhere; a curved tet CAN degenerate at a single GP, and
        // abandoning remains the fail-visible choice there too.
        // Unreachable-by-construction in the normal Newton flow —
        // update() screens the same condition first and returns -1 (step-cut)
        // before any assembly — but fail safe if an accessor is ever called
        // on a degenerate trial state (mirrors formResidAndTangentFinite).
        if (factor < 0.0) {
            opserr << "BezierTet10::formCore - degenerate Jacobian at GP " << gp
                   << " (element " << this->getTag() << ")\n";
            return;
        }

        const Vector &sigma = theMaterial[gp]->getStress();
        for (int I = 0; I < NELD; I++) {
            double s = 0.0;
            for (int k = 0; k < NSTRESS; k++)
                s += B[k][I] * sigma(k);
            fInt(I) += factor * s;
        }

        if (tangFlag && K != 0) {
            // D is GENERALLY UNSYMMETRIC (non-associated flow: DruckerPrager
            // with rho_bar != rho, ManzariDafalias always), so the full BᵀDB
            // must be assembled — an upper-triangle-and-mirror shortcut writes
            // the Bᵀ(I)DB(J) contraction into K(J,I), whose true value is the
            // Dᵀ contraction, silently corrupting the tangent the moment any
            // GP yields (TIMs report item 1). Cache DB once per GP; this is
            // no more work than the old triangular double contraction.
            const Matrix &D = theMaterial[gp]->getTangent();
            double DB[NSTRESS][NELD];
            for (int k = 0; k < NSTRESS; k++)
                for (int Jc = 0; Jc < NELD; Jc++) {
                    double sum = 0.0;
                    for (int l = 0; l < NSTRESS; l++)
                        sum += D(k, l) * B[l][Jc];
                    DB[k][Jc] = sum;
                }
            for (int I = 0; I < NELD; I++)
                for (int Jc = 0; Jc < NELD; Jc++) {
                    double sum = 0.0;
                    for (int k = 0; k < NSTRESS; k++)
                        sum += B[k][I] * DB[k][Jc];
                    (*K)(I, Jc) += factor * sum;
                }
        }
    }
}


// ═══════════════════════════════════════════════════════════════════
//  -geom finite (Ladruno): updated-Lagrangian finite-strain path
// ═══════════════════════════════════════════════════════════════════
//
// The material is a FiniteStrainNDMaterial driven by setTrialF(F): it returns
// the Cauchy stress σ (getStress) and the FULL 4th-order spatial constitutive
// modulus c (getSpatialTangentTensor). Per the LOCKED seam-3 contract the
// element owns the geometric stiffness, so it forms the consistent tangent on
// the current configuration:
//     f = ∫ Bᵀ σ dv ,   K = ∫ (∂Nₐ/∂x_j) a_ijkl (∂N_b/∂x_l) dv ,   dv = J dV
// with a_ijkl = c_ijkl − σ_il δ_jk and spatial gradients
// ∂Nₐ/∂xⱼ = Σ_k (∂Nₐ/∂X_k)(F⁻¹)_kj. F is built from the REFERENCE Bernstein
// gradients (computeJacobian on controlPts); for straight-sided tets detJ_ref
// is constant per element but det F varies (quadratic u → linear F). This is
// the 10-node/4-GP analogue of LadrunoBrick's finite path — F-bar (bbar+finite)
// is a step-2 follow-up, rejected at parse time.

// Inverse of a row-major 3×3 via adjugate/det. Returns det F; fills Finv (row-
// major), zero-filled if |det| underflows (one source of truth for the 9-term
// adjugate, the kind of formula that drifts under copy-paste).
static double
invert3x3(const double F[9], double Finv[9])
{
    double det = F[0]*(F[4]*F[8]-F[5]*F[7]) - F[1]*(F[3]*F[8]-F[5]*F[6])
               + F[2]*(F[3]*F[7]-F[4]*F[6]);
    double id = (det != 0.0) ? 1.0 / det : 0.0;
    Finv[0] =  (F[4]*F[8]-F[5]*F[7])*id;  Finv[1] = -(F[1]*F[8]-F[2]*F[7])*id;
    Finv[2] =  (F[1]*F[5]-F[2]*F[4])*id;  Finv[3] = -(F[3]*F[8]-F[5]*F[6])*id;
    Finv[4] =  (F[0]*F[8]-F[2]*F[6])*id;  Finv[5] = -(F[0]*F[5]-F[2]*F[3])*id;
    Finv[6] =  (F[3]*F[7]-F[4]*F[6])*id;  Finv[7] = -(F[0]*F[7]-F[1]*F[6])*id;
    Finv[8] =  (F[0]*F[4]-F[1]*F[3])*id;
    return det;
}

bool BezierTet10::isFinite(void) const
{
    return theGeom != 0 &&
           theGeom->getStrainMeasure() ==
             SolidTransformation::StrainMeasure::DeformationGradient;
}

double BezierTet10::deformationGradient(const double dN_dX[3][NEN],
                                        double F[9]) const
{
    for (int i = 0; i < 9; i++) F[i] = 0.0;
    F[0] = F[4] = F[8] = 1.0;
    for (int a = 0; a < NEN; a++) {
        const Vector &ua = theNodes[a]->getTrialDisp();   // TOTAL trial disp
        for (int i = 0; i < 3; i++)
            for (int J = 0; J < 3; J++)
                F[3*i + J] += ua(i) * dN_dX[J][a];
    }
    return F[0]*(F[4]*F[8]-F[5]*F[7]) - F[1]*(F[3]*F[8]-F[5]*F[6])
         + F[2]*(F[3]*F[7]-F[4]*F[6]);
}

// F-bar centroid data (bbar + finite). F₀ at the TET centroid — barycentric
// L=(¼,¼,¼), NOT the brick's isoparametric (0,0,0). Returns J₀ = det F₀; fills
// G0[k][b] = ∂N_b/∂x_k|_centroid (from F₀⁻¹) when requested. dSNPO eq 15.5.  // Ladruno
double BezierTet10::centroidFbar(double (*G0)[NEN]) const
{
    double dN0[3][NEN];
    shapeDerivatives(0.25, 0.25, 0.25, dN0);     // tet centroid

    double J0m[3][3];
    double dN0_dX[3][NEN];                        // reference ∂Nₐ/∂X at the centroid
    double detJ0 = computeJacobian(dN0, J0m, dN0_dX);
    if (fabs(detJ0) <= 0.0)
        return 0.0;                               // degenerate; caller's J₀≤0 guard fires

    double F0[9];
    double J0 = deformationGradient(dN0_dX, F0);

    if (G0 != 0) {
        double Fi[9];
        invert3x3(F0, Fi);                        // F₀⁻¹ (row-major)
        for (int b = 0; b < NEN; b++)
            for (int k = 0; k < 3; k++) {
                double s = 0.0;
                for (int m = 0; m < 3; m++)
                    s += dN0_dX[m][b] * Fi[3*m + k];   // ∂N_b/∂x_k at the centroid
                G0[k][b] = s;
            }
    }
    return J0;
}

// F-bar MEAN-DILATATION data — the consistent analogue of centroidFbar with the
// single centroid point replaced by volume averages over the Gauss points:
//   J̄    = (∫ J dV₀)/V₀                = (Σ_gp J_gp dV₀_gp)/(Σ_gp dV₀_gp)
//   Ḡ_kb = (∫ ∂N_b/∂x_k dv)/v          = (Σ_gp g_kb,gp dv_gp)/(Σ_gp dv_gp)
// with dV₀ = w·|detJ_ref| (reference) and dv = J·dV₀ (current). This is the
// classic mean-dilatation (Nagtegaal–Parks–Rice / Simo–Taylor–Pister) variant;
// at small strain it reduces to the element's volume-averaged small-strain bbar.
// Returns J̄ (0.0 on a degenerate/inverted GP so the caller's J̄≤0 guard fires).  // Ladruno
double BezierTet10::fbarMeanDilatation(double (*Gbar)[NEN]) const
{
    double sumJ_dV0 = 0.0, sum_dV0 = 0.0;        // Σ J dV₀ (= v), Σ dV₀ (= V₀)
    double Gacc[3][NEN];
    if (Gbar != 0)
        for (int k = 0; k < 3; k++)
            for (int b = 0; b < NEN; b++) Gacc[k][b] = 0.0;

    for (int gp = 0; gp < NGAUSS; gp++) {
        double dN[3][NEN];
        shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);
        double Jm[3][3], dN_dX[3][NEN];
        double detJ = computeJacobian(dN, Jm, dN_dX);
        if (fabs(detJ) <= 0.0)
            return 0.0;                          // degenerate; caller's J̄≤0 guard fires

        double F[9];
        double Jgp = deformationGradient(dN_dX, F);
        double dV0 = GP4_w[gp] * fabs(detJ);     // reference measure
        sumJ_dV0 += Jgp * dV0;
        sum_dV0  += dV0;

        if (Gbar != 0) {
            if (Jgp <= 0.0)
                return 0.0;                      // can't push the gradient on an inverted GP
            double Fi[9];
            invert3x3(F, Fi);
            double dv = Jgp * dV0;               // current measure
            for (int b = 0; b < NEN; b++)
                for (int k = 0; k < 3; k++) {
                    double s = 0.0;
                    for (int m = 0; m < 3; m++)
                        s += dN_dX[m][b] * Fi[3*m + k];   // ∂N_b/∂x_k at this GP
                    Gacc[k][b] += s * dv;
                }
        }
    }
    if (sum_dV0 <= 0.0)
        return 0.0;
    double Jbar = sumJ_dV0 / sum_dV0;            // = v/V₀
    if (Gbar != 0 && sumJ_dV0 > 0.0)            // v = Σ dv = Σ J dV₀ = sumJ_dV0
        for (int k = 0; k < 3; k++)
            for (int b = 0; b < NEN; b++) Gbar[k][b] = Gacc[k][b] / sumJ_dV0;
    return Jbar;
}

int BezierTet10::updateFinite(void)
{
    static Matrix Fm(3, 3);

    // bbar + finite = F-bar: every GP is driven by F̄ = (Ĵ/J)^(1/3) F so they
    // share one bar dilatation Ĵ (volumetric-locking cure). std uses F. Ĵ is the
    // centroid J₀ (-fbar centroid) or the volume average J̄ (-fbar mean_dilatation).
    const bool useFbar = useBbar;
    double Jhat = 1.0;
    if (useFbar) {
        Jhat = (fbarMode == FBAR_MEAN) ? this->fbarMeanDilatation()
                                       : this->centroidFbar();
        if (Jhat <= 0.0) {
            opserr << "BezierTet10::updateFinite - non-positive F-bar dilatation Ĵ = "
                   << Jhat << " (element " << this->getTag() << ", -fbar "
                   << (fbarMode == FBAR_MEAN ? "mean_dilatation" : "centroid") << ")\n";
            return -1;
        }
    }

    for (int gp = 0; gp < NGAUSS; gp++) {
        double dN[3][NEN];
        shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);

        double J[3][3];
        double dN_dx[3][NEN];                    // reference ∂Nₐ/∂X
        double detJ = computeJacobian(dN, J, dN_dx);
        if (fabs(detJ) <= 0.0) {
            opserr << "BezierTet10::updateFinite - degenerate Jacobian at GP " << gp
                   << " (element " << this->getTag() << ")\n";
            return -1;
        }

        double F[9];
        double Jdet = deformationGradient(dN_dx, F);   // TOTAL F; material derives F_Δ
        if (Jdet <= 0.0) {
            opserr << "BezierTet10::updateFinite - non-positive det F = " << Jdet
                   << " at GP " << gp << " (element " << this->getTag() << ")\n";
            return -1;
        }
        if (useFbar) {
            double s = pow(Jhat / Jdet, 1.0 / 3.0);    // F̄ = (Ĵ/J)^(1/3) F (det F̄ = Ĵ > 0)
            for (int n = 0; n < 9; n++) F[n] *= s;
        }
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                Fm(r, c) = F[3*r + c];

        FiniteStrainNDMaterial *fsm =
            dynamic_cast<FiniteStrainNDMaterial *>(theMaterial[gp]);
        if (fsm == 0) {
            opserr << "BezierTet10::updateFinite - material at GP " << gp
                   << " is not a FiniteStrainNDMaterial (element " << this->getTag()
                   << ")\n";
            return -1;
        }
        if (fsm->setTrialF(Fm) < 0) {
            opserr << "BezierTet10::updateFinite - setTrialF failed at GP " << gp
                   << " (element " << this->getTag() << ", det F<=0?)\n";
            return -1;
        }
    }
    return 0;
}

void BezierTet10::formResidAndTangentFinite(int tangFlag, Vector &fInt, Matrix *K)
{
    fInt.Zero();
    if (tangFlag && K != 0)
        K->Zero();

    // F-bar (bbar + finite): the eq 15.10 tangent gains a (generally UNSYMMETRIC)
    // coupling to an element-wide gradient operator Ĝ. The residual is unchanged —
    // F-bar enters it only through σ̄ (set in updateFinite, eq 15.9). Ĝ is the
    // centroid gradient G₀ (-fbar centroid) or the volume-averaged Ḡ (-fbar
    // mean_dilatation); compute it once.  // Ladruno
    const bool useFbar = useBbar;
    double Ghat[3][NEN];
    if (useFbar && tangFlag && K != 0) {
        if (fbarMode == FBAR_MEAN) this->fbarMeanDilatation(Ghat);
        else                       this->centroidFbar(Ghat);
    }

    for (int gp = 0; gp < NGAUSS; gp++) {
        double dN[3][NEN];
        shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);

        double Jm[3][3];
        double dN_dX[3][NEN];                    // reference ∂Nₐ/∂X
        double detJ = computeJacobian(dN, Jm, dN_dX);
        // Both degeneracy guards ABANDON the pass (return), not skip-GP (continue):
        // a partially-assembled fInt/K is never correct, and for a straight-sided
        // tet detJ_ref is constant so one bad GP means the element is degenerate
        // everywhere. These paths are unreachable-by-construction in the normal
        // Newton flow — update()→updateFinite() screens BOTH conditions first and
        // returns -1 (step-cut) before any assembly — but the guards fail safe if
        // an accessor is ever called on a degenerate trial state.  // Ladruno
        if (fabs(detJ) <= 0.0) {
            opserr << "BezierTet10::formResidAndTangentFinite - degenerate Jacobian "
                      "at GP " << gp << " (element " << this->getTag() << ")\n";
            return;
        }

        double F[9];
        double J = deformationGradient(dN_dX, F);
        if (J <= 0.0) {
            opserr << "BezierTet10::formResidAndTangentFinite - non-positive det F = "
                   << J << " at GP " << gp << " (element " << this->getTag() << ")\n";
            return;
        }

        double Fi[9];
        invert3x3(F, Fi);                         // F⁻¹ (row-major)

        double g[3][NEN];                         // spatial gradients ∂Nₐ/∂x_j
        for (int a = 0; a < NEN; a++)
            for (int j = 0; j < 3; j++) {
                double s = 0.0;
                for (int k = 0; k < 3; k++)
                    s += dN_dX[k][a] * Fi[3*k + j];
                g[j][a] = s;
            }

        double dv = J * fabs(detJ) * GP4_w[gp];   // current-config measure (w·|detJ_ref|·J)

        const Vector &stress = theMaterial[gp]->getStress();   // Cauchy σ (set via setTrialF)
        double sig[3][3];
        sig[0][0]=stress(0); sig[1][1]=stress(1); sig[2][2]=stress(2);
        sig[0][1]=sig[1][0]=stress(3);
        sig[1][2]=sig[2][1]=stress(4);
        sig[2][0]=sig[0][2]=stress(5);

        // residual f_{a,i} = ∫ σ_ij ∂Nₐ/∂x_j dv. NO body force here — it is applied
        // in the global frame in getResistingForce (the small-strain contract).
        for (int a = 0; a < NEN; a++)
            for (int i = 0; i < 3; i++) {
                double fi = 0.0;
                for (int j = 0; j < 3; j++)
                    fi += sig[i][j] * g[j][a];
                fInt(3*a + i) += fi * dv;
            }

        if (tangFlag && K != 0) {
            // FULL 4th-order spatial constitutive modulus c_ijkl. The material's
            // 6×6 getTangent is LOSSY in (k,l) (c = (1/2J)[D:L:B] is not minor-
            // symmetric); the consistent tangent needs the full tensor (see
            // FiniteStrainNDMaterial.h). a_ijkl = c_ijkl − σ_il δ_jk folds in the
            // element-owned geometric (initial-stress) term.  // Ladruno
            FiniteStrainNDMaterial *fsm =
                dynamic_cast<FiniteStrainNDMaterial *>(theMaterial[gp]);
            static double cmat[3][3][3][3];
            if (fsm == 0 || fsm->getSpatialTangentTensor(cmat) != 0) {
                opserr << "BezierTet10::formResidAndTangentFinite - material at GP "
                       << gp << " did not provide a full spatial tangent "
                          "(getSpatialTangentTensor); element " << this->getTag() << "\n";
                return;
            }
            double a4[3][3][3][3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            a4[i][j][k][l] = cmat[i][j][k][l]
                                           - sig[i][l] * (j == k ? 1.0 : 0.0);

            for (int a = 0; a < NEN; a++)
                for (int bn = 0; bn < NEN; bn++)
                    for (int i = 0; i < 3; i++)
                        for (int k = 0; k < 3; k++) {
                            double s = 0.0;
                            for (int j = 0; j < 3; j++)
                                for (int l = 0; l < 3; l++)
                                    s += g[j][a] * a4[i][j][k][l] * g[l][bn];
                            (*K)(3*a + i, 3*bn + k) += s * dv;
                        }

            // F-bar additional (generally UNSYMMETRIC) stiffness, dSNPO eq 15.10:
            //   K_{(a,i)(b,k)} += ∫ (Σ_j g[j][a] q_ij)(Ĝ[k][b] − g[k][b]) dv,
            // with q the matrix form of the eq 15.11 tensor at F=F̄,
            //   q_ij = (1/3) a_ijpp − (2/3) σ̄_ij ,
            // using the SAME a4 = c̄ − σ̄δ modulus as the std term above (NOT the
            // material part c̄ alone — the −(2/3)σ̄ is the spatial initial-stress
            // part, NOT a (1/3)c shortcut). Ĝ is the centroid G₀ or the volume
            // average Ḡ; (Ĝ − g) vanishes when Ĝ = g, collapsing to plain F.  // Ladruno
            if (useFbar) {
                double M[3][3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        M[i][j] = (a4[i][j][0][0] + a4[i][j][1][1] + a4[i][j][2][2]) / 3.0
                                - (2.0 / 3.0) * sig[i][j];

                double Lfac[NEN][3];             // Lfac[a][i] = Σ_j g[j][a] q_ij
                for (int a = 0; a < NEN; a++)
                    for (int i = 0; i < 3; i++) {
                        double s = 0.0;
                        for (int j = 0; j < 3; j++) s += g[j][a] * M[i][j];
                        Lfac[a][i] = s;
                    }

                for (int a = 0; a < NEN; a++)
                    for (int i = 0; i < 3; i++)
                        for (int bn = 0; bn < NEN; bn++)
                            for (int kk = 0; kk < 3; kk++)
                                (*K)(3*a + i, 3*bn + kk)
                                    += Lfac[a][i] * (Ghat[kk][bn] - g[kk][bn]) * dv;
            }
        }
    }
}


const Vector &BezierTet10::getResistingForceIncInertia()
{
    // SNAPSHOT into a function-local buffer before anything else runs:
    // getMass() and the betaK Rayleigh re-entry (getTangentStiff) must not be
    // able to reach the vector being accumulated. Returning `res` also keeps it
    // distinct from P_return, which getResistingForce() returns. Same operations
    // in the same order -- ((f - Q) + M*a) + R -- so results are bit-identical
    // (WP-115; LEDGER_quirks "MUST snapshot the shared static `resid`").  // Ladruno
    static Vector res(NELD);
    res = this->getResistingForce();

    const Matrix &M = this->getMass();
    bool hasMass = false;
    for (int i = 0; i < NELD && !hasMass; i++)
        if (M(i, i) != 0.0) hasMass = true;

    if (hasMass) {
        static Vector a(NELD);   // Ladruno (WP-124 stage 5)
        LadrunoShell::addNodalInertia(res, M, theNodes, NEN, NDOF, false, a);
    }

    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        res += this->getRayleighDampingForces();

    return res;
}


// ═══════════════════════════════════════════════════════════════════
//  SHAPE FUNCTIONS AND DERIVATIVES
// ═══════════════════════════════════════════════════════════════════

void BezierTet10::shapeFunctions(double L1, double L2, double L3,
                                 double N[NEN]) const
{
    double L4 = 1.0 - L1 - L2 - L3;

    N[0] = L1 * L1;        // N1  vertex 1
    N[1] = L2 * L2;        // N2  vertex 2
    N[2] = L3 * L3;        // N3  vertex 3
    N[3] = L4 * L4;        // N4  vertex 4
    N[4] = 2.0 * L1 * L2;  // N5  edge (1-2)
    N[5] = 2.0 * L2 * L3;  // N6  edge (2-3)
    N[6] = 2.0 * L1 * L3;  // N7  edge (1-3)
    N[7] = 2.0 * L1 * L4;  // N8  edge (1-4)
    N[8] = 2.0 * L3 * L4;  // N9  edge (3-4)
    N[9] = 2.0 * L2 * L4;  // N10 edge (2-4)
}


void BezierTet10::shapeDerivatives(double L1, double L2, double L3,
                                   double dN[3][NEN]) const
{
    //  ∂Nₐ/∂L1, ∂Nₐ/∂L2, ∂Nₐ/∂L3 with L4 = 1-L1-L2-L3 (∂L4/∂Lk = -1).
    double L4 = 1.0 - L1 - L2 - L3;

    // ∂N/∂L1
    dN[0][0] =  2.0 * L1;          // L1²
    dN[0][1] =  0.0;               // L2²
    dN[0][2] =  0.0;               // L3²
    dN[0][3] = -2.0 * L4;          // L4²
    dN[0][4] =  2.0 * L2;          // 2 L1 L2
    dN[0][5] =  0.0;               // 2 L2 L3
    dN[0][6] =  2.0 * L3;          // 2 L1 L3
    dN[0][7] =  2.0 * (L4 - L1);   // 2 L1 L4
    dN[0][8] = -2.0 * L3;          // 2 L3 L4
    dN[0][9] = -2.0 * L2;          // 2 L2 L4

    // ∂N/∂L2
    dN[1][0] =  0.0;
    dN[1][1] =  2.0 * L2;
    dN[1][2] =  0.0;
    dN[1][3] = -2.0 * L4;
    dN[1][4] =  2.0 * L1;
    dN[1][5] =  2.0 * L3;
    dN[1][6] =  0.0;
    dN[1][7] = -2.0 * L1;
    dN[1][8] = -2.0 * L3;
    dN[1][9] =  2.0 * (L4 - L2);

    // ∂N/∂L3
    dN[2][0] =  0.0;
    dN[2][1] =  0.0;
    dN[2][2] =  2.0 * L3;
    dN[2][3] = -2.0 * L4;
    dN[2][4] =  0.0;
    dN[2][5] =  2.0 * L2;
    dN[2][6] =  2.0 * L1;
    dN[2][7] = -2.0 * L1;
    dN[2][8] =  2.0 * (L4 - L3);
    dN[2][9] = -2.0 * L2;
}


// ═══════════════════════════════════════════════════════════════════
//  JACOBIAN AND B-MATRIX
// ═══════════════════════════════════════════════════════════════════

double BezierTet10::computeJacobian(const double dN[3][NEN],
                                    double J[3][3],
                                    double dN_dx[3][NEN]) const
{
    //  J[i][j] = ∂x_j/∂L_i = Σₐ ∂Nₐ/∂L_i · Pₐ_j   (control points).
    //  dN/dx = J⁻¹ · dN/dL. Returns SIGNED detJ; callers use |detJ| for the
    //  integration measure (orientation-robust — ADR O11).
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            J[i][j] = 0.0;

    for (int a = 0; a < NEN; a++)
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                J[i][j] += dN[i][a] * controlPts[a][j];

    double detJ =
        J[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
      - J[0][1] * (J[1][0]*J[2][2] - J[1][2]*J[2][0])
      + J[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0]);

    if (fabs(detJ) > 0.0) {
        double inv = 1.0 / detJ;
        double Ji[3][3];
        Ji[0][0] =  (J[1][1]*J[2][2] - J[1][2]*J[2][1]) * inv;
        Ji[0][1] = -(J[0][1]*J[2][2] - J[0][2]*J[2][1]) * inv;
        Ji[0][2] =  (J[0][1]*J[1][2] - J[0][2]*J[1][1]) * inv;
        Ji[1][0] = -(J[1][0]*J[2][2] - J[1][2]*J[2][0]) * inv;
        Ji[1][1] =  (J[0][0]*J[2][2] - J[0][2]*J[2][0]) * inv;
        Ji[1][2] = -(J[0][0]*J[1][2] - J[0][2]*J[1][0]) * inv;
        Ji[2][0] =  (J[1][0]*J[2][1] - J[1][1]*J[2][0]) * inv;
        Ji[2][1] = -(J[0][0]*J[2][1] - J[0][1]*J[2][0]) * inv;
        Ji[2][2] =  (J[0][0]*J[1][1] - J[0][1]*J[1][0]) * inv;

        for (int a = 0; a < NEN; a++)
            for (int i = 0; i < 3; i++)
                dN_dx[i][a] = Ji[i][0]*dN[0][a]
                            + Ji[i][1]*dN[1][a]
                            + Ji[i][2]*dN[2][a];
    }

    return detJ;
}


void BezierTet10::computeBMatrix(const double dN_dx[3][NEN],
                                 double B[NSTRESS][NELD]) const
{
    //  3D strain-displacement (Kadapa Eq. 38). Voigt {xx,yy,zz,xy,yz,zx}.
    //   Bₐ = [ N,x  0    0   ]
    //        [ 0    N,y  0   ]
    //        [ 0    0    N,z ]
    //        [ N,y  N,x  0   ]
    //        [ 0    N,z  N,y ]
    //        [ N,z  0    N,x ]
    for (int i = 0; i < NSTRESS; i++)
        for (int j = 0; j < NELD; j++)
            B[i][j] = 0.0;

    for (int a = 0; a < NEN; a++) {
        double bx = dN_dx[0][a], by = dN_dx[1][a], bz = dN_dx[2][a];
        int cx = 3*a, cy = 3*a + 1, cz = 3*a + 2;

        B[0][cx] = bx;
        B[1][cy] = by;
        B[2][cz] = bz;
        B[3][cx] = by;  B[3][cy] = bx;
        B[4][cy] = bz;  B[4][cz] = by;
        B[5][cx] = bz;  B[5][cz] = bx;
    }
}


void BezierTet10::computeBBarMatrix(const double dN_dx[3][NEN],
                                    const double dN_avg[3][NEN],
                                    double Bbar[NSTRESS][NELD]) const
{
    //  B-bar (Kadapa Eq. 45, full 3D). The three normal-strain rows get the
    //  volumetric (1/3) split with the volume-averaged derivatives B̄; the
    //  three shear rows are unchanged.
    for (int i = 0; i < NSTRESS; i++)
        for (int j = 0; j < NELD; j++)
            Bbar[i][j] = 0.0;

    for (int a = 0; a < NEN; a++) {
        double B1 = dN_dx[0][a], B2 = dN_dx[1][a], B3 = dN_dx[2][a];
        double A1 = dN_avg[0][a], A2 = dN_avg[1][a], A3 = dN_avg[2][a];
        int cx = 3*a, cy = 3*a + 1, cz = 3*a + 2;

        // normal-strain rows (volumetric split)
        Bbar[0][cx] = (A1 + 2.0*B1) / 3.0;
        Bbar[1][cx] = (A1 - B1) / 3.0;
        Bbar[2][cx] = (A1 - B1) / 3.0;
        Bbar[0][cy] = (A2 - B2) / 3.0;
        Bbar[1][cy] = (A2 + 2.0*B2) / 3.0;
        Bbar[2][cy] = (A2 - B2) / 3.0;
        Bbar[0][cz] = (A3 - B3) / 3.0;
        Bbar[1][cz] = (A3 - B3) / 3.0;
        Bbar[2][cz] = (A3 + 2.0*B3) / 3.0;

        // shear rows (unchanged)
        Bbar[3][cx] = B2;  Bbar[3][cy] = B1;
        Bbar[4][cy] = B3;  Bbar[4][cz] = B2;
        Bbar[5][cx] = B3;  Bbar[5][cz] = B1;
    }
}


void BezierTet10::computeVolumeAveragedDerivatives(double dN_avg[3][NEN],
                                                   double &volume) const
{
    //  B̄ₐ = (1/Vₑ) ∫_Ω ∂Nₐ/∂xᵢ dΩ  (4-pt rule, Kadapa Remark 3: full rule).
    for (int i = 0; i < 3; i++)
        for (int a = 0; a < NEN; a++)
            dN_avg[i][a] = 0.0;

    volume = 0.0;

    for (int gp = 0; gp < NGAUSS; gp++) {
        double w = GP4_w[gp];
        double dN[3][NEN];
        shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);

        double J[3][3];
        double dN_dx[3][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);

        double factor = w * fabs(detJ);
        volume += factor;

        for (int a = 0; a < NEN; a++)
            for (int i = 0; i < 3; i++)
                dN_avg[i][a] += factor * dN_dx[i][a];
    }

    if (volume > 0.0) {
        double invVol = 1.0 / volume;
        for (int a = 0; a < NEN; a++)
            for (int i = 0; i < 3; i++)
                dN_avg[i][a] *= invVol;
    }
}


// ═══════════════════════════════════════════════════════════════════
//  LAGRANGE → BÉZIER CONTROL POINT MAPPING
// ═══════════════════════════════════════════════════════════════════

void BezierTet10::computeControlPoints()
{
    //  Kadapa Eq. 11-13, edge-wise:
    //    vertices:  Pᵢ = Xᵢ                                  (i = 1..4)
    //    mid-edge:  P = 2[X_mid - 0.25 Xₐ - 0.25 X_b]        for edge (a,b)
    //  Straight-sided ⇒ P = X for every node.
    double X[NEN][3];
    for (int i = 0; i < NEN; i++) {
        const Vector &crds = theNodes[i]->getCrds();
        X[i][0] = crds(0);
        X[i][1] = crds(1);
        X[i][2] = crds(2);
    }

    // Vertices
    for (int i = 0; i < 4; i++)
        for (int d = 0; d < 3; d++)
            controlPts[i][d] = X[i][d];

    // Mid-edges
    for (int e = 0; e < 6; e++) {
        int a = edgeV[e][0], bnode = edgeV[e][1];
        for (int d = 0; d < 3; d++)
            controlPts[4 + e][d] =
                2.0 * (X[4 + e][d] - 0.25 * X[a][d] - 0.25 * X[bnode][d]);
    }
}


double BezierTet10::computeVolume() const
{
    double vol = 0.0;
    for (int gp = 0; gp < NGAUSS; gp++) {
        double dN[3][NEN];
        shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);
        double J[3][3];
        double dN_dx[3][NEN];
        double detJ = computeJacobian(dN, J, dN_dx);
        vol += GP4_w[gp] * fabs(detJ);
    }
    return vol;
}


// ═══════════════════════════════════════════════════════════════════
//  CHARACTERISTIC LENGTH  (crack-band regularization, e.g. ASDConcrete)
// ═══════════════════════════════════════════════════════════════════
//
//  Crack-band materials (ASDConcrete3D, etc.) call
//  ops_TheActiveElement->getCharacteristicLength() exactly once — on the
//  first setTrialStrain — to regularize the softening branch by element size.
//
//  Element's base default returns the MINIMUM inter-node distance. On a
//  quadratic element that distance is corner-to-mid-edge ≈ ½ the true edge
//  length, which under-estimates the band width and over-softens the response.
//
//  We instead return an element-size equivalent from the integrated volume:
//  the leg of a right tetrahedron (three mutually perpendicular legs) of
//  equal volume,
//
//      lch = cbrt(6 · V),
//
//  the 3D analogue of BezierTri6's sqrt(2·A). It recovers the true edge length
//  for a right-corner tet and is geometry-true for distorted straight-sided
//  elements; curved Bézier edges shift the factor only in the safe
//  (under-estimating) direction.
double BezierTet10::getCharacteristicLength(void)
{
    double V = this->computeVolume();
    if (V <= 0.0)
        return Element::getCharacteristicLength();  // degenerate: fall back
    return cbrt(6.0 * V);
}


// ═══════════════════════════════════════════════════════════════════
//  getInterpolationWeights (Ladruno, ADR 20 §9)
//
//  Quadratic Bernstein shape weights at the barycentric natural coordinate
//  xi = (L1,L2,L3), with L4 = 1-L1-L2-L3 (same node order as the GP/mass
//  interpolation, vertices N1..N4 then mid-edges N5..N10). Lets
//  LadrunoEmbeddedRebar embed a rebar node in this host without re-supplying
//  the weights by hand. N is resized to NEN (=10).
// ═══════════════════════════════════════════════════════════════════
int BezierTet10::getInterpolationWeights(const Vector &xi, Vector &N)
{
    if (xi.Size() < 3) {
        opserr << "BezierTet10::getInterpolationWeights - xi needs 3 barycentric "
                  "coords (L1,L2,L3)\n";
        return -1;
    }
    if (N.Size() != NEN)
        N.resize(NEN);
    double Nloc[NEN];
    this->shapeFunctions(xi(0), xi(1), xi(2), Nloc);
    for (int i = 0; i < NEN; i++)
        N(i) = Nloc[i];
    return 0;
}


// ═══════════════════════════════════════════════════════════════════
//  getInterpolationGradients (Ladruno, ADR 23 §3, Phase 2 UR)
//
//  Cartesian shape-function gradients dN_a/dx_j at the barycentric natural
//  coordinate xi = (L1,L2,L3), evaluated on the REFERENCE control points (the
//  same geometry getInterpolationWeights uses): dN/dx = J^-1 dN/dL via
//  computeJacobian. Used by LadrunoEmbeddedNode's rotation (UR) tie to read the
//  host continuum rotation θ = ½ curl(u) = skew(∇u) at the embedded point. dNdx
//  is resized to NEN(=10) x 3, dNdx(a,j) = ∂N_a/∂x_j. Returns -1 on a degenerate
//  (det J == 0) host.
// ═══════════════════════════════════════════════════════════════════
int BezierTet10::getInterpolationGradients(const Vector &xi, Matrix &dNdx)
{
    if (xi.Size() < 3) {
        opserr << "BezierTet10::getInterpolationGradients - xi needs 3 barycentric "
                  "coords (L1,L2,L3)\n";
        return -1;
    }
    if (dNdx.noRows() != NEN || dNdx.noCols() != 3)
        dNdx.resize(NEN, 3);

    double dN[3][NEN], J[3][3], dN_dx[3][NEN];
    this->shapeDerivatives(xi(0), xi(1), xi(2), dN);
    double detJ = this->computeJacobian(dN, J, dN_dx);
    if (fabs(detJ) <= 0.0) {
        opserr << "BezierTet10::getInterpolationGradients - degenerate Jacobian "
                  "(detJ=" << detJ << ") at the requested natural coordinate\n";
        return -1;
    }
    for (int a = 0; a < NEN; a++)
        for (int j = 0; j < 3; j++)
            dNdx(a, j) = dN_dx[j][a];   // dN_dx[0..2][a] = ∂N_a/∂x_j
    return 0;
}


// ═══════════════════════════════════════════════════════════════════
//  SERIALIZATION (sendSelf / recvSelf)
// ═══════════════════════════════════════════════════════════════════

int BezierTet10::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    // tag + 10 nodes + matClassTag + matDbTag + useBbar + cMass + geomID + fbarMode = 17
    // Ladruno: slot 15 carries the geometry-method id; slot 16 the F-bar variant.
    static ID iData(17);
    iData(0) = this->getTag();
    for (int i = 0; i < NEN; i++)
        iData(i + 1) = connectedExternalNodes(i);
    iData(11) = theMaterial[0]->getClassTag();
    int matDbTag = theMaterial[0]->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        for (int i = 0; i < NGAUSS; i++)
            theMaterial[i]->setDbTag(matDbTag);
    }
    iData(12) = matDbTag;
    iData(13) = useBbar ? 1 : 0;
    iData(14) = cMass ? 1 : 0;
    iData(15) = theGeom->getMethodID();   // Ladruno — geometry method
    iData(16) = fbarMode;                 // Ladruno — F-bar variant

    res += theChannel.sendID(this->getDbTag(), commitTag, iData);

    // rho, pressure, b[0..2] = 5
    static Vector dData(5);
    dData(0) = rho;
    dData(1) = pressure;
    dData(2) = b[0];
    dData(3) = b[1];
    dData(4) = b[2];

    res += theChannel.sendVector(this->getDbTag(), commitTag, dData);

    for (int i = 0; i < NGAUSS; i++)
        res += theMaterial[i]->sendSelf(commitTag, theChannel);

    return res;
}


int BezierTet10::recvSelf(int commitTag, Channel &theChannel,
                          FEM_ObjectBroker &theBroker)
{
    int res = 0;

    static ID iData(17);
    res += theChannel.recvID(this->getDbTag(), commitTag, iData);

    this->setTag(iData(0));
    for (int i = 0; i < NEN; i++)
        connectedExternalNodes(i) = iData(i + 1);
    int matClassTag = iData(11);
    int matDbTag = iData(12);
    useBbar = (iData(13) == 1);
    cMass = (iData(14) == 1);
    fbarMode = iData(16);                 // Ladruno — F-bar variant

    // Ladruno: rebuild the geometry method from the serialized id (Linear
    // fallback for an unknown/zero id) — else a parallel worker silently
    // reverts -geom corot to linear.
    if (theGeom != 0)
        delete theGeom;
    theGeom = SolidTransformation::create(iData(15));
    if (theGeom == 0)
        theGeom = new SolidTransformationLinear();

    static Vector dData(5);
    res += theChannel.recvVector(this->getDbTag(), commitTag, dData);

    rho = dData(0);
    pressure = dData(1);
    b[0] = dData(2);
    b[1] = dData(3);
    b[2] = dData(4);

    if (theMaterial == 0) {
        theMaterial = new NDMaterial*[NGAUSS];
        for (int i = 0; i < NGAUSS; i++)
            theMaterial[i] = 0;
    }

    for (int i = 0; i < NGAUSS; i++) {
        if (theMaterial[i] == 0) {
            theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
            if (theMaterial[i] == 0) {
                opserr << "BezierTet10::recvSelf - material creation failed\n";
                return -1;
            }
        }
        theMaterial[i]->setDbTag(matDbTag);
        res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    // Ladruno (ADR-77 review wave): cMass (a mass-formula branch) is
    // sig-exempt as construction-fixed, but recvSelf just rewrote it -- a
    // guard hit on a live element would serve the pre-recv mass structure.
    massCache.invalidate();
    // Ladruno (WP-124 C6): Ki was formed from the pre-recv material, thickness and
    // formulation; a restore into a LIVE element must not keep serving it.
    LadrunoShell::dropKi(Ki);

    return res;
}


// ═══════════════════════════════════════════════════════════════════
//  PRINT AND DISPLAY
// ═══════════════════════════════════════════════════════════════════

void BezierTet10::Print(OPS_Stream &s, int flag)
{
    if (flag == 2) {
        s << "{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"BezierTet10\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < NEN; i++) {
            s << connectedExternalNodes(i);
            if (i < NEN - 1) s << ", ";
        }
        s << "]}";
        return;
    }

    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nBezierTet10, element id:  " << this->getTag() << endln;
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tmass density:  " << rho << endln;
        s << "\tB-bar formulation: " << (useBbar ? "yes" : "no") << endln;
        s << "\tbody forces:  " << b[0] << " " << b[1] << " " << b[2] << endln;

        s << "\tStress (Gauss):" << endln;
        for (int i = 0; i < NGAUSS; i++)
            theMaterial[i]->Print(s, flag);
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"BezierTet10\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < NEN; i++) {
            s << connectedExternalNodes(i);
            if (i < NEN - 1) s << ", ";
        }
        s << "], ";
        s << "\"bbar\": " << (useBbar ? "true" : "false") << ", ";
        s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
    }
}


int BezierTet10::displaySelf(Renderer &theViewer, int displayMode,
                             float fact, const char **modes, int numMode)
{
    // Draw the 6 edges, each as 2 segments through its mid-edge node.
    static Vector v1(3), v2(3);

    for (int e = 0; e < 6; e++) {
        int seq[3] = {edgeV[e][0], 4 + e, edgeV[e][1]};
        for (int seg = 0; seg < 2; seg++) {
            int n1 = seq[seg];
            int n2 = seq[seg + 1];

            const Vector &c1 = theNodes[n1]->getCrds();
            const Vector &c2 = theNodes[n2]->getCrds();
            for (int i = 0; i < 3; i++) {
                v1(i) = c1(i);
                v2(i) = c2(i);
            }

            if (displayMode >= 0) {
                const Vector &d1 = theNodes[n1]->getDisp();
                const Vector &d2 = theNodes[n2]->getDisp();
                for (int i = 0; i < 3; i++) {
                    v1(i) += fact * d1(i);
                    v2(i) += fact * d2(i);
                }
            }

            theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
        }
    }

    return 0;
}


// ═══════════════════════════════════════════════════════════════════
//  RESPONSE AND RECORDING
// ═══════════════════════════════════════════════════════════════════

Response *BezierTet10::setResponse(const char **argv, int argc,
                                   OPS_Stream &output)
{
    Response *theResponse = 0;

    // ─── MPCO_Ladruno geometry self-declaration (contract Part A) ──
    //  topology=tet (simplex), family=bernstein, paramDomain=bary,
    //  rational=0 ⇒ no controlPointWeights. numCtrl=10, numGP=4.
    if (strcmp(argv[0], "basisInfo") == 0) {
        output.tag("ElementBasis");
        output.attr("topology",    "tet");
        output.attr("family",      "bernstein");
        output.attr("paramDomain", "bary");
        output.attr("rational",    0);
        output.attr("numCtrl",     NEN);      // 10 control points
        output.attr("numGP",       NGAUSS);   // 4 Gauss points (result stations)
        output.attr("orderU",      2);        // total polynomial degree (simplex)
        output.endTag();                      // ElementBasis
        return new ElementResponse(this, 101, ID(1));   // non-null sentinel
    }
    if (strcmp(argv[0], "integrationPoints") == 0)
        return new ElementResponse(this, 102, Matrix(NGAUSS, 3));
    if (strcmp(argv[0], "integrationWeights") == 0)
        return new ElementResponse(this, 103, Vector(NGAUSS));

    output.tag("ElementOutput");
    output.attr("eleType", "BezierTet10");
    output.attr("eleTag", this->getTag());
    for (int i = 0; i < NEN; i++) {
        char buf[16];
        sprintf(buf, "node%d", i + 1);
        output.attr(buf, connectedExternalNodes(i));
    }

    // ─── Element resisting forces ─────────────────────────────
    if (strcmp(argv[0], "force") == 0 ||
        strcmp(argv[0], "forces") == 0) {

        for (int i = 1; i <= NEN; i++) {
            char buf[32];
            sprintf(buf, "P1_%d", i); output.tag("ResponseType", buf);
            sprintf(buf, "P2_%d", i); output.tag("ResponseType", buf);
            sprintf(buf, "P3_%d", i); output.tag("ResponseType", buf);
        }
        theResponse = new ElementResponse(this, 1, Vector(NELD));
    }

    // ─── Stiffness matrix ─────────────────────────────────────
    else if (strcmp(argv[0], "stiff") == 0 ||
             strcmp(argv[0], "stiffness") == 0) {
        theResponse = new ElementResponse(this, 2, Matrix(NELD, NELD));
    }

    // ─── Material / integration point response ────────────────
    else if (strcmp(argv[0], "material") == 0 ||
             strcmp(argv[0], "integrPoint") == 0) {

        if (argc < 2) {
            opserr << "BezierTet10::setResponse - need GP number after 'material'\n";
            output.endTag();
            return 0;
        }

        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= NGAUSS) {
            output.tag("GaussPoint");
            output.attr("number", pointNum);
            output.attr("L1", GP4_L[pointNum-1][0]);
            output.attr("L2", GP4_L[pointNum-1][1]);
            output.attr("L3", GP4_L[pointNum-1][2]);

            double N[NEN];
            shapeFunctions(GP4_L[pointNum-1][0], GP4_L[pointNum-1][1],
                           GP4_L[pointNum-1][2], N);
            double x = 0.0, y = 0.0, z = 0.0;
            for (int a = 0; a < NEN; a++) {
                x += N[a] * controlPts[a][0];
                y += N[a] * controlPts[a][1];
                z += N[a] * controlPts[a][2];
            }
            output.attr("xLoc", x);
            output.attr("yLoc", y);
            output.attr("zLoc", z);

            theResponse = theMaterial[pointNum-1]->setResponse(
                &argv[2], argc-2, output);

            output.endTag();  // GaussPoint
        }
    }

    // ─── Stresses at ALL Gauss points ─────────────────────────
    else if (strcmp(argv[0], "stress") == 0 ||
             strcmp(argv[0], "stresses") == 0) {

        for (int i = 0; i < NGAUSS; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.tag("NdMaterialOutput");
            output.attr("classType", theMaterial[i]->getClassTag());
            output.attr("tag", theMaterial[i]->getTag());
            output.tag("ResponseType", "sigma_xx");
            output.tag("ResponseType", "sigma_yy");
            output.tag("ResponseType", "sigma_zz");
            output.tag("ResponseType", "sigma_xy");
            output.tag("ResponseType", "sigma_yz");
            output.tag("ResponseType", "sigma_zx");
            output.endTag();  // NdMaterialOutput
            output.endTag();  // GaussPoint
        }
        theResponse = new ElementResponse(this, 3, Vector(NSTRESS * NGAUSS));
    }

    // ─── Strains at ALL Gauss points ──────────────────────────
    else if (strcmp(argv[0], "strain") == 0 ||
             strcmp(argv[0], "strains") == 0) {

        for (int i = 0; i < NGAUSS; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.tag("NdMaterialOutput");
            output.tag("ResponseType", "eps_xx");
            output.tag("ResponseType", "eps_yy");
            output.tag("ResponseType", "eps_zz");
            output.tag("ResponseType", "gamma_xy");
            output.tag("ResponseType", "gamma_yz");
            output.tag("ResponseType", "gamma_zx");
            output.endTag();  // NdMaterialOutput
            output.endTag();  // GaussPoint
        }
        theResponse = new ElementResponse(this, 4, Vector(NSTRESS * NGAUSS));
    }

    // ─── Gauss point physical coordinates ─────────────────────
    else if (strcmp(argv[0], "gaussPoint") == 0 ||
             strcmp(argv[0], "gaussCoord") == 0 ||
             strcmp(argv[0], "gpCoord") == 0) {

        for (int i = 0; i < NGAUSS; i++) {
            output.tag("ResponseType", "x");
            output.tag("ResponseType", "y");
            output.tag("ResponseType", "z");
        }
        theResponse = new ElementResponse(this, 5, Vector(3 * NGAUSS));
    }

    // ─── Pressure contribution ────────────────────────────────
    else if (strcmp(argv[0], "pressure") == 0) {
        theResponse = new ElementResponse(this, 6, Vector(NELD));
    }

    // ─── Characteristic length ────────────────────────────────
    // Element-size lch used by crack-band materials (ASDConcrete). Lets the
    // user inspect the value the material regularizes with. See
    // getCharacteristicLength().
    else if (strcmp(argv[0], "charLength") == 0 ||
             strcmp(argv[0], "characteristicLength") == 0) {
        output.tag("ResponseType", "lch");
        theResponse = new ElementResponse(this, 7, Vector(1));
    }

    // Ladruno (WP-124 C2): close ElementOutput FIRST, then fall back to the base
    // vocabulary (globalForce, dampingForce, dynamicForce, inertialForce), which
    // this element never offered -- a recorder on those tokens wrote nothing.
    return LadrunoShell::finishResponse(this, theResponse, argv, argc, output);
}


int BezierTet10::getResponse(int responseID, Information &eleInfo)
{
    switch (responseID) {

    case 1:  // resisting forces
        return eleInfo.setVector(this->getResistingForce());

    case 2:  // tangent stiffness
        return eleInfo.setMatrix(this->getTangentStiff());

    case 3: {
        // Stresses at all Gauss points (committed)
        static Vector stressVec(NSTRESS * NGAUSS);
        for (int i = 0; i < NGAUSS; i++) {
            const Vector &sigma = theMaterial[i]->getStress();
            for (int j = 0; j < NSTRESS; j++)
                stressVec(i * NSTRESS + j) = sigma(j);
        }
        return eleInfo.setVector(stressVec);
    }

    case 4: {
        // Strains at all Gauss points (γ engineering shear)
        static Vector strainVec(NSTRESS * NGAUSS);
        for (int i = 0; i < NGAUSS; i++) {
            const Vector &eps = theMaterial[i]->getStrain();
            for (int j = 0; j < NSTRESS; j++)
                strainVec(i * NSTRESS + j) = eps(j);
        }
        return eleInfo.setVector(strainVec);
    }

    case 5: {
        // Gauss point physical coordinates  x(ξ) = Σ Nₐ(ξ) Pₐ
        static Vector gpCoords(3 * NGAUSS);
        for (int gp = 0; gp < NGAUSS; gp++) {
            double N[NEN];
            shapeFunctions(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], N);
            double x = 0.0, y = 0.0, z = 0.0;
            for (int a = 0; a < NEN; a++) {
                x += N[a] * controlPts[a][0];
                y += N[a] * controlPts[a][1];
                z += N[a] * controlPts[a][2];
            }
            gpCoords(3*gp)     = x;
            gpCoords(3*gp + 1) = y;
            gpCoords(3*gp + 2) = z;
        }
        return eleInfo.setVector(gpCoords);
    }

    case 6: {
        // Pressure contribution (volume hack, +z), 30-vector
        static Vector pressVec(NELD);
        pressVec.Zero();
        if (pressure != 0.0) {
            for (int gp = 0; gp < NGAUSS; gp++) {
                double w = GP4_w[gp];
                double N[NEN];
                shapeFunctions(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], N);
                double dN[3][NEN];
                shapeDerivatives(GP4_L[gp][0], GP4_L[gp][1], GP4_L[gp][2], dN);
                double J[3][3];
                double dN_dx[3][NEN];
                double detJ = computeJacobian(dN, J, dN_dx);
                double factor = w * fabs(detJ) * pressure;
                for (int a = 0; a < NEN; a++)
                    pressVec(3*a + 2) += factor * N[a];
            }
        }
        return eleInfo.setVector(pressVec);
    }

    case 7: {
        // ─── Characteristic length ─────────────────────────────
        // The element-size lch crack-band materials regularize with.
        static Vector lchVec(1);
        lchVec(0) = this->getCharacteristicLength();
        return eleInfo.setVector(lchVec);
    }

    // ─── MPCO_Ladruno geometry probes (contract Part A) ──────────
    case 101:
        // basisInfo sentinel — metadata emitted via the stream already.
        return 0;

    case 102: {
        // GP_PARAM — barycentric coords (L1,L2,L3) of the 4 GPs; L4=1-ΣL.
        static Matrix L(NGAUSS, 3);
        for (int g = 0; g < NGAUSS; g++) {
            L(g, 0) = GP4_L[g][0];
            L(g, 1) = GP4_L[g][1];
            L(g, 2) = GP4_L[g][2];
        }
        return eleInfo.setMatrix(L);
    }

    case 103: {
        // GP_WEIGHT — weights, same GP order (sum = 1/6 = ref tet volume).
        static Vector w(NGAUSS);
        for (int g = 0; g < NGAUSS; g++)
            w(g) = GP4_w[g];
        return eleInfo.setVector(w);
    }

    default:
        // Ladruno (WP-124 C2): the base IDs 111111..444444 (globalForce, ...)
        return this->Element::getResponse(responseID, eleInfo);
    }
}

// ═══════════════════════════════════════════════════════════════════
//  PARAMETER INTERFACE
// ═══════════════════════════════════════════════════════════════════
//
// setParameter is what `parameter` / `addToParameter $ptag element $tag
// <args>` reach. Three routes, mirroring SixNodeTri / LadrunoBrick (and
// BezierTri6, kept in lockstep):
//
//   "pressure"            → register THIS element (responseID 2); the
//                           Parameter drives `pressure` through
//                           updateParameter below.
//   "material" $gp <args> → delegate to the NDMaterial at GP $gp (1-based).
//   anything else         → broadcast to every GP material. This is the
//                           stress-control path: commitStressIncrementXX
//                           etc. register the MATERIALS with the Parameter,
//                           which then dispatches updateParameter straight
//                           to them (no element-level hop needed).
//
// Without this override, Element::setParameter returns -1 and
// addToParameter is a SILENT no-op — staged stress control (apeGmsh
// ops.initial_stress / STKO stressControl) does nothing through this
// element while FourNodeTetrahedron/brick hosts work fine.

int
BezierTet10::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;

    // "pressure" on the element itself (+z volume hack, as the response)
    if (strcmp(argv[0], "pressure") == 0)
        return param.addObject(2, this);

    // a specific Gauss-point material: "material $gp <args>"
    if ((strstr(argv[0], "material") != 0) &&
        (strcmp(argv[0], "materialState") != 0))
        return LadrunoShell::forwardToMaterialPoint(theMaterial, NGAUSS, false, argv, argc, param);

    // otherwise a forall-material parameter — broadcast to every GP
    return LadrunoShell::forwardToMaterials(theMaterial, NGAUSS, argv, argc, param);
}

int
BezierTet10::updateParameter(int parameterID, Information &info)
{
    switch (parameterID) {
    case 2:
        // pressure enters getResistingForce on the fly — no nodal-load
        // recompute needed.
        pressure = info.theDouble;
        return 0;
    default:
        return -1;
    }
}
