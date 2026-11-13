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
//   The 3D sibling of BezierTri6 (ADR 06_bezier_tet10, deferred under 04 D10).
//
// Based on:
//   Kadapa, C. "Novel quadratic Bézier triangular and tetrahedral elements
//   using existing mesh generators: Applications to linear nearly
//   incompressible elastostatics and implicit and explicit elastodynamics."
//   International Journal for Numerical Methods in Engineering,
//   2019; 117(5):543-573. doi:10.1002/nme.5967 (§5)
//
// Features:
//   - Quadratic Bernstein polynomial shape functions (all nonnegative)
//   - Pure displacement and B-bar formulations (Eq. 38 / 45, full 3D)
//   - Consistent and lumped mass matrices (all-positive lumped mass ρVe/10, Eq. 57)
//   - Lagrange-to-Bézier control point mapping (Eq. 11-13, edge-wise)
//   - 4-point tetrahedral quadrature (K / F / B-bar average)
//
// Node order (matches OpenSees TenNodeTetrahedron, jaabell/Larenas shp3d):
//   barycentric (L1,L2,L3,L4), L4 = 1-L1-L2-L3
//     N1=L1²  N2=L2²  N3=L3²  N4=L4²                          (4 vertices)
//     N5=2L1L2 (1-2)  N6=2L2L3 (2-3)  N7=2L1L3 (1-3)
//     N8=2L1L4 (1-4)  N9=2L3L4 (3-4)  N10=2L2L4 (2-4)         (6 mid-edges)
//
// v1 limitation: validated for STRAIGHT-SIDED meshes only (mid-edge nodes at
// edge midpoints ⇒ P = X). Curved elements emit a warning in setDomain.
//
// Integration uses |detJ| (not signed): BᵀDB / NᵀN / J⁻¹ are orientation-
// independent, so K/M/F are correct for either vertex handedness — the hedge
// for the open Gmsh-vs-node-order reconciliation (ADR O11).
//
// Element command:
//   element BezierTet10 $tag $nd1 ... $nd10 $matTag
//                       <-bbar> <-cMass> <-rho $r>
//                       <-bodyForce $b1 $b2 $b3> <-pressure $p>
//                       <-geom linear|corot|finite>
//
//   Geometry method (Ladruno): default linear (small strain). -geom corot adds
//   large-rotation / small-strain corotational kinematics via the shared
//   SolidTransformation layer (std + bbar; pressure unsupported under corot in
//   v1). -geom finite is genuine large-strain updated-Lagrangian: per GP it
//   builds F = I + Σ uₐ⊗∂Nₐ/∂X from the reference Bernstein gradients and drives
//   a FiniteStrainNDMaterial (e.g. nDMaterial LogStrain) via setTrialF(F),
//   assembling ∫Bᵀσ dv (Cauchy σ) + the full a_ijkl = c_ijkl − σ_il δ_jk tangent.
//   -geom finite -bbar = F-bar (dSNPO §15.1): drives the material with
//   F̄ = (Ĵ/J)^(1/3) F to cure near-incompressible volumetric locking; its
//   tangent is GENERALLY UNSYMMETRIC (use FullGeneral). The bar dilatation Ĵ
//   is selected by -fbar: centroid (default; Ĵ=J₀ at the tet centroid, the
//   LadrunoBrick form) or mean_dilatation (Ĵ=J̄=(∫J dV₀)/V₀, which reduces to
//   the element's mean-dilatation small-strain bbar). pressure is unsupported
//   under finite in v1.
//
//   Mass: default is the all-positive lumped mass ρVe/10 (Kadapa Eq. 57) for
//   explicit dynamics; pass -cMass for the consistent mass (implicit/eigen).
//
// Recorder responses:
//   "forces"              → 30 resisting force components
//   "stresses"            → σ_xx,yy,zz,xy,yz,zx at each GP (24 total)
//   "strains"             → ε_xx,yy,zz,γ_xy,yz,zx at each GP (24 total)
//   "gaussPoint"          → physical x,y,z of each GP (12 total)
//   "material" $gp <args> → delegate to NDMaterial at GP
//   "stiffness"           → 30×30 tangent stiffness
//   "charLength"          → element characteristic length lch (crack-band)

#ifndef BezierTet10_h
#define BezierTet10_h

#include <Element.h>
#include <LadrunoMassCache.h>   // Ladruno (ADR-77 G2 ext)
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <DRMHigherOrderNode.h>   // Ladruno (ADR-88)

class Node;
class NDMaterial;
class Response;
class SolidTransformation;   // Ladruno — geometry-method layer (linear/corot)

// Class tag ELE_TAG_BezierTet10 is defined in classTags.h (= 33001, ladruno band).

class BezierTet10 : public Element, public DRMHigherOrderNode
{
  public:
      // Ladruno (ADR-77 G2 ext): escape = -noMassCache
      void setMassCache(bool s) { massCache.setEnabled(s); }

      // Ladruno (ADR-88): H5DRM free-field interpolation at the 6 mid-edge
      // nodes (v1 straight-sided geometry -- weight 0.5/0.5 on the edge's own
      // 2 corner nodes, reusing edgeV, the same table the element already
      // uses to build its own Bernstein basis). See DRMHigherOrderNode.h and
      // ADR-88 (10ter) for why this is exact, not an approximation, under the
      // element's own straight-edge assumption.
      bool getDRMInterpolation(int localNode,
                                std::vector<int>& primaryLocalNodes,
                                std::vector<double>& weights) const override;

    // ─── F-bar variant ids (bbar + -geom finite) — public so the OPS factory
    // can map the -fbar option. CENTROID = single centroid dilatation J₀ (dSNPO
    // eq 15.5, the LadrunoBrick form). MEAN = volume-averaged J̄ = (∫J dV₀)/V₀
    // (reduces to the element's mean-dilatation small-strain bbar).
    static constexpr int FBAR_CENTROID = 0;
    static constexpr int FBAR_MEAN     = 1;

    // ─── Constructors and Destructor ───────────────────────────

    // Full constructor
    BezierTet10(int tag,
                int nd1, int nd2, int nd3, int nd4, int nd5,
                int nd6, int nd7, int nd8, int nd9, int nd10,
                NDMaterial &m, double rho = 0.0,
                double b1 = 0.0, double b2 = 0.0, double b3 = 0.0,
                bool useBbar = false, bool cMass = false,
                double pressure = 0.0,
                int geomMethodID = 0,    // Ladruno — 0 = SolidTransformation::METHOD_LINEAR
                int fbarMode = 0);       // Ladruno — F-bar variant: 0=centroid, 1=mean-dilatation

    // Null constructor (for parallel/database reconstruction)
    BezierTet10();

    // Destructor
    ~BezierTet10();

    // ─── Element Interface ────────────────────────────────────

    const char *getClassType(void) const { return "BezierTet10"; }
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);

    // Element-size characteristic length for crack-band regularization
    // (e.g. ASDConcrete). Overrides Element's min-inter-node-distance default,
    // which on a quadratic element collapses to ~½ the edge length. See .cpp.
    double getCharacteristicLength(void);

    // Ladruno (ADR 20 §9): quadratic Bernstein shape weights at barycentric
    // natural coord xi = (L1,L2,L3), L4 = 1-L1-L2-L3, for embedded-reinforcement
    // coupling. N sized to 10 (NEN).
    int getInterpolationWeights(const Vector &xi, Vector &N);

    // Ladruno (ADR 23 §3, Phase 2 UR): cartesian shape gradients dN_a/dx_j at the
    // barycentric natural coord xi = (L1,L2,L3) (reference control points), for the
    // node-embedding ROTATION tie. dNdx sized to 10x3, dNdx(a,j) = dN_a/dx_j.
    int getInterpolationGradients(const Vector &xi, Matrix &dNdx);

    void setDomain(Domain *theDomain);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad(void);
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    // Parameter interface (stress-control staging / sensitivity)
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    int displaySelf(Renderer &, int mode, float fact,
                    const char **displayModes = 0, int numModes = 0);

  protected:

  private:
    // ─── Constants ────────────────────────────────────────────
    static constexpr int NEN = 10;     // Number of element nodes
    static constexpr int NDOF = 3;     // DOFs per node
    static constexpr int NELD = 30;    // Total element DOFs (10×3)
    static constexpr int NSTRESS = 6;  // Stress/strain components (3D)
    static constexpr int NGAUSS = 4;   // Gauss points (degree-2 tet rule)

    // ─── Computational Methods ────────────────────────────────

    // Shape functions and derivatives at barycentric (L1, L2, L3)
    void shapeFunctions(double L1, double L2, double L3,
                        double N[NEN]) const;

    void shapeDerivatives(double L1, double L2, double L3,
                          double dN[3][NEN]) const;

    // Jacobian, its (signed) determinant, and physical derivatives
    double computeJacobian(const double dN[3][NEN],
                           double J[3][3],
                           double dN_dx[3][NEN]) const;

    // B and B-bar (strain-displacement) matrices
    void computeBMatrix(const double dN_dx[3][NEN],
                        double B[NSTRESS][NELD]) const;

    void computeBBarMatrix(const double dN_dx[3][NEN],
                           const double dN_avg[3][NEN],
                           double Bbar[NSTRESS][NELD]) const;

    // Volume-averaged derivatives for B-bar (Eq. 46)
    void computeVolumeAveragedDerivatives(double dN_avg[3][NEN],
                                          double &volume) const;

    // Map Lagrange node positions to Bézier control points
    void computeControlPoints();

    // Element volume
    double computeVolume() const;

    // ─── Geometry-method (corot) seams ────────────────────────
    // Ladruno: refresh theGeom from current geometry and return the localized
    // (core-frame) 30-dof trial displacement. Identity for -geom linear.
    const Vector &computeLocalDisp(void);

    // Ladruno: build the strain-displacement matrix B (or B̄) at Gauss point gp
    // and return the integration measure w·|detJ|, or a NEGATIVE value if the
    // element is degenerate (detJ == 0). The single guarded B-assembly path
    // shared by update / formCore / getInitialStiff so the degenerate check and
    // the std-vs-bbar choice can never drift between them.
    double formBAtGauss(int gp, const double dN_avg[3][NEN],
                        double B[NSTRESS][NELD]) const;

    // Ladruno: one Gauss pass building the CORE-frame internal force
    // fInt = ∫ Bᵀσ dΩ (always) and, when tangFlag != 0, the core tangent
    // K = ∫ BᵀDB dΩ — NO body force / pressure / Q. Used by BOTH getResistingForce
    // (tangFlag 0) and getTangentStiff (tangFlag 1), so the fCore fed to
    // globalizeStiff is byte-identical to the one globalizeForce rotates, by
    // construction (not by comment). K may be null when tangFlag == 0.
    void formCore(int tangFlag, Vector &fInt, Matrix *K);

    // ─── Geometry-method (finite, -geom finite) seams ─────────────
    // Ladruno: true when theGeom reports a DeformationGradient strain measure
    // (SolidTransformationFinite). The element then runs the updated-Lagrangian
    // path below instead of the small-strain formCore/formBAtGauss path, and
    // skips the localize/globalize seams (identity for finite).
    bool isFinite(void) const;

    // Ladruno: F (row-major [9]) = δ_iJ + Σₐ uₐ[i] ∂Nₐ/∂X_J from the nodal trial
    // displacements and the REFERENCE Bernstein shape gradients dN_dX[J][a]
    // (= computeJacobian's dN_dx, since controlPts are the reference nodes for
    // straight-sided elements). Returns det F. Total F — the material derives F_Δ.
    double deformationGradient(const double dN_dX[3][NEN], double F[9]) const;

    // Ladruno: -geom finite update — per GP build F and drive the material via
    // setTrialF(F). Returns < 0 on a degenerate Jacobian or det F ≤ 0 so the
    // analysis step-cuts instead of assembling a negative-volume contribution.
    int updateFinite(void);

    // Ladruno: -geom finite assembly — ONE Gauss pass building the internal force
    // fInt = ∫ σ_ij ∂Nₐ/∂x_j dv (always; current config, dv = J·|detJ_ref|·w) and,
    // when tangFlag != 0, the consistent tangent K = ∫ (∂Nₐ/∂x_j) a_ijkl (∂N_b/∂x_l)
    // dv with a_ijkl = c_ijkl − σ_il δ_jk (c via getSpatialTangentTensor — the 6×6
    // getTangent is LOSSY in (k,l); see FiniteStrainNDMaterial.h). NO body force /
    // pressure / Q (applied in the global frame in getResistingForce). One helper
    // feeds both force and tangent so f/K share a single Gauss pass.
    //
    // bbar + finite = F-bar (dSNPO §15.1): updateFinite drives the material with
    // F̄ = (J₀/J)^(1/3) F so every GP shares the centroid dilatation J₀ (the
    // volumetric-locking cure); the residual is unchanged (eq 15.9 — only σ̄
    // changes), and the tangent gains the eq 15.10 coupling ∫ Gᵀ q (G₀−G) dv with
    // q_ij = (1/3) a_ijpp − (2/3) σ̄_ij (eq 15.11, GENERALLY UNSYMMETRIC).
    void formResidAndTangentFinite(int tangFlag, Vector &fInt, Matrix *K);

    // Ladruno: F-bar centroid data (bbar + finite). Returns J₀ = det F₀, F₀ the
    // deformation gradient at the tet centroid (barycentric L=(¼,¼,¼); dSNPO
    // eq 15.5). If G0 != 0, also fills the centroid spatial-gradient operator
    // G0[k][b] = ∂N_b/∂x_k|_centroid (from F₀⁻¹) for the eq 15.10 coupling.
    // Returns 0.0 on a degenerate centroid Jacobian so the caller's J₀≤0 guard fires.
    double centroidFbar(double (*G0)[NEN] = 0) const;

    // Ladruno: F-bar MEAN-DILATATION data (the -fbar mean_dilatation variant).
    // Returns J̄ = (∫ J dV₀)/V₀ (reference-volume average of det F over the GPs).
    // If Gbar != 0, fills the volume-averaged spatial gradient operator
    // Gbar[k][b] = (∫ ∂N_b/∂x_k dv)/v (current-volume average) for the eq 15.10
    // coupling. Returns 0.0 on a degenerate/inverted GP so the J̄≤0 guard fires.
    // This is the consistent analogue of centroidFbar with (J₀,G₀)→(J̄,Ḡ).
    double fbarMeanDilatation(double (*Gbar)[NEN] = 0) const;

    // ─── Static Quadrature Data ───────────────────────────────
    // 4-point rule (degree 2) for stiffness / force / B-bar average.
    static const double GP4_L[][3];   // barycentric (L1,L2,L3), L4=1-ΣL
    static const double GP4_w[];      // weights (sum = 1/6 = ref tet volume)

    // 1D Gauss-Legendre (4-pt) on [0,1], collapsed (Duffy) into a 4×4×4
    // (degree-7) tet rule for the consistent mass — the degree-4 NₐNᵦ
    // integrand becomes degree 6 after the (1-a)² Duffy Jacobian (ADR O12).
    static const double GL4_t[];      // nodes on [0,1]
    static const double GL4_w[];      // weights on [0,1]

    // Edge → (vertexA, vertexB), 0-based, for mid-edge nodes 5..10 (idx 4..9)
    static const int edgeV[6][2];

    // ─── Static Return Objects (Petracca/Abell pattern) ──────
    static Matrix K_return;   // 30×30 stiffness return matrix
    static Matrix M_return;   // 30×30 mass return matrix
    static Vector P_return;   // 30×1 force return vector

    // ─── Member Data ──────────────────────────────────────────
    NDMaterial **theMaterial;   // Array of NGAUSS material pointers

    ID connectedExternalNodes;  // 10 node tags
    Node *theNodes[NEN];        // 10 node pointers
    LadrunoMassCache massCache;   // Ladruno (ADR-77 G2 ext): per-instance mass cache,
                                  // guard-checked (rho/thickness/coords); -noMassCache escape


    double controlPts[NEN][3];  // Bézier control points (P = X if straight)

    double rho;                 // Mass density (0 = take from material)
    double pressure;            // Applied "pressure" (volume hack, +z, as Tri6)
    double b[3];                // Body force per unit volume {b1, b2, b3}

    bool useBbar;               // B-bar for near-incompressibility
    bool cMass;                 // true = consistent mass; false = lumped (ρVe/10)
    int  fbarMode;              // Ladruno — F-bar variant (FBAR_CENTROID / FBAR_MEAN)

    Vector Q;                   // 30×1 applied load vector

    double appliedB[3];         // current body force = loadFactor·data·b
    int applyLoad;              // 1 once a SelfWeight load has been added

    Matrix *Ki;                 // cached initial stiffness (lazy, LadrunoShell::cacheKi)

    SolidTransformation *theGeom;  // Ladruno — geometry method (linear/corot)

    static int numInstances;    // One-time print flag (Abell pattern)
};

#endif // BezierTet10_h
