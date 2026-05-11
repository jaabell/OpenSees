// ============================================================================
// Hoek-Brown Model Parameters
// Based on Hoek & Brown (2018) - "The Hoek-Brown failure criterion and GSI - 2018 edition"
//
// The generalized Hoek-Brown criterion is:
//   σ₁ = σ₃ + σci * (mb * σ₃/σci + s)^a
//
// For rock masses, parameters mb, s, a are derived from intact rock parameter mi,
// Geological Strength Index (GSI), and Disturbance factor D as:
//   mb = mi * exp((GSI - 100) / (28 - 14*D))
//   s  = exp((GSI - 100) / (9 - 3*D))
//   a  = 0.5 + (1/6) * (exp(-GSI/15) - exp(-20/3))
//
// For intact rock (GSI=100, D=0): mb = mi, s = 1, a = 0.5
// ============================================================================

// HB_sigci: Unconfined compressive strength of intact rock (σci)
// Units: Stress (e.g., MPa, kPa)
// Typical range: 1-250 MPa depending on rock type
struct HB_sigci_Name { static constexpr const char* name = "HB_sigci";};
using HB_sigci = ModelParameterType<double, HB_sigci_Name>;

// HB_mi: Material constant for intact rock
// Dimensionless parameter related to the frictional characteristics of intact rock
// Typical values: 
//   4-8 for very soft clay-rich rocks
//   8-12 for soft sedimentary rocks
//   12-18 for medium strength rocks
//   18-25 for hard ignite rocks  
//   25-35 for very hard metamorphic and igneous rocks
struct HB_mi_Name { static constexpr const char* name = "HB_mi";};
using HB_mi = ModelParameterType<double, HB_mi_Name>;

// HB_GSI: Geological Strength Index
// Dimensionless index describing rock mass quality (0-100)
// GSI = 100 for intact rock, lower values for more fractured/disturbed rock
struct HB_GSI_Name { static constexpr const char* name = "HB_GSI";};
using HB_GSI = ModelParameterType<double, HB_GSI_Name>;

// HB_D: Disturbance factor
// Dimensionless factor accounting for blast damage and stress relaxation (0-1)
// D = 0 for undisturbed in-situ rock mass
// D = 1 for very disturbed rock (poor blasting)
struct HB_D_Name { static constexpr const char* name = "HB_D";};
using HB_D = ModelParameterType<double, HB_D_Name>;

// HB_mb: Reduced material constant for rock mass
// Derived parameter: mb = mi * exp((GSI - 100) / (28 - 14*D))
// This can be provided directly or computed from mi, GSI, D
struct HB_mb_Name { static constexpr const char* name = "HB_mb";};
using HB_mb = ModelParameterType<double, HB_mb_Name>;

// HB_s: Rock mass parameter s
// Derived parameter: s = exp((GSI - 100) / (9 - 3*D))
// For intact rock (GSI=100): s = 1
// For highly disturbed rock: s → 0
struct HB_s_Name { static constexpr const char* name = "HB_s";};
using HB_s = ModelParameterType<double, HB_s_Name>;

// HB_a: Rock mass parameter a
// Derived parameter: a = 0.5 + (1/6) * (exp(-GSI/15) - exp(-20/3))
// For intact rock or high quality rock mass: a ≈ 0.5
// For very poor quality rock mass: a → 0.65
struct HB_a_Name { static constexpr const char* name = "HB_a";};
using HB_a = ModelParameterType<double, HB_a_Name>;

// HB_mb_psi: Dilation parameter for plastic flow (non-associated flow rule)
// Similar to mb but controls volumetric plastic strain
// Typically mb_psi < mb to prevent unrealistic dilation
// Set mb_psi = 0 for zero dilatancy
struct HB_mb_psi_Name { static constexpr const char* name = "HB_mb_psi";};
using HB_mb_psi = ModelParameterType<double, HB_mb_psi_Name>;

// HB_ds: Relative stress increment for numerical differentiation
// Small perturbation parameter for computing derivatives
// Typical value: 1e-6 to 1e-8 times the stress level
struct HB_ds_Name { static constexpr const char* name = "HB_ds";};
using HB_ds = ModelParameterType<double, HB_ds_Name>;
