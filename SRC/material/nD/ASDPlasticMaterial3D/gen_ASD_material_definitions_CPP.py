#!/usr/bin/python
from itertools import product

# ============================================================================
# Standard Elasticity Models
# ============================================================================
EL = [
    "LinearIsotropic3D_EL",
    "DuncanChang_EL",
]

# ============================================================================
# Standard Yield Functions
# ============================================================================
YF = [
    "VonMises_YF",
    "DruckerPrager_YF",
    "MohrCoulomb_YF",
    "HoekBrown_YF",
    # "TensionCutoff_YF",
]

# ============================================================================
# Standard Plastic Flow Directions
# ============================================================================
PF = [
    "VonMises_PF",
    "DruckerPrager_PF",
    "ConstantDilatancy_PF",
    "MohrCoulomb_PF",
    "HoekBrown_PF"
]

# ============================================================================
# Possible combinations of IV with YFs
# ============================================================================
IV_YF = {}

IV_YF["VonMises_YF"] = [
    "BackStress<TensorLinearHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>",
    "BackStress<ArmstrongFrederickHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>",
]

IV_YF["DruckerPrager_YF"] = [
    "BackStress<TensorLinearHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>",
    "BackStress<ArmstrongFrederickHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>",
]

IV_YF["MohrCoulomb_YF"] = [
    "BackStress<NullHardeningTensorFunction>"
]

IV_YF["HoekBrown_YF"] = [
    "BackStress<NullHardeningTensorFunction>"
]

IV_YF["TensionCutoff_YF"] = [
    "BackStress<NullHardeningTensorFunction>"
]

# ============================================================================
# Possible combinations of IV with PFs
# ============================================================================
IV_PF = {
    "VonMises_PF": [
    "BackStress<NullHardeningTensorFunction>",
    "BackStress<TensorLinearHardeningFunction>",
    "BackStress<ArmstrongFrederickHardeningFunction>",
    ],
    "DruckerPrager_PF":
    [
    "BackStress<TensorLinearHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>",
    "BackStress<ArmstrongFrederickHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>",
    ]
}

IV_PF["ConstantDilatancy_PF"] = IV_PF["VonMises_PF"]

IV_PF["MohrCoulomb_PF"] = [
    "BackStress<NullHardeningTensorFunction>"
]

IV_PF["HoekBrown_PF"] = [
    "BackStress<NullHardeningTensorFunction>"
]


# ============================================================================
# Template for standard models (YF and PF have separate IV template arguments)
# ============================================================================
template = """
createASDPlasticMaterial3D<
        {EL}, 
        {YF}<
            {IV_YF}
            >, 
        {PF}<
            {IV_PF}
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);
"""


# ============================================================================
# STIFF SOIL MODEL DEFINITIONS
# ============================================================================

# Template for Stiff Soil models (YF and PF share the same IV type)
stiffsoil_template = """
createASDPlasticMaterial3D<
        {EL}, 
        {YF}<
            {IV}
            >, 
        {PF}<
            {IV}
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);
"""

# Stiff Soil model configurations
# Each entry is (Elasticity, YieldFunction, PlasticFlow, InternalVariable)
STIFFSOIL_MODELS = [
    # Deviatoric (Shear) mechanism
    {
        "EL": "StiffSoil_EL",
        "YF": "StiffSoilShear_YF",
        "PF": "StiffSoilShear_PF",
        "IV": "EpsQpShear"
    },
    # Cap (Volumetric) mechanism
    {
        "EL": "StiffSoil_EL",
        "YF": "StiffSoilCap_YF",
        "PF": "StiffSoilCap_PF",
        "IV": "CapPressure"
    },
    # Cap mechanism with linear hardening (alternative)
    {
        "EL": "StiffSoil_EL",
        "YF": "StiffSoilCap_YF",
        "PF": "StiffSoilCap_PF",
        "IV": "CapPressureLinear"
    },
]


# ============================================================================
# Generate output files
# ============================================================================
# The combinations are split over N_UNITS compilation units so they compile in
# parallel (each template instantiation is expensive). The split is contiguous,
# so the order of the combinations (and of the `list` output) is unchanged.
# N_UNITS is fixed so CMakeLists.txt / Makefile do not change when combinations
# are added or removed.
N_UNITS = 8

calls = []
for el, yf, pf in product(EL, YF, PF):
    for iv_yf in IV_YF[yf]:
        for iv_pf in IV_PF[pf]:
            calls.append(template.format(EL=el, YF=yf, PF=pf, IV_YF=iv_yf, IV_PF=iv_pf))
for model in STIFFSOIL_MODELS:
    calls.append(stiffsoil_template.format(**model))

SIGNATURE = ("(int instance_tag, const char* yf_type, const char* pf_type, const char* el_type, "
             "const char* iv_type, std::list<NDMaterial*> &instance_pointers, "
             "std::list<model_spec_t> &available_models)")
HEADER = ("// Auto-generated by gen_ASD_material_definitions_CPP.py - do not edit manually\n"
          "// Edit the lists in gen_ASD_material_definitions_CPP.py and re-run it instead.\n\n")

n = len(calls)
for k in range(N_UNITS):
    chunk = calls[k * n // N_UNITS:(k + 1) * n // N_UNITS]
    with open("ASD_material_definitions_%d.cpp" % k, "w") as fid:
        fid.write(HEADER)
        fid.write("#ifdef _EIGEN3\n\n#include \"ASDPlasticMaterial3DFactory.h\"\n\n")
        fid.write("void ASDPlasticMaterial3D_definitions_%d%s\n{\n" % (k, SIGNATURE))
        for c in chunk:
            fid.write(c)
        fid.write("}\n\n#endif // _EIGEN3\n")

with open("ASD_material_definitions.h", "w") as fid:
    fid.write(HEADER)
    fid.write("#ifndef ASD_material_definitions_h\n#define ASD_material_definitions_h\n\n")
    for k in range(N_UNITS):
        fid.write("void ASDPlasticMaterial3D_definitions_%d%s;\n" % (k, SIGNATURE))
    fid.write("\n// Instantiates every combination, in order\n")
    fid.write("inline void ASDPlasticMaterial3D_all_definitions%s\n{\n" % SIGNATURE)
    for k in range(N_UNITS):
        fid.write("    ASDPlasticMaterial3D_definitions_%d(instance_tag, yf_type, pf_type, el_type, iv_type, "
                  "instance_pointers, available_models);\n" % k)
    fid.write("}\n\n#endif // ASD_material_definitions_h\n")

print("Generated %d combinations in %d units (ASD_material_definitions_*.cpp) and ASD_material_definitions.h"
      % (n, N_UNITS))
