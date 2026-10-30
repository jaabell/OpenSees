"""gradientInelasticBeamColumn section assignment with hinge-type integration rules.

The constructor chooses which of the three section objects (end I, interior,
end J) sits at each integration point from the relative locations returned by
the integration rule. It called getSectionLocations with the element length
member, still 0.0 at that point (the length is only known in setDomain).
HingeMidpoint divides by the length, so with L = 0 the locations of a 4 point
rule with equal hinge lengths become

  xi = [ +inf, NaN, NaN, -inf ]

and the sections were assigned as [end J, end I, end I, end I].

Model: one element of length 2.0 along x, HingeMidpoint rule with hinge
lengths lpI = lpJ = 0.1 and sections tagged 1 (end I), 2 (end J), 3 (interior),
elastic with E = 1000 and Iz = 1, 2, 3 respectively, so the bending stiffness
E*Iz of a section identifies it (1000, 2000, 3000).

With the unit length used for the assignment the locations are
xi = [0.05, 0.2691, 0.7309, 0.95]; with the default section length ratios of
0.1 the points are assigned to [end I, interior, interior, end J]:

  EIz per integration point = [1000, 3000, 3000, 2000]

Before the fix the sections reported [2000, 1000, 1000, 1000]. The 2D and 3D
elements are checked. The value is read with the element response
"section <n> stiffness".
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E = 1000.0
IZ = {1: 1.0, 2: 2.0, 3: 3.0}
EXPECTED = [1000.0, 3000.0, 3000.0, 2000.0]


def _section_ei(ndm):
    ops.wipe()
    ndf = 3 if ndm == 2 else 6
    ops.model("basic", "-ndm", ndm, "-ndf", ndf)
    ops.node(1, 0.0, 0.0, *([0.0] if ndm == 3 else []))
    ops.node(2, 2.0, 0.0, *([0.0] if ndm == 3 else []))
    ops.fix(1, *([1] * ndf))
    for tag, iz in IZ.items():
        if ndm == 2:
            ops.section("Elastic", tag, E, 10.0, iz)
        else:
            ops.section("Elastic", tag, E, 10.0, iz, iz, 400.0, 1.0)
    ops.beamIntegration("HingeMidpoint", 1, 1, 0.1, 2, 0.1, 3)
    if ndm == 2:
        ops.geomTransf("Linear", 1)
    else:
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("gradientInelasticBeamColumn", 1, 1, 2, 1, 1, 0.2)

    values = []
    for n in range(1, 5):
        k = ops.eleResponse(1, "section", n, "stiffness")
        dim = int(round(math.sqrt(len(k))))
        values.append(k[dim + 1])  # row 1, column 1 of the section matrix: E*Iz
    return values


@pytest.mark.parametrize("ndm", [2, 3])
def test_hinge_rule_section_assignment(ndm):
    values = _section_ei(ndm)
    assert values == pytest.approx(EXPECTED, rel=1.0e-12)
