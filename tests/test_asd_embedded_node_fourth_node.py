"""ASDEmbeddedNodeElement with a tetrahedron host given from Python.

The optional 4th retained node is read as a string and converted to an
integer. In Python an integer argument was not returned as a string, so
`element('ASDEmbeddedNodeElement', tag, c, n1, n2, n3, n4)` silently dropped
n4 and embedded the node in the triangle n1-n2-n3.

Corner ux = 0, 1, 2, -1 on the unit tetrahedron, embedded node at
(0.2, 0.1, 0.15), barycentric (0.55, 0.2, 0.1, 0.15): ux = 0.25.
The triangle 1-2-3 alone gives 0.4.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


@pytest.fixture(autouse=True)
def clean_domain():
    ops.wipe()
    yield
    ops.wipe()


@pytest.mark.parametrize("n4", [4, "4"], ids=["int", "str"])
def test_tet4_host_uses_fourth_node(n4):
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, xyz in {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1)}.items():
        ops.node(tag, *map(float, xyz))
    ops.node(5, 0.2, 0.1, 0.15)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for tag, ux in {1: 0.0, 2: 1.0, 3: 2.0, 4: -1.0}.items():
        ops.sp(tag, 1, ux)
        ops.sp(tag, 2, 0.0)
        ops.sp(tag, 3, 0.0)
    ops.element("ASDEmbeddedNodeElement", 1, 5, 1, 2, 3, n4, "-K", 1.0e14)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormUnbalance", 1.0e-8, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(5, 1) == pytest.approx(0.25, abs=1.0e-8)
