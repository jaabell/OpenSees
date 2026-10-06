"""Singular LAPACK factorizations and unknown algorithm names report failure.

The singular model numbers a free, unconnected node first, so the first
pivot of the factorization is zero (LAPACK info == 1). That is the case the
old return code -info+1 turned into 0.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def _build_first_pivot_singular(system):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)    # free and unconnected: equation 0 has no stiffness
    ops.node(2, 1.0)
    ops.node(3, 2.0)
    ops.fix(2, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0)
    ops.element("truss", 1, 2, 3, 1.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0)
    ops.system(system)
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")


@pytest.mark.parametrize("system", ["BandGeneral", "FullGeneral", "BandSPD"])
def test_singular_first_pivot_fails_analyze(system):
    try:
        _build_first_pivot_singular(system)
        assert ops.analyze(1) != 0
    finally:
        ops.wipe()


@pytest.mark.parametrize("system", ["BandGeneral", "FullGeneral", "BandSPD"])
def test_nonsingular_control_still_solves(system):
    try:
        _build_first_pivot_singular(system)
        ops.fix(1, 1)
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        assert ops.nodeDisp(3, 1) == pytest.approx(1.0)
    finally:
        ops.wipe()


def test_unknown_algorithm_is_an_error():
    try:
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.algorithm("Linear")
        with pytest.raises(Exception):
            ops.algorithm("NoSuchAlgorithm")
    finally:
        ops.wipe()
