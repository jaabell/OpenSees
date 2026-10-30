"""Tests for the "localAxes" element response of the frame elements.

ElasticBeam2d/3d, DispBeamColumn2d/3d, ForceBeamColumn2d/3d,
MixedBeamColumn2d/3d and GradientInelasticBeamColumn2d/3d answer
"localAxes" with the 9 direction cosines [x(3) y(3) z(3)] of the element
local frame, taken from the coordinate transformation.

Checked for every element, on an oblique member:
  * the response has 9 finite components,
  * the three axes are unit vectors and mutually orthogonal,
  * the local x axis is the normalized chord from node i to node j,
  * in 3D, the local y axis is z x x with z in the plane of x and vecxz,
    as defined by geomTransf; in 2D, z is the out-of-plane unit vector,
  * where the element forwards unknown tokens to the transformation, the
    packed frame matches the existing "xaxis", "yaxis" and "zaxis" responses.

The Linear, PDelta and Corotational transformations are covered.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


E = 2.0e5
G = 8.0e4
A = 100.0
IZ = 1.0e4
IY = 8.0e3
J = 1.5e4
NP = 5      # Lobatto points (at least 3 for the gradient inelastic rule)
LC = 0.1    # gradient inelastic characteristic length

TOL = 1.0e-12

NODE_I_3D = (0.0, 0.0, 0.0)
NODE_J_3D = (1.0, 2.0, 2.0)    # length 3, oblique
VECXZ = (0.0, 0.0, 1.0)
NODE_I_2D = (0.0, 0.0)
NODE_J_2D = (3.0, 4.0)         # length 5, oblique

ELEMENTS = ["elastic", "disp", "force", "mixed", "gradient"]
TRANSFS = ["Linear", "PDelta", "Corotational"]

# GradientInelasticBeamColumn does not forward unknown tokens to the
# transformation, so "xaxis"/"yaxis"/"zaxis" are not available on it.
FORWARDS_TO_TRANSF = {"elastic", "disp", "force", "mixed"}


def _norm(v):
    return math.sqrt(sum(c * c for c in v))


def _unit(v):
    n = _norm(v)
    return [c / n for c in v]


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def _close(a, b):
    return all(abs(x - y) <= TOL for x, y in zip(a, b))


def _add_element(kind, ndm):
    if kind == "elastic":
        if ndm == 2:
            ops.element("elasticBeamColumn", 10, 1, 2, A, E, IZ, 1)
        else:
            ops.element("elasticBeamColumn", 10, 1, 2, A, E, G, J, IY, IZ, 1)
        return
    if ndm == 2:
        ops.section("Elastic", 1, E, A, IZ)
    else:
        ops.section("Elastic", 1, E, A, IZ, IY, G, J)
    ops.beamIntegration("Lobatto", 1, 1, NP)
    if kind == "disp":
        ops.element("dispBeamColumn", 10, 1, 2, 1, 1)
    elif kind == "force":
        ops.element("forceBeamColumn", 10, 1, 2, 1, 1)
    elif kind == "mixed":
        ops.element("mixedBeamColumn", 10, 1, 2, 1, 1)
    else:
        ops.element("gradientInelasticBeamColumn", 10, 1, 2, 1, 1, LC)


def _build(kind, transf, ndm):
    ops.wipe()
    if ndm == 2:
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, *NODE_I_2D)
        ops.node(2, *NODE_J_2D)
        ops.geomTransf(transf, 1)
    else:
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, *NODE_I_3D)
        ops.node(2, *NODE_J_3D)
        ops.geomTransf(transf, 1, *VECXZ)
    _add_element(kind, ndm)


def _check_frame(la):
    assert len(la) == 9, f"localAxes must have 9 components, got {len(la)}"
    assert all(math.isfinite(c) for c in la)
    x, y, z = la[0:3], la[3:6], la[6:9]
    for v in (x, y, z):
        assert abs(_norm(v) - 1.0) <= TOL, f"axis {v} is not a unit vector"
    assert abs(_dot(x, y)) <= TOL
    assert abs(_dot(x, z)) <= TOL
    assert abs(_dot(y, z)) <= TOL
    return x, y, z


@pytest.mark.parametrize("transf", TRANSFS)
@pytest.mark.parametrize("kind", ELEMENTS)
def test_local_axes_3d(kind, transf):
    _build(kind, transf, 3)
    la = ops.eleResponse(10, "localAxes")
    x, y, z = _check_frame(la)

    chord = [NODE_J_3D[i] - NODE_I_3D[i] for i in range(3)]
    x_ref = _unit(chord)
    y_ref = _unit(_cross(VECXZ, x_ref))
    z_ref = _cross(x_ref, y_ref)
    assert _close(x, x_ref), f"x {x} != {x_ref}"
    assert _close(y, y_ref), f"y {y} != {y_ref}"
    assert _close(z, z_ref), f"z {z} != {z_ref}"

    if kind in FORWARDS_TO_TRANSF:
        assert _close(x, ops.eleResponse(10, "xaxis"))
        assert _close(y, ops.eleResponse(10, "yaxis"))
        assert _close(z, ops.eleResponse(10, "zaxis"))


@pytest.mark.parametrize("transf", TRANSFS)
@pytest.mark.parametrize("kind", ELEMENTS)
def test_local_axes_2d(kind, transf):
    _build(kind, transf, 2)
    la = ops.eleResponse(10, "localAxes")
    x, y, z = _check_frame(la)

    c = _unit([NODE_J_2D[0] - NODE_I_2D[0], NODE_J_2D[1] - NODE_I_2D[1]])
    assert _close(x, [c[0], c[1], 0.0])
    assert _close(y, [-c[1], c[0], 0.0])
    assert _close(z, [0.0, 0.0, 1.0])

    if kind in FORWARDS_TO_TRANSF:
        assert _close(x, ops.eleResponse(10, "xaxis"))
        assert _close(y, ops.eleResponse(10, "yaxis"))
        assert _close(z, ops.eleResponse(10, "zaxis"))
