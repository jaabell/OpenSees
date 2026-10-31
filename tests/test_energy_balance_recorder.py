"""EnergyBalance recorder: column layouts, element mass, regions, -ownedNodes.

Six-column layout (with -time a leading time column):

    time  KE  IE  DW  ULW  RES  ERR%

-v2 layout:

    time  KE_ele KE_nod IE DW_ele DW_nod ULW  [channels]  RES ERR%

Checked:
  * element mass: with Truss -rho (element mass only) KE at the first record
    is close to 1/2 M v0^2, i.e. it is read from the element mass matrix and
    not from the damping matrix (many elements return the same static
    scratch matrix from getMass() and getDamp()).
  * -time and its alias -load give 7 columns per row.
  * -v2 without any energy producer gives the 6 split columns plus RES and
    ERR%; the element mass sits in KE_ele and KE_nod is zero.
  * -region adds one six-column block per region; a region holding every
    element and node reproduces the whole-model block.
  * -ownedNodes restricts the nodal terms to the nodes of a region: KE_nod
    equals 1/2 sum m v^2 over the owned nodes only.
"""
import math

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def build_bar(N, L, E, A, rho):
    """Fixed-free chain of N truss elements along x; returns element length."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    le = L / N
    for i in range(N + 1):
        ops.node(i + 1, i * le, 0.0)
    ops.fix(1, 1, 1)
    for i in range(1, N + 1):
        ops.fix(i + 1, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(N):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", rho)
    return le


def transient_analysis():
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


def read_rows(efile):
    rows = []
    with open(efile) as fh:
        for line in fh:
            vals = line.split()
            if vals:
                rows.append([float(x) for x in vals])
    return rows


def test_element_mass_ke(tmp_path):
    N, L, E, A, rho, v0 = 20, 10.0, 100.0, 1.0, 1.0, 2.0
    efile = str(tmp_path / "bar_energy.txt")
    build_bar(N, L, E, A, rho)
    for i in range(1, N + 1):
        ops.setNodeVel(i + 1, 1, v0, "-commit")
    ops.recorder("EnergyBalance", "-file", efile, "-time")
    transient_analysis()
    ops.analyze(2, 1e-4)
    ops.wipe()

    d = np.atleast_2d(np.loadtxt(efile))
    ke0 = d[0, 1]
    ke_expected = 0.5 * rho * A * L * v0 ** 2
    assert ke0 > 0.5 * ke_expected, (ke0, ke_expected)


@pytest.mark.parametrize("flag", ["-time", "-load"])
def test_time_column_layout(tmp_path, flag):
    N, L, E, A, rho = 5, 1.0, 100.0, 1.0, 1.0
    efile = str(tmp_path / "cols.txt")
    le = build_bar(N, L, E, A, rho)
    ops.setNodeVel(N + 1, 1, 1.0, "-commit")
    ops.recorder("EnergyBalance", "-file", efile, flag)
    transient_analysis()
    ops.analyze(5, 0.5 * le / math.sqrt(E / rho))
    ops.wipe()

    rows = read_rows(efile)
    assert len(rows) == 5
    # time + KE IE DW ULW RES ERR%
    assert all(len(r) == 7 for r in rows), sorted({len(r) for r in rows})


def test_v2_layout_no_producers(tmp_path):
    N, L, E, A, rho = 10, 1.0, 100.0, 1.0, 1.0
    efile = str(tmp_path / "v2_plain.txt")
    le = build_bar(N, L, E, A, rho)
    for i in range(1, N + 1):
        ops.setNodeVel(i + 1, 1, 1.0, "-commit")
    ops.recorder("EnergyBalance", "-file", efile, "-time", "-v2")
    transient_analysis()
    ops.analyze(50, 0.2 * le / math.sqrt(E / rho))
    ops.wipe()

    d = np.atleast_2d(np.loadtxt(efile))
    ncol = d.shape[1]
    # time + 6 split columns + RES ERR%; any channel column (declared by an
    # integrator of this model) must be zero for this run
    assert ncol >= 9, ncol
    if ncol > 9:
        assert np.allclose(d[:, 7:ncol - 2], 0.0, atol=1e-14)
    ke_ele, ke_nod = d[0, 1], d[0, 2]
    assert ke_ele > 0.0 and abs(ke_nod) < 1e-12, (ke_ele, ke_nod)


def test_region_block_matches_model(tmp_path):
    N, L, E, A, rho = 6, 1.0, 100.0, 1.0, 1.0
    efile = str(tmp_path / "region.txt")
    le = build_bar(N, L, E, A, rho)
    # add nodal mass so that the region block carries nodal terms as well
    for i in range(1, N + 1):
        ops.mass(i + 1, 0.1, 0.1)
    # region 7: every element (and so every node); region 8: first half
    ops.region(7, "-ele", *range(1, N + 1))
    ops.region(8, "-ele", *range(1, N // 2 + 1))
    for i in range(1, N + 1):
        ops.setNodeVel(i + 1, 1, 0.5 * i / N, "-commit")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(N + 1, 1.0, 0.0)
    ops.recorder("EnergyBalance", "-file", efile, "-time",
                 "-region", 7, "-region", 8, "-precision", 12)
    transient_analysis()
    ops.analyze(40, 0.2 * le / math.sqrt(E / rho))
    ops.wipe()

    d = np.atleast_2d(np.loadtxt(efile))
    # time + model block + two region blocks
    assert d.shape[1] == 1 + 6 * 3, d.shape
    model = d[:, 1:7]
    whole = d[:, 7:13]
    half = d[:, 13:19]
    scale = np.abs(model[:, :4]).max()
    assert scale > 0.0
    # a region holding the whole model reproduces the model block
    np.testing.assert_allclose(whole, model, rtol=0.0, atol=1e-9 * scale)
    # the loaded tip node is outside the first half: no external work there,
    # and its kinetic energy is a proper part of the total
    assert np.abs(half[:, 3]).max() == 0.0
    assert np.abs(model[:, 3]).max() > 0.0
    assert np.all(half[:, 0] < model[:, 0])


def test_owned_nodes_restricts_nodal_terms(tmp_path):
    efile = str(tmp_path / "owned.txt")
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    masses = {1: 1.0, 2: 2.0, 3: 3.0}
    for tag in masses:
        ops.node(tag, float(tag), "-mass", masses[tag])
    ops.node(10, 0.0)
    ops.fix(10, 1)
    ops.uniaxialMaterial("Elastic", 1, 50.0)
    for tag in masses:
        ops.element("zeroLength", tag, 10, tag, "-mat", 1, "-dir", 1)
    vel = {1: 1.0, 2: -0.5, 3: 0.25}
    for tag in masses:
        ops.setNodeVel(tag, 1, vel[tag], "-commit")
    ops.region(5, "-node", 1, 3)
    ops.recorder("EnergyBalance", "-file", efile, "-time", "-v2",
                 "-ownedNodes", 5, "-precision", 15)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    vrows = []
    for _ in range(5):
        assert ops.analyze(1, 1.0e-3) == 0
        vrows.append({tag: ops.nodeVel(tag, 1) for tag in masses})
    ops.wipe()

    d = np.atleast_2d(np.loadtxt(efile))
    assert d.shape[0] == 5
    for k in range(5):
        ke_owned = sum(0.5 * masses[t] * vrows[k][t] ** 2 for t in (1, 3))
        ke_all = sum(0.5 * masses[t] * vrows[k][t] ** 2 for t in masses)
        assert abs(d[k, 2] - ke_owned) < 1e-10 * ke_all, (d[k, 2], ke_owned)
        assert abs(d[k, 2] - ke_all) > 0.1 * ke_all
        assert d[k, 1] == 0.0   # zeroLength carries no element mass
