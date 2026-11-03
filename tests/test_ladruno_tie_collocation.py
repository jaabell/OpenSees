"""Tests for the LadrunoTie collocation tie (static Lagrange and Penalty).

LadrunoTie projects each slave node onto its nearest master facet, evaluates
the master shape functions N_i at the projection and writes one
EQ_Constraint per slave node per DOF, u_s = sum N_i u_i, which the Lagrange
and Penalty handlers enforce.

Mesh: a non-conforming solid-solid interface at x = 1. The master block is
one stdBrick hexahedron, x in [0, 1]; the slave block is two stdBrick
hexahedra stacked in z, x in [1, 2]. The slave interface has 6 nodes; the two
at z = 0.5 sit at the middle of the master facet edges (weights 0.5/0.5).

- Patch test: u_x = gamma z imposed on every node except the tied ones gives
  the same constant shear stress E gamma / 2 (nu = 0) at every Gauss point of
  both blocks, with Lagrange and with Penalty.
- Split column: a uniform stretch of the tied column reproduces the
  conforming two-element column (stress and the free mid-plane displacement).
- Tie rows: under a general load the tied nodes satisfy u_s = sum N_i u_i.
- Refusals: a node on both surfaces, a slave listed twice, a slave off the
  master surface.
- A slave node without mass is accepted: lumped mass on the tied DOFs is not
  required under Lagrange or Penalty.
- Classic Tcl: the LadrunoTie command of the OpenSees Tcl executable gives the
  same patch-test displacement (skipped when no executable is found).
"""
import os
import subprocess
from pathlib import Path

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# a refused LadrunoTie command raises the interpreter error (an unknown
# command would raise AttributeError instead and fail these tests)
OPS_ERROR = getattr(ops, "OpenSeesError", RuntimeError)


MASTER_FACET = [2, 3, 7, 6]
SLAVE_NODES = [11, 14, 15, 18, 19, 22]
# slave : (master node, weight) pairs of the expected tie row
TIE_RELATIONS = {
    11: [(2, 1.0)], 14: [(3, 1.0)], 19: [(6, 1.0)], 22: [(7, 1.0)],
    15: [(2, 0.5), (6, 0.5)], 18: [(3, 0.5), (7, 0.5)],
}

_MASTER = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
           5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_SLAVE = {11: (1, 0, 0.0), 12: (2, 0, 0.0), 13: (2, 1, 0.0), 14: (1, 1, 0.0),
          15: (1, 0, 0.5), 16: (2, 0, 0.5), 17: (2, 1, 0.5), 18: (1, 1, 0.5),
          19: (1, 0, 1.0), 20: (2, 0, 1.0), 21: (2, 1, 1.0), 22: (1, 1, 1.0)}
_ALL = {**_MASTER, **_SLAVE}


def _mesh(E=1000.0, nu=0.0, tie=True):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _ALL.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, 0.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 2, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    ops.element("stdBrick", 3, 15, 16, 17, 18, 19, 20, 21, 22, 1)
    if tie:
        ops.LadrunoTie("-slaveNodes", len(SLAVE_NODES), *SLAVE_NODES,
                       "-masterFacets", 4, 1, *MASTER_FACET)


def _gp_stresses(ele):
    ops.eleResponse(ele, "forces")
    s = list(ops.eleResponse(ele, "stresses"))
    return [s[i * 6:(i + 1) * 6] for i in range(len(s) // 6)]


def _static_solve(handler):
    if handler == "Penalty":
        ops.constraints("Penalty", 1.0e12, 1.0e12)
    else:
        ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


@pytest.mark.parametrize("handler,rtol", [("Lagrange", 1e-8), ("Penalty", 1e-5)])
def test_constant_stress_patch(handler, rtol):
    """u_x = gamma z on the boundary gives one constant stress in both blocks;
    the z = 0.5 slave nodes need the 0.5/0.5 facet weights."""
    E, gamma = 1000.0, 1e-3
    _mesh(E=E, nu=0.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_NODES:
            continue                          # governed by the tie
        ops.fix(t, 0, 1, 1)
        ops.sp(t, 1, gamma * float(z))
    _static_solve(handler)

    allgp = _gp_stresses(1) + _gp_stresses(2) + _gp_stresses(3)
    ref = allgp[0]
    target = E * gamma / 2.0
    atol = rtol * target
    for v in allgp:
        for si, ri in zip(v, ref):
            assert abs(si - ri) <= atol, f"non-constant stress: {v} vs {ref}"
    comps = [abs(c) for c in ref]
    assert abs(max(comps) - target) <= atol, f"stress {ref} != {target}"
    assert sum(1 for c in comps if c > 1e-3 * target) == 1, f"not pure shear: {ref}"


def _monolithic_axial(E, U):
    """Conforming two-hexahedron column x in [0, 2] with u_x = 0 at x = 0 and
    u_x = U at x = 2. Returns (stress at a Gauss point, u_x at x = 1)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nodes, tag = {}, 1
    for xi in (0.0, 1.0, 2.0):
        for (y, z) in [(0, 0), (1, 0), (1, 1), (0, 1)]:
            nodes[tag] = (xi, y, z)
            ops.node(tag, float(xi), float(y), float(z))
            tag += 1
    ops.nDMaterial("ElasticIsotropic", 1, E, 0.0, 0.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 2, 5, 6, 7, 8, 9, 10, 11, 12, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in nodes.items():
        ops.fix(t, 0, 1, 1)
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, 1, U)
    _static_solve("Transformation")
    return _gp_stresses(1)[0], ops.nodeDisp(5, 1)


@pytest.mark.parametrize("handler,rtol", [("Lagrange", 1e-6), ("Penalty", 1e-5)])
def test_split_column_matches_conforming(handler, rtol):
    E, U = 1000.0, 0.02
    ref_stress, ref_mid = _monolithic_axial(E, U)

    _mesh(E=E, nu=0.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if t in SLAVE_NODES:
            continue
        ops.fix(t, 0, 1, 1)
        if abs(x) < 1e-12:
            ops.sp(t, 1, 0.0)
        elif abs(x - 2.0) < 1e-12:
            ops.sp(t, 1, U)
    _static_solve(handler)

    mid = ops.nodeDisp(2, 1)                  # free master node at x = 1
    assert abs(mid - ref_mid) <= rtol * abs(ref_mid), f"mid-plane {mid} != {ref_mid}"
    got = _gp_stresses(1)[0]
    for gi, ri in zip(got, ref_stress):
        assert abs(gi - ri) <= rtol * E * U, f"stress {got} != {ref_stress}"


def test_tie_rows_hold_under_general_load():
    """Clamped master base, transverse and axial loads on the slave end: the
    tied nodes satisfy u_s = sum N_i u_i on every DOF."""
    _mesh(E=1000.0, nu=0.3)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _ALL.items():
        if abs(x) < 1e-12:
            ops.fix(t, 1, 1, 1)
    ops.load(20, 5.0, 1.0, -2.0)
    ops.load(21, 5.0, -1.0, 3.0)
    _static_solve("Lagrange")
    umax = max(abs(ops.nodeDisp(20, d)) for d in (1, 2, 3))
    assert umax > 1e-6, "slave end did not move"
    for s, rel in TIE_RELATIONS.items():
        for d in (1, 2, 3):
            us = ops.nodeDisp(s, d)
            um = sum(w * ops.nodeDisp(m, d) for (m, w) in rel)
            assert abs(us - um) <= 1e-10 * umax, f"tie row violated at node {s} dof {d}"


def _bare_master():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _MASTER.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.0, 0.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)


def test_refuse_non_disjoint_surfaces():
    _bare_master()
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-slaveNodes", 1, 2, "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_duplicate_slave():
    _bare_master()
    ops.node(99, 1.0, 0.5, 0.5)
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-slaveNodes", 2, 99, 99, "-masterFacets", 4, 1, *MASTER_FACET)


def test_refuse_off_surface_slave():
    _bare_master()
    ops.node(99, 1.5, 0.5, 0.5)               # 0.5 off the x = 1 facet plane
    with pytest.raises(OPS_ERROR):
        ops.LadrunoTie("-slaveNodes", 1, 99, "-masterFacets", 4, 1, *MASTER_FACET)


def test_massless_slave_accepted_static():
    """A slave node with no element and no mass, on the master facet: the tie
    is accepted and Lagrange carries the node with the facet (bilinear
    interpolation at (y, z) = (0.25, 0.5))."""
    _bare_master()
    ops.node(99, 1.0, 0.25, 0.5)
    ops.LadrunoTie("-slaveNodes", 1, 99, "-masterFacets", 4, 1, *MASTER_FACET)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    field = {2: 0.0, 3: 0.4, 7: 1.0, 6: 0.2}  # u_x imposed on the facet nodes
    for t in _MASTER:
        if t in field:
            ops.fix(t, 0, 1, 1)
            ops.sp(t, 1, field[t])
        else:
            ops.fix(t, 1, 1, 1)
    _static_solve("Lagrange")
    y, z = 0.25, 0.5
    expect = ((1 - y) * (1 - z) * field[2] + y * (1 - z) * field[3]
              + y * z * field[7] + (1 - y) * z * field[6])
    assert abs(ops.nodeDisp(99, 1) - expect) <= 1e-12
    assert abs(ops.nodeDisp(99, 2)) <= 1e-12
    assert abs(ops.nodeDisp(99, 3)) <= 1e-12


TCL_SCRIPT = """
model basic -ndm 3 -ndf 3
set coords {
 1 0 0 0  2 1 0 0  3 1 1 0  4 0 1 0  5 0 0 1  6 1 0 1  7 1 1 1  8 0 1 1
 11 1 0 0  12 2 0 0  13 2 1 0  14 1 1 0  15 1 0 0.5  16 2 0 0.5  17 2 1 0.5
 18 1 1 0.5  19 1 0 1  20 2 0 1  21 2 1 1  22 1 1 1
}
foreach {t x y z} $coords { node $t $x $y $z }
nDMaterial ElasticIsotropic 1 1000.0 0.0 0.0
element stdBrick 1 1 2 3 4 5 6 7 8 1
element stdBrick 2 11 12 13 14 15 16 17 18 1
element stdBrick 3 15 16 17 18 19 20 21 22 1
LadrunoTie -slaveNodes 6 11 14 15 18 19 22 -masterFacets 4 1 2 3 7 6
timeSeries Constant 1
pattern Plain 1 1 {
  foreach {t x y z} $coords {
    if {[lsearch {11 14 15 18 19 22} $t] >= 0} { continue }
    fix $t 0 1 1
    sp $t 1 [expr 1.0e-3 * $z]
  }
}
constraints Lagrange
numberer Plain
system FullGeneral
test NormDispIncr 1e-12 20
algorithm Linear
integrator LoadControl 1.0
analysis Static
set ok [analyze 1]
puts "RESULT $ok [nodeDisp 15 1] [nodeDisp 18 1] [nodeDisp 22 1]"
"""


def _tcl_executable():
    env = os.environ.get("OPENSEES_EXE")
    if env and Path(env).is_file():
        return Path(env)
    here = Path(ops.__file__).resolve().parent
    for name in ("OpenSees.exe", "OpenSees"):
        exe = here / name
        if exe.is_file():
            return exe
    return None


@pytest.mark.skipif(_tcl_executable() is None,
                    reason="no OpenSees Tcl executable (set OPENSEES_EXE)")
def test_classic_tcl_command(tmp_path):
    script = tmp_path / "ladruno_tie.tcl"
    script.write_text(TCL_SCRIPT)
    res = subprocess.run([str(_tcl_executable()), str(script)],
                         capture_output=True, text=True, cwd=tmp_path,
                         timeout=300, stdin=subprocess.DEVNULL,
                         env=dict(os.environ))
    out = res.stdout + res.stderr
    assert "invalid command name" not in out, out
    line = [x for x in out.splitlines() if x.startswith("RESULT")]
    assert line, out
    ok, u15, u18, u22 = (float(v) for v in line[-1].split()[1:5])
    assert ok == 0.0
    assert abs(u15 - 0.5e-3) <= 1e-12
    assert abs(u18 - 0.5e-3) <= 1e-12
    assert abs(u22 - 1.0e-3) <= 1e-12
