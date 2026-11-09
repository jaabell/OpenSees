"""Mortar friction path state per (slave node, slave facet, master facet).

Each facet pair runs its own return map on its own weighted slip, in the tangent plane of
its own master facet, and owns its committed slip, engagement origin and tangential
multiplier. With one slot per slave node the pairs sharing a node would overwrite each
other's trial state and read back a traction from another facet's tangent plane.

Checks:

(a) Flat non-matching patch (slave facets straddle master facet boundaries), force-controlled
    uniform shear and compression, stick then slip: every slave node slides (Q - cap)/k with
    cap = min(mu P + c, tauMax) (1e-6), for three facet alignments and three cones.
(b) Creased (roof) interface, cohesion only, pressed down in three displacement-driven steps
    so that each flank slips down its slope; the ridge slave nodes are shared by pairs on
    both flanks. The contact force is Fx = 0, Fz = 2 A (p cos(alpha) + c sin(alpha)) with the
    flank pressure p from the per-commit normal update (1e-5); a model with the ridge nodes
    duplicated per flank gives the same forces (1e-6).
(c) The same roof under force control past the cohesion cap in four steps: every step
    converges and the ridge stays on the symmetry plane (1e-9).

Every test fails on the base with an unknown command.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops




# ---------------------------------------------------------------------------------------- (a)
def _grid(xs, ys, z, tag0):
    """Nodes of a structured grid; returns {(i,j): tag}."""
    tags = {}
    t = tag0
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            ops.node(t, float(x), float(y), float(z))
            tags[(i, j)] = t
            t += 1
    return tags


def _quads(tags, ni, nj):
    q = []
    for j in range(nj):
        for i in range(ni):
            q += [tags[(i, j)], tags[(i + 1, j)], tags[(i + 1, j + 1)], tags[(i, j + 1)]]
    return q


def _trib(xs, i):
    """1-D tributary length of grid node i (the row sum of the linear-element D)."""
    return (xs[min(i + 1, len(xs) - 1)] - xs[max(i - 1, 0)]) / 2.0


def _flat_patch(xs_s, ys_s, xs_m, ys_m, Q, P, mu, c, tmax, kx, nsteps):
    """Fixed master grid on z=0, element-less slave grid 1e-5 below it; consistent (tributary-area)
    nodal loads a_K*(Q, 0, -P); x springs a_K*kx to ground (full slip has no static equilibrium
    without them); y fixed. Returns the x displacements of every slave node, or None."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # the master overhangs the slave on every side, so the slid slave stays fully covered
    xs_m, ys_m = [-0.5] + list(xs_m) + [1.5], [-0.5] + list(ys_m) + [1.5]
    mt = _grid(xs_m, ys_m, 0.0, 1)
    for t in mt.values():
        ops.fix(t, 1, 1, 1)
    st = _grid(xs_s, ys_s, -1.0e-5, 1001)
    ops.contactSurface(1, "-master", 4, *_quads(mt, len(xs_m) - 1, len(ys_m) - 1))
    ops.contactSurface(2, "-slave-segments", 4, *_quads(st, len(xs_s) - 1, len(ys_s) - 1))
    opts = ["-mortar", "-epsN", 1.0e5, "-epsT", 1.0e5, "-outward", 0.0, 0.0, 1.0]
    if mu > 0.0:
        opts += ["-mu", mu]
    if c > 0.0:
        opts += ["-cohesion", c]
    if tmax > 0.0:
        opts += ["-tauMax", tmax]
    ops.contact(1, 1, 2, *opts)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    mat = 0
    for (i, j), t in st.items():
        a = _trib(xs_s, i) * _trib(ys_s, j)              # tributary area of a bilinear grid node
        ops.fix(t, 0, 1, 0)
        ops.load(t, a * Q, 0.0, -a * P)
        mat += 1
        ops.uniaxialMaterial("Elastic", mat, a * kx)
        g = 50000 + t
        ops.node(g, *ops.nodeCoord(t))
        ops.fix(g, 1, 1, 1)
        ops.element("zeroLength", 60000 + t, g, t, "-mat", mat, "-dir", 1)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1.0e-9, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            return None
    return [ops.nodeDisp(t, 1) for t in st.values()]


ALIGN = {
    "3x3-on-2x2": ([0, 1 / 3, 2 / 3, 1], [0, 1 / 3, 2 / 3, 1], [0, 0.5, 1], [0, 0.5, 1]),
    "2x2-on-3x3": ([0, 0.5, 1], [0, 0.5, 1], [0, 1 / 3, 2 / 3, 1], [0, 1 / 3, 2 / 3, 1]),
    "graded": ([0, 0.2, 0.55, 1], [0, 0.3, 1], [0, 0.45, 0.8, 1], [0, 0.6, 1]),
}
CONES = {  # name: (mu, cohesion, tauMax) -> cap at P = 100
    "cohesion": (0.0, 40.0, 0.0),
    "coulomb": (0.3, 0.0, 0.0),
    "coulomb+c capped": (0.3, 20.0, 35.0),
}


@pytest.mark.parametrize("align", list(ALIGN))
@pytest.mark.parametrize("cone", list(CONES))
def test_flat_nonmatching_patch_analytic(align, cone):
    """(a) stick then slip at cap = min(mu*P + c, tauMax); after the last step every slave node
    has slid (Q - cap)/kx, whatever the facet alignment."""
    P, kx, nsteps = 100.0, 1.0e5, 4
    mu, c, tmax = CONES[cone]
    cap = min(mu * P + c, tmax) if tmax > 0 else mu * P + c
    Q = 2.0 * cap                                        # stick in step 1 (Q/4 < cap), slip after
    u = _flat_patch(*ALIGN[align], Q, P, mu, c, tmax, kx, nsteps)
    assert u is not None, f"{align}/{cone}: Newton failed"
    expect = (Q - cap) / kx
    for ux in u:
        assert ux == pytest.approx(expect, rel=1e-6), f"{align}/{cone}: {ux} != (Q-cap)/k={expect}"


# ---------------------------------------------------------------------------------------- (b)
ALPHA = math.radians(20.0)
DELTA = 1.0e-6        # initial vertical interference (small: the re-clip geometry stays O(DELTA))
EPS = 1.0e9
COH = 10.0
KD = 1.0e15           # stiff drivers: displacement control through the Plain-replicating handler
                      # (driver compliance f/KD stays ~1e-7 of the imposed motion)


def _roof(split_ridge, steps, force_control=False, pz=0.0, eps=EPS, delta=DELTA, algo="Newton"):
    """Roof z = -tan(alpha)|x|, ridge along y at x=0. Master: 2 flanks x 1 quad, fixed. Slave:
    2 quads in x per flank x 3 in y (non-matching in x and y, symmetric about the ridge),
    lowered by DELTA. split_ridge
    duplicates the ridge slave nodes (left facets use one copy, right facets the other).
    Displacement control (default): every slave node sits on KD springs (x, y, z) to ground and is
    loaded -KD*s in z, so u_z -> -s. Force control: soft springs, a load -pz*lambda per node.
    algo: the solution algorithm (plain Newton unless a gate needs a globalised one).
    Returns (per-step contact force totals on the slave [(Fx, Fz)], per-step max ridge |u_x|, ok)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ta = math.tan(ALPHA)
    xm, ym = [-1.0, 0.0, 1.0], [0.0, 1.0]
    mt, t = {}, 1
    for j, y in enumerate(ym):
        for i, x in enumerate(xm):
            ops.node(t, x, y, -ta * abs(x))
            ops.fix(t, 1, 1, 1)
            mt[(i, j)] = t
            t += 1
    mq = []
    for i in range(2):
        mq += [mt[(i, 0)], mt[(i + 1, 0)], mt[(i + 1, 1)], mt[(i, 1)]]
    xs, ys = [-1.0, -0.45, 0.0, 0.45, 1.0], [0.0, 0.3, 0.75, 1.0]   # symmetric in x
    st, allslave, ridge = {}, [], []
    t = 101
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            for cp in (("L", "R") if (split_ridge and x == 0.0) else ("",)):
                ops.node(t, x, y, -ta * abs(x) - delta)
                st[(i, j, cp)] = t
                allslave.append(t)
                if x == 0.0:
                    ridge.append(t)
                t += 1

    def sn(i, j, side):
        return st.get((i, j, side), st.get((i, j, "")))

    sq = []
    for j in range(len(ys) - 1):
        for i in range(len(xs) - 1):
            side = "L" if xs[i + 1] <= 0.0 else "R"
            sq += [sn(i, j, side), sn(i + 1, j, side), sn(i + 1, j + 1, side), sn(i, j + 1, side)]
    ops.contactSurface(1, "-master", 4, *mq)
    ops.contactSurface(2, "-slave-segments", 4, *sq)
    ops.contact(1, 1, 2, "-mortar", "-epsN", eps, "-epsT", eps, "-cohesion", COH,
                "-outward", 0.0, 0.0, 1.0)
    kd = 1.0e3 if force_control else KD
    ops.uniaxialMaterial("Elastic", 1, kd)
    for s in allslave:
        g = 50000 + s
        ops.node(g, *ops.nodeCoord(s))
        ops.fix(g, 1, 1, 1)
        ops.element("zeroLength", 60000 + s, g, s, "-mat", 1, 1, 1, "-dir", 1, 2, 3)
    total = steps[-1]
    pload = pz if force_control else KD * total
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for s in allslave:
        ops.load(s, 0.0, 0.0, -pload)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    if force_control:
        ops.test("NormUnbalance", 1.0e-7, 60, 0)
    else:
        ops.test("NormDispIncr", 1.0e-15, 60, 0)
    ops.algorithm(algo)
    ops.analysis("Static")
    out, ux_ridge, lam_prev = [], [], 0.0
    for s in steps:
        lam = s / total
        ops.integrator("LoadControl", lam - lam_prev)
        lam_prev = lam
        if ops.analyze(1) != 0:
            return out, ux_ridge, False
        # contact force on the slave = -(applied load + spring force), summed over the nodes
        fx = fz = 0.0
        for n in allslave:
            u = ops.nodeDisp(n)
            fx += kd * u[0]
            fz += kd * u[2] + pload * lam
        out.append((fx, fz))
        ux_ridge.append(max(abs(ops.nodeDisp(n, 1)) for n in ridge))
    return out, ux_ridge, True


S_STEPS = [4.0e-7, 6.0e-7, 8.0e-7]    # tangential slip s*sin(alpha) >= 14 c/epsT: full slip


def _roof_analytic():
    """Per step (Fx, Fz) on the slave: each flank (area A) carries the pressure p along its normal
    and the cohesion cap c up its slope; p = -(lambda_N + epsN*g), g = -(DELTA + s)*cos(alpha) held
    by the drivers, lambda_N <- lambda_N + epsN*g at every commit (the normal Uzawa update)."""
    A, sa, ca = 1.0 / math.cos(ALPHA), math.sin(ALPHA), math.cos(ALPHA)
    lam, out = 0.0, []
    for s in S_STEPS:
        g = -(DELTA + s) * ca
        p = -(lam + EPS * g)
        out.append((0.0, 2.0 * A * (p * ca + COH * sa)))
        lam += EPS * g
    return out


def test_crease_pressed_roof_analytic():
    """(b) the shared-ridge roof reproduces the analytic force at every step (Fx = 0 by symmetry)."""
    got, _, ok = _roof(False, S_STEPS)
    assert ok, f"Newton failed after {len(got)} step(s)"
    for k, ((fx, fz), (ax, az)) in enumerate(zip(got, _roof_analytic())):
        assert abs(fx - ax) < 1e-5 * az, f"step {k + 1}: spurious Fx = {fx:.6f} (Fz {az:.3f})"
        assert fz == pytest.approx(az, rel=1e-5), f"step {k + 1}: Fz {fz:.6f} != analytic {az:.6f}"


def test_crease_shared_ridge_equals_split_ridge():
    """(b) sharing the ridge slave nodes between the flanks does not change the contact force."""
    shared, _, ok1 = _roof(False, S_STEPS)
    split, _, ok2 = _roof(True, S_STEPS)
    assert ok1 and ok2, f"Newton failed (shared ok={ok1}, split ok={ok2})"
    for k, ((fx1, fz1), (fx2, fz2)) in enumerate(zip(shared, split)):
        assert abs(fx1 - fx2) < 1e-6 * abs(fz2) and abs(fz1 - fz2) < 1e-6 * abs(fz2), (
            f"step {k + 1}: shared ({fx1:.6f}, {fz1:.6f}) != split ({fx2:.6f}, {fz2:.6f})")


# ---------------------------------------------------------------------------------------- (c)
def test_crease_force_control_multistep_cohesion():
    """(c) cohesion only, force control, 4 steps past the cap: every step converges and the ridge
    stays on the symmetry plane (u_x = 0). NewtonLineSearch: the symmetric tangent omits the
    variation of D and M of the thin cross-flank pairs, and plain Newton can cycle on their
    active set; the check is on the path state, not on the Newton rate."""
    got, ux, ok = _roof(False, [0.25, 0.5, 0.75, 1.0], force_control=True, pz=40.0,
                        eps=1.0e6, delta=1.0e-3, algo="NewtonLineSearch")
    assert ok, f"force-controlled cohesion roof failed after {len(got)} converged step(s)"
    assert max(ux) < 1e-9, f"ridge left the symmetry plane: max |u_x| = {max(ux):.3e}"
