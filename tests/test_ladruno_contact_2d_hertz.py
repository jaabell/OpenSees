"""2D Hertz line contact: a rigid cylinder pressed into an elastic half-plane (node-to-segment).

Hertz (Johnson 1985, Section 5.2) for a cylinder of radius R on a half-plane, per unit thickness:

    b^2 = 4 P' R / (pi E*),   p_max = 2 P' / (pi b),   1/E* = (1 - nu^2)/E (rigid cylinder).

A 2D half-plane has a log-divergent surface displacement, so the indentation that produces a
given P' depends on the domain size. The model is driven by an interference delta, P' is measured
(the sum of the contact forces) and the half-width and the peak pressure are compared with Hertz at
the measured P'.

Model: unit-thickness plane-strain quads (E = 1e4, nu = 0.3) on [-2, 2] x [-2, 0], graded away
from a uniform fine zone around the contact, fixed at the base and on rollers at the sides. The
top row is the slave node set; the rigid cylinder is a fixed chain of 2-node master segments
(-outward 0 -1), offset by half a slave spacing so that no slave projects on a master vertex, with
a facet length of the order of delta. A first pass on a provisional mesh measures b, the second
refines around it (2b about 16 elements).

Checks: the summed contact force balances the base reaction (1e-3); ladrunoContactForce equals
kn times the geometric penetration on every active pair (1e-9); the contact half-width a (the
outermost penetrating node) is within one element of b; the peak nodal pressure is within 10
percent of p_max.

Reference: Johnson, K.L. (1985). Contact Mechanics. Cambridge University Press.
doi:10.1017/CBO9781139171731

The file fails on the base with an unknown command.
"""
import math

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


R = 1.0
DELTA = 5.0e-3
E = 1.0e4
NU = 0.3
W2 = 2.0
H = 2.0
KFAC = 1.0e4        # kn = KFAC E dx
RATIO = 1.15
FINE = 2.5


def _hertz2d(P, R, E, nu):
    Es = E / (1.0 - nu * nu)
    b = math.sqrt(4.0 * P * R / (math.pi * Es))
    return b, 2.0 * P / (math.pi * b)


def _graded(x_fine, n_fine, x_far, ratio):
    h = x_fine / n_fine
    xs = [i * h for i in range(n_fine + 1)]
    x, hh = x_fine, h
    while x < x_far - 1e-12:
        hh *= ratio
        x = min(x + hh, x_far)
        xs.append(x)
    return xs


def _run(b_hint, nb):
    dx = b_hint / nb
    n_fine = max(int(round(FINE * b_hint / dx)), 4)
    x_fine = n_fine * dx
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    xh = _graded(x_fine, n_fine, W2, RATIO)
    xs = [-v for v in reversed(xh[1:])] + xh
    yh = _graded(x_fine, n_fine, H, RATIO)
    ys = [-v for v in reversed(yh)]
    nx, ny = len(xs) - 1, len(ys) - 1
    nid, t_ = {}, 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(t_, xs[i], ys[j])
            nid[(i, j)] = t_
            t_ += 1
    base = [nid[(i, 0)] for i in range(nx + 1)]
    for t in base:
        ops.fix(t, 1, 1)
    for j in range(1, ny + 1):
        for i in (0, nx):
            ops.fix(nid[(i, j)], 1, 0)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    eid = 1
    for j in range(ny):
        for i in range(nx):
            ops.element("quad", eid, nid[(i, j)], nid[(i + 1, j)],
                        nid[(i + 1, j + 1)], nid[(i, j + 1)], 1.0, "PlaneStrain", 1)
            eid += 1
    # rigid cylinder: facets of length hm ~ delta (the projection on a curved
    # master is lost beyond about two facet lengths of penetration), offset by
    # dx/2, vertices lowered by the chord sagitta
    q = max(1, int(round(DELTA / dx)))
    hm = q * dx
    sag = hm * hm / (8.0 * R)

    def yarc(x):
        return (R - DELTA) - math.sqrt(max(R * R - x * x, 0.0))
    arc_half = FINE * b_hint
    k0 = int(math.ceil((arc_half - 0.5 * dx) / hm))
    axs = [(-k0 + k) * hm + 0.5 * dx for k in range(2 * k0 + 1)]
    ays = [yarc(x) - sag for x in axs]
    mn, mt = [], 800001
    for x, y in zip(axs, ays):
        ops.node(mt, x, y)
        ops.fix(mt, 1, 1)
        mn.append(mt)
        mt += 1
    chain = []
    for i in range(len(mn) - 1):
        chain += [mn[i], mn[i + 1]]
    ops.contactSurface(10, "-master", 2, *chain)
    lim = axs[-1] - 2.0 * hm
    slaves, sx, sw = [], {}, {}
    for i in range(nx + 1):
        if abs(xs[i]) > lim:
            continue
        t = nid[(i, ny)]
        slaves.append(t)
        sx[t] = xs[i]
        sw[t] = 0.5 * (xs[i + 1] - xs[i - 1])
    ops.contactSurface(20, "-slave", *slaves)
    kn = KFAC * E * dx
    ops.contact(1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, -1.0)
    assert min(min(abs(sx[t] - a) for a in axs) for t in slaves) > 0.25 * dx
    ops.constraints("LadrunoContact")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1.0e-10, 200, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    rc = 0
    for _ in range(2):
        rc = ops.analyze(1)
        if rc != 0:
            break

    def chord(xq):
        k = min(max(int(math.floor((xq - axs[0]) / hm)), 0), len(axs) - 2)
        s = (ays[k + 1] - ays[k]) / hm
        return ays[k] + (xq - axs[k]) * s, s
    ops.reactions()
    Rbase = sum(ops.nodeReaction(t, 2) for t in base)
    act, qerr, nq = [], 0.0, 0
    for t in slaves:
        ux, uy = ops.nodeDisp(t, 1), ops.nodeDisp(t, 2)
        yc, sl = chord(sx[t] + ux)
        pen = (uy - yc) / math.sqrt(1.0 + sl * sl)
        fq = ops.ladrunoContactForce(t)
        if pen > 0.0:
            act.append((sx[t], fq, sw[t]))
            if fq > 0.0:
                qerr = max(qerr, abs(fq - kn * pen) / (kn * pen))
                nq += 1
    Pm = sum(f for _, f, _ in act)
    out = dict(rc=rc, P=Pm, Rbase=Rbase, dx=dx, nq=nq, qerr=qerr)
    if Pm > 0.0:
        b, pmax = _hertz2d(Pm, R, E, NU)
        out.update(b=b, pmax=pmax, a=max(abs(x) for x, _, _ in act),
                   p0=max(f / w for _, f, w in act))
    return out


def test_hertz_line_contact():
    first = _run(b_hint=0.05 * R, nb=6)
    assert first["rc"] == 0 and first["P"] > 0.0, first
    r = _run(b_hint=first["b"], nb=8)
    assert r["rc"] == 0, r
    assert abs(r["P"] - r["Rbase"]) < 1.0e-3 * r["P"], r
    assert r["nq"] > 0 and r["qerr"] < 1.0e-9, r
    assert 2.0 * r["b"] / r["dx"] > 8.0, r
    assert r["a"] <= r["b"] + r["dx"] and r["a"] >= r["b"] - r["dx"], r
    assert abs(r["p0"] / r["pmax"] - 1.0) < 0.10, r
