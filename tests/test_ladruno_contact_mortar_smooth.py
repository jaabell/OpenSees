"""Smoothed (C1) mortar contact law: contact ... -mortar -augment never -smoothN g0 <-smoothT r>.

The penalty law p = min(0, epsN gbar) has a kink at gbar = 0 and the cohesive friction cone
jumps from c to 0 at lift-off, which Newton has to carry through every active-set change.
-smoothN replaces the kink by a C1 quadratic onset over |gbar| < g0 and fades the friction
traction in over the same band (a C1 smoothstep); -smoothT rounds the stick/slip corner of
the return map over |rho - cap| < r cap. With x = -epsN gbar and S = epsN g0:
P = 0 (x <= -S), (x + S)^2 / (4 S) (|x| < S), x (x >= S).

Checks:

- A stiff stdBrick block on springs pressed to a pressure inside the band matches the
  smoothed law (1e-6), beyond the band the unsmoothed law, and pulled off past the band the
  contact carries nothing (master reaction exactly 0).
- A cohesive block under a lateral pull, pressed then pulled off past separation in 12 steps,
  converges in at most 8 iterations per step, carries exactly nothing once open, and the
  lateral pull rides the springs alone (1e-6).
- The creased roof with shared multi-pair ridge nodes and a small band still matches the
  analytic descent (2e-4) in at most 6 iterations per step.
- As g0 -> 0 the state converges monotonically to the unsmoothed one (below 1e-6 at
  g0 = 1e-6).
- Refusals, each for its named reason: without -mortar, with -tie, with an augmenting
  contact (commit or request), -smoothT without friction, out-of-range values.

Every test fails on the base with an unknown command.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

import test_ladruno_contact_mortar_controls as R158


KN = 1.0e6       # epsN
KD = 1.0e3       # spring per top node (4 top nodes)
A = 1.0          # block footprint (unit square)


def _block(flags=(), coh=0.0, mu=0.0, delta=0.0):
    """A stiff unit brick resting on a fixed, larger master quad (mortar, -augment never), exactly
    (delta = 0) or overlapping it by delta, each top node on (x, y, z) springs KD. Returns (top nodes, bottom nodes, master nodes)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = []
    for t, (x, y) in enumerate([(-0.5, -0.5), (1.5, -0.5), (1.5, 1.5), (-0.5, 1.5)], start=1):
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
        m.append(t)
    bot, top = [11, 12, 13, 14], [15, 16, 17, 18]
    for t, (x, y) in zip(bot, [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), -delta)
        ops.node(t + 4, float(x), float(y), 1.0 - delta)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e11, 0.2)
    ops.element("stdBrick", 1, *bot, *top, 1)
    ops.contactSurface(1, "-master", 4, *m)
    ops.contactSurface(2, "-slave-segments", 4, *bot)
    opts = ["-mortar", "-epsN", KN, "-epsT", KN, "-outward", 0.0, 0.0, 1.0, "-augment", "never"]
    if mu > 0.0:
        opts += ["-mu", mu]
    if coh > 0.0:
        opts += ["-cohesion", coh]
    ops.contact(1, 1, 2, *opts, *flags)
    ops.uniaxialMaterial("Elastic", 1, KD)
    for t in top:
        ops.node(100 + t, *ops.nodeCoord(t))
        ops.fix(100 + t, 1, 1, 1)
        ops.element("zeroLength", 200 + t, 100 + t, t, "-mat", 1, 1, 1, "-dir", 1, 2, 3)
    _solver()
    return top, bot, m


def _solver():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1.0e-7, 40, 0)
    ops.algorithm("Newton")


def _law(x, S):
    if S <= 0.0:
        return max(x, 0.0)
    if x <= -S:
        return 0.0
    if x >= S:
        return x
    return (x + S) ** 2 / (4.0 * S)


def _w_analytic(F, g0):
    """Downward displacement w of the rigid block: F = A*P(KN*w) + 4*KD*w (bisection)."""
    S = KN * g0
    lo, hi = -10.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if A * _law(KN * mid, S) + 4 * KD * mid - F > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _apply(top, fz, fx=0.0, steps=1):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, fx / 4.0, 0.0, -fz / 4.0)
    ops.integrator("LoadControl", 1.0 / steps)
    ops.analysis("Static")
    its = []
    for _ in range(steps):
        if ops.analyze(1) != 0:
            return its, False
        its.append(ops.testIter())
    return its, True


def _w(bot):
    return -sum(ops.nodeDisp(n, 3) for n in bot) / len(bot)


# ------------------------------------------------------------------------------------- (law)
@pytest.mark.parametrize("F", [0.5, 2.0, 50.0, -0.5])
def test_block_equilibrium_matches_smoothed_law(F):
    """In the band (F = 0.5, 2: x/S = -0.1..0.8), beyond it (F = 50: the unsmoothed ramp), and pulled
    off past it (F = -0.5: w = F/(4 KD), the contact carries nothing)."""
    g0 = 1.0e-6                                    # S = 1.0 (pressure units)
    top, bot, m = _block(("-smoothN", g0))
    _, ok = _apply(top, F)
    assert ok
    w = _w(bot)
    wa = _w_analytic(F, g0)
    assert w == pytest.approx(wa, rel=1e-6, abs=1e-14), f"w {w:.6e} != smoothed law {wa:.6e}"
    if F >= 50.0:
        assert wa == pytest.approx(_w_analytic(F, 0.0), rel=1e-12)   # outside the band: unsmoothed
    ops.reactions()
    rz = sum(ops.nodeReaction(t, 3) for t in m)
    if F < 0.0:
        assert KN * w < -KN * g0                   # open beyond the band
        assert rz == 0.0, f"tension through an open contact: master Rz = {rz}"


# --------------------------------------------------------------------------------------- (a)
def test_block_liftoff_newton_bounded():
    """Cohesive bond (c = 2) under a lateral pull, pressed then pulled off past separation, force
    control, 12 steps: every step converges in <= 8 iterations (measured 3, 3, 3, 3, 4, 4, 6, 2, 1,
    1, 1, 1), the master carries exactly nothing once the block is open beyond the band, and the
    lateral pull then rides the springs alone. A rigid block flips all four nodes at once, so the
    unsmoothed law converges here too: the check is that the smoothed law and its tangent get
    through separation."""
    g0, coh, H = 1.0e-4, 2.0, 0.5
    top, bot, m = _block(("-consistanttan", "-smoothN", g0, "-smoothT", 0.1), coh=coh)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, H / 4.0, 0.0, -1.0 / 4.0)      # pattern 1: lateral H + press 1 (x lambda)
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in top:
        ops.load(t, 0.0, 0.0, 1.0 / 4.0)           # pull up: net vertical goes -1 -> +1.4
    ops.integrator("LoadControl", 2.4 / 12)
    its = []
    for _ in range(12):
        assert ops.analyze(1) == 0, f"Newton failed through lift-off after {its}"
        its.append(ops.testIter())
    assert max(its) <= 8, f"iterations per step {its} (want <= 8)"
    w = _w(bot)
    assert KN * w < -KN * g0, f"block not open beyond the band (w = {w})"
    ops.reactions()
    rx = sum(ops.nodeReaction(t, 1) for t in m)
    rz = sum(ops.nodeReaction(t, 3) for t in m)
    assert rx == 0.0 and rz == 0.0, f"open contact carries Rx {rx}, Rz {rz}"
    ux = sum(ops.nodeDisp(t, 1) for t in top) / 4.0
    assert ux == pytest.approx(H / (4 * KD), rel=1e-6)


# --------------------------------------------------------------------------------------- (b)
def test_creased_roof_analytic_with_band():
    """The solid creased roof of tests/test_ladruno_contact_mortar_controls.py (mu = 0.2, shared ridge, full slip) with -consistanttan and a band far
    inside the penetration: the analytic descent still holds (2e-4) in <= 6 iterations per step."""
    its, ok, surf, springs = R158._roof158(flags=("-consistanttan", "-smoothN", 1.0e-8), **R158.SOLID)
    assert ok, f"Newton failed after {its}"
    assert max(its) <= 6, f"iterations per step {its}"
    kd = 1.0e3
    _, fz = R158._descent_analytic(len(surf), len(springs), R158.SOLID["pz"], kd, 1.0e6,
                                   R158.SOLID["mu"], R158.SOLID["delta"], 1.0)
    got = len(surf) * R158.SOLID["pz"] + sum(kd * ops.nodeDisp(n, 3) for n in springs)
    assert got == pytest.approx(fz, rel=2e-4)


# --------------------------------------------------------------------------------------- (c)
def test_band_to_zero_converges_to_the_unsmoothed_law():
    """Closed case (the block starts 3e-6 into the foundation -- engaged from the first evaluation
    under both laws -- and is pressed by F = 2 with a lateral pull and a cohesive bond; penetration
    ~5e-6): the smoothed state -> the unsmoothed state as g0 -> 0, monotone, and identical to the
    tolerance once the band no longer reaches the penetration (g0 = 1e-6)."""
    def state(flags):
        top, _, _ = _block(("-consistanttan",) + flags, coh=2.0, delta=3.0e-6)
        _, ok = _apply(top, 2.0, fx=0.5)
        assert ok, f"{flags}: Newton failed"
        return [v for t in top for v in ops.nodeDisp(t)]
    ref = state(())
    scale = max(abs(v) for v in ref)
    errs = []
    for g0 in (3.0e-5, 1.0e-5, 1.0e-6):
        got = state(("-smoothN", g0))
        errs.append(max(abs(a - b) for a, b in zip(got, ref)) / scale)
    assert errs[0] > errs[1] > errs[2], f"not monotone: {errs}"
    assert errs[0] > 1e-2 and errs[2] < 1e-6, f"errors {errs}"


# ------------------------------------------------------------------------------------- (ref)
def _facets():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], start=1):
        ops.node(t, float(x), float(y), 0.0)
        ops.node(t + 10, float(x), float(y), 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 11, 12, 13, 14)


def _contact_fails(opts):
    _facets()
    try:
        ops.contact(1, 1, 2, *opts)
    except Exception:
        return True
    return False


def _refused(capfd, why, opts):
    """Refused for the named reason (a build without the options also raises, with "unexpected
    token")."""
    _facets()
    capfd.readouterr()
    with pytest.raises(ops.OpenSeesError):
        ops.contact(1, 1, 2, *opts)
    err = capfd.readouterr()
    text = err.out + err.err
    assert why in text, text


@pytest.mark.parametrize("opts,why", [
    (("auto", "-smoothN", 1e-4), "-smoothN/-smoothT are -mortar options"),                  # NTS
    (("-mortar", "-epsN", 1e6, "-tie", "-augment", "never", "-smoothN", 1e-4),
     "-smoothN/-smoothT do not apply to -tie"),
    (("-mortar", "-epsN", 1e6, "-smoothN", 1e-4),                                           # commit
     "pure-penalty laws; add -augment never"),
    (("-mortar", "-epsN", 1e6, "-augment", "request", "-smoothN", 1e-4),
     "pure-penalty laws; add -augment never"),
    (("-mortar", "-epsN", 1e6, "-augment", "never", "-smoothT", 0.1),
     "this contact has no friction"),
    (("-mortar", "-epsN", 1e6, "-augment", "never", "-smoothN", 0.0), "-smoothN - need a band g0 > 0"),
    (("-mortar", "-epsN", 1e6, "-augment", "never", "-smoothN", -1e-4), "-smoothN - need a band g0 > 0"),
    (("-mortar", "-epsN", 1e6, "-mu", 0.3, "-augment", "never", "-smoothT", 1.0),
     "-smoothT - need a ratio 0 < r < 1"),
    (("-mortar", "-epsN", 1e6, "-mu", 0.3, "-augment", "never", "-smoothT", 0.0),
     "-smoothT - need a ratio 0 < r < 1"),
])
def test_refusals(capfd, opts, why):
    _refused(capfd, why, opts)


def test_accepts_the_recipe():
    assert not _contact_fails(("-mortar", "-epsN", 1e6, "-cohesion", 1.0, "-augment", "never",
                               "-smoothN", 1e-4, "-smoothT", 0.1))
