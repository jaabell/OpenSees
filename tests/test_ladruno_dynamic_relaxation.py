"""LadrunoDynamicRelaxation quasi-static integrator.

The integrator relaxes a model to static equilibrium with an explicit
leap-frog march on a diagonal fictitious mass M* (Gershgorin row sums of the
element tangents scaled by dt^2 / (4 f^2)) and kinetic damping. It runs under
`analysis Transient` with `algorithm Linear` and `system Diagonal`.

Checks on the von Mises truss (two corotTruss bars, E A = 1e4, half span 1.0,
rise 0.1, no mass, snap-through limit load about 3.80) under a constant
apex load P:
  - P = 3.0: the rest state equals the LoadControl equilibrium (1e-3), with
    kinetic and with viscous damping;
  - P = 1.15 x 3.80: the truss snaps through to the far stable branch
    (uy about -0.217), with -recompute 50 and with the default rebuild of M*
    at kinetic-energy peaks; LoadControl cannot pass the limit point;
  - the documented argument forms are accepted.
Stability margin, `ladrunoDR stabilityMargin` = (omega_max dt / 2)^2:
  - stdBrick block (2 x 2 x 2, ElasticIsotropic and J2Plasticity) pushed to
    a converged static state by Newton, then held by DR with no further
    load: residual, witness displacement and reaction do not drift over the
    hold (1e-9 relative), the reaction equals the Newton one to 1e-9
    (elastic) or 1e-8 (J2Plasticity, whose re-evaluation of the committed
    state offsets it by about 5e-9 in a static re-solve too), and the
    margin is f^2 = 0.25 for the elastic block, in [0.25, 1] for the
    plastic one;
  - the margin equals f^2 for -massSafety 1, 0.75, 0.5 and 0.25; the
    parser clamps and defaults -massSafety and -marginEvery;
  - multiplying E of a truss by r during the march moves the margin to
    f^2 r, with and without -noAutoRefresh, and it exceeds 1 exactly for
    r > 1/f^2 = 4; -marginEvery 0 reports -2 and a unity mass -1.
Failure handling:
  - a step that fails after newStep() and one update() (Newton calls update()
    twice, which the integrator refuses) and is then retried gives the same
    state, bit for bit, as an uninterrupted march of a Steel01 truss chain
    driven by sp, for kinetic damping, kinetic damping with -recompute 50 on
    the failed step, viscous damping, and two consecutive failures;
  - a NaN acceleration stops the step with a non-finite message and the
    time does not advance.
Classic Tcl: `ladrunoDR` answers residualNorm, kineticEnergy and
stabilityMargin and the relaxed state equals LoadControl (run only when an
OpenSees Tcl executable is found).
"""
import math
import os
import subprocess
from pathlib import Path

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


DR = "LadrunoDynamicRelaxation"
H = 0.10
APEX = 2
LIMIT = 3.80
DEFAULT_SAFETY = 0.5


# ---------------------------------------------------------------- truss
def _build_geom():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, -1.0, 0.0)
    ops.node(APEX, 0.0, H)
    ops.node(3, 1.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(APEX, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)
    ops.element("corotTruss", 1, 1, APEX, 1.0, 1)
    ops.element("corotTruss", 2, APEX, 3, 1.0, 1)


def _loadcontrol_uy(target, dlam=0.05):
    _build_geom()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(APEX, 0.0, -1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", dlam)
    ops.analysis("Static")
    while ops.getTime() < target - 1e-9:
        if ops.analyze(1) != 0:
            break
    return ops.getTime(), ops.nodeDisp(APEX, 2)


def _relax(P, integrator_args, chunk=50, dt=1.0, max_steps=400000,
           tol_vel=1e-6, tol_du=1e-8):
    _build_geom()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(APEX, 0.0, -P)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e-6, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    prev = ops.nodeDisp(APEX, 2)
    done = 0
    while done < max_steps:
        if ops.analyze(chunk, dt) != 0:
            return False, ops.nodeDisp(APEX, 2), done
        done += chunk
        vmax = abs(ops.nodeVel(APEX, 2))
        u = ops.nodeDisp(APEX, 2)
        du = abs(u - prev)
        prev = u
        if vmax < tol_vel and du < tol_du:
            return True, u, done
    return False, prev, done


def test_subcritical_matches_loadcontrol():
    _lam, ref = _loadcontrol_uy(3.0)
    conv, uy, steps = _relax(3.0, (DR,))
    assert conv, (steps, uy)
    assert abs(uy - ref) < 1e-3, (uy, ref)


def test_viscous_subcritical_matches_loadcontrol():
    _lam, ref = _loadcontrol_uy(3.0)
    conv, uy, steps = _relax(3.0, (DR, "-damping", "viscous"))
    assert conv, (steps, uy)
    assert abs(uy - ref) < 1e-3, (uy, ref)


def test_loadcontrol_cannot_pass_the_limit():
    lam, _uy = _loadcontrol_uy(1.15 * LIMIT)
    assert lam < LIMIT


@pytest.mark.parametrize("args", [(DR, "-recompute", 50), (DR,)],
                         ids=["recompute", "autoRefresh"])
def test_snap_through_to_far_branch(args):
    conv, uy, steps = _relax(1.15 * LIMIT, args, max_steps=800000)
    assert conv, (steps, uy)
    assert uy < -2.0 * H, uy
    assert abs(uy - (-0.217)) < 0.03, uy


@pytest.mark.parametrize("args", [
    (DR,),
    (DR, "-mass", "gershgorin"),
    (DR, "-mass", "unity"),
    (DR, "-mass", "lumped", 1.0),
    (DR, "-mass", "lumped", "-noAutoRefresh"),
    (DR, "-dt", 1.0, "-recompute", 100),
    (DR, "-damping", "viscous"),
    (DR, "-damping", "viscous", 0.7),
    (DR, "-damping", "viscous", "-massSafety", 0.5),
    (DR, "-noAutoRefresh"),
    (DR, "-massSafety", 0.25),
])
def test_parser_forms(args):
    _build_geom()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(APEX, 0.0, -1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e-6, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(*args)
    ops.analysis("Transient")
    assert ops.analyze(5, 1.0) == 0, args
    assert math.isfinite(ops.nodeDisp(APEX, 2))


def test_query_guards():
    _build_geom()
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(APEX, 0.0, -1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e-6, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(DR)
    ops.analysis("Transient")
    assert ops.analyze(1, 1.0) == 0
    assert ops.ladrunoDR("residualNorm") > 0.0
    assert ops.ladrunoDR("kineticEnergy") >= 0.0
    with pytest.raises(Exception):
        ops.ladrunoDR("notASubcommand")
    with pytest.raises(Exception):
        ops.ladrunoDR()
    ops.integrator("Newmark", 0.5, 0.25)
    with pytest.raises(Exception):
        ops.ladrunoDR("residualNorm")


# ------------------------------------------------- stability margin, block
E_BLOCK, NU = 2.0e5, 0.3
K_BULK = E_BLOCK / (3.0 * (1 - 2 * NU))
G_SHEAR = E_BLOCK / (2.0 * (1 + NU))
SIG0 = 250.0
U_ELASTIC = -2.0e-4
U_PLASTIC = -5.0e-3


def _build_block(plastic, n=2):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    if plastic:
        ops.nDMaterial("J2Plasticity", 1, K_BULK, G_SHEAR, SIG0, SIG0, 0.0, 200.0)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E_BLOCK, NU)
    h = 1.0 / n
    nid, tag = {}, 1
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                nid[(i, j, k)] = tag
                ops.node(tag, i * h, j * h, k * h)
                tag += 1
    eid = 1
    for k in range(n):
        for j in range(n):
            for i in range(n):
                c = [nid[(i, j, k)], nid[(i + 1, j, k)],
                     nid[(i + 1, j + 1, k)], nid[(i, j + 1, k)],
                     nid[(i, j, k + 1)], nid[(i + 1, j, k + 1)],
                     nid[(i + 1, j + 1, k + 1)], nid[(i, j + 1, k + 1)]]
                ops.element("stdBrick", eid, *c, 1)
                eid += 1
    for (i, j, k), t in nid.items():
        if k == 0:
            ops.fix(t, 1, 1, 1)
    tops = [t for (i, j, k), t in nid.items() if k == n]
    witness = nid[(n, n, n // 2)]
    return tops, witness


def _push(tops, u_target, nsteps=10):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in tops:
        ops.sp(t, 3, u_target)
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-10, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
    ops.reactions()
    return sum(ops.nodeReaction(t, 3) for t in tops)


def _dr_after_push(args):
    ops.loadConst("-time", ops.getTime())
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e30, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(*args)
    ops.analysis("Transient")


def _hold(args, nsteps, chunk=500):
    _dr_after_push(args)
    r_first = r_last = None
    done = 0
    while done < nsteps:
        if ops.analyze(chunk, 1.0) != 0:
            return False, r_first, r_last, ops.ladrunoDR("stabilityMargin")
        done += chunk
        r_last = ops.ladrunoDR("residualNorm")
        if r_first is None:
            r_first = r_last
    return True, r_first, r_last, ops.ladrunoDR("stabilityMargin")


@pytest.mark.parametrize("plastic", [False, True], ids=["elastic", "plastic"])
def test_hold_on_equilibrium_stays_put(plastic):
    tops, witness = _build_block(plastic)
    u = U_PLASTIC if plastic else U_ELASTIC
    R0 = _push(tops, u)
    if plastic:
        assert abs(R0) < 0.9 * E_BLOCK * abs(u)   # the block has yielded
    w0 = [ops.nodeDisp(witness, d) for d in (1, 2, 3)]
    scale = abs(R0)

    def reaction():
        ops.reactions()
        return sum(ops.nodeReaction(t, 3) for t in tops)

    _dr_after_push((DR,))
    assert ops.analyze(500, 1.0) == 0
    R_first = reaction()
    assert ops.analyze(2500, 1.0) == 0
    R_last = reaction()
    r_last = ops.ladrunoDR("residualNorm")
    margin = ops.ladrunoDR("stabilityMargin")

    assert r_last < 1e-9 * scale, (r_last, scale)
    w1 = [ops.nodeDisp(witness, d) for d in (1, 2, 3)]
    assert max(abs(a - b) for a, b in zip(w0, w1)) < 1e-9 * abs(u)
    # The march must not move the reaction over the hold.
    assert abs(R_last - R_first) < 1e-9 * scale, (R_first, R_last)
    # Against the Newton state: elastic exact to 1e-9. J2Plasticity
    # re-evaluates its committed state with an offset of about 5e-9 of the
    # reaction at the first evaluation after the commit; a static Newton
    # re-solve with no load change shows the same offset (measured
    # 1.30e-6 static vs 1.50e-6 DR on |R0| = 304.4), so it is not produced
    # by the march.
    tol0 = 1e-8 if plastic else 1e-9
    assert abs(R_first - R0) < tol0 * scale, (R0, R_first)
    assert DEFAULT_SAFETY ** 2 - 1e-9 <= margin <= 1.0, margin
    if not plastic:
        assert margin == pytest.approx(DEFAULT_SAFETY ** 2, rel=1e-6)


@pytest.mark.parametrize("f", [1.0, 0.75, 0.5, 0.25])
def test_margin_equals_safety_squared(f):
    tops, _w = _build_block(False)
    _push(tops, U_ELASTIC)
    ok, _r0, _r1, margin = _hold((DR, "-massSafety", f), 100, chunk=50)
    assert ok
    assert margin == pytest.approx(f * f, rel=1e-9)


@pytest.mark.parametrize("args,expect", [
    ((DR,), DEFAULT_SAFETY ** 2),
    ((DR, "-massSafety", 0.25), 0.0625),
    ((DR, "-massSafety", 1.0), 1.0),
    ((DR, "-massSafety", 4.0), 1.0),
    ((DR, "-massSafety", 0.0), DEFAULT_SAFETY ** 2),
    ((DR, "-massSafety", -1.0), DEFAULT_SAFETY ** 2),
    ((DR, "-mass", "gershgorin", "-massSafety", 0.5), 0.25),
    ((DR, "-massSafety", 0.5, "-damping", "viscous"), 0.25),
    ((DR, "-marginEvery", 100), DEFAULT_SAFETY ** 2),
    ((DR, "-marginEvery", 0), -2.0),
    ((DR, "-marginEvery", -5), -2.0),
    ((DR, "-mass", "unity", "-massSafety", 0.25), -1.0),
])
def test_margin_parser(args, expect):
    tops, _w = _build_block(False)
    _push(tops, U_ELASTIC)
    _dr_after_push(args)
    assert ops.analyze(20, 1.0) == 0, args
    assert ops.ladrunoDR("stabilityMargin") == pytest.approx(expect, rel=1e-9)


def _stiffen_and_read(factor, extra_args, steps=2000, chunk=50):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 10.0, 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1e30, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(DR, "-marginEvery", 50, *extra_args)
    ops.analysis("Transient")
    ops.analyze(50, 1.0)
    before = ops.ladrunoDR("stabilityMargin")
    ops.parameter(1, "element", 1, "E")
    ops.updateParameter(1, 1000.0 * factor)
    done = 0
    while done < steps:
        if ops.analyze(chunk, 1.0) != 0:
            break
        done += chunk
    return before, ops.ladrunoDR("stabilityMargin")


@pytest.mark.parametrize("factor", [1.5, 2.0, 4.0, 6.0, 8.0])
@pytest.mark.parametrize("refresh", ["auto", "noAutoRefresh"])
def test_margin_tracks_a_known_stiffening(factor, refresh):
    extra = [] if refresh == "auto" else ["-noAutoRefresh"]
    before, after = _stiffen_and_read(factor, extra)
    assert before == pytest.approx(DEFAULT_SAFETY ** 2, rel=1e-9)
    assert after == pytest.approx(DEFAULT_SAFETY ** 2 * factor, rel=1e-6)
    assert (after > 1.0 + 1e-6) == (factor > 1.0 / DEFAULT_SAFETY ** 2)


def test_margin_independent_of_refresh_policy():
    _b, with_refresh = _stiffen_and_read(6.0, [])
    _b, without = _stiffen_and_read(6.0, ["-noAutoRefresh"])
    assert with_refresh == pytest.approx(without, rel=1e-12)
    assert without > 1.0


def test_probe_off_reports_not_measured():
    _b, after = _stiffen_and_read(6.0, ["-marginEvery", 0])
    assert after == -2.0


# ------------------------------------------------- failed and retried step
RATE = 1.0e-5
N_STEPS = 400
FAIL_AT = 150


def _build_chain():
    # Steel01 chain driven by sp at the free end: u = -t (the DR clock);
    # 4e-3 total over two bars of length 1 is above the yield strain 1.25e-3
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i, x in enumerate((0.0, 1.0, 2.0), 1):
        ops.node(i, x, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(3, 0, 1)
    ops.uniaxialMaterial("Steel01", 1, 250.0, 2.0e5, 0.01)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.element("truss", 2, 2, 3, 1.5, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(3, 1, -1.0)
    return [1, 2, 3]


def _dr_chain(options):
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1.0e-12, 1, 0)
    ops.algorithm("Linear")
    ops.integrator(DR, "-dt", RATE, *options)
    ops.analysis("Transient")


def _state(nodes):
    return [ops.getTime()] + [u for nd in nodes for u in ops.nodeDisp(nd)]


def _fail_once():
    # Newton calls update() a second time in the step, which DR refuses,
    # after newStep() and one update() have advanced the march
    ops.test("NormUnbalance", 1.0e-300, 2, 0)
    ops.algorithm("Newton")
    rc = ops.analyze(1, RATE)
    ops.test("NormUnbalance", 1.0e-12, 1, 0)
    ops.algorithm("Linear")
    return rc


OPTS = {
    "kinetic": ["-damping", "kinetic"],
    "kinetic_recompute": ["-damping", "kinetic", "-recompute", 50],
    "viscous": ["-damping", "viscous", 1.0],
}


def _uninterrupted(opts):
    nodes = _build_chain()
    _dr_chain(opts)
    assert ops.analyze(N_STEPS, RATE) == 0
    return _state(nodes)


@pytest.mark.parametrize("name", sorted(OPTS))
def test_failed_then_retried_equals_uninterrupted(name):
    ref = _uninterrupted(OPTS[name])
    nodes = _build_chain()
    _dr_chain(OPTS[name])
    assert ops.analyze(FAIL_AT - 1, RATE) == 0
    before = _state(nodes)
    assert _fail_once() < 0
    assert _state(nodes) == before
    assert ops.analyze(N_STEPS - (FAIL_AT - 1), RATE) == 0
    got = _state(nodes)
    assert got == ref, max(abs(a - b) for a, b in zip(got, ref))


def test_two_failures_then_retry():
    ref = _uninterrupted(OPTS["kinetic_recompute"])
    nodes = _build_chain()
    _dr_chain(OPTS["kinetic_recompute"])
    assert ops.analyze(FAIL_AT - 1, RATE) == 0
    assert _fail_once() < 0
    assert _fail_once() < 0
    assert ops.analyze(N_STEPS - (FAIL_AT - 1), RATE) == 0
    assert _state(nodes) == ref


def test_nan_acceleration_stops_the_step(capfd):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.setNodeDisp(2, 1, float("nan"), "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(DR)
    ops.analysis("Transient")
    capfd.readouterr()
    rc = ops.analyze(1, 1.0e-3)
    err = capfd.readouterr().err
    assert rc != 0
    assert "non-finite" in err, err
    assert ops.getTime() == pytest.approx(0.0, abs=1e-15)


# --------------------------------------------------------- classic Tcl
TCL_SCRIPT = """
proc build_arch {} {
    wipe
    model basic -ndm 2 -ndf 2
    node 1 -1.0 0.0
    node 2  0.0 0.1
    node 3  1.0 0.0
    fix 1 1 1
    fix 3 1 1
    fix 2 1 0
    uniaxialMaterial Elastic 1 1.0e4
    element corotTruss 1 1 2 1.0 1
    element corotTruss 2 2 3 1.0 1
}
build_arch
timeSeries Linear 1
pattern Plain 1 1 { load 2 0.0 -1.0 }
system BandGeneral
numberer Plain
constraints Plain
test NormDispIncr 1.0e-10 50 0
algorithm Newton
integrator LoadControl 0.05
analysis Static
set guard [catch {ladrunoDR residualNorm} msg]
while {[getTime] < 3.0 - 1.0e-9} {
    if {[analyze 1] != 0} { break }
}
set refUy [nodeDisp 2 2]

build_arch
timeSeries Constant 1
pattern Plain 1 1 { load 2 0.0 -3.0 }
constraints Plain
numberer Plain
system Diagonal
test NormUnbalance 1.0e30 1 0
algorithm Linear
integrator LadrunoDynamicRelaxation
analysis Transient
analyze 1 1.0
set r0 [ladrunoDR residualNorm]
analyze 2000 1.0
set r1 [ladrunoDR residualNorm]
set ke [ladrunoDR kineticEnergy]
set sm [ladrunoDR stabilityMargin]
set bad [catch {ladrunoDR notASubcommand} msg]
puts "RESULT $guard $r0 $r1 $ke $sm $bad [nodeDisp 2 2] $refUy"
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
    script = tmp_path / "dynamic_relaxation.tcl"
    script.write_text(TCL_SCRIPT)
    res = subprocess.run([str(_tcl_executable()), str(script)],
                         capture_output=True, text=True, cwd=tmp_path,
                         timeout=300, stdin=subprocess.DEVNULL,
                         env=dict(os.environ))
    out = res.stdout + res.stderr
    assert "invalid command name" not in out, out
    line = [x for x in out.splitlines() if x.startswith("RESULT")]
    assert line, out
    guard, r0, r1, ke, sm, bad, uy, ref = (float(v) for v in line[-1].split()[1:9])
    assert guard == 1.0
    assert r0 > 1.0
    assert r1 < 1.0e-9 * r0
    assert ke >= 0.0
    # the corotational tangent stiffens as the arch deflects, so the margin
    # moves above f^2 but stays inside the stability limit
    assert DEFAULT_SAFETY ** 2 - 1e-9 <= sm <= 1.0
    assert bad == 1.0
    assert abs(uy - ref) < 1.0e-6
