"""responseSpectrumPerMode: per-mode response-spectrum resultants in one call.

For every requested direction d and mode k the command imposes the modal
displacement of `responseSpectrumAnalysis d ... -mode k`,
u_k = phi_k s_k Gamma_k Sa(T_k) / lambda_k, on the elements a requested
quantity needs, and returns, as the increment over the committed state and
flat [dir][mode][component]:

    reactions        sum of the reactions of every supported node
    reaction:<name>  the same over a named set of supported nodes
    ele:<name>       sum of eleForce of an element group at its cut nodes
    inertia:<name>   sum of the modal inertia forces M phi s Gamma Sa

Reference: the existing per-mode loop, `responseSpectrumAnalysis d -mode k`
(which commits u_k from a zero state) followed by `reactions`, nodeReaction,
eleForce and lambda_k * nodeMass * nodeDisp, on a two-storey 2D elastic
frame. Each component must agree to 1e-10 of the sum of the absolute values
of the terms that make it up (the round-off floor of a sum that cancels).
-about adds r x F about a point to the moments. The domain is left at its
committed state. The tests fail on a build without responseSpectrumPerMode
(unknown command).
"""
import math

import numpy as np
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


TN = [0.02, 0.1, 0.3, 0.6, 1.0, 3.0]
SA = [0.40, 1.00, 1.20, 0.90, 0.60, 0.20]
SPEC = ['-Tn', *TN, '-Sa', *SA]
NMODES = 3   # 12 equations: ARPACK needs nev + 8 below that
ABOUT = (3.0, -1.0)
RTOL = 1e-10

SUPPORTS = [1, 2]
STOREY_NODES = [3, 4, 5, 6]
FIRST_COLUMNS = [1, 2]


def _frame():
    """two-storey one-bay elastic frame; rotational masses keep M of full
    rank, which the ARPACK subspace (nev + 8 = 11 of 12 equations) needs"""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    crd = {1: (0.0, 0.0), 2: (6.0, 0.0), 3: (0.0, 3.0), 4: (6.0, 3.0),
           5: (0.0, 6.0), 6: (6.0, 6.0)}
    for n, (x, y) in crd.items():
        ops.node(n, x, y)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.geomTransf('Linear', 1)
    E = 2.5e7
    col = (0.16, 2.1e-3)
    beam = (0.12, 1.6e-3)
    conn = [(1, 1, 3, col), (2, 2, 4, col), (3, 3, 5, col), (4, 4, 6, col),
            (5, 3, 4, beam), (6, 5, 6, beam)]
    for e, i, j, (A, Iz) in conn:
        ops.element('elasticBeamColumn', e, i, j, A, E, Iz, 1)
    for n, m in zip(STOREY_NODES, (12.0, 12.0, 9.0, 9.0)):
        ops.mass(n, m, 0.8 * m, 0.05 * m)
    lam = ops.eigen(NMODES)
    ops.modalProperties()
    return np.asarray(lam)


def _moment(f, x, p):
    """add r x F about p to the in-plane moment of a 2D force block"""
    f = np.array(f, dtype=float)
    if p is not None:
        f[2] += (x[0] - p[0]) * f[1] - (x[1] - p[1]) * f[0]
    return f


def _reference(lam, dirs, modes):
    """the existing per-mode loop, values and term scales [d][k][c]"""
    keys = ['reactions', 'reaction:left', 'ele:storey1', 'inertia:all',
            'inertia:roof']
    out = {k: np.zeros((len(dirs), len(modes), 3)) for k in keys}
    scl = {k: np.zeros((len(dirs), len(modes), 3)) for k in keys}

    def put(key, d, a, terms):
        terms = np.asarray(terms).reshape(-1, 3)
        out[key][d, a] = terms.sum(axis=0)
        scl[key][d, a] = np.abs(terms).sum(axis=0)

    for di, d in enumerate(dirs):
        for a, k in enumerate(modes):
            ops.responseSpectrumAnalysis(d, *SPEC, '-mode', k)
            ops.reactions()
            # reactions: the element blocks at the supports
            blocks = []
            for e in FIRST_COLUMNS:
                f = np.array(ops.eleForce(e))
                n = ops.eleNodes(e)[0]
                blocks.append((n, f[0:3]))
            put('reactions', di, a,
                [_moment(ops.nodeReaction(n), ops.nodeCoord(n), ABOUT)
                 for n in SUPPORTS])
            scl['reactions'][di, a] = np.abs(
                [_moment(b, ops.nodeCoord(n), ABOUT) for n, b in blocks]).sum(axis=0)
            put('reaction:left', di, a, [ops.nodeReaction(1)])
            put('ele:storey1', di, a,
                [np.array(ops.eleForce(e))[0:3] for e in FIRST_COLUMNS])
            terms = []
            for n in STOREY_NODES:
                m = np.array(ops.nodeMass(n))
                u = np.array(ops.nodeDisp(n))
                terms.append(lam[k - 1] * m * u)
            put('inertia:all', di, a, terms)
            put('inertia:roof', di, a,
                [_moment(t, ops.nodeCoord(n), ABOUT)
                 for t, n in zip(terms[2:], STOREY_NODES[2:])])
    return out, scl


def _per_mode(dirs, modes=None):
    extra = ['-modes', *modes] if modes else []
    out = ops.responseSpectrumPerMode(
        *SPEC, '-dir', *dirs, *extra,
        '-reactions', '-about', *ABOUT,
        '-reactionSet', 'left', 1,
        '-eleGroup', 'storey1', '-ele', *FIRST_COLUMNS, '-cut', *SUPPORTS,
        '-inertia', 'all', *STOREY_NODES,
        '-inertia', 'roof', 5, 6, '-about', *ABOUT)
    nd, nm, nc = (int(v) for v in out['shape'])
    return {k: (np.array(v).reshape(nd, nm, nc) if ':' in k or k == 'reactions'
                else np.array(v)) for k, v in out.items()}


@pytest.mark.parametrize("modes", [None, [3, 1]])
def test_per_mode_equals_existing_loop(modes):
    lam = _frame()
    dirs = [1, 2]
    got = _per_mode(dirs, modes)
    use = modes if modes else list(range(1, NMODES + 1))
    assert got['shape'].tolist() == [2, len(use), 3]
    assert got['modes'].tolist() == use
    assert got['numSupports'].tolist() == [2]
    # the domain is left at its committed (zero) state
    for n in STOREY_NODES:
        assert np.allclose(ops.nodeDisp(n), 0.0, atol=0.0)

    ref, scl = _reference(lam, dirs, use)
    for key in ref:
        diff = np.abs(got[key] - ref[key])
        tol = RTOL * np.maximum(scl[key], 1e-300) + 1e-14
        assert np.all(diff <= tol), (key, float(np.max(diff)))
        assert np.max(np.abs(ref[key])) > 0.0

    # per-mode scalars: period, Sa and Gamma of modalProperties
    w = np.sqrt(lam[np.array(use) - 1])
    assert np.allclose(got['omega'], w, rtol=1e-12)
    assert np.allclose(got['period'], 2 * math.pi / w, rtol=1e-12)
    assert np.allclose(got['Sa'][:len(use)], np.interp(2 * math.pi / w, TN, SA),
                       rtol=1e-12)


def test_base_shear_equals_inertia():
    """equilibrium per mode: K u_k = lambda_k M u_k, so the horizontal base
    reaction equals minus the sum of the modal inertia forces"""
    lam = _frame()
    got = _per_mode([1])
    shear = got['reactions'][0, :, 0]
    inertia = got['inertia:all'][0, :, 0]
    assert np.allclose(shear, -inertia, rtol=1e-8, atol=1e-10 * np.max(np.abs(inertia)))


def test_refusals():
    _frame()
    with pytest.raises(Exception):      # nothing requested
        ops.responseSpectrumPerMode(*SPEC, '-dir', 1)
    with pytest.raises(Exception):      # direction out of range
        ops.responseSpectrumPerMode(*SPEC, '-dir', 4, '-reactions')
    with pytest.raises(Exception):      # mode out of range
        ops.responseSpectrumPerMode(*SPEC, '-dir', 1, '-modes', 9, '-reactions')
    with pytest.raises(Exception):      # not a supported node
        ops.responseSpectrumPerMode(*SPEC, '-dir', 1, '-reactionSet', 'x', 3)
    with pytest.raises(Exception):      # unknown option
        ops.responseSpectrumPerMode(*SPEC, '-dir', 1, '-reactions', '-bogus')
    with pytest.raises(Exception):      # no spectrum
        ops.responseSpectrumPerMode('-dir', 1, '-reactions')
