"""Self-tests of the .ladruno reader (ladruno_format.py, ladruno_basis.py).

Pure Python (h5py + numpy), no OpenSees run. Checked against the synthetic
conformant file of ladruno_synthetic.py:
  - the validator accepts a conformant file and rejects targeted corruptions;
  - the reader and the canonical normalizers return the values written;
  - the basis reconstruction maps quad Gauss points to the right global coords;
  - the round-trip oracle reproduces GLOBAL_GP_COORDS from GP_PARAM.

The same reader runs on real recorder output in test_ladruno_recorder_parity.py.
"""

from __future__ import annotations

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")  # the reader needs h5py
import ladruno_basis as lb
import ladruno_format as lf
import ladruno_synthetic


@pytest.fixture()
def synth(tmp_path):
    p = str(tmp_path / "synthetic.ladruno")
    ladruno_synthetic.build(p)
    return p


# ---- L2: validator ----------------------------------------------------------


def test_validator_accepts_conformant_file(synth):
    problems = lf.validate(synth)
    assert problems == [], f"synthetic file should be conformant, got: {problems}"


def test_validator_catches_bad_generator(synth, tmp_path):
    import h5py

    with h5py.File(synth, "a") as f:
        f["INFO"].attrs["GENERATOR"] = np.bytes_(b"STKO")
    problems = lf.validate(synth)
    assert any("GENERATOR" in p for p in problems)


def test_validator_catches_column_count_mismatch(synth):
    import h5py

    with h5py.File(synth, "a") as f:
        # break the NUM_COLUMNS == sum(MULT*NUM_COMP) invariant
        f["MODEL_STAGE[0]/RESULTS/ON_ELEMENTS/stress/156-FourNodeQuad[0:0]"].attrs[
            "NUM_COLUMNS"
        ] = 99
    problems = lf.validate(synth)
    assert any("NUM_COLUMNS" in p for p in problems)


def test_validator_catches_gp_param_shape(synth):
    import h5py

    with h5py.File(synth, "a") as f:
        g = f["MODEL_STAGE[0]/MODEL/ELEMENTS/156-FourNodeQuad[0]"]
        del g["QUADRATURE/GP_PARAM"]
        g["QUADRATURE"].create_dataset("GP_PARAM", data=np.zeros((3, 2)))  # wrong nGP
    problems = lf.validate(synth)
    assert any("GP_PARAM" in p for p in problems)


# ---- reader + normalizers ---------------------------------------------------


def test_reader_info_and_stage(synth):
    with lf.LadrunoReader(synth) as r:
        assert r.info()["GENERATOR"] == "Ladruno"
        assert r.stages() == ["MODEL_STAGE[0]"]
        assert r.stage_attrs("MODEL_STAGE[0]")["KIND"] == "static"


def test_normalize_nodal(synth):
    with lf.LadrunoReader(synth) as r:
        nd = lf.normalize_nodal(r, "MODEL_STAGE[0]")
    # node 1, Ux, step 0 -> first entry of (k+1)*1e-3*arange(1..10) = 1e-3*1
    assert nd[(1, "DISPLACEMENT:Ux", 0)] == pytest.approx(1e-3)
    assert nd[(5, "DISPLACEMENT:Uy", 1)] == pytest.approx(2e-3 * 10)


def test_normalize_element_quad_stress(synth):
    with lf.LadrunoReader(synth) as r:
        ed = lf.normalize_element(r, "MODEL_STAGE[0]")
    # quad 10, gp 2, no section/fiber, sigma_yy, step 0.
    # row data = arange(12); gp2 block starts at col 6 -> [6,7,8]=[xx,yy,xy]; yy=7
    assert ed[(10, 2, -1, -1, "stress:sigma_yy", 0)] == pytest.approx(7.0)


def test_normalize_element_fiber_multiplicity(synth):
    with lf.LadrunoReader(synth) as r:
        ed = lf.normalize_element(r, "MODEL_STAGE[0]")
    # beam 20, gp0 has 3 fibers (multiplicity 3) -> fiber ids 0,1,2 at cols 0,1,2
    # data = arange(6)*10 -> fiber1 of gp0 = col1 = 10
    assert ed[(20, 0, 1, 1, "section.fiber.stress:sigma", 0)] == pytest.approx(10.0)
    # gp1 fiber2 -> col 5 = 50
    assert ed[(20, 1, 1, 2, "section.fiber.stress:sigma", 0)] == pytest.approx(50.0)


# ---- L3: basis reconstruction ----------------------------------------------


def test_reconstruct_quad_gauss_points(synth):
    """Quad on the unit square: GP at +/-1/sqrt3 -> global 0.5 +/- 0.5/sqrt3."""
    with lf.LadrunoReader(synth) as r:
        g = r.element_groups("MODEL_STAGE[0]")["156-FourNodeQuad[0]"]
        nid, coords = r.nodes("MODEL_STAGE[0]")
    tag2row = {int(t): i for i, t in enumerate(nid)}
    conn = g["connectivity"][0, 1:]  # node tags of element 10
    Xc = coords[[tag2row[int(t)] for t in conn]]

    half = 0.5
    d = 0.5 / np.sqrt(3.0)
    expected = np.array(
        [[half - d, half - d], [half + d, half - d],
         [half + d, half + d], [half - d, half + d]]
    )
    for j in range(g["num_gp"]):
        xg = lb.reconstruct_global(
            Xc, g["family"], g["topology"], g["order"], g["param_domain"], g["gp_param"][j]
        )
        assert np.allclose(xg, expected[j], atol=1e-12), f"gp {j}: {xg} != {expected[j]}"


def test_reconstruct_bernstein_line_endpoints():
    """Sanity: a linear Bernstein line still interpolates its endpoints."""
    Xc = np.array([[0.0, 0.0], [2.0, 0.0]])
    x0 = lb.reconstruct_global(Xc, "bernstein", "line", (1,), "[0,1]", np.array([0.0]))
    x1 = lb.reconstruct_global(Xc, "bernstein", "line", (1,), "[0,1]", np.array([1.0]))
    assert np.allclose(x0, [0.0, 0.0])
    assert np.allclose(x1, [2.0, 0.0])


# ---- shape_functions (the GLOBAL_GP_COORDS oracle basis) -------------------


@pytest.mark.parametrize(
    "topo,nn,xi",
    [
        ("line", 2, [0.3]),
        ("tri", 3, [0.2, 0.5]),
        ("tri", 6, [0.2, 0.5]),
        ("quad", 4, [0.3, -0.4]),
        ("quad", 9, [0.3, -0.4]),
        ("tet", 4, [0.2, 0.3, 0.1]),
        ("tet", 10, [0.2, 0.3, 0.1]),
        ("hex", 8, [0.1, -0.2, 0.3]),
        ("hex", 20, [0.1, -0.2, 0.3]),
    ],
)
def test_shape_functions_partition_of_unity(topo, nn, xi):
    """Every supported basis must form a partition of unity (sum N_i == 1) -- an
    independent correctness check on the oracle basis, basis-by-basis."""
    N = lb.shape_functions(topo, nn, np.array(xi))
    assert N.shape == (nn,)
    assert np.isclose(N.sum(), 1.0, atol=1e-12), f"{topo}{nn}: sum={N.sum()}"


# ---- round-trip oracle -------------------------------------------------------


def test_round_trip_oracle_passes(synth):
    """The synthetic fixture carries GLOBAL_GP_COORDS for the quad, beam, and tri
    (bary) groups; the oracle must reconstruct them all to <=1e-12."""
    problems = lf.round_trip_oracle(synth)
    assert problems == [], f"oracle should be clean, got: {problems}"


def test_round_trip_oracle_catches_corruption(synth):
    """Perturb one GLOBAL_GP_COORDS entry -> the oracle must flag a mismatch."""
    import h5py

    with h5py.File(synth, "a") as f:
        g = f["MODEL_STAGE[0]/MODEL/ELEMENTS/156-FourNodeQuad[0]/GLOBAL_GP_COORDS"]
        g[0, 0] += 0.5  # move a Gauss point off its reconstructed location
    problems = lf.round_trip_oracle(synth)
    assert any("round-trip" in p and "FourNodeQuad" in p for p in problems), problems


def test_validator_partition_manifest(synth):
    """A ".part-N.ladruno" file declares INFO/PARTITIONED + PARTITION_ID +
    NUM_PARTITIONS; the validator accepts a consistent manifest and rejects a
    PARTITION_ID outside [0, NUM_PARTITIONS) (the contiguous-from-0 contract)."""
    import h5py

    with h5py.File(synth, "a") as f:
        f["INFO"].attrs["PARTITIONED"] = 1
        f["INFO"].attrs["PARTITION_ID"] = 2
        f["INFO"].attrs["NUM_PARTITIONS"] = 4
    assert lf.validate(synth) == []

    with h5py.File(synth, "a") as f:
        f["INFO"].attrs["PARTITION_ID"] = 4  # out of [0, 4)
    assert any("PARTITION_ID" in p for p in lf.validate(synth))


def test_validator_envelopes(synth):
    """ENVELOPES/<family>/<name>/{ID,MIN,MAX,ABSMAX,ARG_STEP}: validator
    accepts a consistent envelope (ABSMAX==max(|MIN|,|MAX|), MIN<=MAX) and the
    reader.envelopes() returns it; a broken ABSMAX is rejected."""
    import h5py

    mn = np.array([[-2.0, -1.0], [0.0, -3.0]])
    mx = np.array([[1.0, 4.0], [5.0, 2.0]])
    with h5py.File(synth, "a") as f:
        env = (f["MODEL_STAGE[0]/RESULTS"].create_group("ENVELOPES")
               .create_group("ON_NODES").create_group("DISPLACEMENT"))
        env.create_dataset("ID", data=np.array([[1], [2]], dtype=np.int64))
        env.create_dataset("MIN", data=mn)
        env.create_dataset("MAX", data=mx)
        env.create_dataset("ABSMAX", data=np.maximum(np.abs(mn), np.abs(mx)))
        env.create_dataset("ARG_STEP", data=np.zeros((2, 2), dtype=np.int32))
    assert lf.validate(synth) == []

    with lf.LadrunoReader(synth) as r:
        env = r.envelopes("MODEL_STAGE[0]")
        assert "ON_NODES" in env and "DISPLACEMENT" in env["ON_NODES"]

    with h5py.File(synth, "a") as f:
        f["MODEL_STAGE[0]/RESULTS/ENVELOPES/ON_NODES/DISPLACEMENT/ABSMAX"][0, 0] = 99.0
    assert any("ABSMAX" in p for p in lf.validate(synth))


def test_validator_partition_reduction(synth):
    """PARTITION_REDUCTION: optional (older files lack it); when present
    NONE|SUM|UNSUPPORTED; in a PARTITIONED file an envelope must be NONE, because a
    per-partition extreme of a partial result cannot be recombined."""
    import h5py

    with h5py.File(synth, "a") as f:
        on_nodes = f["MODEL_STAGE[0]/RESULTS/ON_NODES"]
        first = on_nodes[next(iter(on_nodes))]
        first.attrs["PARTITION_REDUCTION"] = "SUM"
    assert lf.validate(synth) == []                       # valid value, serial file

    with h5py.File(synth, "a") as f:
        on_nodes = f["MODEL_STAGE[0]/RESULTS/ON_NODES"]
        on_nodes[next(iter(on_nodes))].attrs["PARTITION_REDUCTION"] = "MAX"
    assert any("PARTITION_REDUCTION" in p for p in lf.validate(synth))
    with h5py.File(synth, "a") as f:
        on_nodes = f["MODEL_STAGE[0]/RESULTS/ON_NODES"]
        on_nodes[next(iter(on_nodes))].attrs["PARTITION_REDUCTION"] = "NONE"

    mn = np.array([[0.0, 5.0]])
    mx = np.array([[0.0, 10.0]])
    with h5py.File(synth, "a") as f:
        env = (f["MODEL_STAGE[0]/RESULTS"].create_group("ENVELOPES")
               .create_group("ON_NODES").create_group("REACTION_FORCE"))
        env.create_dataset("ID", data=np.array([[1]], dtype=np.int64))
        env.create_dataset("MIN", data=mn)
        env.create_dataset("MAX", data=mx)
        env.create_dataset("ABSMAX", data=np.maximum(np.abs(mn), np.abs(mx)))
        env.create_dataset("ARG_STEP", data=np.zeros((1, 2), dtype=np.int32))
        env.attrs["PARTITION_REDUCTION"] = "SUM"
    assert lf.validate(synth) == []                       # serial: a SUM envelope is fine

    with h5py.File(synth, "a") as f:
        f["INFO"].attrs["PARTITIONED"] = 1
        f["INFO"].attrs["PARTITION_ID"] = 0
        f["INFO"].attrs["NUM_PARTITIONS"] = 2
    assert any("unrecoverable" in p for p in lf.validate(synth))   # partitioned: refused


def test_validator_accepts_bary_simplex_ndir(synth):
    """The tri group uses PARAM_DOMAIN='bary' with NDIR=2 decoupled from ORDER=(1,)
    -- validate() must accept it."""
    with lf.LadrunoReader(synth) as r:
        g = r.element_groups("MODEL_STAGE[0]")["33-Tri31[0]"]
    assert g["param_domain"] == "bary"
    assert g["ndir"] == 2
    assert g["order"] == (1,)
    assert lf.validate(synth) == []


def _add_element_envelope(f, ncol_block=(4, 3), name="stress/204-Quad[201:0:0]"):
    """Build a nested element ENVELOPE leaf (ENVELOPES/ON_ELEMENTS/<display>/<bucket>)
    with a structured per-column COLUMN_MAP, mirroring writeElementColumnMap +
    EnvelopeSink output. Returns (mult, ncomp, num_cols)."""
    import h5py

    mult, ncomp = ncol_block
    num_cols = mult * ncomp
    nelem = 2
    rng = np.arange(nelem * num_cols, dtype=float).reshape(nelem, num_cols)
    mn = -rng - 1.0
    mx = rng + 2.0
    # leaf nests one level under ON_ELEMENTS -> create the <display> parent first
    on_e = f["MODEL_STAGE[0]/RESULTS"].create_group("ENVELOPES").create_group("ON_ELEMENTS")
    leaf = on_e.create_group(name)  # h5py creates the intermediate "<display>" group
    leaf.attrs["COMPONENTS"] = ""   # element result-level COMPONENTS is empty
    leaf.attrs["NUM_COLUMNS"] = num_cols
    leaf.create_dataset("ID", data=np.array([[10], [11]], dtype=np.int64))
    leaf.create_dataset("MIN", data=mn)
    leaf.create_dataset("MAX", data=mx)
    leaf.create_dataset("ABSMAX", data=np.maximum(np.abs(mn), np.abs(mx)))
    leaf.create_dataset("ARG_STEP", data=np.zeros((nelem, num_cols), dtype=np.int32))
    cm = leaf.create_group("COLUMN_MAP")
    # one block: GAUSS_ID increments across MULTIPLICITY, each with ncomp components
    cm.create_dataset("LEVELS", data=np.array([[lf.LEVEL_GAUSS]], dtype=np.int32))
    cm.create_dataset("GAUSS_ID", data=np.array([[1]], dtype=np.int32))
    cm.create_dataset("SECTION_TAG", data=np.array([[-1]], dtype=np.int32))
    cm.create_dataset("FIBER_ID", data=np.array([[-1]], dtype=np.int32))
    cm.create_dataset("NUM_COMP", data=np.array([[ncomp]], dtype=np.int32))
    cm.create_dataset("MULTIPLICITY", data=np.array([[mult]], dtype=np.int32))
    cm.attrs["COMP_NAMES"] = ",".join(f"s{j}" for j in range(ncomp))
    return mult, ncomp, num_cols


def test_reader_nested_element_envelope_with_column_map(synth):
    """Element envelopes nest under ON_ELEMENTS/<display>/<bucket> and carry the
    structured per-column COLUMN_MAP (item (1)). The reader must find the leaf and
    decode the gauss/section/fiber columns; validate() must accept it."""
    import h5py

    with h5py.File(synth, "a") as f:
        mult, ncomp, num_cols = _add_element_envelope(f)

    assert lf.validate(synth) == []

    with lf.LadrunoReader(synth) as r:
        env = r.envelopes("MODEL_STAGE[0]")
        assert "ON_ELEMENTS" in env
        leaf = env["ON_ELEMENTS"]["stress/204-Quad[201:0:0]"]
        assert leaf["min"].shape == (2, num_cols)
        cols = leaf["columns"]
        assert len(cols) == num_cols
        # GAUSS_ID increments across the MULTIPLICITY block (1..mult); each GP has
        # ncomp components named s0.. and no section/fiber level (-1).
        assert [c["gauss_id"] for c in cols[:ncomp]] == [1] * ncomp
        assert [c["gauss_id"] for c in cols[ncomp:2 * ncomp]] == [2] * ncomp
        assert {c["section_tag"] for c in cols} == {-1}
        assert {c["fiber_id"] for c in cols} == {-1}
        assert cols[0]["comp"] == "s0"


def test_validator_rejects_bad_element_envelope_column_map(synth):
    """A COLUMN_MAP whose sum(MULT*NUM_COMP) does not span the envelope columns is
    rejected (the element ENVELOPES COLUMN_MAP is validated like the time-series
    ON_ELEMENTS path)."""
    import h5py

    with h5py.File(synth, "a") as f:
        _add_element_envelope(f)
        # corrupt: claim 5 comps/block while the data is 3 -> column total mismatch
        cm = f["MODEL_STAGE[0]/RESULTS/ENVELOPES/ON_ELEMENTS/stress/204-Quad[201:0:0]/COLUMN_MAP"]
        cm["NUM_COMP"][...] = 5
    assert any("column total" in p or "NUM_COMP" in p for p in lf.validate(synth))
