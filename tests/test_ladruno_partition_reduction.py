"""recorder ladruno: PARTITION_REDUCTION of each result group.

Every result group says how a reader combines it across `.part-N` files: NONE
(consistent: kinematics, element results), SUM (additive partials: reactions,
unbalanced loads), UNSUPPORTED (not recoverable by a sum: energyBalance). In a
PARTITIONED run the recorder refuses `-envelope` of a SUM/UNSUPPORTED result,
because a per-partition extreme of a partial cannot be recombined.

Single process here. The partitioned leg runs in a SUBPROCESS with
PMI_SIZE=2/PMI_RANK=0, which is how the recorder detects an
interpreter-per-rank run (the launcher environment must be set before
`import opensees`).
"""
import os
import subprocess
import sys
import textwrap

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

h5py = pytest.importorskip("h5py")

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


MODEL = textwrap.dedent("""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 1.0e6)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.0, 1.0); ops.fix(2, 1, 0)
    ops.node(3, 1.0, 1.0); ops.fix(3, 1, 0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.element("Truss", 2, 1, 3, 1.0, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -10.0); ops.load(3, 0.0, -20.0)
""")
RUN = textwrap.dedent("""
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.5); ops.analysis("Static")
    ops.analyze(2)
    ops.wipe()
""")


def _scalar(v):
    """The recorder writes attributes as rank-1 arrays (same writer as COMPONENTS)."""
    if v is None:
        return None
    try:
        v = v.flat[0]
    except AttributeError:
        pass
    return v.decode() if isinstance(v, bytes) else v


def _attrs(path, family="ON_NODES", envelopes=False):
    """{result name: PARTITION_REDUCTION} for one family of one file."""
    out = {}
    with h5py.File(path, "r") as f:
        stage = [k for k in f if k.startswith("MODEL_STAGE")][0]
        base = f"{stage}/RESULTS/" + ("ENVELOPES/" if envelopes else "") + family
        if base not in f:
            return out
        for name, g in f[base].items():
            v = _scalar(g.attrs.get("PARTITION_REDUCTION"))
            out[name] = None if v is None else str(v)
    return out


def test_streaming_result_groups_carry_partition_reduction(tmp_path):
    out = str(tmp_path / "pr.ladruno")
    exec(MODEL, {"ops": ops})
    ops.recorder("ladruno", out, "-N", "displacement", "reactionForce", "unbalancedForce")
    exec(RUN, {"ops": ops})
    a = _attrs(out)
    assert a.get("DISPLACEMENT") == "NONE", a
    assert a.get("REACTION_FORCE") == "SUM", a
    assert a.get("UNBALANCED_FORCE") == "SUM", a


def test_energy_balance_is_unsupported(tmp_path):
    out = str(tmp_path / "pr_energy.ladruno")
    exec(MODEL, {"ops": ops})
    ops.recorder("ladruno", out, "-N", "displacement", "-G", "energy")
    exec(RUN, {"ops": ops})
    dom = _attrs(out, family="ON_DOMAIN")
    assert dom, "no ON_DOMAIN energy group was written"
    assert set(dom.values()) == {"UNSUPPORTED"}, dom


def test_serial_envelope_of_reaction_is_kept(tmp_path):
    """One partition only: a SUM envelope is exact, so it is recorded (and flagged)."""
    out = str(tmp_path / "pr_env.ladruno")
    exec(MODEL, {"ops": ops})
    ops.recorder("ladruno", out, "-N", "displacement", "reactionForce", "-envelope")
    exec(RUN, {"ops": ops})
    env = _attrs(out, envelopes=True)
    assert env.get("REACTION_FORCE") == "SUM", env
    assert env.get("DISPLACEMENT") == "NONE", env


def test_partitioned_envelope_of_reaction_is_refused(tmp_path):
    """Rank 0 of a 2-rank interpreter-per-rank run: the displacement envelope is kept,
    the reaction envelope is refused with a warning (a per-partition extreme of a partial
    cannot be recombined)."""
    out = str(tmp_path / "pr_part.ladruno")
    mod = ops.__name__
    here = os.path.dirname(os.path.abspath(ops.__file__))
    moddir = os.path.dirname(here) if "." in mod else here
    script = "\n".join([
        "import os, sys",
        f"sys.path.insert(0, {moddir!r})",
        f"os.add_dll_directory({here!r})" if hasattr(os, "add_dll_directory") else "",
        f"import {mod} as ops" if "." not in mod
        else f"import {mod.rsplit('.', 1)[0]}.opensees as ops",
        MODEL,
        f"ops.recorder('ladruno', {out!r}, '-N', 'displacement', 'reactionForce', '-envelope')",
        RUN,
    ])
    env = dict(os.environ, PMI_SIZE="2", PMI_RANK="0")
    r = subprocess.run([sys.executable, "-c", script], env=env, capture_output=True,
                       text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
    part = str(tmp_path / "pr_part.part-0.ladruno")
    assert os.path.exists(part), os.listdir(tmp_path)          # partitioned naming engaged
    with h5py.File(part, "r") as f:
        assert int(_scalar(f["INFO"].attrs["PARTITIONED"])) == 1
    env_attrs = _attrs(part, envelopes=True)
    assert "REACTION_FORCE" not in env_attrs, env_attrs          # refused
    assert env_attrs.get("DISPLACEMENT") == "NONE", env_attrs    # consistent: kept
    assert "is NOT recorded in this partitioned run" in (r.stdout + r.stderr)
