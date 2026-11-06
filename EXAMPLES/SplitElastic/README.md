# Paired elastic exterior element

`SplitElastic` uses one physical mesh for two half-property elastic branches.
Ordinary nodes (2 DOFs in 2-D, 3 in 3-D) share both branches. Exterior nodes
(4 or 6 DOFs) carry two independent translation sets. Specify full parent
`E nu rho`; the element supplies half stiffness and mass for each branch.
Fix the first branch at an outer 3-D node with `fix $n 1 1 1 0 0 0` (2-D:
`fix $n 1 1 0 0`). The second branch remains free. No automatic fixities or
interface constraints are added.

```
element SplitElastic tag n1 n2 n3 n4 E nu rho -coupling eta -couplingFactors f1 f2 f3 f4
# 3-D uses eight nodes and eight factors; dimensions follow the model NDM.
```

The relative connector contributed at a paired node is
`c_i = eta * f_i * integral(N_i dV)`, without a branch half factor. Element
contributions sum to the nodal tributary-volume coefficient. It is passive
and exactly zero for common branch motion. Ordinary Rayleigh/region damping
is independent of this connector. Elastic operators are cached once.

Default mass is **row-summed lumped** in both dimensions. `-cMass` selects
consistent mass. `Diagonal`/`MPIDiagonal` sum each element matrix's columns
onto the diagonal only with `-lumped`; without it they keep the diagonal of a
consistent mass, which loses mass. The generic solver parsers (OpenSeesPy)
accept `-lumped` upstream, but the traditional Tcl `system Diagonal|MPIDiagonal`
command ignored it until the fix on this branch (also on
`fix/tcl-diagonal-lumped`). The connector also restricts explicit time
steps: its isolated rigid differential mode requires
`dt < rho/(2*eta*f_i)`; the full mesh can require a smaller step.

Translational `UniformExcitation` in directions `1..ndm` applies the same
physical ground acceleration to both branches. Extra node DOFs are branch
translations, not rotations; generic six-DOF visualization/rotation features
must not be assumed compatible. Body/pressure element loads and nonlinear
materials are not implemented. A nonlinear version needs independent branch
material histories and independent validation of its absorption behavior.

Full documentation is maintained in the matching
`OpenSeesDocumentation` feature branch, `SplitElastic.rst`.

## Verification

Build `OpenSees`, `OpenSeesPy`, and optionally `OpenSeesMP`. This CMake build
names the Python library `OpenSeesPy.so`, but its import name is `opensees`.
A symlink in the ignored build directory makes it directly importable:

```sh
ln -s OpenSeesPy.so build/Release/opensees.so
PYTHONPATH=build/Release python3 -m unittest discover -s tests -p test_split_elastic.py -v
```

The tests require NumPy. They compare assembled stiffness, mass and damping
with independent stock quads/bricks and separate viscous connectors, including
distorted geometry, shared nodes, plane stress/strain, and nonzero Poisson's
ratio. They also check lumping, affine patches, branch stresses, connector
passivity/common-mode transparency, region overrides, database restoration,
malformed inputs, translational uniform excitation and the explicit damping
limit. Source-driven 2-D and 3-D histories are compared under Newmark and
ExplicitDifference. The explicit 3-D stock reference uses independently
assembled row-summed nodal mass, avoiding the consistent-diagonal trap.

## Matched MPI example

```sh
python3 EXAMPLES/SplitElastic/compare_split_elastic.py \
  --binary build/Release/OpenSeesMP --ranks 4 --steps 4000 --repeats 2 \
  --output /tmp/split-elastic-comparison
```

All decks, histories, final interior states, logs and `summary.json` go to the
specified external directory. The script compares the legacy two-solid-mesh
construction plus nodal dashpots with the combined element, using one binary,
matching variable elastic properties, geometric partitions, source, time step,
integration, damping and outputs. Interior quads remain the same in both.
Two repetitions run in legacy/combined/combined/legacy order. Ranks bind to
cores; MPI startup and model creation are excluded from the analysis timer.
Peak resident memory is the sum of each process's `/proc/self/status` VmHWM,
including shared library pages counted per process. It is not a measurement of
unique physical RAM or the simultaneous aggregate peak.

The example uses a small rectangle and simple column partitions, which are
not independently load-balanced for each implementation. It validates MPI
behavior and measures this case's throughput; it does not reproduce the full
paper geometry or establish a comparison with PML. Use `--ranks 1` and the
serial `OpenSees` binary for a serial comparison.

`column3d.tcl` supplies a small 3-D mixed-DOF MPI smoke example. From separate
external scratch directories, run the same file with serial `OpenSees` and
`mpirun --bind-to core -np 4 OpenSeesMP`. Each run emits `column3d-1.txt` through
`column3d-4.txt`; matching histories verify shared paired-node assembly between
ranks. The output files belong in scratch, not in this source repository.
