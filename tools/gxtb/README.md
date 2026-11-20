# Periodic g-xTB in CP2K

`FORCE_EVAL / METHOD GXTB` evaluates a neutral, three-dimensionally periodic
cell in the CP2K process. CP2K constructs the
Born-von Karman supercell from `REPEATS na nb nc`, passes its lattice and
atoms to a Gamma-point electronic solver, and returns the energy per input
cell, the mean force on equivalent images and the pressure-positive virial.
One MPI rank calls the provider and broadcasts the results. The CP2K changes
and adapters were written against public interfaces.

## Open-source CPU provider

The build uses Paderborn's public `save_tblite` snapshot
`35e7942b60edd89bb407ab3da5768d3410af83f5` for the periodic electronic
solver. Its original g-xTB dispersion dependency was
`grimme-lab/dftd@99d64ee83832fda03df2c7dd7b7fd9f4c9ebb098`. The source for
that revision could not be verified. The electronic basis, integrals,
Hamiltonian, exchange, diagonalization and SCF iteration come from this
electronic solver. [Electronic solver archive](https://github.com/DCM-Uni-Paderborn/Periodic-g-xTB-Benchmarks/tree/866691a93b36ee582e49d896dfda5b91524dce12/validation/accelerated_exchange/root_a_partial_bridge/final_35e7942_20260717/provider/source)

CP2K implements g-xTB 2.0.1 D4Srev in `src/gxtb_revd4.cpp` as part of
`libcp2k`. It computes coordination numbers, pair-specific C6 interpolation
on the imaginary frequency grid, screened damping, periodic pair and
Axilrod-Teller-Muto (ATM) sums, and analytic charge, coordinate and strain
derivatives. The electronic solver receives CP2K's dispersion functions
when it enters an energy/force evaluation. Every SCF cycle
uses the current Mulliken charges for the two-body energy and its charge
potential. The three-body term uses neutral reference C6 at the current
geometry. The dispersion uses the 2.0.1 program's hard cutoffs of 30, 50 and
25 bohr for coordination, two-body and three-body terms, respectively.
Energy derivatives at cutoff crossings can depend on the displacement step.
Parameter tables have a reproducible numerical import from the public
g-xTB and D4 data distributed with axtb 0.2.0. The earlier Rust
implementation is available as an optional numerical reference.
[Parameter provenance and regeneration](revd4/PARAMETERS.md)

Build CP2K from this branch with its usual CMake configuration. The plugin
requires Linux, CMake 3.22 or newer, C, C++17 and Fortran compilers with
OpenMP, BLAS/LAPACK, curl, GNU tar, patch, flock and
coreutils. Run from `tools/gxtb` with your CP2K executable on PATH.

```sh
bash revd4/build.sh "$PWD/build-gxtb"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cp2k.ssmp -i revd4/nacl.inp
```

The example selects `METHOD GXTB`, `LIBRARY ./build-gxtb/install/lib/libgxtb_save.so`
and `REPEATS 1 1 1`. Its expected distorted NaCl energy is
-622.7008860203 Ha. The provider sets the SCF accuracy to 0.01, uses zero
electronic temperature and starts from EEQBC charges. It takes its exchange
settings from the electronic solver. For runs elsewhere, set `LIBRARY` to
the absolute path. To make this the linked electronic solver, configure
CP2K with `-DCP2K_GXTB_PROVIDER_LIBRARY=/absolute/path/to/libgxtb_save.so`
and rebuild. This enables `METHOD GXTB` without a `LIBRARY` keyword. The
provider library and its compiler and BLAS runtimes must be accessible at
run time. The separately built CLI `build-gxtb/install/bin/tblite` compiles
the same CP2K dispersion source directly. Neither build requires Rust.
The build script verifies its source downloads by SHA-256.

The dispersion currently searches basis atom pairs in quadratic time and
stores interactions within a fixed real-space cutoff. Once the geometry
is prepared, each SCF charge update visits the unique interacting atom
pairs. For `M` orbitals per input cell
and `m` replicated cells the dense electronic eigensolve costs O((mM)^3)
time and O((mM)^2) memory per SCF iteration. These are algorithmic estimates.

### Qualification

The scripts below take a NumPy/ASE Python environment and a new directory
for each CP2K run. The direct comparison uses a separate process running
the same provider. Finite differences use independent CP2K energy runs.

```sh
prefix="$PWD/build-gxtb/install"
python3 revd4/check.py /path/to/libcp2k.so
python3 check_periodic.py "$(command -v cp2k.ssmp)" checks-open-nacl --engine save-tblite --library "$prefix/lib/libgxtb_save.so" --reference-cli "$prefix/bin/tblite" --all-derivatives --force-step 0.0001 --strain-step 0.00002
python3 check_periodic.py "$(command -v cp2k.ssmp)" checks-open-si --engine save-tblite --library "$prefix/lib/libgxtb_save.so" --reference-cli "$prefix/bin/tblite" --all-derivatives --structure silicon
```

Use `--linked-provider` in place of `--library PATH` to exercise the linked
provider. To compare both implementations term by term, build the Rust
oracle separately (tested with Rust 1.91.1).

```sh
cargo build --release --locked --manifest-path revd4/Cargo.toml --target-dir "$PWD/build-gxtb/oracle"
python3 revd4/check.py /path/to/libcp2k.so --reference "$PWD/build-gxtb/oracle/release/libcp2k_gxtb_revd4.so"
```

| Structure | Displacement (angstrom) | Strain step | Maximum force error (eV/angstrom) | Maximum stress error (eV/angstrom^3) |
| --- | ---: | ---: | ---: | ---: |
| Distorted NaCl, 2 atoms | 1e-4 | 2e-5 | 3.14e-5 | 9.10e-8 |
| Distorted diamond Si, 8 atoms | 2.5e-4 | 1e-5 | 1.36e-5 | 4.41e-8 |

The table covers all 3N forces and six strains. For both cells CP2K and
the standalone program agree within 1.2e-7 eV in energy, 1.6e-8 eV/angstrom
in force and 3.4e-9 eV/angstrom^3 in stress. `REPEATS 2 1 1` and a
separately supplied doubled structure agree in energy per cell to
1.5e-11 eV. The NaCl force difference at the coarser 1e-3 angstrom
step is 0.030234 eV/angstrom, exceeding the 2e-4 threshold. A separate run
with the dispersion energy, potential and derivatives all set to zero
still gives 0.030134 eV/angstrom at that step. Its cause remains open.

At the archived molecular water density the two-body and ATM terms give
-4.42494068950e-5 and +3.03371243557e-7 Ha, respectively, within
6.36e-17 Ha of the individual reference terms. On a skew periodic cell,
finite differences of the native dispersion give maximum errors of 5.43e-14
Ha/e in the charge potential, 4.03e-15 Ha/bohr in the gradient and 1.28e-13
Ha in the strain derivative. Against the Rust oracle, 206 three-atom
geometries covering Z=1..103 and two skew periodic cells (including Ar
self images) agree to 6.0e-19 Ha in per-atom two-body energies, 3.9e-21 Ha
in per-atom ATM energies, 8.7e-19 Ha/bohr in gradient components and
3.1e-18 Ha in strain components. These compare implementations of the
same model. Finite differences test the consistency of its derivatives.
With `--archive /path/to/axtb/crates/gxtb/tests/reference_v2/reference`, the
runner compares 19 molecular cases to archived Grimme program results.
Using the printed charges (six decimal places), the largest differences
are 1.64e-8 Ha for the two-body term and 1.57e-14 Ha for ATM.
[Molecular reference coordinates, SCF density and energies](https://github.com/atomsforhumanity/axtb/tree/65d0f64a9a8211f73b0738ff908b7a414975f50b/crates/gxtb/tests/reference_v2/reference/water)

For the separate distorted `water.xyz` in this directory the **full**
molecular g-xTB energy is -76.430905244116 Ha in both the standalone solver
and Grimme's executable. CP2K at 250 angstrom gives -76.430901049380 Ha,
a residual of +4.19e-6 Ha. Standard molecular GFN2 on the same geometry
agrees with CP2K's stock tblite 0.8.0 CLI within 8.9e-16 Ha in energy
and 9.1e-11 Ha/bohr in gradient at `--acc 0.0001 --etemp 300 --guess sad`.
These checks cover a CPU explicit Gamma supercell. Primitive-cell Bloch
sampling and GPU execution are outside this provider's tested path.
A CP2K build linked with its stock tblite held a Quickstep/GFN2 and a GXTB
force environment in the same process. Running water/GFN2, NaCl/GXTB and
water/GFN2 in succession reproduced both reference energies within 2e-15 Ha.
Counting host calls for a NaCl single point gave 16 two-body energy/potential
evaluations, one ATM evaluation, one two-body gradient evaluation and one
geometry creation/destruction pair.
The two-rank MPI NaCl run passed direct energy, force, stress, displacement,
shear and explicit-supercell comparisons.

## Atomli CPU provider

Build CP2K on this branch with its usual CMake configuration. The CP2K build
needs no electronic provider installed. The tested Atomli adapter requires
Linux, a C compiler, CPython development headers and the shared Python library.
Run these commands from `tools/gxtb` with the CP2K executable on your PATH.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install "atomli-tb==0.1.18" ase numpy cffi setuptools
python build_atomli.py libgxtb_atomli.so
export PYTHONPATH="$(python -c 'import sysconfig; print(sysconfig.get_path("purelib"))')${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONHOME="$(python -c 'import sys; print(sys.base_prefix)')"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cp2k.ssmp -i native.inp
```

`native.inp` contains a distorted two-atom NaCl cell. Its expected energy with
this backend is approximately -622.522431269 Ha. A library pathname containing
a slash is resolved relative to the run directory or can be an absolute path.
The adapter is tied to the Python installation used to build it.

Atomli needs fixed exchange radii in angstrom with `0 < Rd < Rc` and a
shortest supercell lattice translation greater than `2 Rc`. The radii and
supercell size affect the electronic model and require convergence checks.
The adapter uses CPU g-xTB and a fresh electronic solve at each geometry.
Charge, spin, shifted k-point meshes and symmetry reduction have no controls
in this interface. For a fixed basis per atom, a dense diagonalization costs
O((mM)^3) time and O((mM)^2) storage per SCF iteration, where M is the input
cell orbital count and m is the product of the three repeats. Iteration count,
exchange and integral construction also affect the elapsed time.

### Atomli end-to-end checks

Pass executable paths and new output directories to the scripts. Each check
writes its input decks, structures, CP2K outputs and JSON results. The default
socket run starts the independent ASE i-PI client. Supplying `--library`
exercises the in-process force evaluator.

```sh
python check_periodic.py "$(command -v cp2k.ssmp)" checks-nacl --engine atomli --library ./libgxtb_atomli.so --all-derivatives --check-loader-failure
python check_periodic.py "$(command -v cp2k.ssmp)" checks-si --engine atomli --library ./libgxtb_atomli.so --structure silicon --all-derivatives
python check_periodic.py "$(command -v cp2k.psmp)" checks-mpi --engine atomli --library ./libgxtb_atomli.so --ranks 2
python check_periodic.py "$(command -v cp2k.psmp)" checks-socket --engine atomli --ranks 2
python check_molecular_limit.py "$(command -v cp2k.ssmp)" ./libgxtb_atomli.so /path/to/grimme/xtb checks-water
```

The following are maximum absolute component differences from local CPU runs.
The direct reference uses the same periodic backend. The finite differences
use energies from separate CP2K processes.

| Structure | Direct energy (eV/cell) | Direct forces (eV/angstrom) | Direct stress (eV/angstrom^3) | Force differences (eV/angstrom) | Strain differences (eV/angstrom^3) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Distorted NaCl, 2 atoms | 1.17e-8 | 4.48e-10 | 5.81e-10 | 6.62e-9 | 3.93e-10 |
| Distorted diamond Si, 8 atoms | 3.66e-8 | 1.24e-8 | 6.29e-8 | 3.74e-7 | 4.86e-8 |

Both cells were checked for all 3N atomic force and six independent strain
components. Two-rank MPI runs of the native and socket paths also passed the
NaCl direct, displacement, shear and supercell checks. A deliberately missing
CFFI module caused CP2K to abort before printing an energy. The existing
Quickstep H2 regtest passed both of its registered matchers.
`REPEATS 2 1 1` reproduced a separately supplied supercell to the
printed precision in energy and forces, and within 7e-15 eV/angstrom^3 in
stress. Displacement/strain steps were 1e-3 angstrom/2e-4 for NaCl and
2.5e-4 angstrom/1e-5 for Si. For Si, the original 1e-3 angstrom displacement
gave a maximum force error of 2.86e-4 eV/angstrom; a 5e-5 strain gave a
maximum stress error of 1.07e-4 eV/angstrom^3. A second Si strain check at
5e-6 agrees with the 1e-5 result. The thresholds for both structures are
2e-4 eV/angstrom and 2e-5 eV/angstrom^3.

The independent water control uses Grimme's modified xtb 6.7.1 with
`--gxtb --grad --acc 0.01`. The isolated Atomli energy differs from
-76.430905244116 Ha by -6.23e-7 Ha. With a periodic 80 angstrom cell,
changing the radii from 4/6 to 24/36 angstrom leaves a +0.02478353 Ha
difference. At 250 angstrom and radii 70/105 angstrom the maximum force
component difference from the official result is 9.41e-4 Ha/bohr. CP2K and
direct periodic Atomli agree within 6e-10 Ha for every cell and radius tested.

In a separate release comparison, the 250 angstrom periodic-minus-molecular
water energy was -1.34e-8 Ha with Atomli 0.1.9 and +0.024784575 Ha with
0.1.11, 0.1.12 and 0.1.18. Atomli 0.1.11 documents setting the periodic
bond-order correction to zero. Zeroing only the corresponding model
coefficients in a temporary copy of Grimme's molecular executable reproduced
the water periodic version jump to 4.9e-10 Ha. [Release notes](https://github.com/atomag/atomli-assets/releases/tag/v0.1.11)

## Quickstep bond-order experiment

`bond_order_gamma.f90` independently evaluates the closed-shell Mayer
bond-order contribution in a real explicit supercell. It takes the
spin-summed AO density, overlap, AO-to-atom map and symmetric atom-pair
weights. It returns the energy and its derivatives with respect to both
AO matrices at fixed pair weights. For a bond between atoms A and B the
energy is the pair weight times the sum of (PS)_ij (PS)_ji over i on A
and j on B.

`check_bond_order.py` runs Quickstep with `METHOD XTB`, `GFN_TYPE TBLITE`,
and `METHOD GFN2` and reads its converged `P_CSR_WRITE` and `S_CSR_WRITE`
matrices. It compiles the Fortran operator separately and applies it to those
fixed densities. The H-H fixtures use the g-xTB bond-correlation parameters
and minimum-image distances. This is an operator test. A full self-consistent
g-xTB calculation also needs the g-xTB basis, one-electron Hamiltonian,
the other exchange and response terms, and screened dispersion.

```sh
FC=gfortran python check_bond_order.py /path/to/cp2k.ssmp checks-bond-order
```

The executable must include tblite. The runner uses NumPy, ASE and a
Fortran compiler and runs seven single-rank calculations in fresh directories.
The isolated H2, 250 angstrom periodic H2, boundary-crossing H2 and a doubled
cell all gave -0.0098632492 Ha per molecule. Disabling the term in the
independent molecular reference raised the H2 energy by +0.0098632492 Ha.
On the distorted four-H chain, finite differences of all independent
entries of P and S agreed with the analytic matrix responses within
3e-14 Ha. On the full 3x1x1 mesh the runner transforms the complex matrices
into a three-cell supercell with the phase for each AO's lattice image.
The bond-order energy per primitive cell agreed with an independently
converged explicit Gamma supercell to 1.3e-13 Ha. The k-space density
and overlap response contractions agreed with central differences within
1e-12 Ha. The supercell transform here is a dense numerical reference.
Nuclear forces, stress and self-consistent changes in the electronic
density require additional model terms and integration.

## Other implementations

The official development repository publishes Grimme's molecular executable.
Paderborn's Dynamics of Condensed Matter group publishes a separate CP2K
implementation with primitive-cell Bloch sampling, Brillouin-zone-coupled
exchange and symmetry reduction. The Atomli wheel used here declares the
PolyForm Noncommercial 1.0.0 license. Install it separately for uses its
terms permit.

- [Grimme g-xTB distribution and method reference](https://github.com/grimme-lab/g-xtb)
- [Atomli binary distribution](https://github.com/atomag/atomli-assets/releases/tag/v0.1.18)
- [Paderborn implementation](https://github.com/DCM-Uni-Paderborn/cp2k/tree/g-xTB-pbc)
- [Paderborn validation data and source provenance](https://github.com/DCM-Uni-Paderborn/Periodic-g-xTB-Benchmarks)
