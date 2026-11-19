# Periodic g-xTB in CP2K

`FORCE_EVAL / METHOD GXTB` evaluates a neutral, three-dimensionally periodic
cell through a shared library in the CP2K process. CP2K constructs the
Born-von Karman supercell from `REPEATS na nb nc`, passes its lattice and
atoms to a Gamma-point electronic solver, and returns the energy per input
cell, the mean force on equivalent images and the pressure-positive virial.
One MPI rank calls the provider and broadcasts the results. The CP2K changes
and adapters were written against public interfaces. Electronic structure and
parameters come from the separately installed provider.

**Physical limitation in the tested backend.** For the water geometry in this
directory, the Atomli 0.1.18 periodic energy at 250 angstrom exceeds the
official molecular g-xTB energy by 0.02478395 Ha. The discrepancy is also
present when evaluating periodic Atomli directly. This branch is experimental.

## Run a single point

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

## End-to-end checks

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

## Other implementations

`save_tblite.c` is an optional C adapter to the periodic save_tblite API.
Compile against a build exporting `tblite_new_gxtb_calculator`; omit
`EXCHANGE_RADII` for this provider. The adapter compiled against the snapshot
`35e7942b60edd89bb407ab3da5768d3410af83f5` headers. Its runtime has not
been qualified here. Building that snapshot stopped at its pinned
`grimme-lab/dftd` dependency `99d64ee83832fda03df2c7dd7b7fd9f4c9ebb098`.

```sh
cc -std=c11 -shared -fPIC save_tblite.c $(pkg-config --cflags --libs tblite) -o libgxtb_tblite.so
```

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
