# g-xTB 2.0.1 parameters

`src/gxtb_revd4_data.inc` contains numerical model parameters for H through
Lr. The g-xTB radii, charges and interpolation parameters originate with
T. Froitzheim, M. Mueller, A. Hansen and S. Grimme. The
imaginary-frequency polarizabilities originate with DFT-D4. The numerical
source for this import is Colton Hicks's axtb release 0.2.0, commit
`65d0f64a9a8211f73b0738ff908b7a414975f50b`. Its g-xTB tables are distributed
under LGPL-3.0-or-later. The upstream notice records the original works,
licenses and the authors' permission for the axtb distribution.
[Upstream notice and license texts](data-licenses/) [g-xTB method](https://doi.org/10.26434/chemrxiv-2025-bjxvt)
[Release data](https://github.com/atomsforhumanity/axtb/tree/65d0f64a9a8211f73b0738ff908b7a414975f50b/crates/gxtb/src/data)

To regenerate the table, check out this revision of axtb, then run the
following from the CP2K source root.

```sh
python3 tools/gxtb/revd4/import_parameters.py /path/to/axtb src/gxtb_revd4_data.inc
python3 tools/gxtb/revd4/import_electronic.py /path/to/axtb src/gxtb_electronic_data.inc
python3 tools/gxtb/revd4/import_basis.py /path/to/axtb src/gxtb_basis_data.inc
```

The importer checks the SHA-256 of each source and reads literal numerical
arrays. It subtracts the charge-scaled secondary-system polarizability,
applies the reference-charge normalization and folds the square root of
each Casimir-Polder quadrature weight into the reference spectrum. At run
time CP2K interpolates these spectra at each atom's periodic coordination
number for each partner element. Their dot product gives the geometry part
of C6. The energy, SCF potential, periodic image sums and derivatives
are evaluated by the CP2K implementation.

The electronic tables provide shell interactions, spin constants, exchange,
bond correlation, repulsion and core increments. The basis import reads
the numerical primitives and charge-dependent contraction coefficients
from the release's `crates/gxtb/param/basisq`, along with the numeric
Hamiltonian, projector and EEQBC parameters. The shared spherical Gaussian
primitives and diatomic scaling are taken from the standard CP2K tblite
dependency at build time. The CP2K algorithms are independent implementations.
[q-vSZP basis reference](https://doi.org/10.1063/5.0133026)
