#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Build an optional CP2K g-xTB provider with the current Python environment."""

import argparse
import inspect
import sysconfig
import tempfile
from importlib.metadata import version
from pathlib import Path

from cffi import FFI
from atomli_model import atomli_model


def register():
    from _cp2k_atomli import ffi
    from atomli import Atoms
    from ase.data import chemical_symbols
    from ase.units import Bohr, Hartree
    import numpy as np

    @ffi.def_extern(error=0)
    def atomli_evaluate(n, z, xyz, cell, radii, energy, gradient, strain, message, capacity):
        try:
            def doubles(pointer, count):
                return np.frombuffer(ffi.buffer(pointer, count * 8), dtype=np.float64)

            rd, rc = doubles(radii, 2) * Bohr
            atoms = Atoms(symbols=[chemical_symbols[z[i]] for i in range(n)],
                          positions=doubles(xyz, n * 3).reshape(n, 3) * Bohr,
                          cell=doubles(cell, 9).reshape(3, 3) * Bohr, pbc=True)
            atoms.calc = atomli_model((rd, rc))
            energy[0] = atoms.get_potential_energy() / Hartree
            doubles(gradient, n * 3)[:] = -np.asarray(atoms.get_forces()).ravel() * Bohr / Hartree
            stress = np.asarray(atoms.get_stress())
            tensor = stress[[0, 5, 4, 5, 1, 3, 4, 3, 2]]
            volume = abs(np.linalg.det(doubles(cell, 9).reshape(3, 3))) * Bohr**3
            doubles(strain, 9)[:] = tensor * volume / Hartree
            return 1
        except Exception as error:
            text = str(error).encode()[:capacity - 1]
            ffi.memmove(message, text + b"\0", len(text) + 1)
            return -1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library", type=Path)
    args = parser.parse_args()
    ffi = FFI()
    ffi.embedding_api("""
        int atomli_evaluate(int, const int *, const double *, const double *,
                           const double *, double *, double *, double *, char *, size_t);
    """)
    label = version("atomli-tb")
    here = Path(__file__).resolve().parent
    soname = sysconfig.get_config_var("INSTSONAME") or sysconfig.get_config_var("LDLIBRARY")
    ffi.set_source("_cp2k_atomli", (here / "atomli_adapter.c").read_text(),
                   include_dirs=[str(here.parents[1] / "src")], libraries=["dl"],
                   define_macros=[("LIBPYTHON_NAME", f'"{soname}"'),
                                  ("ATOMLI_VERSION", f'"{label}"')],
                   # Python extensions and CFFI callbacks live for the interpreter lifetime.
                   extra_link_args=["-Wl,-z,nodelete", "-Wl,--disable-new-dtags",
                                    "-Wl,-rpath," + sysconfig.get_config_var("LIBDIR")])
    ffi.embedding_init_code("import site\nsite.addsitedir(" + repr(sysconfig.get_path("purelib"))
                            + ")\n" + inspect.getsource(atomli_model)
                            + "\n" + inspect.getsource(register) + "\nregister()\n")
    with tempfile.TemporaryDirectory() as temp:
        ffi.compile(tmpdir=temp, target=str(args.library.resolve()))


if __name__ == "__main__":
    main()
