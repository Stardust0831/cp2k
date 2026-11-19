#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Evaluate periodic g-xTB forces for a CP2K i-PI server."""

import argparse
from math import prod

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.socketio import SocketClient
from ase.io import read


class PeriodicGXTB(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def __init__(self, engine, repeat, radii=None):
        super().__init__()
        self.repeat = tuple(repeat)
        self.engine = engine
        if len(self.repeat) != 3 or min(self.repeat) < 1:
            raise ValueError("Three positive supercell factors are required")
        if engine == "save-tblite":
            from tblite.ase import TBLite

            if radii is not None:
                raise ValueError("Fixed exchange radii apply only to Atomli")
            self.model = TBLite(method="g-xTB", guess="eeqbc", verbosity=0)
        elif engine == "atomli":
            from atomli_model import atomli_model

            self.model = atomli_model(radii)
        else:
            raise ValueError(f"Unknown g-xTB engine {engine}")

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        if not self.atoms.pbc.all() or self.atoms.cell.volume <= 0:
            raise ValueError("A nonsingular three-dimensional periodic cell is required")
        cell = self.atoms.repeat(self.repeat)
        if self.engine == "atomli":
            from atomli import Atoms

            cell = Atoms(
                symbols=cell.get_chemical_symbols(),
                positions=cell.positions,
                cell=cell.cell.array,
                pbc=True,
            )
        cell.calc = self.model
        count = prod(self.repeat)
        energy = cell.get_potential_energy()
        forces = np.asarray(cell.get_forces()).reshape(count, len(self.atoms), 3)
        stress = np.asarray(cell.get_stress()).reshape(6)
        if not (np.isfinite(energy) and np.isfinite(forces).all() and np.isfinite(stress).all()):
            raise ValueError("The periodic g-xTB evaluation returned a non-finite result")
        # ASE tiles in complete copies of the input atom list. Moving a basis
        # atom moves each image, so its per-cell force is the image average.
        self.results.update(energy=energy / count, forces=forces.mean(axis=0), stress=stress)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geometry", help="XYZ file used by CP2K, in the same atom order")
    parser.add_argument("--engine", required=True, choices=("save-tblite", "atomli"))
    parser.add_argument("--repeat", type=int, nargs=3, default=(1, 1, 1), metavar=("NA", "NB", "NC"))
    parser.add_argument("--radii", type=float, nargs=2, metavar=("RD", "RC"),
                        help="Fixed Atomli exchange radii in angstrom")
    parser.add_argument("--socket", required=True, help="CP2K MOTION/DRIVER/HOST with PREFIX ipi")
    args = parser.parse_args()
    atoms = read(args.geometry, index=0)
    atoms.pbc = True
    atoms.calc = PeriodicGXTB(args.engine, args.repeat, args.radii)
    client = SocketClient(unixsocket=args.socket)
    client.run(atoms, use_stress=True)


if __name__ == "__main__":
    main()
