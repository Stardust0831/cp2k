#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Compare the neutral water molecular limit with the Grimme g-xTB binary."""

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
from ase.io import read, write
from ase.units import Bohr, Hartree
from atomli import Atoms
from atomli.calculators.xtb import XTB

from check_periodic import HERE, run_cp2k
from periodic_client import PeriodicGXTB


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("cp2k", "library", "xtb", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True)
    water = read(HERE / "water.xyz")
    write(out / "water.xyz", water, format="xyz")
    with (out / "xtb.log").open("w") as log:
        subprocess.run([args.xtb.resolve(), "water.xyz", "--gxtb", "--grad", "--acc", "0.01"],
                       cwd=out, stdout=log, stderr=subprocess.STDOUT, check=True)
    values = [line for line in (out / "water.engrad").read_text().splitlines()
              if line.strip() and not line.startswith("#")]
    n = int(values[0])
    assert n == len(water)
    e0 = float(values[1])
    f0 = -np.array(values[2:2 + 3*n], float).reshape(n, 3)

    mol = Atoms(symbols=water.get_chemical_symbols(), positions=water.positions)
    mol.calc = XTB(method="gxtb", device="cpu", scf={"energy_threshold_ev": 1e-9,
                   "density_threshold": 1e-9, "max_iterations": 250})
    data = {"xtb_energy_Ha": e0, "isolated_atomli_minus_xtb_Ha": mol.get_potential_energy()/Hartree - e0,
            "isolated_force_error_Ha_per_bohr": float(abs(np.asarray(mol.get_forces())*Bohr/Hartree - f0).max()),
            "periodic": []}
    for length, rd, rc in [(20, 4, 6), (40, 4, 6), (80, 4, 6), (80, 8, 12),
                           (80, 16, 24), (80, 24, 36), (160, 40, 60), (250, 70, 105)]:
        water.set_cell(np.eye(3) * length)
        water.pbc = True
        energy, force, stress = run_cp2k(out/f"L{length}-Rc{rc}", water, [args.cp2k.resolve()],
                                         "atomli", (1, 1, 1), (rd, rc), args.library.resolve())
        water.calc = PeriodicGXTB("atomli", (1, 1, 1), (rd, rc))
        data["periodic"].append({
            "cell_A": length, "radii_A": [rd, rc], "energy_Ha": energy/Hartree,
            "minus_xtb_Ha": energy/Hartree - e0,
            "force_error_Ha_per_bohr": float(abs(force*Bohr/Hartree - f0).max()),
            "minus_direct_periodic_Ha": (energy - water.get_potential_energy())/Hartree,
        })
    print(json.dumps(data, indent=2))
    (out / "comparison.json").write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
