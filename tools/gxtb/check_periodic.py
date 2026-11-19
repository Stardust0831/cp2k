#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Check periodic CP2K energies, forces and stress against a direct backend and differences."""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
import uuid
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import write
from ase.units import Bohr, Hartree, bar

from periodic_client import PeriodicGXTB


HERE = Path(__file__).resolve().parent


def run_cp2k(directory, atoms, command, engine, repeat, radii, library):
    directory.mkdir()
    socket = "gxtb-" + uuid.uuid4().hex
    address = Path("/tmp", "ipi_" + socket)
    template = (HERE / "periodic.inp").read_text()
    cell = "\n".join(
        f"      {axis} " + " ".join(f"{v:.15g}" for v in row)
        for axis, row in zip("ABC", atoms.cell)
    )
    deck = re.sub(r"      A [^\n]+\n      B [^\n]+\n      C [^\n]+", cell, template)
    deck = deck.replace("water.xyz", "structure.xyz")
    if library:
        deck = re.sub(r"&MOTION\n.*?&END MOTION\n", "", deck, flags=re.DOTALL)
        options = ["  METHOD GXTB", "  &GXTB", f"    LIBRARY {library}",
                   "    REPEATS " + " ".join(map(str, repeat))]
        if radii:
            options.append("    EXCHANGE_RADII [angstrom] " + " ".join(map(str, radii)))
        deck = deck.replace("  METHOD IPI", "\n".join([*options, "  &END GXTB"]))
    else:
        deck = deck.replace("HOST gxtb", "HOST " + socket)
    (directory / "run.inp").write_text(deck)
    write(directory / "structure.xyz", atoms, format="xyz")
    with (directory / "cp2k.log").open("w") as log:
        process = subprocess.Popen([*command, "-i", "run.inp"], cwd=directory, stdout=log, stderr=log)
        try:
            if library is None:
                deadline = time.monotonic() + 30
                while not address.is_socket():
                    if process.poll() is not None or time.monotonic() > deadline:
                        raise RuntimeError(f"CP2K did not open its socket; see {directory}")
                    time.sleep(0.02)
                args = [sys.executable, str(HERE / "periodic_client.py"), "structure.xyz",
                        "--engine", engine, "--socket", socket, "--repeat", *map(str, repeat)]
                if radii:
                    args.extend(["--radii", *map(str, radii)])
                with (directory / "client.log").open("w") as clientlog:
                    subprocess.run(args, cwd=directory, stdout=clientlog, stderr=clientlog,
                                   timeout=180, check=True)
            if process.wait(timeout=180) != 0:
                raise RuntimeError(f"CP2K failed; see {directory}")
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
            address.unlink(missing_ok=True)
    output = (directory / "cp2k.log").read_text()
    if "PROGRAM ENDED AT" not in output:
        raise RuntimeError(f"Missing CP2K completion marker in {directory}")
    energy = float(re.search(r"ENERGY\| Total FORCE_EVAL .*?\[hartree\]\s+(\S+)", output)[1])
    forces = np.array(re.findall(r"FORCES\|\s+\d+\s+(\S+)\s+(\S+)\s+(\S+) ", output), float)
    lines = output.split("Analytical stress tensor [bar]", 1)[1]
    stress = np.array(re.findall(r"STRESS\| +[xyz] +([-+\d.][\d.E+-]*) +(\S+) +(\S+)", lines)[:3], float)
    # CP2K prints pressure-positive stress. ASE uses tension-positive stress.
    if forces.shape != (len(atoms), 3) or stress.shape != (3, 3):
        raise ValueError(f"Incomplete CP2K force or stress output in {directory}")
    return energy * Hartree, forces * Hartree / Bohr, -stress * bar


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path, help="New directory for full inputs and logs")
    parser.add_argument("--engine", required=True, choices=("save-tblite", "atomli"))
    parser.add_argument("--library", type=Path, help="Run METHOD GXTB with this shared library (default: i-PI)")
    parser.add_argument("--check-loader-failure", action="store_true",
                        help="Simulate a missing dependency in the Atomli CFFI adapter")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--structure", choices=("nacl", "silicon"), default="nacl")
    parser.add_argument("--all-derivatives", action="store_true",
                        help="Check every force component and six independent strain derivatives")
    parser.add_argument("--force-step", type=float, help="Displacement in angstrom")
    parser.add_argument("--strain-step", type=float, help="Dimensionless homogeneous strain")
    args = parser.parse_args()
    binary = args.binary.resolve()
    out = args.output.resolve()
    out.mkdir(parents=True)
    if args.ranks < 1:
        parser.error("--ranks must be positive")
    if args.force_step is None:
        args.force_step = 0.00025 if args.structure == "silicon" else 0.001
    if args.strain_step is None:
        args.strain_step = 0.00001 if args.structure == "silicon" else 0.0002
    if min(args.force_step, args.strain_step) <= 0:
        parser.error("Finite-difference steps must be positive")
    command = ([args.mpiexec, "-n", str(args.ranks)] if args.ranks > 1 else []) + [binary]
    library = args.library.resolve(strict=True) if args.library else None
    if args.check_loader_failure and (not library or args.engine != "atomli"):
        parser.error("--check-loader-failure requires the Atomli shared library")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    if args.structure == "nacl":
        atoms = Atoms("NaCl", scaled_positions=[[0.05, 0.04, 0.03], [0.47, 0.51, 0.49]],
                      cell=[[7.2, 0, 0], [0.36, 7.4, 0], [0.45, 0.18, 7.0]], pbc=True)
        radii = (2.0, 3.3)
    else:
        atoms = bulk("Si", "diamond", a=5.43, cubic=True)
        atoms.cell[1, 0] += 0.04
        atoms.positions[1] += [0.016, -0.023, 0.011]
        radii = (2.1, 2.6)
    if args.engine != "atomli":
        radii = None
    report = {}

    def evaluate(label, structure, repeat=(1, 1, 1)):
        return run_cp2k(out / label, structure, command, args.engine, repeat, radii, library)

    def check(name, error, tolerance):
        report[name] = {"max_abs_error": float(error), "tolerance": tolerance}
        (out / "checks.json").write_text(json.dumps(report, indent=2) + "\n")

    energy, forces, stress = evaluate("reference", atoms)
    if args.check_loader_failure:
        case = out / "failed-python-init"
        case.mkdir()
        for name in ("run.inp", "structure.xyz"):
            shutil.copyfile(out / "reference" / name, case / name)
        (case / "_cffi_backend.py").write_text("raise ImportError('missing test dependency')\n")
        environment = dict(os.environ)
        environment["PYTHONPATH"] = os.pathsep.join([str(case), environment.get("PYTHONPATH", "")])
        with (case / "cp2k.log").open("w") as log:
            failed = subprocess.run([*command, "-i", "run.inp"], cwd=case, env=environment,
                                    stdout=log, stderr=subprocess.STDOUT, timeout=30)
        output = (case / "cp2k.log").read_text(errors="replace")
        if (failed.returncode == 0 or "Atomli Python initialization failed" not in output
                or "ENERGY| Total FORCE_EVAL" in output):
            raise AssertionError(f"Provider initialization failure was not propagated in {case}")
    direct = atoms.copy()
    direct.calc = PeriodicGXTB(args.engine, (1, 1, 1), radii)
    check("energy_eV", abs(energy - direct.get_potential_energy()), 1e-5)
    check("forces_eV_per_A", abs(forces - direct.get_forces()).max(), 1e-5)
    check("stress_eV_per_A3", abs(stress - direct.get_stress(voigt=False)).max(), 1e-6)

    def derivative(label, perturb, step):
        plus, minus = atoms.copy(), atoms.copy()
        perturb(plus, step)
        perturb(minus, -step)
        return (evaluate(label + "_plus", plus)[0] - evaluate(label + "_minus", minus)[0]) / (2 * step)

    force_errors = []
    for atom, axis in (np.ndindex(len(atoms), 3) if args.all_derivatives else [(1, 0)]):
        def displace(structure, delta):
            structure.positions[atom, axis] += delta
        diff = derivative(f"atom{atom}_{axis}", displace, args.force_step)
        force_errors.append(abs(forces[atom, axis] + diff))
    check("fd_force_eV_per_A", max(force_errors), 2e-4)

    stress_errors = []
    for row, col in (zip(*np.tril_indices(3)) if args.all_derivatives else [(1, 0)]):
        def deform(structure, amount):
            transform = np.eye(3)
            transform[row, col] += amount
            structure.set_cell(structure.cell @ transform.T, scale_atoms=True)
        diff = derivative(f"strain{row}_{col}", deform, args.strain_step)
        stress_errors.append(abs(stress[row, col] - diff / atoms.get_volume()))
    check("fd_stress_eV_per_A3", max(stress_errors), 2e-5)
    repeat = (2, 1, 1)
    e_fold, f_fold, s_fold = evaluate("folded", atoms, repeat)
    e_super, f_super, s_super = evaluate("explicit", atoms.repeat(repeat))
    check("supercell_energy_eV", abs(e_fold - e_super / 2), 1e-5)
    check("supercell_forces_eV_per_A",
          abs(f_fold - f_super.reshape(2, len(atoms), 3).mean(axis=0)).max(), 1e-5)
    check("supercell_stress_eV_per_A3", abs(s_fold - s_super).max(), 1e-6)
    print(json.dumps(report, indent=2))
    failures = [name for name, result in report.items()
                if not result["max_abs_error"] <= result["tolerance"]]
    if failures:
        raise SystemExit("Failed checks: " + ", ".join(failures))


if __name__ == "__main__":
    main()
