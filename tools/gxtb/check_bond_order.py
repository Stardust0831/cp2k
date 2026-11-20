"""Exercise the Gamma bond-order contraction with matrices from Quickstep GFN2."""
import argparse
import ctypes
import json
import math
import os
from pathlib import Path
import shlex
import subprocess

import numpy as np
from ase.geometry import find_mic
from ase.units import Bohr


def build_kernel(directory):
    library = directory / "libbondorder.so"
    source = Path(__file__).with_name("bond_order_gamma.f90").resolve()
    subprocess.run([*shlex.split(os.getenv("FC", "gfortran")), "-O2", "-fPIC", "-shared",
                    str(source), "-o", str(library)], cwd=directory, check=True)
    routine = ctypes.CDLL(str(library)).gxtb_gamma_bond_order
    real_matrix = np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags="F_CONTIGUOUS")
    indices = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS")
    routine.argtypes = [ctypes.c_int, ctypes.c_int, indices, real_matrix, real_matrix, real_matrix,
                        ctypes.POINTER(ctypes.c_double), real_matrix, real_matrix]
    routine.restype = None

    def evaluate(owner, weights, density, overlap):
        n, m = density.shape
        assert n == m == len(owner) and overlap.shape == (n, n)
        nat = len(weights)
        assert weights.shape == (nat, nat) and 1 <= min(owner) <= max(owner) <= nat
        fock = np.empty((n, n), order="F")
        overlap_response = np.empty_like(fock)
        energy = ctypes.c_double()
        routine(n, nat, np.ascontiguousarray(owner, dtype=np.int32), np.asfortranarray(weights),
                np.asfortranarray(density), np.asfortranarray(overlap),
                ctypes.byref(energy), fock, overlap_response)
        return energy.value, fock, overlap_response

    return evaluate


def read_matrix(path, n):
    entries = np.loadtxt(path, ndmin=2)
    rows, cols = entries[:, :2].astype(int).T - 1
    assert entries.shape[1] in (3, 4) and min(rows) >= 0 and max(rows) < n
    assert min(cols) >= 0 and max(cols) < n
    values = entries[:, 2] if entries.shape[1] == 3 else entries[:, 2] + 1j*entries[:, 3]
    matrix = np.zeros((n, n), dtype=values.dtype)
    matrix[rows, cols] = values
    np.testing.assert_allclose(matrix, matrix.T.conj(), atol=1e-11, rtol=0)
    return matrix


def quickstep(executable, directory, coords, cell, periodic, mesh=1):
    directory.mkdir()
    xyz = "\n".join(f"      H {x:.12f} {y:.12f} {z:.12f}" for x, y, z in coords)
    lattice = "\n".join(f"      {abc} {x} {y} {z}" for abc, (x, y, z) in zip("ABC", cell))
    kpoints = "" if mesh == 1 else f"""
    &KPOINTS
      SCHEME MONKHORST-PACK {mesh} 1 1
      FULL_GRID TRUE
      SYMMETRY FALSE
      WAVEFUNCTIONS COMPLEX
      VERBOSE TRUE
    &END
"""
    input_text = f"""&GLOBAL
  PROJECT probe
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
{kpoints}
    &QS
      METHOD XTB
      EPS_DEFAULT 1.0E-12
      &XTB
        GFN_TYPE TBLITE
        SCC_MIXER TBLITE
        &TBLITE
          METHOD GFN2
          ACCURACY 0.01
        &END
      &END
    &END
    &SCF
      SCF_GUESS MOPAC
      EPS_SCF 1.0E-11
      MAX_SCF 200
    &END
    &PRINT
      &S_CSR_WRITE ON
      &END
      &P_CSR_WRITE ON
      &END
    &END
  &END
  &SUBSYS
    &CELL
{lattice}
      PERIODIC {'XYZ' if periodic else 'NONE'}
    &END
    &COORD
{xyz}
    &END
  &END
&END
"""
    (directory / "cp2k.inp").write_text(input_text)
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    subprocess.run([executable, "-i", "cp2k.inp", "-o", "cp2k.log"], cwd=directory, env=env,
                   check=True, timeout=120)
    log = (directory / "cp2k.log").read_text()
    assert "SCF run converged" in log and "ENERGY| Total FORCE_EVAL" in log
    n = len(coords)  # One 1s AO per H in the Quickstep GFN2 fixture.
    if mesh == 1:
        p, s = [read_matrix(next(directory.glob(f"probe-{key}_SPIN_1-*.csr")), n) for key in ("P", "S")]
        np.testing.assert_allclose(np.trace(p @ s), n, atol=1e-8, rtol=0)
        return p, s
    blocks = []
    for key in ("P", "S"):
        blocks.append(np.stack([read_matrix(next(directory.glob(f"probe-{key}_SPIN_1_K_{k}-*.csr")), n)
                                for k in range(1, mesh+1)]))
    p, s = blocks
    np.testing.assert_allclose(np.trace(p @ s, axis1=1, axis2=2).mean(), n, atol=1e-8, rtol=0)
    return p, s


def hydrogen_kernel(coords, cell, periodic):
    # Molecular g-xTB H-H bond-order scale, switching radius in bohr and
    # damping width in bohr. Minimum-image extension for this Gamma test.
    i, j = np.triu_indices(len(coords), k=1)
    dr = coords[i] - coords[j]
    distance = find_mic(dr, cell)[1] if periodic else np.linalg.norm(dr, axis=1)
    weight = np.array([-0.0098632492 * 0.5 * math.erfc((r/Bohr - 14.3839259001) / 2.1823)
                       for r in distance])
    pairs = np.zeros((len(coords), len(coords)))
    pairs[i, j] = weight
    pairs[j, i] = weight
    return pairs


def check_responses(evaluate, owner, weights, p, s):
    energy, fock, ds = evaluate(owner, weights, p, s)
    np.testing.assert_allclose(0.5 * np.sum(fock * p), energy, atol=1e-12, rtol=0)
    np.testing.assert_allclose(0.5 * np.sum(ds * s), energy, atol=1e-12, rtol=0)
    errors = []
    step = 1e-4
    for vary_density, analytic in ((True, fock), (False, ds)):
        max_error = 0.
        for i, j in zip(*np.triu_indices(len(p))):
            d = np.zeros_like(p)
            d[i, j] = d[j, i] = step
            if vary_density:
                plus, minus = evaluate(owner, weights, p+d, s), evaluate(owner, weights, p-d, s)
            else:
                plus, minus = evaluate(owner, weights, p, s+d), evaluate(owner, weights, p, s-d)
            observed = (plus[0] - minus[0]) / (2*step)
            expected = analytic[i, j] * (1 if i == j else 2)
            max_error = max(max_error, abs(observed - expected))
        assert max_error < 1e-10
        errors.append(max_error)
    return dict(energy_Ha=energy, max_dP_error=errors[0], max_dS_error=errors[1])


def check_bloch_folding(executable, directory, evaluate, coords, cell):
    # Three equally weighted k points, ordered as in CP2K's full odd MP grid.
    wavevectors = np.array([-1., 0., 1.]) / 3
    nk = len(wavevectors)
    n = len(coords)
    pk, sk = quickstep(executable, directory/"bloch", coords, cell, True, mesh=nk)
    images = np.arange(nk)
    positions = np.vstack([coords + i*cell[0] for i in images])
    supercell = cell.copy()
    supercell[0] *= nk
    p, s = quickstep(executable, directory/"three_cells", positions, supercell, True)
    weights = hydrogen_kernel(positions, supercell, True)
    owner = np.arange(1, nk*n + 1)
    # Quickstep's real-space neighbor list uses centered fractional
    # coordinates. Attach each AO's corresponding integer-cell phase.
    offsets = np.floor(np.linalg.solve(cell.T, coords.T)[0] + 0.5).astype(int)
    phase = np.exp(2j*np.pi*(images[:, None, None] + offsets[None, :, None])*wavevectors) / np.sqrt(nk)
    basis = np.einsum("aik,ij->aikj", phase, np.eye(n)).reshape(nk*n, nk, n)

    def to_supercell(blocks):
        dense = np.einsum("aki,kij,bkj->ab", basis, blocks, basis.conj())
        assert np.max(abs(dense.imag)) < 1e-11
        return dense.real

    pp, ss = to_supercell(pk), to_supercell(sk)
    np.testing.assert_allclose(ss, s, atol=2e-11, rtol=0)
    np.testing.assert_allclose(pp, p, atol=1e-8, rtol=0)
    e, f, ds = evaluate(owner, weights, pp, ss)
    direct = evaluate(owner, weights, p, s)[0]
    np.testing.assert_allclose(e/nk, direct/nk, atol=1e-10, rtol=0)
    # Hermitian variations paired by time reversal give real supercell
    # matrices. Project the response back to each k and test dE/cell.
    rng = np.random.default_rng(42)
    a = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    direction = a + a.T.conj()
    variations = np.stack([direction, direction.real, direction.conj()])
    zero = np.zeros_like(variations)
    errors = []
    step = 1e-5
    for dp, d_overlap, response in ((variations, zero, f), (zero, variations, ds)):
        energies = []
        for sign in (1, -1):
            perturbed_p = to_supercell(pk + sign*step*dp)
            perturbed_s = to_supercell(sk + sign*step*d_overlap)
            energies.append(evaluate(owner, weights, perturbed_p, perturbed_s)[0])
        observed = (energies[0] - energies[1]) / (2*step*nk)
        derivative = np.einsum("aki,ab,bkj->kij", basis.conj(), response, basis)
        expected = np.vdot(derivative, variations).real / nk
        assert abs(observed - expected) < 2e-10
        errors.append(abs(observed - expected))
    return dict(energy_per_cell_Ha=e/nk, direct_supercell_per_cell_Ha=direct/nk,
                delta_Ha_per_cell=(e-direct)/nk, max_delta_P=float(np.max(abs(pp-p))),
                max_delta_S=float(np.max(abs(ss-s))), dPk_error=errors[0], dSk_error=errors[1])


def main(executable, directory):
    directory.mkdir(parents=True, exist_ok=False)
    evaluate = build_kernel(directory)
    box = np.diag([24., 24., 24.])
    dimer = np.array([[2., 2., 2.], [2.74, 2., 2.]])
    wrap = np.array([[23.63, 2., 2.], [.37, 2., 2.]])
    chain = np.array([[1.0, 2., 2.], [1.79, 2., 2.], [3.18, 2.1, 2.], [4.13, 2., 2.]])
    cases = [("isolated", dimer, box, False), ("large", dimer, np.eye(3)*250., True),
             ("wrapped", wrap, box, True),
             ("repeated", np.vstack([dimer, dimer + [24., 0., 0.]]), np.diag([48., 24., 24.]), True),
             ("distorted_chain", chain, np.diag([6.7, 20., 20.]), True)]
    results = {}
    for name, coords, cell, periodic in cases:
        p, s = quickstep(executable, directory/name, coords, cell, periodic)
        weights = hydrogen_kernel(coords, cell, periodic)
        owner = np.arange(1, len(p)+1)
        results[name] = check_responses(evaluate, owner, weights, p, s)
        if name == "distorted_chain":
            # An independent pair sum with multiple AOs per atom checks the
            # general AO mapping on the same nonsymmetric PS input.
            grouped = np.array([1, 1, 2, 2])
            pair = np.array([[0., -0.02], [-0.02, 0.]])
            grouped_energy = evaluate(grouped, pair, p, s)[0]
            ps = p @ s
            assert np.max(abs(ps - ps.T)) > 1e-3
            expected = -0.02 * np.trace(ps[:2, 2:] @ ps[2:, :2])
            np.testing.assert_allclose(grouped_energy, expected, atol=1e-14)
    for name in ("isolated", "large", "wrapped", "repeated"):
        expected = -0.0098632492 * (2 if name == "repeated" else 1)
        np.testing.assert_allclose(results[name]["energy_Ha"], expected, atol=1e-10, rtol=0)
    results["bloch_3x1x1"] = check_bloch_folding(executable, directory, evaluate, chain, np.diag([6.7, 20., 20.]))
    (directory/"results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cp2k", type=Path)
    parser.add_argument("out", type=Path, help="new directory for inputs and results")
    args = parser.parse_args()
    main(args.cp2k.resolve(), args.out.resolve())
