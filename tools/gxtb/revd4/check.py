"""Validate the D4Srev C ABI against molecular reference data and periodic differences."""
import argparse
import ctypes as ct
import json
from pathlib import Path

import numpy as np


class Dispersion:
    def __init__(self, library):
        self.lib = ct.CDLL(str(library.resolve()))
        reals = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")
        ints = np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS")
        self.lib.cp2k_revd4_new.argtypes = [ct.c_int, ints, reals, reals, ints, ct.c_char_p, ct.c_size_t]
        self.lib.cp2k_revd4_new.restype = ct.c_void_p
        for name, signature in {
            "free": [], "scf": [reals]*3, "atm": [reals]*3, "derivatives": [reals]*3,
        }.items():
            routine = getattr(self.lib, "cp2k_revd4_" + name)
            routine.argtypes = [ct.c_void_p, *signature]
            routine.restype = None

    def __call__(self, z, xyz, cell, pbc, q):
        n = len(z)
        z = np.ascontiguousarray(z, dtype=np.int32)
        pbc = np.ascontiguousarray(pbc, dtype=np.int32)
        xyz, cell, q = [np.ascontiguousarray(value, dtype=float) for value in (xyz, cell, q)]
        msg = ct.create_string_buffer(256)
        handle = self.lib.cp2k_revd4_new(n, z, xyz, cell, pbc, msg, len(msg))
        if not handle:
            raise ValueError(msg.value.decode())
        try:
            e2, e3, v = np.empty((3,n))
            g2, g3 = np.empty((2,n,3))
            s2, s3 = np.empty((2,3,3))
            self.lib.cp2k_revd4_scf(handle, q, e2, v)
            self.lib.cp2k_revd4_atm(handle, e3, g3, s3)
            self.lib.cp2k_revd4_derivatives(handle, q, g2, s2)
            return np.array([sum(e2),sum(e3)]), v, g2+g3, s2+s3
        finally:
            self.lib.cp2k_revd4_free(handle)


def validate(calc):
    # Grimme 2.0.1 program, last SCF density, atomic units; source in README.
    water = np.array([[0,0,.119262],[0,.763239,-.477047],[0,-.763239,-.477047]]) / .529177210903
    q = [-.13335217008675593-.4956216773371711,.31448692371196296,.3144869237119633]
    reference = np.array([-4.42494068950694e-5,3.0337124355412764e-7])
    e, *_ = calc([8,1,1], water, np.zeros((3,3)), [0,0,0], q)
    np.testing.assert_allclose(e, reference, atol=2e-13, rtol=0.)
    residual = np.max(abs(e-reference))
    large, *_ = calc([8,1,1], water, np.eye(3)*250/.529177210903, [1,1,1], q)
    np.testing.assert_allclose(e, large, atol=1e-16, rtol=0.)
    # A distorted, partially charged cell samples cross-boundary pairs,
    # self images and neutral-C6 ATM triples with repeated atom indices.
    lattice = np.array([[14.1,0,0],[2.2,13.8,0],[.9,-1.5,14.6]])
    coords = np.array([[1.2,2.1,2.3],[3.0,2.2,2.6],[.8,3.6,2.1]])
    charge = np.array([-.53,.29,.24])
    numbers = np.array([8,1,1])
    periodic = [1,1,1]

    def evaluate(x=coords, box=lattice, q=charge):
        return calc(numbers,x,box,periodic,q)

    energy, v, g, s = evaluate()
    step = 2e-5
    errors = [0.,0.,0.]
    for i in range(len(numbers)):
        dq = np.zeros(len(numbers))
        dq[i] = step
        de = (evaluate(q=charge+dq)[0] - evaluate(q=charge-dq)[0])/(2*step)
        np.testing.assert_allclose(de[1], 0., atol=1e-16)
        errors[0] = max(errors[0],abs(de[0]-v[i]))
        for a in range(3):
            dx = np.zeros_like(coords)
            dx[i,a] = step
            fd = (sum(evaluate(coords+dx)[0])-sum(evaluate(coords-dx)[0]))/(2*step)
            errors[1] = max(errors[1],abs(fd-g[i,a]))
    for a in range(3):
        for b in range(3):
            strain = np.zeros((3,3))
            strain[a,b] = step
            plus, minus = np.eye(3)+strain, np.eye(3)-strain
            fd = (sum(evaluate(coords@plus,lattice@plus)[0])
                  -sum(evaluate(coords@minus,lattice@minus)[0]))/(2*step)
            errors[2] = max(errors[2],abs(fd-s[a,b]))
    assert max(errors) < 2e-10, errors
    wrapped = coords.copy()
    wrapped[1] += lattice[0]-2*lattice[1]
    for actual, original in zip(evaluate(wrapped), (energy,v,g,s)):
        np.testing.assert_allclose(actual, original, atol=1e-13, rtol=0.)
    box = lattice.copy()
    box[0] *= 2
    e_super, v_super, g_super, s_super = calc(np.tile(numbers,2), np.vstack([coords,coords+lattice[0]]),
                                              box, periodic, np.tile(charge,2))
    np.testing.assert_allclose(e_super/2, energy, atol=1e-13, rtol=0.)
    np.testing.assert_allclose(v_super, np.tile(v,2), atol=1e-13, rtol=0.)
    np.testing.assert_allclose(g_super, np.tile(g,(2,1)), atol=1e-13, rtol=0.)
    np.testing.assert_allclose(s_super/2, s, atol=1e-13, rtol=0.)
    one, _, one_g, one_s = calc([18], [[0,0,0]], lattice, periodic, [0.])
    assert one[0] < 0. and one[1] != 0. and np.linalg.norm(one_s)>0.
    np.testing.assert_allclose(one_g, 0., atol=1e-14)
    return dict(water_energy_Ha=e.tolist(), water_reference_max_error_Ha=residual,
                cell_energy_Ha=energy.tolist(), max_dEdq_error_Ha_per_e=errors[0],
                max_gradient_error_Ha_per_bohr=errors[1], max_strain_error_Ha=errors[2],
                supercell_component_error_Ha=float(max(abs(e_super/2-energy))),
                argon_self_image_energies_Ha=one.tolist())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library", type=Path)
    args = parser.parse_args()
    print(json.dumps(validate(Dispersion(args.library)), indent=2))
