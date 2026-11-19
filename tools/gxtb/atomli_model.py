# SPDX-License-Identifier: GPL-2.0-or-later
def atomli_model(radii):
    from atomli.calculators.xtb import XTB

    if radii is None or not 0 < radii[0] < radii[1]:
        raise ValueError("Atomli requires 0 < Rd < Rc in angstrom")
    return XTB(method="gxtb", device="cpu", gamma_exchange_radii=tuple(radii),
               scf={"energy_threshold_ev": 1e-9, "density_threshold": 1e-9,
                    "max_iterations": 250})
