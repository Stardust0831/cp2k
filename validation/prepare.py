"""Prepare small lifecycle fixtures and copy the unchanged solid benchmarks."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("source", type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parent
source = args.source.resolve()
inputs = root / "inputs"
inputs.mkdir(exist_ok=True)
upstream = source / "tests/QS/regtest-cusolver"
real = (upstream / "Si8-generalized.inp").read_text()
# The additional full-spectrum/block-layout cases use the same diagnostic
# accuracy as the upstream complex generalized fixture. Keep checks enabled.
real = real.replace("EPS_CHECK_DIAG 1.0E-14", "EPS_CHECK_DIAG 1.0E-12")
complex_input = (upstream / "Si8-generalized-complex.inp").read_text()
fixtures = [
    ("real-standard", real.replace("DIRECT_GENERALIZED_DIAGONALIZATION T", "DIRECT_GENERALIZED_DIAGONALIZATION F"), 20, 64, "cp_fm_diag_cusolver"),
    ("real-subset", real, 20, 64, "cp_fm_general_cusolver"),
    ("real-partial", real, 1, 8, "cp_fm_general_cusolver"),
    ("real-full", real, -1, 64, "cp_fm_general_cusolver"),
    ("complex-subset", complex_input, 20, 64, "cp_cfm_general_cusolver"),
    ("complex-full", complex_input, -1, 64, "cp_cfm_general_cusolver"),
]
manifest = {"fixtures": [], "files": {}}
for name, text, added_mos, block_size, backend in fixtures:
    # Preserve 104 basis functions while avoiding the poorly conditioned
    # high-energy subspace of DZVP-GTH in the additional full-spectrum cases.
    if added_mos == -1:
        text = text.replace("BASIS_SET_FILE_NAME GTH_BASIS_SETS", "BASIS_SET_FILE_NAME BASIS_MOLOPT")
        text = text.replace("BASIS_SET DZVP-GTH", "BASIS_SET DZVP-MOLOPT-SR-GTH")
    text = text.replace("PROJECT Si8-generalized-complex", "PROJECT " + name)
    text = text.replace("PROJECT Si8-generalized", "PROJECT " + name)
    text = text.replace("RUN_TYPE ENERGY", "RUN_TYPE ENERGY_FORCE")
    text = text.replace("ADDED_MOS 20 20", f"ADDED_MOS {added_mos} {added_mos}")
    text = text.replace("&END GLOBAL", "  BLACS_GRID SQUARE\n  &FM\n"
                        f"    NROW_BLOCKS {block_size}\n    NCOL_BLOCKS {block_size}\n"
                        "    FORCE_BLOCK_SIZE T\n  &END FM\n&END GLOBAL")
    text = text.replace("../sample_xyz/SI_8.xyz", "SI_8.xyz")
    (inputs / (name + ".inp")).write_text(text)
    manifest["fixtures"].append(dict(name=name, backend=backend, atoms=8,
                                     nmo=104 if added_mos == -1 else 16 + added_mos,
                                     matrix_size=104, block_size=block_size,
                                     basis="DZVP-MOLOPT-SR-GTH" if added_mos == -1 else "DZVP-GTH"))
shutil.copy2(source / "tests/QS/sample_xyz/SI_8.xyz", inputs / "SI_8.xyz")
old = root.parents[1] / "cp2k-elpa-dbcsr-retest-20261105/raw/inputs"
for family in ("ag-gamma", "diamond-gamma"):
    target = inputs / family
    target.mkdir(exist_ok=True)
    previous = old / (family + "-cusolvermp-direct")
    for name in ("input.inp", "BASIS_MOLOPT", "POTENTIAL", "provenance.json"):
        shutil.copy2(previous / name, target / name)
for path in sorted(inputs.rglob("*")):
    if path.is_file():
        manifest["files"][str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
manifest["provenance"] = {
    "small_fixture_source": str(upstream),
    "small_fixture_changes": "RUN_TYPE ENERGY_FORCE, explicit square FM blocks, requested MO counts, EPS_CHECK_DIAG 1e-12, full-spectrum cases use DZVP-MOLOPT-SR-GTH (same 104 functions); same Si8 geometry and five SCF steps",
    "solid_inputs": str(old),
    "solid_input_changes": "none; exact original Multiwfn-derived direct-generalized input, basis and potential",
}
(root / "fixtures.json").write_text(json.dumps(manifest, indent=2) + "\n")
(root / "TEST_DIRS").write_text("QS/regtest-cusolver cusolvermp\n")
print(json.dumps(manifest, indent=2))
