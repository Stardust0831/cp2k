"""Compare actual energies/forces, solver paths, binaries and rank mappings."""
import hashlib
import json
from pathlib import Path
import re

root = Path(__file__).resolve().parent
results = root / "results"
request = json.loads((root / "request.json").read_text())
fixtures = json.loads((root / "fixtures.json").read_text())
number = r"[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?"


def timers(text):
    answer = {}
    for line in text.splitlines():
        match = re.match(r"\s+(\w+)\s+(\d+)\s+", line)
        if match and len(line.split()) == 7:
            try:
                answer[match[1]] = dict(calls=int(match[2]), max_seconds=float(line.split()[-1]))
            except ValueError:
                pass
    return answer


def mapping(path):
    rows = [json.loads(p.read_text()) for p in sorted(path.glob("rank-*.json"))]
    assert len(rows) == 4 and {r["rank"] for r in rows} == set(range(4)), path
    assert len({r["host"] for r in rows}) == 1 and len({r["gpu"] for r in rows}) == 4
    assert all(len(r["cores"]) == 4 for r in rows)
    assert len({tuple(core) for r in rows for core in r["cores"]}) == 16
    return [{k: r[k] for k in ("rank", "host", "gpu", "cores", "job")} for r in rows]


def observables(path):
    text = path.read_text(errors="replace")
    energies = re.findall(r"ENERGY\| Total FORCE_EVAL[^\n]*?(" + number + r")\s*$", text, re.M)
    forces = re.findall(r"^\s*FORCES\|\s+(\d+)\s+(" + number + r")\s+(" + number +
                        r")\s+(" + number + r")\s+" + number + r"\s*$", text, re.M)
    assert len(energies) == 1, (path, energies)
    assert "SCF run converged" in text and "PROGRAM ENDED" in text, path
    return dict(energy_ha=float(energies[0]), forces=[float(x.replace("D", "E"))
                for row in forces for x in row[1:]], atoms=len(forces), timers=timers(text),
                scf_steps=[int(x) for x in re.findall(r"SCF run converged in\s+(\d+) steps", text)])


def compare(a, b, energy_tol, force_tol):
    assert len(a["forces"]) == len(b["forces"])
    delta_e = abs(a["energy_ha"] - b["energy_ha"])
    delta_f = max(abs(x - y) for x, y in zip(a["forces"], b["forces"]))
    assert delta_e <= energy_tol and delta_f <= force_tol, (delta_e, delta_f)
    return dict(energy_delta_ha=delta_e, max_force_component_delta_ha_bohr=delta_f)


variants = {}
common_mapping = None
for variant in request["variants"]:
    name = variant["name"]
    stage = results / name
    report = (results / (name + "-regtest.log")).read_text()
    for kind in ("FAILED", "WRONG"):
        assert re.search(r"Number of\s+" + kind + r"\s+tests\s+0\b", report), (name, kind)
    lifecycle_log = (results / (name + "-lifecycle.log")).read_text()
    assert "CP2K_CUSOLVER_LIFECYCLE_COMPLETE" in lifecycle_log, name
    assert "diagonalization failed" not in lifecycle_log, (name, "orthogonality diagnostic")
    lifecycle = stage / "lifecycle"
    rows = {}
    for path in lifecycle.glob("observables-g*.jsonl"):
        for line in path.read_text().splitlines():
            row = json.loads(line)
            key = f"{row['phase']}/{row['case']}/{row['group']}/{row['step']}"
            assert key not in rows and len(row["forces"]) == 24
            rows[key] = row
    assert len(rows) == 18, (name, len(rows))
    path_evidence = {}
    for fixture in fixtures["fixtures"]:
        output = lifecycle / ("world-" + fixture["name"] + "-g0.out")
        text = output.read_text(errors="replace")
        assert "diagonalization failed" not in text, (name, output, "orthogonality diagnostic")
        backend = fixture["backend"]
        measured = timers(text)
        assert measured[backend]["calls"] >= 2, (name, output, measured)
        nmos = [int(x) for x in re.findall(r"Number of molecular orbitals:\s+(\d+)", text)]
        assert nmos and set(nmos) == {fixture["nmo"]}, (fixture["name"], nmos)
        path_evidence[fixture["name"]] = dict(backend=backend, timer=measured[backend], nmos=nmos)
    for group, case in ((0, "real-subset"), (1, "complex-subset")):
        assert rows[f"split/{case}/{group}/0"]["ranks"] == 2
    for fixture in fixtures["fixtures"]:
        first = rows[f"world/{fixture['name']}/0/0"]["energy_ha"]
        changed = rows[f"world/{fixture['name']}/0/1"]["energy_ha"]
        assert abs(first - changed) > 1e-9, (name, fixture["name"], "geometry did not affect energy")

    current_mapping = mapping(lifecycle)
    if common_mapping is None:
        common_mapping = current_mapping
    assert current_mapping == common_mapping
    solids = {}
    for family, atoms in (("ag-gamma", 500), ("diamond-gamma", 1000)):
        folder = stage / family
        assert mapping(folder) == common_mapping
        for filename in ("input.inp", "BASIS_MOLOPT", "POTENTIAL"):
            assert (folder / filename).read_bytes() == (root / "inputs" / family / filename).read_bytes()
        row = observables(folder / "cp2k.out")
        assert row["atoms"] == atoms
        assert row["timers"]["cp_fm_general_cusolver"]["calls"] > 1
        solids[family] = row
    counters = []
    if request.get("count_context_calls", True):
        counters = [json.loads(p.read_text()) for p in sorted(lifecycle.glob("cusolver-context-counts-rank-*.json"))]
        assert len(counters) == 4 and {r["rank"] for r in counters} == set(range(4))
        assert all(r["created"] > 0 and r["created"] == r["destroyed"] for r in counters), (name, counters)
    variants[name] = dict(lifecycle=rows, lifecycle_path_evidence=path_evidence,
                         context_counts=counters, solids=solids)

comparisons = {}
base = variants["baseline"]
for name, variant in variants.items():
    if name == "baseline":
        continue
    assert variant["lifecycle"].keys() == base["lifecycle"].keys()
    comparisons[name] = dict(lifecycle={key: compare(base["lifecycle"][key], value, 1e-8, 1e-7)
        for key, value in variant["lifecycle"].items()}, solids={})
    if request.get("count_context_calls", True) and name == "reuse":
        assert all(candidate["created"] < old["created"] for old, candidate in
                   zip(base["context_counts"], variant["context_counts"])), "No measured GPU-context reuse"
    for family, row in variant["solids"].items():
        baseline = base["solids"][family]
        assert baseline["scf_steps"] == row["scf_steps"], (name, family, "SCF step count changed")
        comparison = compare(baseline, row, 1e-6, 1e-5)
        for key in ("CP2K", "scf_env_do_scf", "eigensolver_generalized", "cp_fm_general_cusolver"):
            comparison[key] = dict(baseline_seconds=baseline["timers"][key]["max_seconds"],
                                   candidate_seconds=row["timers"][key]["max_seconds"],
                                   speedup=baseline["timers"][key]["max_seconds"] / row["timers"][key]["max_seconds"])
        comparisons[name]["solids"][family] = comparison

result = dict(source_base=request["source_base"], image_sha256=request["image_sha256"],
              rank_mapping=common_mapping, comparisons=comparisons, variants=variants,
              accepted=True)
(results / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(dict(accepted=True, comparisons=comparisons), indent=2))
