# Independent cuSOLVERMp optimization validation

Both patches are independently compared with the same unmodified baseline on four V100 GPUs on SAI (4 MPI ranks x 4 OpenMP threads; 16 distinct physical cores of an AMD Ryzen Threadripper PRO 5995WX). All three builds run in one allocation with identical GPU/core mapping and the full retained dependency environment. Scheduler-default memory is used.

Final Slurm job: `2234980`, `16v100n66`, COMPLETED / 0:0; 3.081 allocated GPU-hours including application builds and validation.

[Raw inputs and outputs](evidence.tar.gz) · [Archive SHA256](evidence.tar.gz.sha256) · [Numerical audit](audit-result.json) · [Validation harness](validation/README.md) · [Scheduler accounting](accounting-all.psv)

## Results

Times below are wall time from process launch to return, including MPI/container startup and cleanup. Parentheses give baseline/candidate ratios. Each large input is run once per build; these are observations for this allocation, not statistical speedup guarantees.

| Case | Matrix / requested vectors | SCF steps (all builds) | Baseline (s) | Reuse (s, ratio) | Transfers (s, ratio) |
|---|---:|---:|---:|---:|---:|
| Ag500 Gamma | 12500 / 3250 | 10 | 281.566 | 268.549 (1.048x) | 261.282 (1.078x) |
| Diamond1000 Gamma | 13000 / 2000 | 8 | 167.437 | 153.205 (1.093x) | 158.144 (1.059x) |

The transfer patch reduces the eigenvector device-to-host payload by 74.0% for Ag and 84.6% for diamond, calculated from the requested vector counts. This is not a measurement of total PCIe traffic. It still computes the full spectrum and retains the same matrix/workspace allocations; input matrix uploads are unchanged. The observed whole-process timing differences cannot be attributed entirely to reduced copies.

| Build | Successful handle creations per MPI rank | Destructions |
|---|---|---|
| Baseline | 81, 81, 81, 81 | Equal to creations on every rank |
| Execution-resource reuse | 8, 8, 8, 8 | Equal to creations on every rank |
| Requested result transfers | 81, 81, 81, 81 | Equal to creations on every rank |

## Correctness coverage

- All four original `QS/regtest-cusolver` regressions pass for every build, with their original inputs, references and tolerances.
- Eighteen energy/force evaluations per build exercise standard real, generalized real/complex, partial/full spectra, partial blocks and process columns with no requested vectors. They include atom displacement, independent two-rank subcommunicators, environment destruction and world-communicator recreation.
- Every lifecycle fixture executes the intended GPU solver. The added fixtures pass without orthogonality warnings at `EPS_CHECK_DIAG 1e-12`. Every handle creation has a matching destruction.
- Ag500 and diamond1000 ENERGY_FORCE calculations converge with identical SCF step counts. Their Multiwfn-derived inputs, basis and pseudopotential files are byte-identical across builds and unchanged from the earlier matched benchmark.

| Patch | Lifecycle max delta E (Ha) | Lifecycle max delta force (Ha/bohr) | Solids max delta E (Ha) | Solids max delta force (Ha/bohr) |
|---|---:|---:|---:|---:|
| Execution-resource reuse | 2.49e-14 | 1.03e-13 | 1.82e-12 | 1e-13 |
| Requested result transfers | 2.49e-14 | 1.21e-13 | 3.64e-12 | 1e-13 |

Acceptance limits: lifecycle delta E <= 1e-8 Ha and maximum force-component delta <= 1e-7 Ha/bohr; solids <= 1e-6 Ha and <= 1e-5 Ha/bohr. Tests use the NCCL cuSOLVERMp build; the legacy CAL runtime was not exercised.

## Provenance and reproduction

- Baseline: `5d4c07bd2d03f4106c5e72172c7a20d5e697dcdd`.
- Execution-resource reuse: `d774c4b95ffe603ba2e883521cf8e2789799bc44`.
- Requested result transfers: `00ccf6dba2862c45b315676ef1da2dd550d2d58a`.
- Container SHA256: `790326ae501f7db678a6cb5d41f4aa72cbaf02fdbe1aa3110369b68eea320c88`.
- Evidence archive SHA256: `ea73d3ed73540bf51418ebc2a0f5e8e202bd360f32dff09c6353d14e865558ee`.
- `evidence.tar.gz` contains the exact inputs, raw outputs, application build logs, CMake caches, source/binary hashes, rank/GPU/core mapping, handle counters and the complete validation harness. Binaries remain on SAI; their retained copies were hash-verified during collection.
- `collection/manifest.json` within the archive hashes every retained file. `audit.py` reproduces the numerical/path/resource checks from those raw files. The audit was independently regenerated from an extracted copy and saved as `audit-result.json`; its acceptance/comparisons match the archived `results/audit.log`. In the original harness, the audit-stage command receipt overwrote `results/audit.json`, so that archived filename contains command timing and exit status rather than numerical results. The original evidence archive is preserved unchanged.
- The harness builds each patch independently, restoring the pinned baseline sources between builds. It reuses existing third-party dependencies. The selected base includes the direct-generalized/DBCSR correctness fix; the next upstream commit only changes unrelated DFTD4/TBLITE dependency interfaces.
- To reproduce on SAI, prepare a fresh run directory with the archived harness and inputs, export the complete pinned baseline under `source/`, set `request.json` paths for the existing container/dependency environment, and regenerate matching `job.sbatch` paths. Run `run_validation.py` in the documented four-GPU allocation; do not reuse the collected `results/` or `artifacts/` directories. `validation/README.md` describes the contract; the standalone `lifecycle.c` driver and prepared fixtures can also be used with another matching libcp2k build.

## Additional full-spectrum fixture selection

Initial baseline-only attempts exposed two harness startup issues (MPI thread level and a parent directory), followed by orthogonality sensitivity in the additional full-spectrum DZVP-GTH Si8 fixture. The records were retained. A separate baseline diagnostic measured the Gamma overlap condition number as 4.894e6 for DZVP-GTH and 7.651e2 for DZVP-MOLOPT-SR-GTH; the latter keeps the same 104 basis functions and passed real and complex full-spectrum solves without warnings. Only the two added full-spectrum fixtures use MOLOPT-SR. The final 1e-12/no-warning criterion was retained; original upstream regressions and the two large benchmark inputs were not changed.
