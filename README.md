# CP2K host-transfer experiments on four V100 GPUs

Both changes passed the numerical checks, but neither demonstrated a meaningful end-to-end speedup on these workloads. The cuSOLVERMp backend itself became slightly slower despite small changes in the total timer. ELPA reduced the intended Q traffic exactly, while whole-calculation times remained within the variation of repeated runs.

The cuSOLVERMp change should remain experimental on this evidence. The ELPA change removes unnecessary transfers with a small, direct change, but these results do not support a performance claim for the full calculation.

## Three paired converged runs

Hardware per calculation: SAI, 4 NVIDIA V100 GPUs, 16 physical CPU cores, 4 MPI ranks × 4 OpenMP threads, one GPU per rank. Slurm allocated 32 hardware threads and its default 60 GiB host memory; no memory size was requested. Each backend/system comparison used one allocation and an unchanged GPU/core mapping for all six runs. Different comparisons ran on different nodes. Order was baseline/candidate, candidate/baseline, baseline/candidate.

The primary metric is CP2K's final **TOTAL MAXIMUM** timer, in seconds. Speedup is baseline median divided by candidate median. Inputs are the unchanged Multiwfn-derived Ag500 and diamond1000 Γ-point ENERGY_FORCE fixtures, with the same calculation settings and default block size 64. The ELPA path includes GPU DBCSR transformations.

| Change | Case | Baseline runs 1 / 2 / 3 (s) | Candidate runs 1 / 2 / 3 (s) | Median baseline → candidate (s) | Speedup |
|---|---|---:|---:|---:|---:|
| cuSOLVERMp host pinning | Ag500 Γ | 252.432 / 246.354 / 247.244 | 246.082 / 247.086 / 246.523 | 247.244 → 246.523 | 1.0029× |
| cuSOLVERMp host pinning | Diamond1000 Γ | 134.455 / 134.510 / 133.918 | 133.826 / 134.169 / 133.964 | 134.455 → 133.964 | 1.0037× |
| ELPA requested Q | Ag500 Γ | 312.260 / 313.585 / 313.319 | 312.430 / 312.866 / 313.120 | 313.319 → 312.866 | 1.0014× |
| ELPA requested Q | Diamond1000 Γ | 217.419 / 216.531 / 218.918 | 217.950 / 218.899 / 215.929 | 217.419 → 217.950 | 0.9976× |

All 24 runs converged in the expected 10 SCF steps for Ag and 8 for diamond. All 500/1000 atom force vectors were checked; maximum within-comparison energy and force-component differences were **3.638e-12 Ha** and **1.000e-13 Ha/bohr**, below the fixed 1e-6 / 1e-5 limits.

| Change | Case | Diagonalization-stage median (s) | Backend median (s) | Whole-process median (s) |
|---|---|---:|---:|---:|
| cuSOLVERMp host pinning | Ag500 Γ | 78.337 → 78.517 | 58.803 → 59.174 | 251.508 → 251.292 |
| cuSOLVERMp host pinning | Diamond1000 Γ | 66.603 → 66.635 | 51.279 → 51.462 | 138.104 → 138.159 |
| ELPA requested Q | Ag500 Γ | 137.907 → 136.976 | 85.847 → 85.465 | 317.420 → 317.447 |
| ELPA requested Q | Diamond1000 Γ | 142.015 → 142.196 | 73.994 → 74.171 | 221.356 → 222.477 |

Backend timers have different boundaries: cuSOLVERMp includes the generalized reduction, whereas ELPA solves the standard problem. These experiments compare each backend against itself, not those two backend timers against each other.

## Changes and correctness

**cuSOLVERMp:** register the existing real/complex generalized A, B and requested Z host buffers for one solve, then unregister after the final stream synchronization. If any registration is unavailable or overlaps existing registrations, undo the registrations created by this call before using the original path. No caller-owned registration is released. This rollback is necessary: a CUDA transfer spanning a partly registered range fails with `invalid argument`. Zero local output columns are skipped. Device allocation, mathematical work, interfaces and the standard-solver path are retained.

The final source also completed an independent converged Ag500 ENERGY_FORCE run with API tracing. All **120 matrix registrations** and **120 unregisters** succeeded; the accumulated registered size was **28,250,000,000 bytes** over ten solves across four ranks. Thus the result was not caused by falling back to pageable transfers. That instrumented run is excluded from the timing table. Its energy/force differences were 3.638e-12 Ha / 1.000e-13 Ha/bohr.

**ELPA:** move the existing local requested-column calculation before the Q transfers and use it in the synchronous/asynchronous H2D and D2H copy sizes in `trans_ev_band_to_full`. Allocation, registration, interfaces and layout are retained. The operation already computes only these columns.

API interception on an uneven 193×193 matrix and a 2×2 process grid confirmed the exact requested transfer sizes, including zero-column ranks. For real `nev=16`, each direction dropped from 297,992 to 24,704 bytes; `nev=95` used 146,680 bytes. Complex `nev=95` used 293,360 instead of 595,984 bytes. A full-spectrum request retained the original sizes. Both repeated eigensolves passed analytic eigenvalue, residual and orthogonality checks. Instrumentation ran separately from performance measurements.

For the large fixtures, this Q segment alone requires 325 MB rather than 1,250 MB per direction for Ag (N=12500, nev=3250), and 208 MB rather than 1,352 MB for diamond (N=13000, nev=2000): 74.0% and 84.6% fewer bytes. These are not reductions in all transfers or in total runtime.

Additional validation:

- Four unchanged `QS/regtest-cusolver` cases passed for baseline and final logic.
- Eighteen libcp2k energy/force evaluations per variant covered repeated real/complex solves, changed geometry, partial/full spectra, split communicators and recreation. Maximum v2 comparison differences were 3.91e-14 Ha and 1.61e-13 Ha/bohr.
- CUDA registration checks exercised repeated use, pre-registered buffers, partial overlap, aliasing, zero local output columns and real H2D/D2H copies after rollback.
- ELPA passed 24 mode combinations, each with two solves: real/complex double precision, nev=16/95/193, and blocking/nonblocking collectives, for both libraries. The fixed spectral/residual/orthogonality threshold was 1e-10 with explicit finiteness checks.
- All 10 cases from `QS/regtest-elpa` and `QS/regtest-elpa-2` passed with each library, including complex k-point and redistribution cases.
- CP2K's official precommit check passed. No test reference or tolerance was changed.

## Build provenance

CP2K baseline: `178c0111101a982eae8d237ab62dc58b3811ecac`, including the earlier resource-reuse and requested-eigenpair changes, with the previously enabled dependency features including DLAF. Both versions recompiled the wrapper using the original compiler/options and relinked against the same 1502 unaffected objects. The rebuilt baseline library was byte-identical to the earlier accepted library. Loader paths and library hashes were verified for every comparison.

Performance build: `083280fced5bb24c51ccdc3f60f9b5f3f06d8c70`. Final fork source: `4df30133c6d5dd2e79ad525fc782809c4095c391`; the final change is only official C formatting, with identical C tokens (`format-only.json`). The independent full Ag500 registration trace used the compiled final source.

ELPA baseline: release 2026.02.001, git `2f09197d8342243f431025686ed59f7bbf2cc6fa`; candidate `34e538ac6baab9047ef04ede215cad359885f62d`. Both were built with GCC 13.3.0 and identical options (including GPU streams). The site's existing `cuda-aware-cpu`, `mpi-communicator` and `gpu-backtransform-host-staging` patches were retained identically. Only the band-to-full source differed; `config.h` differed only in the build timestamp. The installed library was not modified. This validation therefore describes the site's established host-staging MPI configuration.

## Evidence and reproduction

`evidence.tar.gz` includes frozen input files, all final calculation logs, atom-wise force data, timer records, rank/GPU/core mappings, GPU/host memory samples, Slurm accounting, source patches, library/source hashes and the build/test scripts. `manifest.sha256.json` covers every bundled file. Existing dependencies/images are identified by the recorded SAI paths and hashes; the harness assumes that environment rather than downloading or rebuilding its dependency stack.

The authoritative performance outputs are `results/benchmarks/`, final pinning checks are `results/v2-functional/` and `results/final-pin-trace/`, and the ELPA transfer verification is `results/q-trace/`. Older harness variants are retained for the audit trail; the accepted comparisons use `run_shard.py`, `rank-v2.py`, `container-v2.py` and `audit_shards.py`. Check the recorded sources and library paths before rerunning.

Performance jobs: Ag/cuSOLVERMp 2274443; diamond/cuSOLVERMp 2274444; Ag/ELPA 2274445; diamond/ELPA 2274446. Final-source Ag validation/registration trace: 2274754. These tests used the `cp2k-ci` service account and the 16V100 partition. The source changes and this evidence remain in Stardust0831's forks.
