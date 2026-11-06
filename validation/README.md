# cuSOLVERMp end-to-end validation

External end-to-end validation for the two independent cuSOLVERMp changes.
The completed results and provenance are in [the report](../README.md);
[evidence.tar.gz](../evidence.tar.gz) preserves the exact executed harness and raw files.

The three independent variants are the unmodified base, GPU-context reuse only,
and reduced result transfers only. They run sequentially in one allocation on
the same four V100 GPUs, 4 MPI ranks x 4 OpenMP threads, with no explicit memory
request. Reuse the full previously validated dependency set; never silently
disable optional libraries to configure a new baseline.

## Bounded matrix

1. All existing `QS/regtest-cusolver` tests, preserving references/tolerances.
2. One libcp2k process exercises real standard, real generalized and complex
   generalized solvers with two geometry evaluations per environment. It
   creates/destroys six environments on four ranks, then two independent
   two-rank subcommunicators, frees them and recreates a four-rank environment.
   This checks repeated calls, changed matrices, communicator lifetime, cleanup,
   and recreation. Scientific matrices/workspaces are not intended to be cached.
3. Si8 subsets use 36 of 104 eigenvectors with blocks of 64 (some ranks own zero
   retained result columns); a 17-vector case uses blocks of 8 (partial block);
   real and complex full-spectrum cases request all 104 eigenvectors. H2O alone
   is below the cuSOLVER threshold and is not proof of standard GPU execution.
   The standard path is an unchanged control for the transfer PR, whose scope is
   generalized real and complex solves only.
4. Ag500 Gamma and diamond1000 Gamma each run once per build, using exact inputs,
   basis files and pseudopotentials from the accepted direct-generalized retest.
   Both must converge with unchanged step counts. Compare final energy and every
   force component as well as total/SCF/whole-eigensolver/backend timers.

Small lifecycle comparisons require energy delta <= 1e-8 Ha and max force delta
<= 1e-7 Ha/bohr, in addition to existing regression tolerances. Solids retain
the prior 1e-6 Ha / 1e-5 Ha/bohr acceptance thresholds. Timer ratios are measured
results, not required pass thresholds or pure eigensolver comparisons.

The additional lifecycle fixtures use `EPS_CHECK_DIAG 1e-12`, matching the
upstream complex generalized fixture, and must produce no orthogonality warning.
Full-spectrum fixtures use DZVP-MOLOPT-SR-GTH, retaining the same 104 basis
functions. A baseline-only diagnostic measured the Gamma overlap condition
number as 4.894e6 with DZVP-GTH and 7.651e2 with MOLOPT-SR; the original full
spectrum produced orthogonality warnings/aborts, while real and complex
MOLOPT-SR full-spectrum solves passed without warnings. No further tolerance
relaxation was used. Original upstream regression inputs and references remain
unchanged; all variants use identical new fixtures.

## Preparation contract

- Stage the chosen complete source under `source/`, with exact SHA in request.
- Supply `request.json` with the same image/dependency fields as the prior
  direct-DBCSR validation, new campaign path and `variants` in this form:
  `[{"name":"baseline","source_sha":"..."},
  {"name":"reuse","source_sha":"...","patch":"reuse.patch","patch_sha256":"..."},
  {"name":"transfers","source_sha":"...","patch":"transfers.patch","patch_sha256":"..."}]`.
  Each build gets its exact source revision in `source/REVISION`.
- The prepared fixtures are included under `inputs/`, with hashes in
  `fixtures.json`; use them directly. `prepare.py` records their original
  derivation and expects the earlier local benchmark directory if rerun.
- Copy `fypp.py` and its package from the prior source-validation bundle if the
  target container has no fypp executable. `configure.py` preserves dependency
  options from the dependency snapshot.
- Use the prior clean Slurm/module environment with campaign paths updated,
  four ranks and four physical cores per rank. Call `run_validation.py` once.
- Each variant verifies `ldd` resolves the new `build/src/libcp2k.so` before
  executing anything, then saves binary/cache/source hashes and build logs.
- A tiny external interposer counts successful cuSOLVERMp handle creation and
  destruction during the lifecycle driver only. Counts must balance for every
  rank, and the reuse variant must create fewer handles than baseline. It is
  never loaded into the timed solid benchmarks. Set `count_context_calls=false`
  in the request only if explicitly recording that this extra check is omitted.
- After completion retain `results/`, `artifacts/`, the exact patches, source
  hashes, fixtures, scheduler `sacct` and GPU accounting. Run `audit.py` again
  locally against collected evidence before reporting acceptance.

## Recheck the archived numerical results locally

No GPUs are needed for this audit. From the evidence branch root, run:

```sh
sha256sum -c evidence.tar.gz.sha256
mkdir audit-copy
tar -xzf evidence.tar.gz -C audit-copy
python3 audit-copy/audit.py
```

The archive has 326 hashed files plus `collection/manifest.json`. It contains
all three builds' original regression, lifecycle and solid outputs. Run the
verification on a separate copy because `audit.py` writes `results/audit.json`.
The published `audit-result.json` was independently regenerated this way and
its acceptance and comparison data match the archived `results/audit.log`.
The as-executed harness writes its command receipt over `results/audit.json`;
the raw audit is retained in `results/audit.log`.

For a new GPU run, use a fresh directory and export the full pinned CP2K source
into `source/` without a `.git` directory. Copy the prepared harness and inputs,
adapt `request.json`, container/dependency paths and `job.sbatch` to that run,
and retain the four-GPU/rank/core mapping. Do not copy collected `results/` or
`artifacts/` into the new run. The archive records only selected source/data
files for provenance, not a complete CP2K source tree. The SAI dependency paths
are environment-specific; another cluster must provide equivalent dependencies.

