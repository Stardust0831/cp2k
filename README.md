# SAI multi-GPU cuSOLVERMp tester

This tester runs selected upstream CP2K regression cases inside the full daily
SAI Apptainer image on **4 V100 GPUs, with 4 MPI ranks × 4 OpenMP threads**.
Each rank uses a distinct GPU. The suite contains nine cases, including six
that require evidence of actual cuSOLVERMp execution, covering real and complex
generalized diagonalization, tail blocks, and matrix block sizes.

The original references and tolerances are retained. Adapted ELPA inputs change
only the eigensolver selection and remove ELPA-specific options. The other three
cases exercise small-matrix fallback. This is a focused backend suite.

- [Latest dashboard report](https://raw.githubusercontent.com/Stardust0831/cp2k/sai-cusolvermp-reports/report.txt)
- [GitHub Actions workflow](https://github.com/Stardust0831/cp2k/actions/workflows/sai-cusolvermp.yml)
- [Tester implementation](https://github.com/Stardust0831/cp2k/tree/feat/sai-cusolvermp-actions/tools/regtesting/sai)

## Latest completed run

| Item | Value |
| --- | --- |
| Action | [41450889239](https://github.com/Stardust0831/cp2k/actions/runs/41450889239) |
| Tested CP2K commit | [`f824184648872c2b8c1d85ed961b78e6ec3d6a14`](https://github.com/cp2k/cp2k/commit/f824184648872c2b8c1d85ed961b78e6ec3d6a14) |
| Result | 9/9 cases passed; six cuSOLVERMp paths verified |
| Slurm | Job 2132061, COMPLETED / 0:0, 2m25s |
| Allocation | 16V100, one node, four GPUs; no explicit memory request |

[Download the raw files](https://github.com/Stardust0831/cp2k/actions/runs/41450889239/artifacts/13143644116)
for original and adapted inputs, the input diff, full CP2K outputs, rank/GPU
records, image identity, and Slurm accounting. GitHub artifact downloads require
sign-in and are retained for 90 days; copies are also retained on SAI.

The report identifies the source of the tested binary, independently of the
workflow commit. The image SHA-256 is recorded in the report and raw files.

## Rollout status

The initial self-hosted run and dashboard report are validated. Daily scheduling
and automatic report refresh are not enabled yet. They must be in place before
submitting the dashboard registration upstream. This branch provides the public
report endpoint for the fork's draft registration PR.
