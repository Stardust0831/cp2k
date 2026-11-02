# SAI cuSOLVERMp tester

The fork-only `SAI cuSOLVERMp (4 V100)` workflow runs on a SAI self-hosted
runner as `stardust`. It submits one Slurm job to `16V100`: four V100 GPUs,
four MPI ranks and four OpenMP threads per rank, with no explicit memory request.
The GitHub runner coordinates the job on the login host; calculations and the
upstream regression harness run inside the existing full daily CP2K SIF.

The suite and container launcher are reused from the pinned
`AI4SAI/sai-hpc-software` tester commit. Nine upstream inputs cover real and
complex generalized diagonalization, tail blocks, redistribution, block sizes,
and small-matrix fallback. Six cases require positive cuSOLVERMp timer calls.
All upstream numerical references and tolerances are retained. Original and
adapted inputs, their diff, raw outputs, rank/GPU assignments, and Slurm
accounting are retained in the uploaded artifact.

By default the workflow selects the latest published development/16V100 image.
An explicit image manifest can be provided with `workflow_dispatch`.
The image is hashed and tests are read from its exact source commit, not from
the workflow branch. The adapter records this independent test result without
changing the image's publication manifest or build recipe identity.

The workflow runs on pushes to `feat/sai-cusolvermp-actions` and manual dispatch.
It does not run pull-request code. A runner with the `sai-cp2k-v100` label and
access to the SAI software project is required. For the initial run, GitHub
connectivity uses a loopback-only SSH reverse forward to the owner's local
HTTP proxy on port 61314. An always-available connection and runner are needed
before scheduling daily jobs.
