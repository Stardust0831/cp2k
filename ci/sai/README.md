# SAI multi-GPU CI target

SAI provides a free four-V100 testing environment. The `sai-cusolvermp` target uses cp2k-ci's
existing SSH remote runner for daily and pull-request testing. An organization member can request a
run with `/cp2kci sai-cusolvermp`. The target also runs automatically for changes matching its
`trigger_path`.

The relay fetches the exact commit SHA supplied by cp2k-ci, including a PR's merge commit. SAI uses
that commit from its daily source cache when available; otherwise the relay uploads a compressed Git
archive of the build and test sources. SAI does not need outbound GitHub access. The trusted,
separately installed launcher submits one Slurm job with four MPI ranks and four OpenMP threads per
rank, assigning a distinct GPU to each rank. It does not request an explicit memory limit.

The job reuses the accepted daily container and its CMake feature settings. Headers, static
libraries, and build metadata omitted from the runtime image are restored from that image's accepted
auxiliary cache using the publisher's existing validation and relocation helper. Only CP2K is
rebuilt. The job runs the candidate's unchanged `QS/regtest-cusolver`, `QS/regtest-elpa`,
`QS/regtest-elpa-2`, and `QS/regtest-elpa-qr` suites with their original references and tolerances.
Acceptance additionally checks the built revision, enabled features, dynamic library resolution,
real and complex cuSOLVERMp calls, ELPA calls, rank/GPU mapping, and final Slurm status.

## Installation

Install a reviewed copy of this directory on both the relay and SAI. The launcher and batch scripts
are deployment code: they are never taken from the candidate checkout. Python 3.12, Git, and OpenSSH
are required. The SAI batch script uses the site's existing Apptainer, compiler, and MPI modules.

On the relay, install `runner.py` as `~/ci/sai/runner.py`, make it executable, and place this
`config.json` alongside it, substituting deployment paths:

```json
{
  "work_root": "/home/cp2k/sai-ci-state",
  "ssh_host": "sai-cp2k-build",
  "remote_command": "/org/org-open-sai/cp2k-ci/ci/sai/runner.py"
}
```

Configure `sai-cp2k-build` in the relay's SSH configuration to reach the SAI service account, with a
dedicated identity and a pinned host key. On SAI, install this directory at the `remote_command`
location and use:

```json
{
  "work_root": "/org/org-open-sai/cp2k-ci/sai-ci-runs",
  "image_catalog": "/org/org-open-sai/cp2k-ci/cp2k-cusolvermp-tester/delivery/containers/software/cp2k/development",
  "auxiliary_cache": "/org/org-open-sai/proj-sai-hpc-software/cache/cp2k-auxiliary",
  "source_repository": "/org/org-open-sai/cp2k-ci/cp2k-cusolvermp-tester/delivery/cache/repositories/cp2k"
}
```

`source_repository` is optional and must point to a trusted site cache. The launcher only reads
commit objects; it does not change that repository or reuse its working tree. The service account
also needs read access to the accepted auxiliary cache and the publisher's controller directory
recorded in the image manifest. Both are mounted read-only in the test container.

Restrict the relay key on SAI with the authorized-key options
`restrict,command="/org/org-open-sai/cp2k-ci/ci/sai/runner.py --receive"`. The forced command
accepts only the configured CI request syntax. Restrict the relay's inbound CI key with
`restrict,command="/home/cp2k/ci/sai/runner.py"`, substituting its installed path.

Finally, the upstream CI administrator configures the `sai-cp2k-ci` SSH alias, identity, and pinned
relay host key in cp2k-ci's existing `ssh-config` Secret. The relay command is
`ci/sai/runner.py BRANCH COMMIT_SHA`. Credentials and local `config.json` files are not committed.
The CI-side `cpu: 0.1` and `nodepool: default-pool` reserve only the lightweight SSH bridge; the
batch script requests the actual SAI GPU allocation.

## Reports and validation

The SSH output supplies `CommitSHA`, `Summary`, and `Status` for the normal cp2k-ci GitHub check and
report upload. Raw source, configuration, build logs, inputs, outputs, rank records, and accounting
remain in the reported `RawFiles` directory on SAI. The existing daily-image dashboard tester
remains separate. On interruption, the receiver cancels its recorded Slurm allocation.

Run the protocol checks with:

```sh
python3 -m unittest discover -s ci/sai -p 'test_*.py' -v
```

Before enabling the upstream target, invoke the deployed relay with `master` and a full CP2K commit
SHA. Verify a real build and four-GPU test, the reported SHA, solver coverage, and completed Slurm
accounting. An upstream CI invocation is a separate deployment check after its SSH credential is
configured.
