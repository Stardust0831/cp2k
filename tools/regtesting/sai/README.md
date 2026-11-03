# SAI four-GPU cuSOLVERMp tester

This fork's Action uses a small relay VPS to test the complete CP2K image
published by the existing daily SAI build. The workflow checks hourly for a new
accepted `development/16V100` image. An unchanged image and controller version
allocate no GPUs; `workflow_dispatch` with `force=true` requests a retest.
Scheduled workflows must be present on the fork's default branch.

The relay runs on `DO2-cp2k`. A dedicated SSH key can invoke only `inspect`,
`latest`, `start`, `status`, and `fetch` on SAI. The server selects the image,
inputs, resources, and job script. The client cannot upload a script or request
a shell. An Ubuntu Apptainer submission container exposes the Slurm clients and
the tester's private state directory. Slurm then runs the regressions inside the
full published CP2K SIF on four V100 GPUs, with four MPI ranks and four OpenMP
threads per rank. There is no explicit Slurm memory request.

The existing controller is reused from
[`AI4SAI/sai-hpc-software@f6eab22`](https://github.com/AI4SAI/sai-hpc-software/tree/f6eab22f4580968afa74dc11fe28033af7c59cf6/controller).
It selects nine upstream cases from `regtest-cusolver`, `regtest-elpa`,
`regtest-elpa-2`, and `regtest-elpa-qr` at the exact image source revision.
ELPA-specific settings are adapted for cuSOLVERMp; numerical references and
tolerances are retained. Acceptance requires numerical checks, six positive
cuSOLVERMp timer paths, four distinct rank/GPU assignments, the full image's
feature flags, and successful final Slurm accounting.

## Results and recovery

Each Action uploads the original and adapted inputs, their diff, raw output,
rank/GPU records, image/source identity, and Slurm accounting. Artifacts are kept
for 90 days on GitHub; raw runs also remain under the SAI state directory.
A separate GitHub-hosted job updates `sai-cusolvermp-reports`, including failed
tests and infrastructure errors. The VPS runner has only a read-only GitHub
token. The tested CP2K commit, not the workflow commit, appears in the report.

* [Public dashboard report](https://raw.githubusercontent.com/Stardust0831/cp2k/sai-cusolvermp-reports/report.txt)
* [Actions](https://github.com/Stardust0831/cp2k/actions/workflows/sai-cusolvermp.yml)

An interrupted Action resumes the recorded allocation instead of submitting a
second job. Submission uses a held job so a lost response can be recovered by
its unique name and comment before releasing it. Ambiguous submissions stop for
inspection. A manual hold after release is preserved. The Action waits at most
45 minutes and does not cancel an allocation when its own connection expires.

If the VPS is offline before an Action starts, no fresh report can be published;
the dashboard's report-age indication remains necessary. The tester observes
accepted daily images, so build/publication failures remain the responsibility
of the existing build workflow.

## Ordinary-user deployment

No root privileges are used. SAI requires an existing Apptainer installation,
Slurm project access, and a writable `~/.ssh/authorized_keys`. The relay requires
Linux Landlock ABI 4 or newer, user systemd with seccomp support, and permission
to enable user lingering. These are checked on the actual hosts.

SAI root: `/org/org-open-sai/proj-sai-hpc-software/cp2k-cusolvermp-tester`.
Relay root: `/home/cp2k/ci`.

1. On SAI, create private `code`, `state/runtime-tests`, `scratch`, and `cache`
   directories. Export the pinned repository's `controller/` to
   `code/controller/`. Install `server.py` and `deployment/entry.py` in `code/`.
2. Download Ubuntu Base 24.04.5 from
   `https://cdimage.ubuntu.com/ubuntu-base/releases/24.04/release/ubuntu-base-24.04.5-base-amd64.tar.gz`.
   `deployment/build-image.py` verifies its SHA-256 and builds `ubuntu-submit.sif`
   using the matching Ubuntu 24.04 host Python/git runtime. Build scratch lives
   in the project directory. Keep `ubuntu-submit.sif.sha256`.
3. Create `code/version.json` with the controller commit, SHA-256 of `entry.py`
   and `server.py`, and submission-image SHA-256. Update it when deploying trusted
   controller changes. The workflow never installs or changes this code.
4. Generate a dedicated relay SSH key under `ci/ssh/`. On SAI append its public
   key with `from="RELAY_IP",restrict,command="/usr/bin/python3
   /org/org-open-sai/proj-sai-hpc-software/cp2k-cusolvermp-tester/code/entry.py"`
   as one line. Preserve existing keys and a backup. Pin the already-verified
   SAI host key in `ci/ssh/known_hosts`; configure alias `sai-cp2k-tester` in
   `ci/ssh/config` with that identity, strict host checking, and no forwarding.
5. Install the GitHub runner under `ci/runner/`, register it to this fork with
   label `sai-cp2k-vps`, and set its work directory to `ci/work/`. Registration
   uses a short-lived token; no credential belongs in the repository.
6. Install `deployment/sandbox.py` under `ci/deployment/`. Install
   `deployment/sai-cp2k-runner.service` in `~/.config/systemd/user/`. Create
   `ci/work/home` and `ci/tmp`, run `loginctl enable-linger`, then
   `systemctl --user enable --now sai-cp2k-runner.service`.

The service applies Landlock to the runner and its children: writable runner,
work, and temporary directories; read-only dedicated SSH credentials and
deployment code; no access to personal home files. TCP connections are limited
to ports 443 and 12022, and TCP listeners are denied. Directory names remain
visible because the runner checks its parent directories. The systemd seccomp policy
blocks new Unix sockets, including the user-service control socket, and prevents
mount/ptrace operations. .NET diagnostics are disabled. The runner fails closed
if Landlock is unavailable.

The allowed TCP ports are not an IP allowlist. This is an ordinary-user sandbox,
not a separate kernel or UID. On SAI, the fixed command protocol is the entry
boundary; the container alone would not make an arbitrary shell safe. Only the
fixed, trusted controller may submit jobs, and the dedicated key cannot access
the user's normal SSH shell or SFTP service.

Before enabling scheduling, verify allowed operations and rejection of shell,
SCP/SFTP, path traversal, forwarding, personal-file access, deployment writes,
Unix sockets, and disallowed TCP operations. Then run a real Action and confirm
its artifact and public report agree with Slurm accounting. Run the local
recovery checks with:

```sh
python3 -m unittest discover -s tools/regtesting/sai/tests -v
```
