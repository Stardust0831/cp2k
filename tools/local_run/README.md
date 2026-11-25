# Local MPI/OpenMP launcher

`cp2k_local.py` runs one Open MPI application on a Linux workstation. Two invocations by the same
user coordinate physical-core reservations, including when their allowed CPU masks contain different
SMT siblings of the same core. Each rank receives one different physical core per OpenMP thread, and each
thread uses one singleton CPU place. No root privileges, scheduler, service, or toolchain setup changes
are needed.

## Run CP2K

Requirements: Python 3.9 or newer with `os.pidfd_open` and `signal.pidfd_send_signal`, Linux 5.3 or
newer with readable `/proc` and CPU sysfs, Open MPI 4 or 5, and a matching `cp2k.psmp` build. The
OpenMP runtime must support explicit numeric `OMP_PLACES` (GNU libgomp and LLVM libomp). Use a
separate working directory for every calculation to protect output and restart files.

```shell
python3 /path/to/cp2k/tools/local_run/cp2k_local.py \
  --ranks 2 --threads 4 --reserve-cores 2 -- cp2k.psmp -i input.inp -o output.out
```

The command may be run concurrently in another terminal. On a 16-core workstation, two requests
for 2 ranks with 3 threads each fit while leaving 2 whole cores unused by both calculations. Ranks and
threads are resolved independently, so CLI and environment values can be mixed:

```shell
# Terminal A, in calculation-a/: ranks from CLI, threads from environment
OMP_NUM_THREADS=3 python3 /path/to/cp2k/tools/local_run/cp2k_local.py \
  -n 2 --reserve-cores 2 -- cp2k.psmp -i input.inp -o output.out

# Terminal B, in calculation-b/: ranks from environment, threads from CLI
MPI_RANKS=2 python3 /path/to/cp2k/tools/local_run/cp2k_local.py \
  -t 3 --reserve-cores 2 -- cp2k.psmp -i input.inp -o output.out
```

| Setting | First choice | Otherwise, if set | Default |
| --- | --- | --- | --- |
| MPI ranks | `-n` / `--ranks` | `MPI_RANKS` | 1 |
| OpenMP threads per rank | `-t` / `--threads` | `OMP_NUM_THREADS` | 1 |

`MPI_RANKS` follows the spelling used by CP2K benchmarks; MPI has no universal environment variable
for the launch size. Process-local `OMPI_COMM_WORLD_SIZE` and scheduler variables are not launch-size
defaults. With neither CLI options nor these two environment variables, the launcher uses 1 rank and
1 thread. Each selected value must be a positive decimal integer: empty, nonnumeric, zero, negative
values and nested OpenMP lists such as `OMP_NUM_THREADS=4,2` are rejected. An explicit CLI value
overrides its environment variable without validating the unused value, so
`OMP_NUM_THREADS=4,2 python3 .../cp2k_local.py -t 3 -- cp2k.psmp ...` uses three threads per rank.

Use `--mpirun /path/to/mpirun` to select an installation. `--dry-run` prints the current plan and
releases its temporary reservation; it does not promise the same CPUs will remain available later.
Plans on stderr show online CPUs, the caller's allowed CPUs, all physical sibling groups, system
reserve groups, currently available groups, the selected CPUs per rank, and the actual command.
CPU numbers are Linux OS CPU IDs, not hwloc logical indexes. Selection is deterministic first-fit;
this is not a NUMA, cache, GPU, or heterogeneous-core optimizer.

## Reservation and placement

The launcher reads `thread_siblings_list` for every online CPU under `/sys/devices/system/cpu`.
Only online CPUs in the caller's `sched_getaffinity` mask can be selected. It uses one hardware
thread per physical core, but reserves the **whole sibling group**, including siblings outside that
mask. It does not infer the number of cores from `nproc` or instantaneous CPU utilization.

`--reserve-cores N` (default **1**, explicit **0** allowed) leaves the first N machine-wide sibling
groups, in CPU-number order, out of participating calculations. These groups need not belong to the
caller's allowed mask. All concurrent invocations share these reserve locks: a second call with a
smaller N cannot consume the first call's reserve. A larger N fails if it conflicts with an already
running calculation; it never evicts it. Remaining allocation is atomic and fails immediately with
the requested/available physical-core counts if there is insufficient capacity.

The pool uses a fixed `/tmp/cp2k-local-UID` directory, owned by the user with mode 0700, and checked
0600 lock files. Shared locks protect the system reserve; exclusive locks protect calculation cores;
a short allocation lock prevents partial concurrent claims. Lock files are never unlinked, and their
existence does **not** mean a reservation is active. Changing TMPDIR does not create a second pool.

Open MPI starts locally with explicit `--bind-to none`; a per-rank wrapper then applies the selected
Linux CPU mask **before** executing CP2K and initializing OpenMP. It sets `OMP_NUM_THREADS`,
singleton `OMP_PLACES`, `OMP_PROC_BIND=close`, `OMP_DYNAMIC=FALSE`, `OMP_MAX_ACTIVE_LEVELS=1` and
`OMP_THREAD_LIMIT`. This deliberate per-run policy replaces inherited OpenMP placement settings;
`GOMP_CPU_AFFINITY`, `KMP_AFFINITY` and `KMP_HW_SUBSET` are removed from the application environment.
MKL, OpenBLAS and BLIS thread limits are set to 1 to avoid independent nested BLAS teams. An OpenMP
BLAS build can instead follow the OpenMP limits; inspect the actual linked libraries when comparing
performance. No settings in the calling shell or other programs are modified.

## Cancellation and scope

Ctrl-C, SIGTERM and SIGHUP to the launcher cancel the run and return 128 plus the signal number.
A separate supervisor holds reservations until all its MPI descendants exit, adopting orphaned
children with Linux's child-subreaper facility. After two seconds of cancellation it sends SIGKILL
to remaining descendants. A pipe also detects an abrupt launcher death, including SIGKILL, and
starts cleanup. Reservations remain held while a child cannot be reaped. Normal runs return the
MPI launcher's exit code; setup failures return 2.

Cancel the **launcher**, not its internal supervisor. Directly killing the supervisor with SIGKILL,
programs deliberately escaping process supervision, changing CPU topology during a run, and deleting
or replacing the pool directory while participants run are unsupported. Do not remove pool files to
cancel a calculation.

This is cooperation among the same user's participating launchers on the same node and shared `/tmp`
mount, not system-wide resource isolation. Other users, ordinary `mpirun` commands, and other software
can use the same cores, including the system reserve. Separate containers or mount namespaces with
different `/tmp` pools do not coordinate. Virtual machines use the topology exposed by their guest
kernel; that does not establish exclusive ownership of the host's physical CPUs.

Only one local MPI application and one level of OpenMP parallelism are supported. Scheduler jobs
and nested MPI launches are rejected. There are no arbitrary MPI placement options, MPMD, remote
hosts, MPI oversubscription, or automatic tuning. Use the scheduler/site launcher for cluster jobs.
Libraries that subsequently change affinity or thread counts themselves are outside this policy.

## Checks

Small regression tests cover SMT sibling masks, offline CPUs, atomic competing reservations, shared
system headroom, capacity rollback, pool ownership/permissions and cleanup of a detached orphan that
ignores SIGTERM, plus independent CLI/environment/default selection and invalid-value handling,
without requiring MPI:

```shell
python3 -m unittest discover -s tools/local_run -v
```

For actual MPI/OpenMP affinity, compile the included probe with the same MPI and OpenMP toolchain as
the application, then run two copies concurrently in different directories:

```shell
mpicc -O2 -fopenmp tools/local_run/affinity_probe.c -o /tmp/cp2k-affinity-probe
# Run this command in each probe directory; rank-N.json records every computing thread.
python3 /path/to/cp2k/tools/local_run/cp2k_local.py \
  -n 2 -t 3 --reserve-cores 2 -- /tmp/cp2k-affinity-probe
```

Each thread's recorded mask must be a singleton in its rank's selected CPU list. Map all recorded
CPUs back to the printed sibling groups: rank groups and the two calculations must be disjoint,
and neither may intersect `system_cores`. The probe keeps each run alive for five seconds after
computing so concurrency is observable. `CP2K_PROBE_HOLD` changes that interval for cancellation
checks. A passing probe establishes affinity, not CP2K numerical correctness or a performance gain.

The initial runtime validation used Open MPI 5.0.8 and GNU libgomp on Linux/WSL with 16 exposed
cores and two siblings per core, including concurrent runs restricted to complementary sibling
masks. Open MPI 4, LLVM libomp and bare-metal host CPU ownership have not yet been validated.
