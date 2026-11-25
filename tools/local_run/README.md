# Local MPI/OpenMP launcher

Add the executable script's directory to PATH for this shell, then run two calculations in separate
directories. They share a same-user physical-core pool:

```shell
export PATH=/path/to/cp2k/tools/local_run:$PATH
cp2k_local.py -h                 # --help is equivalent
# Terminal A: CLI ranks, environment threads
OMP_NUM_THREADS=3 cp2k_local.py -n 2 -- cp2k.psmp -i input.inp
# Terminal B: environment ranks, CLI threads
MPI_RANKS=2 cp2k_local.py -t 3 -- cp2k.psmp -i input.inp
# An absolute path also works, without changing PATH:
/path/to/cp2k/tools/local_run/cp2k_local.py -n 2 -t 3 -- cp2k.psmp -i input.inp
```

Requirements: Linux, Python 3.9+, **Open MPI 5**, its matching `hwloc-calc`, and an OpenMP runtime
that respects the rank's affinity mask (`OMP_PLACES=threads`; tested with GNU libgomp).
`--mpirun /path/to/mpirun` selects MPI. The launcher looks for `hwloc-calc` beside it, then in PATH;
both must see the same hwloc topology. Open MPI 4 and LLVM/Intel OpenMP are untested and
unsupported.

Ranks and threads resolve independently: `-n`/`--ranks` > `MPI_RANKS` > **1**, and `-t`/`--threads`
\> `OMP_NUM_THREADS` > **1**. `MPI_RANKS` is CP2K benchmark spelling, not a universal MPI variable.
Only the selected value is checked: it must be one positive decimal integer, not an empty value or
OpenMP list such as `4,2`. A valid CLI value overrides even an invalid environment value.
Process-local MPI sizes and scheduler variables are not defaults.

The launcher reads online CPUs and complete SMT sibling groups from sysfs, intersects them with
`sched_getaffinity`, and selects one allowed hardware thread per physical core. It reserves the
whole sibling group, so two callers seeing complementary siblings cannot take the same core. **One
entire system core is left unused by default**; `--reserve-cores N` changes this, including explicit
0\. The first N machine-wide groups are shared system reservations: a concurrent caller with a
smaller reserve cannot take them; a larger reserve conflicting with a live calculation fails.
Insufficient capacity fails immediately. Plans print topology, allowed/available groups and each
rank's OS CPUs.

`hwloc-calc` converts OS CPU IDs to **logical PU indexes** for an Open MPI rankfile. The printed
native command is equivalent to:

```shell
OMP_NUM_THREADS=3 OMP_PLACES=threads OMP_PROC_BIND=close OMP_DYNAMIC=FALSE \
  mpirun --host localhost:2 -n 2 --map-by rankfile:file=rankfile:hwtcpus \
  --bind-to hwthread cp2k.psmp -i input.inp
```

The actual rankfile is an inherited anonymous file (`/proc/self/fd/N`); its contents and command are
printed for review. There is no per-rank Python wrapper: the launcher **execs mpirun**, which starts
the application directly. Open MPI binds each rank to its selected hardware threads; OpenMP binds
computing threads to singleton places within that mask. The four shown OpenMP values are set for
this run; BLAS settings are left to the user. Use single-level OpenMP with no conflicting vendor
affinity or thread-limit overrides. `--dry-run` prints a temporary reservation and releases it.

The private `/tmp/cp2k-local-UID` pool uses a short allocation mutex and kernel `flock`
reservations. `mpirun` inherits the held descriptors, so no Python monitor, PID database or cleanup
service is needed. Normal exit and **one Ctrl-C/SIGINT or SIGTERM** use MPI's cleanup and release
the locks. The application's/MPI's exit status is returned directly. **SIGKILL, repeated forced
interrupts, detached descendants and MPI crashes are outside this guarantee**: check for surviving
ranks before starting another run after such a failure. Do not remove pool files while jobs run;
persistent empty files are not active reservations. Only cooperating invocations of the same user
sharing `/tmp` coordinate; other software/users remain free to use these CPUs, including the system
reserve.

This is for one local application, outside schedulers. Multi-node jobs, CPU hotplug, different
container pools and NUMA/GPU tuning are outside scope. WSL uses guest-exposed topology, not a
promise of exclusive Windows host cores. For another MPI build, first verify inherited locks and
placement:

```shell
python3 -m unittest discover -s tools/local_run -v
mpicc -O2 -fopenmp tools/local_run/affinity_probe.c -o /tmp/affinity-probe
# Run in two directories concurrently; inspect rank-N.json and the printed plans.
OMP_NUM_THREADS=3 cp2k_local.py -n 2 -- /tmp/affinity-probe
```

The probe computes on every thread and holds for five seconds (`CP2K_PROBE_HOLD` changes this). Each
thread must have one CPU, and both calculations' physical groups must be disjoint and exclude the
system reserve. Initial native placement/lock checks use Open MPI 5.0.8 and hwloc 2.12.1.
