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

For a single calculation, standard Open MPI 5 also provides explicit placement. With `mpirun` and
`cp2k.psmp` in PATH, these local commands use two ranks:

```shell
# Pure MPI: one OpenMP thread per rank, one physical core per rank.
OMP_NUM_THREADS=1 OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE \
  mpirun --host localhost:2 -n 2 --map-by core --bind-to core --report-bindings \
  cp2k.psmp -i input.inp
# Hybrid: three OpenMP threads and three physical cores per rank.
OMP_NUM_THREADS=3 OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE \
  mpirun --host localhost:2 -n 2 --map-by slot:PE=3 --bind-to core --report-bindings \
  cp2k.psmp -i input.inp
```

`--bind-to core` and `OMP_PLACES=cores` allow a core's SMT siblings; they do not pin each thread to
one logical CPU. Both commands were checked with the probe below on Open MPI 5.0.8: pure MPI threads
had OS CPU masks `{0,1}` and `{2,3}`; hybrid rank 0 had `{0,1}`, `{2,3}`, `{4,5}`, and rank 1 had
`{6,7}`, `{8,9}`, `{10,11}`. `--report-bindings` reports MPI placement; the probe checks actual
computing-thread masks. Exact singleton masks require hardware-thread placement, such as the
launcher-generated rankfile below. Choosing hwthreads alone does not ensure different physical
cores. **Independent bare `mpirun` jobs do not coordinate**: two simultaneous hybrid probes used the
same six cores. Use the launcher to coordinate participating jobs of the same user.

The launcher reads online CPUs and complete SMT sibling groups from sysfs, intersects them with
`sched_getaffinity`, and selects one allowed hardware thread per physical core. It reserves the
whole sibling group, so two callers seeing complementary siblings cannot take the same core. All
free allowed physical cores are available; choose the calculation's size with ranks and threads.
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
reservations. `mpirun` inherits the held descriptors. A lock is associated with an **open file
description** and automatically releases when its last referring descriptor closes, including on
SIGKILL. Lock release does **not** prove that computing ranks have exited. Empty lock files remain
and are reused; their existence does not indicate a held lock. Do not delete pool files while jobs
or launchers may use them. Only cooperating invocations of the same user sharing `/tmp` coordinate;
other software remains free to use these CPUs.

Use one Ctrl-C for a foreground job, or save the MPI PID separately for each background job:

```shell
cp2k_local.py -n 2 -t 3 -- cp2k.psmp -i input.inp >run.log 2>&1 &
mpi_pid=$!
# Later, inspect and cancel only this job; TERM is also kill/pkill's default signal.
ps -o pid,ppid,sid,comm,args -p "$mpi_pid"
kill -TERM "$mpi_pid"
wait "$mpi_pid"
```

Exec preserves this PID, but process names can change: locally the MPI process had `comm=prterun`
and executable `prte`, while application ranks were named `affinity_probe`. Matching the launcher or
`mpirun` name can miss the process; broad name matching can kill other jobs. Prefer the saved PID
and verify its identity before signalling. TERM lets MPI attempt cleanup; its exit status is
returned directly. Record/check this job's rank PIDs too, especially after forced termination;
neither waiting for the MPI PID nor successfully reacquiring cores guarantees all ranks are gone.

Isolated 2-rank × 3-thread probes on Open MPI 5.0.8 gave these results five seconds after signalling
one exact PID (no name-based killing):

| Target and signal | MPI exit           | Remaining ranks                                | Same-core dry-run |
| ----------------- | ------------------ | ---------------------------------------------- | ----------------- |
| MPI PID, TERM     | 1                  | Two exited zombies, reaped by the test harness | Succeeded         |
| MPI PID, KILL     | Killed by signal 9 | Both alive                                     | Succeeded         |
| Rank 0, TERM      | 143                | None                                           | Succeeded         |
| Rank 0, KILL      | 137                | None                                           | Succeeded         |

The live ranks after KILL were explicitly terminated by their recorded PIDs before another
calculation was started. All four cases then ran successfully on the same cores, reusing the same
empty lock-file inodes. This is observed MPI behavior, not a promise of descendant cleanup for other
builds, applications, detached processes or repeated forced interrupts. There is no launcher process
monitor or orphan collector.

This is for one local application, outside schedulers. Multi-node jobs, CPU hotplug, different
container pools and NUMA/GPU tuning are outside scope. WSL uses guest-exposed topology, not a
promise of exclusive Windows host cores. On shared clusters, use the site's existing Slurm
allocation and MPI launch policy for CPU ownership across users; this launcher provides no
cross-user coordination and does not require Slurm locally. For another MPI build, first verify
inherited locks and placement:

```shell
python3 -m unittest discover -s tools/local_run -v
mpicc -O2 -fopenmp tools/local_run/affinity_probe.c -o /tmp/affinity-probe
# Run in two directories concurrently; inspect rank-N.json and the printed plans.
OMP_NUM_THREADS=3 cp2k_local.py -n 2 -- /tmp/affinity-probe
```

The probe computes on every thread and holds for five seconds (`CP2K_PROBE_HOLD` changes this). Each
thread launched through `cp2k_local.py` must have one CPU, and both calculations' physical groups
must be disjoint. To check the standard commands above, substitute `/tmp/affinity-probe` for
`cp2k.psmp -i input.inp` and inspect `rank-N.json` in separate working directories. Local tests use
Linux/WSL, Open MPI 5.0.8, hwloc 2.12.1 and GNU libgomp; the Conda MPI needed
`OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm OMPI_MCA_coll=^ucc` for reliable startup. Those transport
settings are specific to the test environment, not launcher defaults.
