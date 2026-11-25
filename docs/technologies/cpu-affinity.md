# CPU affinity

CPU affinity controls where CP2K's MPI ranks and OpenMP threads run and can affect performance. These examples use `cp2k.psmp` with Open MPI 5 on a single Linux node. Use the MPI installation that matches your CP2K build and cores available to your job.

## MPI and OpenMP runs

For two MPI ranks with one thread each, allow two physical cores:

```shell
OMP_NUM_THREADS=1 \
  mpirun -np 2 --map-by core --bind-to core --report-bindings \
  cp2k.psmp -i input.inp -o output.out
```

For two MPI ranks with three OpenMP threads each, allow six physical cores:

```shell
OMP_NUM_THREADS=3 OMP_PLACES=cores OMP_PROC_BIND=close \
  mpirun -np 2 --map-by slot:PE=3 --bind-to core --report-bindings \
  cp2k.psmp -i input.inp -o output.out
```

Open MPI does not infer the OpenMP thread count. `PE=3` assigns three processing elements to each rank; Open MPI counts physical cores by default. Keep ranks times threads per rank within the physical cores available to the job. `OMP_PLACES=cores` and `OMP_PROC_BIND=close` bind the OpenMP threads to successive core places within each rank's allocation.

## NUMA and L3 cache placement

Inspect the topology with `lstopo` (from hwloc) or `lscpu -e=CPU,CORE,SOCKET,NODE,CACHE`. NUMA domains describe memory locality; L3 cache groups identify cores sharing a cache, and their boundaries need not coincide. Keeping a rank's threads within a domain can help locality, but cache capacity, memory bandwidth, and the CP2K calculation also affect the best placement. CPU binding alone does not determine memory page placement; allocation and first-touch behavior matter too.

### Example: AMD Ryzen 9 9950X3D2 Dual Edition

For this example, assume two CCDs, each with its own L3 cache and eight physical cores, and that the job has all 16 physical cores available. Confirm this topology with `lstopo` before using the command:

```shell
OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close \
  mpirun -np 4 --map-by ppr:2:l3cache:PE=4 --bind-to core --report-bindings \
  cp2k.psmp -i input.inp -o output.out
```

`ppr:2:l3cache` places two ranks in each L3 group. Each rank gets four cores for its four threads: two ranks times four threads use eight cores per L3 group and 16 cores in total. This illustrates the specified topology and has not been benchmarked on this processor; compare timings for a representative CP2K input before choosing a placement.

## Verify the bindings

`--report-bindings` prints each rank's actual binding; Open MPI 5 also provides `--display bindings`. Check for disjoint core sets of the requested size. In the L3 example, each pair of ranks should occupy the eight cores sharing one L3 cache. Set `OMP_DISPLAY_AFFINITY=TRUE` to inspect OpenMP thread bindings as well. On SMT machines, a core place can include multiple sibling hardware threads, so a thread is bound to a physical core rather than necessarily one logical CPU.

Independent `mpirun` invocations do not coordinate their CPU assignments and can select the same cores. On shared clusters, use the scheduler's allocation and the site's launch instructions.

See the [MPI/ScaLAPACK guidance](libraries.md#mpi-and-scalapack-required-for-mpi-parallel-builds), Open MPI's [mapping and binding options](https://docs.open-mpi.org/en/v5.0.11/man-openmpi/man1/mpirun.1.html#mapping-ranking-and-binding-options), and its [memory affinity guidance](https://docs.open-mpi.org/en/v5.0.11/tuning-apps/affinity.html#memory-affinity).
