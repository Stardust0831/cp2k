CPU / ELPA GPU / cuSOLVERMp benchmark data for CP2K issue #6371

Original inputs and complete output logs for nine converged ENERGY_FORCE calculations: Ag Gamma (500 atoms), Ag 2x2x2 k points (192 atoms), and diamond Gamma (1000 atoms), each using CPU ScaLAPACK, ELPA GPU and cuSOLVERMp.

[Download the data archive](cp2k-elpa-cusolvermp-20261025.zip) · [Archive SHA-256](cp2k-elpa-cusolvermp-20261025.zip.sha256)

The archive includes a detailed README, result tables, exact inputs and full logs, Multiwfn generation records, CPU/GPU runtime metadata, and per-file SHA-256 checksums.

CPU: 16 physical cores of an AMD Ryzen Threadripper PRO 5995WX with default ScaLAPACK. GPU runs: the same CPU allocation plus four V100-SXM2 GPUs, 16 GiB each. All configurations use 4 MPI ranks x 4 OpenMP threads. CP2K source: cfe5e5d87aaf302953f14eb43587bef2411a8fb2.

| Case | CPU (s) | ELPA GPU (s) | cuSOLVERMp (s) | cuSOLVERMp speedup vs CPU | vs ELPA GPU |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ag, 500 atoms, Gamma | 1606.878 | 587.373 | 268.342 | 5.99x | 2.19x |
| Ag, 192 atoms, 2x2x2 k points | 823.453 | 430.174 | 161.720 | 5.09x | 2.66x |
| Diamond, 1000 atoms, Gamma | 1245.299 | 441.717 | 149.788 | 8.31x | 2.95x |

ELPA uses NVIDIA GPU kernels; cuSOLVERMp uses direct generalized diagonalization. The timings cover the complete calculation. Original files are preserved byte for byte.

Related discussion: https://github.com/cp2k/cp2k/issues/6371
