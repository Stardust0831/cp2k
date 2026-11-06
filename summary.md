## SAI cuSOLVERMp — 4 V100

Status: **OK**

CP2K: `bc428196533ec1afec9844324d6ed8d98daaa661`

Full daily SIF; 4 MPI ranks × 4 OpenMP threads; no explicit memory request.

Slurm `2235226`: COMPLETED, 00:02:23, `16v100n62`.

9/9 cases passed; 8 numerical checks; 6 cuSOLVERMp execution paths verified.

| Input | Backend coverage |
| --- | --- |
| `QS/regtest-cusolver/H2O-6.inp` | `fallback` |
| `QS/regtest-cusolver/Si8-generalized.inp` | `cp_fm_general_cusolver` |
| `QS/regtest-cusolver/Si8-generalized-complex.inp` | `cp_cfm_general_cusolver` |
| `QS/regtest-elpa/cplx_kp.inp` | `fallback` |
| `QS/regtest-elpa/tailblock_hang.inp` | `cp_fm_diag_cusolver` |
| `QS/regtest-elpa/H13-tailblock.inp` | `cp_fm_diag_cusolver` |
| `QS/regtest-elpa-2/H2O-elpa-redistribute.inp` | `fallback` |
| `QS/regtest-elpa-qr/si128-elpa.inp` | `cp_fm_diag_cusolver` |
| `QS/regtest-elpa-qr/si128-elpa-qr.inp` | `cp_fm_diag_cusolver` |

[Action](https://github.com/Stardust0831/cp2k/actions/runs/42152801944)

Raw artifacts include original/adapted inputs, outputs, backend timers, GPU/rank mapping, image and source identity, and final Slurm accounting.
