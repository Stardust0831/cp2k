"""Reuse the full installed dependency set without rebuilding third-party code."""
import json
import os
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parent
request = json.loads((root / "request.json").read_text())
prefix = Path(request["dependency_prefix"])
cache = {}
for line in (prefix / "share/sai/CMakeCache.txt").read_text().splitlines():
    if not line or line.startswith(("#", "//")) or ":" not in line or "=" not in line:
        continue
    key, value = line.split("=", 1)
    name, kind = key.split(":", 1)
    if kind != "INTERNAL":
        cache[name] = (kind, value)

names = {
    "CMAKE_BUILD_TYPE", "CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER", "CMAKE_Fortran_COMPILER",
    "CMAKE_C_FLAGS", "CMAKE_CXX_FLAGS", "CMAKE_Fortran_FLAGS", "CMAKE_CUDA_ARCHITECTURES",
    "CMAKE_PREFIX_PATH", "CMAKE_LIBRARY_PATH", "DBCSR_DIR", "Libint2_DIR", "libxsmm_DIR",
    "HDF5_ROOT", "HDF5_PREFER_PARALLEL", "HDF5_USE_STATIC_LIBRARIES", "ZLIB_ROOT",
    "BUILD_SHARED_LIBS",
}
names.update(k for k in cache if k.startswith(("CP2K_USE_", "CP2K_ENABLE_")))
names.update(k for k in cache if k.startswith("CP2K_") and
             (k.endswith("_ROOT") or k in ("CP2K_WITH_GPU", "CP2K_CUSOLVERMP_USE_NCCL",
                                          "CP2K_BLAS_VENDOR", "CP2K_BLAS_INTERFACE")))
names.update(k for k, (_, v) in cache.items() if k.endswith("_DIR") and
             v.startswith((str(prefix) + "/", "/opt/devtools/")))
args = ["cmake", "-S", str(root / "source"), "-B", str(root / "build")]
for name in sorted(names & cache.keys()):
    kind, value = cache[name]
    args.append(f"-D{name}:{'STRING' if kind == 'UNINITIALIZED' else kind}={value}")
args += ["-DCMAKE_INSTALL_PREFIX=" + str(root / "install"),
         "-DCMAKE_PROJECT_cp2k_INCLUDE=" + str(root / "prerequisites.cmake"),
         "-DFYPP_EXECUTABLE=" + str(root / "fypp.py")]
receipt = json.loads((prefix / "share/sai/spack-dependencies.json").read_text())
assert receipt["dependencies_installed"]
for name in ("libint", "libxsmm", "hdf5", "plumed"):
    os.environ["CP2K_" + name.upper() + "_ROOT"] = receipt["roots"][name]
(root / "results/configure-command.json").write_text(json.dumps(args, indent=2) + "\n")
subprocess.run(args, check=True)
