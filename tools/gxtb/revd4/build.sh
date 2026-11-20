#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-or-later
set -euo pipefail

if [[ $# != 1 ]]; then
    echo "Usage: $0 BUILD_DIRECTORY" >&2
    exit 2
fi
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$1"
work=$(cd "$1" && pwd)
mkdir -p "$work/downloads"
cd "$work/downloads"
exec 9>"$work/build.lock"
flock 9

fetch() {
    local file=$1 url=$2 digest=$3
    if [[ -f "$file" ]] && sha256sum --check --status <<<"$digest  $file"; then
        return
    fi
    curl --fail --location --retry 3 --output "$file.part" "$url"
    sha256sum --check <<<"$digest  $file.part"
    mv -f "$file.part" "$file"
}

archive=https://raw.githubusercontent.com/DCM-Uni-Paderborn/Periodic-g-xTB-Benchmarks/866691a93b36ee582e49d896dfda5b91524dce12
stock=https://raw.githubusercontent.com/tblite/tblite/fa8a4416e8fe093d0075bc10ac875494c2a449a9/src/tblite/disp
fetch save_tblite.tar.gz "$archive/validation/accelerated_exchange/root_a_partial_bridge/final_35e7942_20260717/provider/source/save_tblite-35e7942-source.tar.gz" \
    2535519767302bc851c30512852e4ca031fadc0f97d0ece860e983effedbfd28
fetch tblite-0.8.0.tar.xz https://www.cp2k.org/static/downloads/tblite-0.8.0.tar.xz \
    2d00fd612b3fcbb69365de53577080457ceae221cd151be52e7a9977b4832ad3
fetch multicharge.tar.gz https://codeload.github.com/thfroitzheim/multicharge/tar.gz/8abf1566a13c9c343d42f27ba9ab7d63b2da2454 \
    66585e4a055568b0ecf7a0542cd2ad01982e31625da62d74a2050980f3367c18
fetch d4.f90 "$stock/d4.f90" 0a16308335e8eb62d1487c8725705ea45dd1c5c8c8eab4425656a4fff63d3f53
fetch cache.f90 "$stock/cache.f90" c55d909c7c2ed7431e9bc1edd55f1ba95093e7b6dc4b105f7fbf513ea3770412

# A different patch or source pin gets a separate, fully staged source tree.
key=$(cat "$here/build.sh" "$here/save-tblite.patch" | sha256sum | cut -c1-16)
src="$work/sources-$key"
if [[ ! -d "$src" ]]; then
    stage=$(mktemp -d "$work/staging.XXXXXXXX")
    trap 'rm -rf -- "$stage"' EXIT
    mkdir -p "$stage"/{provider,stock,multicharge}
    tar -xf save_tblite.tar.gz -C "$stage/provider"
    tar -xf tblite-0.8.0.tar.xz --strip-components=1 -C "$stage/stock"
    tar -xf multicharge.tar.gz --strip-components=1 -C "$stage/multicharge"
    cp d4.f90 cache.f90 "$stage/provider/src/tblite/disp/"
    patch --batch --fuzz=0 -d "$stage/provider" -p1 -i "$here/save-tblite.patch"
    mv "$stage" "$src"
    trap - EXIT
fi
cmake -E copy_if_different "$here/tblite_revd4.f90" "$src/provider/src/tblite/disp/revd4.f90"

cargo build --release --locked --manifest-path "$here/Cargo.toml" --target-dir "$work/cargo-target"
build="$work/cmake-$key"
cmake -S "$here" -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$work/install" -DREVD4_LIBRARY="$work/cargo-target/release/libcp2k_gxtb_revd4.so" \
    -DPROVIDER_SOURCE="$src/provider" -DSTOCK_SOURCE="$src/stock" -DMULTICHARGE_SOURCE="$src/multicharge"
cmake --build "$build" --target gxtb_save tblite-exe --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-4}"
cmake --install "$build" --component gxtb_runtime
