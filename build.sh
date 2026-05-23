#!/bin/bash
# Build the WebAssembly bindings under src/ and drop the resulting
# .js + .wasm pair into dist/. The output is intended to be vendored into
# https://github.com/MathieuDutSik/mathieudutsik.github.io under wasm/ via
# the site repo's pull-wasm.sh.
#
# Requires:
#   * emscripten (emcc) on PATH      — `brew install emscripten`
#   * boost headers                  — `brew install boost`
#   * eigen headers                  — `brew install eigen`
#   * a checkout of polyhedral_common as a sibling of this repo's parent
#     dir, i.e. at $(dirname $(dirname $PWD))/GITmathieu/polyhedral_common
#     (override with the POLYHEDRAL_COMMON env var)

set -euo pipefail

cd "$(dirname "$0")"

# --- Paths -------------------------------------------------------------------

POLYHEDRAL_COMMON_DEFAULT="$(cd ../../GITmathieu/polyhedral_common && pwd)"
POLYHEDRAL_COMMON="${POLYHEDRAL_COMMON:-$POLYHEDRAL_COMMON_DEFAULT}"

if [ ! -d "$POLYHEDRAL_COMMON/src_copos" ]; then
    echo "polyhedral_common not found at $POLYHEDRAL_COMMON" >&2
    echo "set POLYHEDRAL_COMMON env var to its absolute path" >&2
    exit 1
fi

BCPP="$POLYHEDRAL_COMMON/basic_common_cpp"
BOOST_INC="${BOOST_INC:-/opt/homebrew/opt/boost/include}"
EIGEN_INC="${EIGEN_INC:-/opt/homebrew/include/eigen3}"

for d in "$BCPP" "$BOOST_INC" "$EIGEN_INC"; do
    if [ ! -d "$d" ]; then
        echo "Required include dir missing: $d" >&2
        exit 1
    fi
done

mkdir -p dist

# --- Compiler flags ----------------------------------------------------------

INCLUDES=(
    -I"$BCPP/src_basic"
    -I"$BCPP/src_number"
    -I"$BCPP/src_matrix"
    -I"$BCPP/src_comb"
    -I"$BCPP/src_graph"
    -I"$BCPP/sparse-map/include/tsl"
    -I"$BCPP/robin-map/include/tsl"
    -I"$BCPP/hopscotch-map/include/tsl"
    -I"$POLYHEDRAL_COMMON/src_copos"
    -I"$POLYHEDRAL_COMMON/src_poly"
    -I"$POLYHEDRAL_COMMON/src_latt"
    -I"$POLYHEDRAL_COMMON/src_isotropy"
    -I"$POLYHEDRAL_COMMON/src_group"
    -I"$POLYHEDRAL_COMMON/src_short"
    -I"$POLYHEDRAL_COMMON/src_sparse_solver"
    -I"$POLYHEDRAL_COMMON/permutalib/src"
    -I"$BOOST_INC"
    -I"$EIGEN_INC"
)

CXXFLAGS=(
    -std=c++20
    -O3
    -flto
    -Wall
    -Wno-deprecated-declarations
    -D_LIBCPP_ENABLE_CXX20_REMOVED_TYPE_TRAITS
    -DWASM_PLATFORM
)

# Emscripten-specific link flags.
#   --bind                          register embind bindings
#   MODULARIZE=1 + EXPORT_NAME      classic factory: `createCoposModule()`
#                                   returns a Promise<Module>. Works as a
#                                   plain <script> tag — no ES-module setup
#                                   required, no bundler.
#   ALLOW_MEMORY_GROWTH=1           heap can grow past the initial 16 MB
#                                   (some matrices balloon transiently).
#   DISABLE_EXCEPTION_CATCHING=0    the C++ wrapper uses try/catch to convert
#                                   all errors into the result's `error`
#                                   field, so exceptions must be caught.
#   ENVIRONMENT=web,worker,node     supported runtimes; node only used for
#                                   the local smoke test below.
LINKFLAGS=(
    --bind
    -s MODULARIZE=1
    -s EXPORT_NAME=createCoposModule
    -s ALLOW_MEMORY_GROWTH=1
    -s DISABLE_EXCEPTION_CATCHING=0
    -s ENVIRONMENT=web,worker,node
)

# --- Build -------------------------------------------------------------------

SRC=src/copos_bindings.cpp
OUT_JS=dist/copos.js
OUT_WASM=dist/copos.wasm

echo "==> Building $SRC -> $OUT_JS, $OUT_WASM"
emcc "${CXXFLAGS[@]}" "${INCLUDES[@]}" "${LINKFLAGS[@]}" "$SRC" -o "$OUT_JS"

echo "==> Artifact sizes:"
ls -lh dist/copos.js dist/copos.wasm

# --- Smoke test under node ---------------------------------------------------

echo "==> Running node smoke test"
node --no-warnings -e '
const create = require("./dist/copos.js");
create().then(m => {
    const r1 = m.testCopositivity(3, ["1","0","0","0","1","0","0","0","1"]);
    console.log("I_3 :", r1.isCopositive, JSON.stringify({nature: r1.nature, error: r1.error}));
    if (r1.error || !r1.isCopositive) process.exit(1);

    const r2 = m.testCopositivity(3, ["-1","0","0","0","1","0","0","0","1"]);
    const witness = [];
    for (let i = 0; i < r2.witness.size(); ++i) witness.push(r2.witness.get(i));
    console.log("Diag(-1,1,1) :", r2.isCopositive, JSON.stringify({nature: r2.nature, witness, error: r2.error}));
    if (r2.error || r2.isCopositive) process.exit(1);

    // Bad input — should populate r.error.
    const r3 = m.testCopositivity(3, ["1-","0","0","0","1","0","0","0","1"]);
    console.log("Bad parse :", JSON.stringify({error: r3.error}));
    if (!r3.error) { console.error("expected an error"); process.exit(1); }

    // Non-symmetric.
    const r4 = m.testCopositivity(2, ["1","1","0","1"]);
    console.log("Non-symmetric :", JSON.stringify({error: r4.error}));
    if (!r4.error) { console.error("expected an error"); process.exit(1); }

    console.log("smoke test OK");
});
'
