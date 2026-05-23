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

mkdir -p dist build

# --- nauty (built once under emscripten) ------------------------------------
# SHORT_/LATT_ functionality pulls in GRAPH_traces.h which depends on
# libnauty. nauty's configure runs feature probes that crash under wasm-ld;
# we disable them with --host=wasm32-unknown-emscripten + --disable-* below.
# Mirrors the recipe in polyhedral_common/src_export_wasm/build.sh.
NAUTY_TAG="${NAUTY_TAG:-2.9.3}"
NAUTY_REPO="${NAUTY_REPO:-https://github.com/MathieuDutSik/nauty}"
NAUTY_SRC="$(pwd)/build/nauty-src"
NAUTY_INSTALL="$(pwd)/build/nauty-install"
NAUTY_LIB="$NAUTY_INSTALL/lib/libnauty.a"
NAUTY_INC="$NAUTY_INSTALL/include"

ensure_wasm_nauty() {
    if [ -f "$NAUTY_LIB" ]; then
        echo "==> nauty (wasm) already built: $NAUTY_LIB"
        return
    fi
    if [ ! -d "$NAUTY_SRC" ]; then
        echo "==> Cloning nauty $NAUTY_TAG from $NAUTY_REPO"
        git clone --depth 1 --branch "$NAUTY_TAG" "$NAUTY_REPO" "$NAUTY_SRC"
    fi
    echo "==> Building nauty under emscripten in $NAUTY_SRC"
    (
        cd "$NAUTY_SRC"
        emconfigure ./configure \
            --host=wasm32-unknown-emscripten \
            --prefix="$NAUTY_INSTALL" \
            --disable-popcnt \
            --disable-clz \
            --disable-tls
        emmake make
        emmake make install
    )
    if [ ! -f "$NAUTY_LIB" ]; then
        echo "ERROR: $NAUTY_LIB not produced by emmake install" >&2
        exit 1
    fi
}

ensure_wasm_nauty

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
    -I"$POLYHEDRAL_COMMON/src_perfect"
    -I"$POLYHEDRAL_COMMON/src_lorentzian"
    -I"$NAUTY_INC"
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
    -s EXPORT_NAME=createPolyhedralModule
    -s ALLOW_MEMORY_GROWTH=1
    -s DISABLE_EXCEPTION_CATCHING=0
    -s ENVIRONMENT=web,worker,node
)

# --- Build -------------------------------------------------------------------

SRC=src/polyhedral_bindings.cpp
OUT_JS=dist/polyhedral.js
OUT_WASM=dist/polyhedral.wasm

echo "==> Building $SRC -> $OUT_JS, $OUT_WASM"
emcc "${CXXFLAGS[@]}" "${INCLUDES[@]}" "${LINKFLAGS[@]}" "$SRC" "$NAUTY_LIB" -o "$OUT_JS"

echo "==> Artifact sizes:"
ls -lh "$OUT_JS" "$OUT_WASM"

# --- Smoke test under node ---------------------------------------------------

echo "==> Running node smoke test"
node --no-warnings -e '
const create = require("./dist/polyhedral.js");
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

    // Strict-copositivity: identity is strictly copositive.
    const r5 = m.testStrictCopositivity(3, ["1","0","0","0","1","0","0","0","1"]);
    console.log("strict I_3 :", r5.isCopositive, JSON.stringify({nature: r5.nature, error: r5.error}));
    if (r5.error || !r5.isCopositive) process.exit(1);

    // The zero matrix is copositive but NOT strictly copositive
    // (any non-zero v gives 0, which is not > 0).
    const r6 = m.testStrictCopositivity(2, ["0","0","0","0"]);
    console.log("strict 0 :", r6.isCopositive, JSON.stringify({nature: r6.nature, error: r6.error}));
    if (r6.error || r6.isCopositive) process.exit(1);

    // Copositive factorization: I_2 = 1*(1,0)(1,0)^T + 1*(0,1)(0,1)^T.
    const r7 = m.testCopositiveFactorization(2, ["1","0","0","1"]);
    console.log("fact I_2 :", JSON.stringify({hasFactorization: r7.hasFactorization,
        nBlocks: r7.nBlocks, error: r7.error}));
    if (r7.error || !r7.hasFactorization) process.exit(1);

    // diag(-1, 1) is not even copositive, so not completely positive.
    const r8 = m.testCopositiveFactorization(2, ["-1","0","0","1"]);
    console.log("fact diag(-1,1) :", JSON.stringify({hasFactorization: r8.hasFactorization,
        certificateLen: r8.certificateFlat ? r8.certificateFlat.size() : 0, error: r8.error}));
    if (r8.error || r8.hasFactorization) process.exit(1);

    // SHV realizability: the LLL kernel inside requires a SQUARE input
    // (rows == cols == n). The natural input is a basis of n vectors in
    // dimension n, here I_2 — realizable as the shortest vectors of Z^2.
    const r9 = m.testShortestVectorsRealizability(2, 2, ["1","0","0","1"]);
    console.log("shv I_2 :", JSON.stringify({realizable: r9.realizable,
        gramLen: r9.gramFlat ? r9.gramFlat.size() : 0, error: r9.error}));
    if (r9.error) process.exit(1);

    // SHV automorphism group: same input.
    const r10 = m.shortestVectorsAutomorphismGroup(2, 2, ["1","0","0","1"]);
    console.log("shv-aut I_2 :", JSON.stringify({nGenerators: r10.nGenerators,
        flatLen: r10.generatorsFlat ? r10.generatorsFlat.size() : 0, error: r10.error}));
    if (r10.error) process.exit(1);

    // Gram canonicalisation: I_2 should canonicalise to itself.
    const r11 = m.gramCanonicalForm(2, ["1","0","0","1"]);
    const canonical = [];
    if (r11.canonicalFlat) {
        for (let i = 0; i < r11.canonicalFlat.size(); ++i) canonical.push(r11.canonicalFlat.get(i));
    }
    console.log("canon I_2 :", JSON.stringify({dim: r11.dim, canonical, error: r11.error}));
    if (r11.error || r11.dim !== 2) process.exit(1);

    console.log("smoke test OK");
});
'
