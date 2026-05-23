# wasm_polyhedral

WebAssembly bindings exposing selected entry points of
[polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common) to
JavaScript. Built with Emscripten + embind.

The artifacts are vendored into the interactive Tools page on
<https://mathieudutsik.github.io/Tools.html>.

## Layout

```
src/                  C++ files registering EMSCRIPTEN_BINDINGS(...)
build.sh              Compile every src/*.cpp into dist/<name>.js + .wasm
dist/                 Built artifacts. Committed so consumers don't need a
                      C++ toolchain to deploy the site.
```

## Bindings

The current module exposes one combined binary (`polyhedral.{js,wasm}`)
holding every WASM entry point. Each entry registers via
`emscripten::function(...)` in `src/polyhedral_bindings.cpp`:

| JS function name                       | Polyhedral_common routine                |
|----------------------------------------|------------------------------------------|
| `testCopositivity`                     | `TestCopositivity`                       |
| `testStrictCopositivity`               | `TestStrictCopositivity`                 |
| `testCopositiveFactorization`          | `TestingAttemptStrictPositivity`         |
| `testShortestVectorsRealizability`     | `SHORT_TestRealizabilityShortestFamily`  |
| `shortestVectorsAutomorphismGroup`     | `SHORT_GetStabilizer`                    |
| `gramCanonicalForm`                    | `ComputeCanonicalForm`                   |

## Building locally

Requires Emscripten + Boost + Eigen.

```bash
brew install emscripten boost eigen
./build.sh
```

By default the script expects `polyhedral_common` to live at
`../../GITmathieu/polyhedral_common`. Override with:

```bash
POLYHEDRAL_COMMON=/path/to/polyhedral_common ./build.sh
```

The script runs a node-based smoke test at the end. If it prints
`smoke test OK`, the artifacts in `dist/` are ready to be vendored.

## Consuming from the website

In the static-site repo (`mathieudutsik.github.io`):

```bash
./pull-wasm.sh         # copies dist/polyhedral.{js,wasm} into wasm/
```

then commit the result. The page loads `/wasm/polyhedral.js` and calls:

```js
const m = await createPolyhedralModule();
const r = m.testCopositivity(3, ["1","0","0","0","1","0","0","0","1"]);
// r.isCopositive : bool
// r.nature       : string
// r.witness      : VectorString  (use r.witness.size() / r.witness.get(i))
```

## Adding a new binding

1. Add the entry-point function inside `src/polyhedral_bindings.cpp` and
   register it in the `EMSCRIPTEN_BINDINGS(polyhedral_module)` block.
2. Extend the node smoke test at the bottom of `build.sh` to exercise it.
3. Bump the table above and (if the binding needs a new include path)
   the `INCLUDES` array in `build.sh`.
