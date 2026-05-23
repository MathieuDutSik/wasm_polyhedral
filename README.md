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

| Binding | Entry symbol       | Source                   |
|---------|--------------------|--------------------------|
| `copos` | `testCopositivity` | `src/copos_bindings.cpp` |

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
./pull-wasm.sh         # copies dist/copos.{js,wasm} into wasm/
```

then commit the result. The page loads `/wasm/copos.js` and calls:

```js
const m = await createCoposModule();
const r = m.testCopositivity(3, ["1","0","0","0","1","0","0","0","1"]);
// r.isCopositive : bool
// r.nature       : string
// r.witness      : VectorString  (use r.witness.size() / r.witness.get(i))
```

## Adding a new binding

1. Drop a new `src/<name>_bindings.cpp` next to `copos_bindings.cpp`.
2. Extend `build.sh` to compile it (single-source build for now —
   intentional, since each binding becomes one cohesive `.wasm` module).
3. Bump the table above and add a corresponding `pull-wasm.sh` line in
   the site repo.
