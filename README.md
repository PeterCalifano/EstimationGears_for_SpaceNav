# EstimationGears for Space Navigation

EstimationGears is a dual MATLAB and C++ library of reusable estimation building blocks for spacecraft navigation. MATLAB is the mature algorithm implementation; the C++20/CUDA tree provides a native library, examples, installation exports, and optional Python and MATLAB wrappers.

## Requirements

- MATLAB R2023b or newer for the MATLAB library. The current local acceptance environment is R2024b.
- CMake 3.15 or newer, a C++20 compiler, Eigen 3.4, and Python 3.12 or newer.
- CUDA 12.6 or newer only when the native CUDA feature is enabled.
- Catch2 for native tests; CMake may fetch it when explicitly permitted.

Initialize the declared source/runtime checkouts when MATLAB functionality is required:

```bash
git submodule update --init --recursive
```

Ordinary native CMake builds do not add or modify the MATLAB-only checkouts under `lib/`.

## MATLAB

Run the repository setup from MATLAB:

```matlab
SetupPaths_EstimationGears
```

MATLAB unit tests live under `tests/matlab/`. Mission-specific filter tailoring lives under `matlab/+filter_tailoring/`; shared filter APIs stay outside that package.

## Native CPU build

```bash
./build_lib.sh -N
```

The default build is `RelWithDebInfo` with tests enabled. A portable release build is:

```bash
./build_lib.sh -N -t release -D CPU_ENABLE_NATIVE_TUNING=OFF
```

The canonical CUDA option is project-qualified:

```bash
./build_lib.sh -N \
  -D EstimationGears_for_SpaceNav_ENABLE_CUDA=ON \
  -D CMAKE_CUDA_ARCHITECTURES=89
```

The historical top-level `ENABLE_CUDA` spelling remains a compatibility alias.

## Wrappers

The checkout already declares `lib/wrap`. Normal configure and build commands never move that checkout or change its gitlink.

```bash
./build_lib.sh -N -p --gtwrap-root lib/wrap
./build_lib.sh -N -m --gtwrap-root /path/to/fixed/wrap
```

Python wheels co-locate declared project-owned shared runtime libraries. The generated `_wrapper_build.py` file is checkout-only metadata and is never installed or included in a wheel.

The pinned `lib/wrap` revision supports Python generation but predates the MATLAB scalar/error-path fixes required by the wrapper regression suite. See `doc/wrappers.md` for the required MATLAB gtwrap baseline; changing the gitlink is an explicit maintenance operation.

## Installation and consumption

```bash
cmake -S . -B /tmp/estimation-gears-build -GNinja \
  -DCMAKE_INSTALL_PREFIX=/tmp/estimation-gears-install \
  -DEstimationGears_for_SpaceNav_ENABLE_CUDA=OFF
cmake --build /tmp/estimation-gears-build --target install
```

A consumer links the exported target:

```cmake
find_package(EstimationGears_for_SpaceNav REQUIRED CONFIG)
target_link_libraries(
    my_target
    PRIVATE EstimationGears_for_SpaceNav::EstimationGears_for_SpaceNav)
```

See `examples/template_consumer_project/` for a disposable installed-package consumer.

## Documentation and CI

Build Doxygen documentation locally with:

```bash
cmake --preset docs
cmake --build --preset docs
```

The Linux CPU workflow builds and tests the native library. The Pages workflow builds and uploads the documentation artifact, while its deployment job is intentionally disabled with an always-false condition.

Facility code is derived from signed `cpp_cuda_template_project` tag `v2.0.1` and then tailored to this repository. Normative development rules and the tailoring ledger are in `AGENTS.md`.
