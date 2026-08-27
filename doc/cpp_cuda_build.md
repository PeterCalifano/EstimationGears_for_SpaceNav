# C++ and CUDA Build Guide

## Native build contract

The native library uses C++20, CMake 3.15 or newer, and Eigen 3.4. CUDA is optional and requires toolkit 12.6 or newer. Build directories must be outside the source tree.

```bash
cmake -S . -B /tmp/estimation-gears-cpu -GNinja \
  -DEstimationGears_for_SpaceNav_ENABLE_CUDA=OFF \
  -DCPU_ENABLE_NATIVE_TUNING=OFF
cmake --build /tmp/estimation-gears-cpu --parallel
ctest --test-dir /tmp/estimation-gears-cpu --output-on-failure
```

`EstimationGears_for_SpaceNav_ENABLE_CUDA` is canonical. At the top level only, `ENABLE_CUDA` is migrated as a compatibility alias.

## Main options

| Option | Default | Purpose |
|---|---:|---|
| `BUILD_SHARED_LIBS` | `ON` | Select shared or static project library |
| `ENABLE_TESTS` | `ON` | Register Catch2 and pytest tests |
| `ENABLE_PYTHON_TESTS` | `ON` | Discover `test*.py` through CTest |
| `EstimationGears_for_SpaceNav_BUILD_PROGRAMS` | `ON` | Build native programs at top level |
| `EstimationGears_for_SpaceNav_BUILD_EXAMPLES` | `ON` | Build examples at top level |
| `CPU_ENABLE_NATIVE_TUNING` | `ON` | Enable host-specific optimized CPU flags |
| `ENABLE_TBB` | `OFF` | Link oneTBB |
| `ENABLE_OPENGL` | `OFF` | Link OpenGL |
| `SANITIZE_BUILD` | `OFF` | Enable selected sanitizers |

OptiX, TensorRT, and ZeroMQ are not EstimationGears features.

## CUDA

```bash
cmake -S . -B /tmp/estimation-gears-cuda -GNinja \
  -DEstimationGears_for_SpaceNav_ENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES="89;120"
cmake --build /tmp/estimation-gears-cuda --parallel
ctest --test-dir /tmp/estimation-gears-cuda --output-on-failure
```

Set `CUDA_ARCHITECTURES` or `CMAKE_CUDA_ARCHITECTURES` explicitly for reproducible cross-machine builds. The CUDA test fixture initializes the runtime once and gates dependent tests.

## Nested and checkout behavior

The root build never recursively adds arbitrary `lib/*` directories. MATLAB dependencies remain source/runtime checkouts, and gtwrap is resolved only by the wrapper facility. Checkout updates and submodule initialization are explicit maintenance actions.

## Installation

Both build-tree and install-tree consumers use `EstimationGears_for_SpaceNav::EstimationGears_for_SpaceNav`. The installed package resolves Eigen and any enabled exported dependency through `find_dependency`.
