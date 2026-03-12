# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EstimationGears_for_SpaceNav is a dual-language (MATLAB + C++) library of general-purpose estimation algorithm building blocks for spacecraft navigation. It provides filter implementations (EKF, UKF, SRIF, batch least squares), shared filter infrastructure, and evaluation utilities.

## Build Commands (C++)

The C++ part uses CMake (≥3.15) with C++17. Dependencies: Eigen3 (required), GTSAM CMake tools, Catch2 (fetched if not found), Doxygen (optional).

```bash
# Full configure + build + test (default: RelWithDebInfo)
./build_lib.sh

# Debug build with Ninja
./build_lib.sh -t debug -N

# Build only (skip configure)
./build_lib.sh -r

# Skip tests
./build_lib.sh --skip-tests

# Clean rebuild
./build_lib.sh --clean

# With Python/MATLAB wrappers
./build_lib.sh -p   # Python
./build_lib.sh -m   # MATLAB

# Pass extra CMake defines
./build_lib.sh -DENABLE_OMP=ON
```

Tests run automatically via CTest after build (forced on for Release). To run tests manually:
```bash
ctest --test-dir build --output-on-failure
```

The legacy VSCode helper `.vscode/cmake_run.sh` does a simpler configure+build into `build/`.

## MATLAB

Run `SetupPaths_EstimationGears.m` from the repo root to set up MATLAB paths (adds `matlab/`, `simulink/`, `lib/`, `tests/` to path).

MATLAB tests live in `tests/matlab/` with subdirs per module (ekf_modules, jacobians, shared_filter_models, srif_modules, uncertainty_propagation, test_helpers).

## Architecture

### Submodules (in `lib/`)
- **MathCore_for_SpaceNav** — core math library (algebra, Givens rotations, etc.)
- **SimulationGears_for_SpaceNav** — simulation infrastructure
- **UnitTesting4MATLAB** — MATLAB test framework

### MATLAB modules (`matlab/`)
- **`+filter_templates_impl/`** — Standard filter interface functions (package namespace). All filters use the same function signatures for modularity: `computeDynFcn`, `computeDynMatrix`, `computeMeasPred`, `computeObsMatrix`, `computeMeasResiduals`, `computeProcessNoiseCov`. Tailoring is done via config structs, not by changing filter mechanization code. This function-based approach (vs OOP) maintains MATLAB Coder/Simulink compatibility.
- **`ekf_modules/`** — EKF implementations: full-covariance sliding window, UD factorized variant (Agee-Turner), Givens rotation EKF, and component functions for MSCKF-style operations
- **`sigma_points_filters_modules/`** — Square-root UKF and adaptive SR-USKF
- **`srif_modules/`** — Square Root Information Filter via Givens rotations
- **`batch_least_squares_modules/`** — Weighted/recursive/total LS solvers, LOESS, nonlinear LS
- **`sharedFiltersModules/`** — Shared building blocks: dynamics models, observation models, process noise, integrators, adaptive modules, measurement buffering
- **`filters_eval_utils/`** — Filter consistency checks (NES, NME), Mahalanobis distance, estimation error evaluation
- **`datastructs/`** — Enum definitions for covariance models

### C++ (`src/`)
The C++ side is a template/scaffold using CMake conventions from GTSAM. Source files in `src/` are compiled into a shared library. Tests in `tests/` use Catch2 with `catch_discover_tests()`.

## Key Conventions

- The CMake variable `project_name` is set externally (by `build_lib.sh` or parent CMake) and flows through `${project_name}` in CMakeLists.txt
- Build output goes to `build/` (gitignored via absence from tracking, present in `.gitignore` patterns)
- The `cmake/` directory contains reusable CMake modules (FindMKL, HandleCUDA, HandleOpenMP, etc.) adapted from a shared template
- Devcontainer support via `configure_devcontainer.sh` (Ubuntu/Debian base, optional CUDA/ROS)
