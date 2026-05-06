# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EstimationGears_for_SpaceNav is a dual-language (MATLAB + C++) library of general-purpose estimation algorithm building blocks for spacecraft navigation. It provides filter implementations (EKF, UKF, SRIF, batch least squares), shared filter infrastructure, and evaluation utilities. The MATLAB side is the actively maintained and developed part; the C++ side is a scaffold/template.

## Build Commands (C++)

The C++ part uses CMake (>=3.15) with C++17. Dependencies: Eigen3 (required), GTSAM CMake tools, Catch2 (fetched if not found), Doxygen (optional).

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

## MATLAB Setup and Tests

Run `SetupPaths_EstimationGears.m` from the repo root to set up MATLAB paths (adds `matlab/`, `simulink/`, `lib/`, `tests/` to path).

MATLAB tests live in `tests/matlab/` with subdirs per module (ekf_modules, jacobians, shared_filter_models, srif_modules, uncertainty_propagation, test_helpers). Tests use `matlab.unittest.TestCase` with `verifyEqual()`/`assertDifference()` assertions. Some tests load reference data from `.mat` files and validate against known solutions (e.g., Tapley ch5.6.4 for SRIF). MEX equivalence testing is also supported where applicable.

## Architecture

### Submodules (in `lib/`)

- **MathCore_for_SpaceNav** — core math library (algebra, Givens rotations, etc.)
- **SimulationGears_for_SpaceNav** — simulation infrastructure
- **UnitTesting4MATLAB** — MATLAB test framework (provides `assertDifference()` and other helpers in `lib/UnitTesting4MATLAB/utils`)
- **wrap** — GTSAM wrap tool for Python/MATLAB C++ bindings

### MATLAB Module Map (`matlab/`)

- **`+filter_templates_impl/`** — Standard filter interface functions (MATLAB package namespace). All filters call into this package for modularity. The 6 core interface functions are: `computeDynFcn`, `computeDynMatrix`, `computeMeasPred`, `computeObsMatrix`, `computeMeasResiduals`, `computeProcessNoiseCov`. Also includes `ManageMeasLatency` for backward propagation to measurement timestamps. Tailoring is done via config structs, not by changing filter mechanization code. This function-based approach (vs OOP) maintains MATLAB Coder/Simulink compatibility.
- **`ekf_modules/`** — EKF implementations (see detailed breakdown below)
- **`sigma_points_filters_modules/`** — Square-root UKF (`SR_UKF_kernel_SGNexe3`) and adaptive SR-USKF (`adaptiveSRUSKF_ObsUp`); common subfolder has `computeFactorProcessNoiseCov`
- **`srif_modules/`** — Square Root Information Filter via Givens rotations (`GivensRotSRIF`). Information form avoids explicit covariance inversion. Reference: Tapley 2004 ch5, Mourikis MSCKF 2007
- **`batch_least_squares_modules/`** — `WeightedLS`, `RecursiveWeightedLS`, `solveNonlinLS` (Gauss-Newton), `SolveTLS` (total LS), `Regress1DLpNorm`, `POLRLSstep`, `LOESS`
- **`sharedFiltersModules/`** — Shared building blocks (see detailed breakdown below)
- **`filters_eval_utils/`** — `computeEstimError` (additive + multiplicative/quaternion errors), `filterNMEtest` (Normalized Mean Error), `filterNEStest` (Normalized Estimation Error Squared), `evalFilterConsistency` (plots + statistics), `EvalRE` (relative error)
- **`datastructs/`** — Enum definitions (`EnumManCovModel`: MAG_DIR_THR, HERA_GNC, MAG_DIR_DIRECT, GATES)
- **`utils/`** — `GenCubeVertices` for landmark generation; `.legacy/` has visualization and statistics helpers

### EKF Modules Detail (`ekf_modules/`)

Three sub-architectures:

1. **`full-covariance-sliding/`** — Main active filter. Entry point: `EKF_SlideWindow_step` orchestrates time update -> measurement update -> state management. Subfolder `modules/` has `EKF_SlideWindow_FullCov_TimeUp` (STM-based propagation), `EKF_SlideWindow_FullCov_ObsUp` (multi-measurement fusion: LIDAR centroiding, feature tracking, AI/CRA biases), `EKF_SlideWindow_StateManagementStep` (window sliding/marginalization), `EKF_SlideWindow_AdaptivityManagementStep`. Enum `EnumMeasDelayManagementMode`: NONE, BACKWARD_PROP, ADJUST_DELTASTATE.

2. **`UD-variant/`** — UD decomposition EKF (P = U*D*U^T). Key files: `UDdecomposition`, `UDCov_TimeUp` (WMGS orthogonalization), `UDobsUpdate_ModAgeeTurner`, `UDrank1Up_AgeeTurner`/`UDrank1Up_ModAgeeTurner`, `EKF_UDcov_ObsUpDT`/`EKF_UDcov_TimeUpDT`.

3. **`components/`** — Reusable EKF building blocks for MSCKF-style ops: `AugmentStateWithNewCameraPose`, `MarginalizeSlidingWindowPose`, `UpdateFilterStateBuffers`, `UpdateStateOrdering`, `UpdateFullStateCovariance`, `UpdateGlobalQuat`, `GivensRotEKF`, `ApplyManoeuvreDeltaV`, centroiding/tracking measurement model evaluators.

### Shared Filter Modules Detail (`sharedFiltersModules/`)

Organized into subdirectories:

- **`dynamicsModels/`** — `EvalFilterDynOrbit` (main orbit dynamics evaluator), `EvalFilterDynOrbit_FixedEph`. Sub-dirs:
  - `RHSmodels/`: `evalRHS_DynLEO`, `evalRHS_DynFOGM`, `evalRHS_VariationalEqs`, `evalRHS_ContinuousTimeLinCov`, `evalRHS_RotatingFrame`
  - `RHSmodels/acceleration_components/`: `evalRHS_ExponentialAtmDrag`, `evalRHS_ZonalHarmonics20`, `CheckForEclipseMainSphereBody`, `evalAtmExpDensity`
  - `JacobianModels/`: `evalJAC_DynLEO`, `evalJAC_DynFOGM`, `evalJAC_InertialPosVelDyn`, `evalJAC_InertialMainBodyGrav`, `evalJAC_SRPwithBias`
  - `JacobianModels/jacobian_components/`: `evalJAC_AtmExpDrag`, `evalJAC_3rdBodyGrav`, `evalJAC_ZonalHarmonics20`
  - `attitude/`: `evalRHS_QuatKin`, `BuildQuatOmegaMatrix`, `ComputeAngVelFromIMU`
  - `STMmodels/`: `getDiscreteTimeSTM`
- **`observationModels/`** — `AnalyticalCoBMeasModel` (center-of-brightness with pinhole camera), `ComputeCamRelPoses`, `Pixel2LoS_NoDistorsion`, `normalizedProjectIDP`/`pinholeProjectIDP`/`pinholeProjectSymHP` (projection models), IDP<->EP transforms. Subfolder `Jacobians/` has corresponding measurement Jacobians.
- **`processNoise/`** — `GetDiscreteQforPosVelSNC` (SNC for pos/vel), `evalProcessNoiseResidualAccel`, `computeProcessNoiseCovGMresAccel` (Gauss-Markov), `evalMappedProcessNoiseFOGM`, `ComputeManoeuvreInputNoise`
- **`integratorsModules/`** — `IntegratorStepRK4` (fixed-step RK4), `IntegratorStepRK8`, `ADIntegratorStepRK45` (adaptive RK45), `ADPropagationFcn`
- **`adaptiveModules/`** — `AdaptMeasCov` (R adaptation via forgetting factor), `AdaptProcessCov` (Q adaptation), `AdaptQCovASNC`, `AdaptRQcovs` (joint R+Q)

### C++ (`src/`)

Template/scaffold using GTSAM CMake conventions. Not actively developed. Source modules: `utils/`, `template_src/`, `template_src_kernels/`, `wrapped_impl/`, `bin/`. Tests in `tests/` use Catch2 with `catch_discover_tests()`.

## Key Conventions

### MATLAB Code Patterns

- All filter functions are marked `%#codegen` for MATLAB Coder compatibility
- Memory management uses `coder.nullcopy()` for uninitialized arrays, `coder.const()` for compile-time constants, `coder.target('MATLAB')` for conditional MATLAB-only paths
- Configuration is passed via nested structs, not classes. Key struct families:
  - `strDynParams` — dynamics parameters (gravity, atmosphere, ephemerides, unmodeled acceleration stats)
  - `strFilterConstConfig` — immutable filter config (state size, index mappings, measurement vector layout); uses `coder.mustBeConst`
  - `strFilterMutabConfig` — mutable filter state (window counter, tracking mode flags, sliding mode)
  - `strMeasModelParams` — measurement model data (DCMs, camera params, pose buffers, timestamps)
  - `strMeasBus` — measurement bus (validity flags, timetags, feature keypoints, LIDAR data)
  - `strStatesIdx` — state vector index mapping (posVelIdx, unmodelAccIdx, AImeasBiasIdx, CRAmeasBiasIdx)
- State vector layout: `[pos(3) | vel(3) | unmodeled_accel(3) | AI_bias(3) | CRA_bias(3) | window_poses(3xN)]`, accessed via `strStatesIdx`
- Standard function signatures use Hungarian-ish prefixes: `d` (double), `str` (struct), `ui8`/`ui16`/`ui32` (unsigned int), `b` (boolean), `i8` (int8)
- Functions use `arguments` blocks for input validation where applicable
- Numerical noise trimming pattern: `dMatrix(abs(dMatrix) < eps) = 0`

### CMake / C++

- The CMake variable `project_name` is set externally (by `build_lib.sh` or parent CMake) and flows through `${project_name}` in CMakeLists.txt
- Build output goes to `build/` (gitignored)
- The `cmake/` directory contains reusable CMake modules (FindMKL, HandleCUDA, HandleOpenMP, etc.) adapted from a shared template
- Devcontainer support via `configure_devcontainer.sh` (Ubuntu/Debian base, optional CUDA/ROS)
