# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Read `AGENTS.md` first. It is the authoritative source for language conventions, template tailoring, Git safety, testing ownership, and staged-review requirements. This file records repository facts and commands only.

## Project Overview

EstimationGears_for_SpaceNav is a dual-language (MATLAB + C++) library of general-purpose estimation algorithm building blocks for spacecraft navigation. It provides filter implementations (EKF, UKF, SRIF, batch least squares), shared filter infrastructure, and evaluation utilities. The MATLAB side is the actively maintained and developed part; the C++ side is a scaffold/template.

## Build Commands (C++)

The native library uses CMake >=3.15 and C++20. Eigen 3.4 is required; Catch2 and pytest own runtime tests, Doxygen is optional, and CUDA >=12.6 is optional. Native CMake configuration does not recursively add the MATLAB-only checkouts under `lib/`.

```bash
# Full configure + build + test (default: RelWithDebInfo)
./build_lib.sh -N

# Debug build with Ninja
./build_lib.sh -N -t debug

# Build only (skip configure)
./build_lib.sh -r

# Skip tests
./build_lib.sh --skip-tests

# Clean rebuild
./build_lib.sh -N --clean

# With Python/MATLAB wrappers
./build_lib.sh -N -p --gtwrap-root lib/wrap   # Python
./build_lib.sh -N -m --gtwrap-root /path/to/fixed/wrap  # MATLAB; see doc/wrappers.md

# Pass extra CMake defines
./build_lib.sh -N -D ENABLE_TBB=ON

# Optional CUDA (canonical project-qualified option)
./build_lib.sh -N -D EstimationGears_for_SpaceNav_ENABLE_CUDA=ON

# Documentation
cmake --preset docs
cmake --build --preset docs
```

Tests run automatically unless `--skip-tests` is given. To run them manually:

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

- **`+filter_tailoring/`** — Mission-specific tailoring layer. This package is intended to contain all and only the functions that a user is expected to edit manually to adapt the generic filter implementations: `BuildArchitectureTemplate`, `BuildInputStructsTemplate`, `ComputeMeasPred`, `ComputeObsMatrix`, `ComputeMeasResiduals`, and `ComputeProcessNoiseCov`. Shared/public filter entrypoints such as `ComputeDynFcn`, `ComputeDynMatrix`, `PropagateDyn`, and `ManageMeasLatency` live outside this package.
- **`ekf_modules/`** — EKF implementations (see detailed breakdown below)
- **`sigma_points_filters_modules/`** — Adaptive square-root UKF observation update (`SR_UKF_Adaptive_ObsUp`) plus the UD-form observation-update wrapper (`SRUSKF_UDcov_ObsUpDT`); common subfolder has `ComputeFactorProcessNoiseCov`
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
- **`integratorsModules/`** — `IntegratorStepRK4` (fixed-step RK4), `IntegratorStepRK8`, `ADIntegratorStepRK45` (adaptive RK45), `ADPropagationFcn`, `PropagateDyn`
- **`adaptiveModules/`** — `AdaptMeasCov` (R adaptation via forgetting factor), `AdaptProcessCov` (Q adaptation), `AdaptQCovASNC`, `AdaptRQcovs` (joint R+Q)

### C++ (`src/`)

C++20/CUDA scaffold with install/export support and optional Python/MATLAB gtwrap bindings. Source modules are `utils/`, `template_src/`, `template_src_kernels/`, `wrapped_impl/`, and `bin/`. Tests use Catch2, pytest, and wrapper-specific MATLAB regression cases.

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

### Native facilities

- CMake project/package: `EstimationGears_for_SpaceNav`; exported target: `EstimationGears_for_SpaceNav::EstimationGears_for_SpaceNav`
- Canonical CUDA option: `EstimationGears_for_SpaceNav_ENABLE_CUDA`; top-level `ENABLE_CUDA` remains a compatibility alias
- Python package: `EstimationGears_for_SpaceNav`; wrapper namespace: `estimation_gears`
- Build output defaults to `build/`; fresh acceptance builds and consumers use disposable out-of-tree directories
- `configure_devcontainer.sh` configures Ubuntu/Debian images with optional CUDA and generic ROS tooling; the project itself has no ROS overlay
- `build_linux.yml` owns portable CPU CI. `docs_pages.yml` builds documentation, while its `deploy` job is manually disabled with `if: ${{ false }}`
