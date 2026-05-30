# Observability analysis implementation plan

## Purpose

Implement a lean, reusable observability-analysis layer for EstimationGears_for_SpaceNav. The package must not become a model library and must not duplicate general math utilities already suitable for MathCore_for_SpaceNav / MathCore_for_ComputerVision. EstimationGears should contain only:

- observability-analysis algorithms;
- light data containers for linearized systems, factor Jacobians, gauges, and analysis reports;
- C++ and MATLAB callers/wrappers that accept externally supplied dynamics, measurement, Jacobian, and gauge providers;
- examples and tests using small synthetic systems or existing repository functions.

All Lie-group primitives, STL/Eigen adapters, CUDA kernels, CUDA memory abstractions, and general linear algebra utilities should live in MathCore. Dynamics, sensor simulation, ray tracing, and scenario generation should remain in SimulationGears or be supplied by user callbacks.

## Source alignment

The design is motivated by the observability-consistency issues in EKF-SLAM/VIO, contact-inertial estimation, invariant EKF formulations, and AstroSLAM-like relative-dynamics visual SLAM. The implementation should evaluate the linearized model actually used by an estimator and compare its nullspace against expected physical gauges.

Repository constraints to respect:

- EstimationGears has MATLAB as the active side and C++ as a scaffold/template. Keep MATLAB code function-based and struct-configured.
- Existing MATLAB filter templates use standardized functions such as `computeDynFcn`, `computeDynMatrix`, `computeMeasPred`, `computeObsMatrix`, `computeMeasResiduals`, and `computeProcessNoiseCov`.
- Existing MATLAB conventions use structs such as `strDynParams`, `strFilterConstConfig`, `strFilterMutabConfig`, `strMeasModelParams`, `strMeasBus`, and `strStatesIdx`.
- Existing MATLAB variable conventions use prefixes such as `d`, `str`, `ui8`, `ui16`, `ui32`, `b`, and `id`.
- Existing MATLAB test conventions use `matlab.unittest.TestCase` under `tests/matlab/`.
- C++ should use the current CMake infrastructure, Eigen, Catch2 tests, and namespace `estimation_gears`/`estimation_gears::observability` for new code.
- MathCore currently uses C++20, `#pragma once`, Eigen 3.4+, header-only public modules under `src/mathcore/`, namespace `mathcore`, and Catch2 tests named `tests/test_<module>.cpp`.

## Cross-repository responsibility split

### EstimationGears_for_SpaceNav

Owns:

- `J` / `H Phi` observability matrix assembly from supplied matrix sequences;
- factor-graph Fisher information / whitened Jacobian analysis from supplied residual Jacobian blocks;
- SVD-based numerical rank, nullity, singular-value reporting, nullspace basis extraction, condition estimates;
- expected-gauge leakage metrics and principal-angle comparisons;
- block observability / contribution scores from supplied state-block index maps;
- wrappers that call existing EstimationGears MATLAB filter-template functions;
- MATLAB provider normalization for function handles, symbolic toolbox objects, and optional CasADi entities;
- examples and tests with synthetic, minimal, or existing repository models.

Does not own:

- SO(3), SE(3), Sim(3), or generic Lie-group exp/log/adjoint implementations;
- CUDA kernels or GPU matrix abstractions;
- generic STL/Eigen adapters;
- spacecraft dynamics models, visual projection models, SPICE, shape-model, ray-tracing, or scene-generation logic;
- full VIO, SLAM, AstroSLAM, or InEKF estimators.

### MathCore_for_SpaceNav / MathCore_for_ComputerVision

Owns or should receive separate Codex tasks for:

- STL-to-Eigen adapters using spans/maps for contiguous storage;
- dense/sparse matrix utility traits usable by EstimationGears;
- SO(3) and SE(3) operations, with scalar templating where practical;
- CUDA-enabled matrix reduction or normal-matrix primitives, if GPU acceleration is required;
- CPU/GPU consistency tests for general math operations.

### SimulationGears_for_SpaceNav

Owns or supplies through callbacks:

- orbital dynamics and variational equations;
- camera, LIDAR, ray-tracing, and landmark visibility models;
- synthetic scenario generation;
- optional reference datasets for examples.

## Minimal C++ architecture

Do not create a model hierarchy. Use simple data structs and free functions.

Proposed directory:

```text
src/observability/
  CMakeLists.txt
  observability_types.h
  observability_analysis.h
  observability_analysis.cpp          # optional; keep header-only if practical
  observability_assembly.h
  observability_assembly.cpp          # optional
  observability_report.h
```

Add `src/observability` as a subdirectory from `src/CMakeLists.txt`. Keep the module independent from any real estimator or sensor model.

### C++ public API sketch

```cpp
namespace estimation_gears::observability {

struct RankOptions {
    bool bWhiten = true;
    bool bScaleColumns = true;
    double dRankTolerance = -1.0;   // negative means automatic tolerance
    double dGaugeTolerance = 1e-9;
};

struct LinearizedSequence {
    std::vector<Eigen::MatrixXd> caTransitionMatrices;   // F_k
    std::vector<Eigen::MatrixXd> caObservationMatrices;  // H_k
    std::vector<Eigen::MatrixXd> caObservationCovariances; // optional R_k
};

struct FactorJacobianSet {
    std::vector<Eigen::MatrixXd> caJacobians;            // J_i
    std::vector<Eigen::MatrixXd> caCovariances;          // optional R_i
    std::vector<std::string> caFactorLabels;
};

struct GaugeSet {
    Eigen::MatrixXd dExpectedGaugeBasis;                 // columns are expected null directions
    std::vector<std::string> caGaugeLabels;
};

struct StateBlockLayout {
    std::vector<std::string> caBlockLabels;
    std::vector<Eigen::Index> caBlockStartIndices;
    std::vector<Eigen::Index> caBlockSizes;
};

struct ObservabilityReport {
    Eigen::VectorXd dSingularValues;
    Eigen::MatrixXd dNullspaceBasis;
    Eigen::VectorXd dGaugeLeakage;
    Eigen::VectorXd dGaugePrincipalAnglesRad;
    Eigen::VectorXd dBlockScores;
    Eigen::Index iRank = 0;
    Eigen::Index iNullity = 0;
    double dConditionNumberObservable = 0.0;
    std::vector<std::string> caWarnings;
};

Eigen::MatrixXd assemble_ltv_observability_matrix(
    const LinearizedSequence& strLinearizedSequence,
    const RankOptions& strOptions);

Eigen::MatrixXd assemble_factor_whitened_jacobian(
    const FactorJacobianSet& strFactorJacobians,
    const RankOptions& strOptions);

ObservabilityReport analyze_jacobian(
    const Eigen::Ref<const Eigen::MatrixXd>& dWhitenedJacobian,
    const GaugeSet* pExpectedGauges,
    const StateBlockLayout* pStateBlockLayout,
    const RankOptions& strOptions);

ObservabilityReport analyze_ltv_sequence(
    const LinearizedSequence& strLinearizedSequence,
    const GaugeSet* pExpectedGauges,
    const StateBlockLayout* pStateBlockLayout,
    const RankOptions& strOptions);

ObservabilityReport analyze_factor_jacobians(
    const FactorJacobianSet& strFactorJacobians,
    const GaugeSet* pExpectedGauges,
    const StateBlockLayout* pStateBlockLayout,
    const RankOptions& strOptions);

} // namespace estimation_gears::observability
```

Sparse variants may be added after the dense API is validated. Do not start with a broad inheritance hierarchy.

### Algorithms in EstimationGears

1. Whiten each measurement or factor block using Cholesky/LDLT of the supplied covariance, or accept pre-whitened Jacobians.
2. Optionally column-scale the assembled Jacobian before SVD, retaining the scaling for diagnostics.
3. Compute rank from SVD of the whitened/scaled Jacobian, not from a raw Hessian rank.
4. Extract computed nullspace basis from right singular vectors associated with small singular values.
5. Compute expected-gauge leakage:

```text
leakage = norm(Jw * N_expected, Frobenius) / max(norm(Jw, Frobenius), eps)
```

6. Compute principal angles between the computed nullspace and expected gauge basis.
7. Compute block scores by column-block norms of the whitened Jacobian or by diagonal-block traces of `Jw' * Jw`.
8. Report warnings when expected gauges leak, rank tolerance is ambiguous, covariances are singular, or dimensions are inconsistent.

## MATLAB architecture

Use a small function-based module, not classes.

Proposed directory:

```text
matlab/observability_analysis/
  ComputeObservabilityReport.m
  AssembleLTVObservabilityMatrix.m
  AssembleFactorJacobianMatrix.m
  AnalyzeWhitenedJacobian.m
  ComputeGaugeLeakage.m
  ComputePrincipalAngles.m
  ComputeBlockObservabilityScores.m
  NormalizeJacobianProvider.m
  EvalJacobianProvider.m
  DetectCasadiAvailability.m
  README.md
```

Use standard EstimationGears-style structs and prefixes:

- `strObsProblem`
- `strObsConfig`
- `strJacobianProvider`
- `strGaugeSet`
- `strStateBlockLayout`
- `strObsReport`

### MATLAB top-level API

```matlab
strObsReport = ComputeObservabilityReport(strObsProblem, strObsConfig)
```

Minimal `strObsProblem` fields:

```matlab
strObsProblem.strJacobianProvider
strObsProblem.strGaugeSet.dExpectedGaugeBasis
strObsProblem.strGaugeSet.caGaugeLabels
strObsProblem.strStateBlockLayout
strObsProblem.dStateNominal
strObsProblem.dTimes
strObsProblem.strParams
strObsProblem.strMeasBus
```

Minimal `strObsConfig` fields:

```matlab
strObsConfig.eAnalysisMode          % "ltv", "factor_jacobian", "prewhitened_jacobian"
strObsConfig.bWhiten
strObsConfig.bScaleColumns
strObsConfig.dRankTolerance         % negative or [] means automatic tolerance
strObsConfig.bComputeGaugeMetrics
strObsConfig.bComputeBlockScores
strObsConfig.bUseSparseAssembly
```

### MATLAB Jacobian providers

Normalize all providers into a single evaluator signature:

```matlab
[dJ, dCov, strMeta] = EvalJacobianProvider(strJacobianProvider, ui32Id, strEvalContext)
```

Supported provider types:

1. Matrix sequence provider:
   - `strJacobianProvider.eType = "matrix_sequence"`
   - fields contain `caTransitionMatrices`, `caObservationMatrices`, `caFactorJacobians`, `caCovariances`.

2. Function-handle provider:
   - `strJacobianProvider.eType = "function_handle"`
   - required handle examples:
     ```matlab
     [dF, dQ, strMeta] = fcnTransitionJacobian(ui32Step, dStateNominal, strParams)
     [dH, dR, strMeta] = fcnObservationJacobian(ui32Step, dStateNominal, strMeasBus, strParams)
     [dJ, dR, strMeta] = fcnFactorJacobian(ui32Factor, dStateNominal, strMeasBus, strParams)
     ```

3. Symbolic toolbox provider:
   - `strJacobianProvider.eType = "symbolic"`
   - accept `sym` or `symfun` Jacobian expressions plus variable ordering fields.
   - convert once with `matlabFunction` and cache the resulting function handle in the provider struct when running in MATLAB.

4. Optional CasADi provider:
   - `strJacobianProvider.eType = "casadi"`
   - branch only when CasADi is detected on the MATLAB path.
   - no hard dependency and no failing tests when CasADi is unavailable.
   - use `DetectCasadiAvailability()` and skip optional tests with `assumeTrue(testCase, bHasCasadi)`.

Symbolic and CasADi branches are MATLAB-only. Core numeric functions should remain compatible with MATLAB Coder where feasible, but the provider-normalization layer can be non-codegen.

## Wrappers and callers

Create thin callers that adapt existing repository functions into the standardized provider form. Examples:

```text
matlab/observability_analysis/wrappers/
  MakeProviderFromFilterTemplate.m
  MakeProviderFromMatrixSequence.m
  MakeProviderFromFactorBlocks.m
```

`MakeProviderFromFilterTemplate` should call existing filter-template functions when supplied in config structs. It should not change EKF mechanization code.

## Stage plan and tests

### Stage 0 — Repository reconnaissance and interface lock

Deliverables:

- confirm MathCore, EstimationGears, and SimulationGears repository paths/submodule names;
- inspect current CMake, MATLAB setup, naming conventions, and test layout;
- open or draft separate MathCore TODOs for Lie-group utilities, adapters, and optional CUDA primitives;
- freeze the EstimationGears observability API names before implementation.

Tests/checks:

- no code yet;
- record exact inspected files in the implementation PR description;
- verify that the plan avoids model-specific VIO/SLAM/AstroSLAM classes in EstimationGears.

### Stage 1 — Core C++ dense algorithms

Deliverables:

- add `src/observability` module;
- implement dense LTV observability matrix assembly;
- implement dense factor-Jacobian whitening/stacking;
- implement SVD rank/nullspace diagnostics;
- implement gauge leakage, principal angles, and block scores;
- add module to CMake.

Tests:

- `tests/test_observability_rank.cpp`: diagonal matrices with known rank/nullity;
- `tests/test_observability_ltv.cpp`: two-step LTV toy system with analytically known observability matrix;
- `tests/test_observability_gauge.cpp`: synthetic Jacobian with exact nullspace and expected-gauge leakage below tolerance;
- `tests/test_observability_blocks.cpp`: block-score dimensions and monotonicity under added rows.

Acceptance:

- `./build_lib.sh` passes with tests enabled;
- no new dependency beyond Eigen and Catch2;
- all algorithm tests use fixed numeric fixtures, not random-only checks.

### Stage 2 — MathCore harmonization and sparse/STL/CUDA boundary

Deliverables in EstimationGears:

- replace any local generic math helper with MathCore calls when available;
- add compile-time seams for future sparse and CUDA backends without implementing kernels locally;
- document which functions are intentionally delegated to MathCore.

Parallel MathCore tasks for Codex, if missing:

- add `mathcore/linalg/eigen_adapters.h` for `std::span`, `std::vector`, and contiguous buffer mapping;
- add `mathcore/lie/so3.h` and `mathcore/lie/se3.h` if not already present;
- add optional CUDA primitive wrappers only in MathCore.

Tests:

- EstimationGears CPU dense tests still pass;
- if sparse support is added, compare sparse and dense assembly on the same toy sequence;
- if CUDA support is available, compare CPU/GPU normal-matrix or reduction output within tolerance, but keep this optional and guarded by `ENABLE_CUDA`.

Acceptance:

- EstimationGears contains no Lie-group implementation except calls into MathCore;
- EstimationGears contains no CUDA kernels;
- CMake works with `ENABLE_CUDA=OFF` by default.

### Stage 3 — MATLAB numeric interface

Deliverables:

- implement MATLAB core functions under `matlab/observability_analysis`;
- implement matrix-sequence and function-handle providers;
- provide one minimal example using matrices only;
- provide one wrapper example calling existing filter-template-style functions.

Tests:

- `tests/matlab/observability_analysis/testComputeObservabilityReportMatrixSequence.m`;
- `tests/matlab/observability_analysis/testComputeObservabilityReportFunctionHandle.m`;
- `tests/matlab/observability_analysis/testGaugeLeakage.m`;
- `tests/matlab/observability_analysis/testPrincipalAngles.m`;
- `tests/matlab/observability_analysis/testBlockScores.m`.

Acceptance:

- tests run with `runtests('tests/matlab/observability_analysis', 'IncludeSubfolders', true)`;
- core numeric MATLAB functions use `arguments` blocks and EstimationGears variable prefixes;
- provider-based tests avoid global state and external datasets.

### Stage 4 — MATLAB symbolic and optional CasADi providers

Deliverables:

- implement Symbolic Toolbox provider normalization using `matlabFunction`;
- implement `DetectCasadiAvailability.m`;
- implement optional CasADi provider path without introducing a hard dependency;
- document limitations of code generation for symbolic/CasADi branches.

Tests:

- `tests/matlab/observability_analysis/testSymbolicProvider.m`: simple nonlinear measurement with symbolic Jacobian compared against analytic numeric matrix;
- `tests/matlab/observability_analysis/testCasadiProviderOptional.m`: skipped unless CasADi is available on path;
- negative test for missing provider fields.

Acceptance:

- test suite passes when Symbolic Toolbox is present;
- CasADi test is skipped, not failed, when CasADi is absent;
- no CasADi import or path manipulation happens during package load.

### Stage 5 — Canonical observability examples without model-library expansion

Deliverables:

- add compact examples using supplied matrices/callbacks rather than full estimator models:
  - VIO/MSCKF-style gauge test: global translation plus yaw expected nullspace;
  - contact-inertial gauge test: absolute position plus yaw expected nullspace;
  - factor-graph/RelDyn-style test: compare visual-only and visual-plus-dynamics synthetic Jacobian conditioning.

Tests:

- `tests/matlab/observability_analysis/testCanonicalVIOGauge.m` checks expected 4-DOF gauge leakage;
- `tests/matlab/observability_analysis/testCanonicalContactGauge.m` checks expected 4-DOF gauge leakage;
- `tests/matlab/observability_analysis/testRelDynInformationIncrease.m` checks that adding a synthetic dynamics factor does not reduce rank and improves selected singular values;
- matching C++ dense tests may be added with fixed matrices saved in source as constants.

Acceptance:

- examples do not implement full VIO, SLAM, InEKF, or AstroSLAM;
- expected-gauge vectors are explicit inputs, not hidden model methods;
- RelDyn example uses toy Jacobians or callbacks and documents that true orbital dynamics should come from SimulationGears or user code.

### Stage 6 — Wrappers, gtwrap boundary, and caller prototypes

Deliverables:

- expose C++ report structs and core functions to wrappers only if useful for current MATLAB/Python workflows;
- keep MATLAB-native implementation available for function-handle, symbolic, and CasADi providers;
- add examples showing how a future estimator calls the observability analyzer with its own `F`, `H`, `J`, `R`, and gauge basis.

Tests:

- C++ wrapper build smoke test when wrapper generation is enabled;
- MATLAB wrapper smoke test with matrix-sequence provider;
- no wrapper test should require GTSAM, CasADi, or CUDA.

Acceptance:

- wrapper interfaces remain algorithmic, not estimator-specific;
- no model class hierarchy appears in generated wrapper APIs.

### Stage 7 — Final review, documentation, and gap closure

Deliverables:

- `matlab/observability_analysis/README.md` with examples, provider schema, rank tolerance, whitening, gauge leakage, and interpretation of nullspace results;
- Doxygen comments for C++ public structs/functions;
- update top-level or module documentation with cross-repo responsibility split;
- final implementation checklist in the PR description.

Review checklist:

- no duplicated Lie-group code in EstimationGears;
- no CUDA kernels in EstimationGears;
- no VIO/SLAM/AstroSLAM full model implementation in EstimationGears;
- C++ tests pass with `ENABLE_CUDA=OFF`;
- optional CUDA tests are guarded by `ENABLE_CUDA`;
- MATLAB tests pass without CasADi;
- optional CasADi tests skip when CasADi is unavailable;
- symbolic tests are either guarded or documented according to toolbox availability;
- docs state that rank is numerical and trajectory/model dependent;
- docs state that matching expected gauges is more important than rank alone;
- examples clearly mark synthetic matrices versus physical model callbacks;
- coding style matches the repository conventions.

## Codex implementation notes

1. Start with Stage 0 and write a short repository reconnaissance note in the PR description before coding.
2. Implement the dense algorithm path first. Do not start by adding abstract base classes or a model registry.
3. Use fixed small test matrices with analytically known rank/nullspace.
4. Keep all general math improvements in MathCore; reference the MathCore task/PR from the EstimationGears PR.
5. Treat MATLAB symbolic and CasADi as provider-normalization paths, not as core dependencies.
6. Keep the first PR small enough to review: Stage 1 plus Stage 3 matrix/function-handle path is the preferred first merge target.
7. Perform Stage 7 before asking for final review.
