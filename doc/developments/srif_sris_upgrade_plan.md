# SRIF/SRIS Upgrade Plan

## Summary

- [ ] Keep the mission-tailoring seam unchanged in `matlab/+filter_tailoring`.
- [ ] Bring SRIF to the same architecture level as the EKF sliding-window stack.
- [ ] Preserve codegen compatibility at every stage.
- [ ] Keep SRIF/SRIS native in square-root information form.
- [ ] Avoid helper sprawl and duplicate orchestration.

## Public Contract To Preserve

- [ ] Keep these tailoring entrypoints unchanged:
  - [ ] `BuildArchitectureTemplate`
  - [ ] `BuildInputStructsTemplate`
  - [ ] `ComputeMeasPred`
  - [ ] `ComputeObsMatrix`
  - [ ] `ComputeMeasResiduals`
  - [ ] `ComputeProcessNoiseCov`
- [ ] Keep these shared/public runtime entrypoints outside `+filter_tailoring`:
  - [ ] `ComputeDynFcn`
  - [ ] `ComputeDynMatrix`
  - [ ] `PropagateDyn`
  - [ ] `ManageMeasLatency`

## Public SRIF/SRIS Entry Points To Add

- [ ] `SRIF_SlideWindow_step`
- [ ] `SRIF_SlideWindow_TimeUp`
- [ ] `SRIF_SlideWindow_ObsUp`
- [ ] `SRIF_SlideWindow_StateManagementStep`
- [ ] `SRIS_SlideWindow_Smooth`
- [ ] Later: `SRIS_FullHistory_Smooth`

## Const-Config And Interface Changes

- [ ] Add a const-config backend selector in `strFilterConstConfig` so EKF and SRIF are chosen at architecture-build time, not by runtime-polymorphic dispatch.
- [ ] Keep smoothing selection separate from forward-filter backend selection.
- [ ] Keep adaptivity outside the SRIF/SRIS estimator kernels; use the standalone external adaptivity layer as the source of truth.
- [ ] Add only the SRIF runtime fields that are strictly required for native information-form execution and verification outputs.

## Stage 0: Land The Runtime/Template Prerequisite Commit

### Implementation

- [ ] Commit the current runtime/template/process-noise cleanup before SRIF work starts.
- [ ] Include the process-noise rename `ComputeLinearizedMappedQcov -> ComputeTrapzMappedProcessNoiseCov` in that prerequisite slice.
- [ ] Keep the current sigma-point/template/runtime changes separate from SRIF work.

### Tests

- [ ] Run the current runtime/template MATLAB test slice.
- [ ] Run the mapped process-noise regression tests after the rename.
- [ ] Run the sigma-point/template tests already added in this repo.

### Exit Criteria

- [ ] The runtime/template cleanup is committed independently and becomes the base for SRIF/SRIS work.

## Stage 1: Freeze Architecture Rules In Docs

### Implementation

- [ ] Finalize this roadmap file.
- [ ] Finalize `doc/developments/srif_sris_architecture_instructions.md`.
- [ ] Lock the ownership rule that generic Givens rotation application belongs in MathCore.
- [ ] Lock the rule that wrappers are allowed only when they improve clarity or specialization and remain efficient.

### Tests

- [ ] No code tests required.
- [ ] Review both docs against the current EKF architecture and the current SRIF low-level implementation before coding starts.

### Exit Criteria

- [ ] The implementation order, ownership rules, helper rules, and codegen requirements are documented and stable.

## Stage 2: Extract Shared Givens Rotation Application To MathCore

### Implementation

- [ ] Move generic Givens rotation application into MathCore as `ApplyGivensRot`.
- [ ] Place the primitive in the MathCore Givens rotations area, not in EstimationGears.
- [ ] Refactor EstimationGears SRIF/EKF callers to use the shared primitive.
- [ ] Keep thin SRIF/EKF wrappers only where they improve readability or specialization.

### Tests

- [ ] Add MathCore MATLAB tests for `ApplyGivensRot`.
- [ ] Add equivalence tests against current row/column Givens helpers in MathCore.
- [ ] Add codegen coverage for the shared primitive.
- [ ] Compare givens rotations against matlab built-in utilities

### Exit Criteria

- [ ] EstimationGears no longer owns the generic Givens-application primitive.

## Stage 3: Repair And Harden The Low-Level SRIF Wrapper

### Implementation

- [ ] Repair `GivensRotSRIF` on top of `ApplyGivensRot`.
- [ ] Keep `GivensRotSRIF` as an SRIF algorithm wrapper, not the owner of the generic rotation primitive. Rename it appropriately to explain its functionality in the style of the repo. Update test cases as well. Make it codegen friendly and efficient, i.e. "architecture" choices must be coder const,
- [ ] Clean the interface, validation, and output contracts.
- [ ] Make covariance reconstruction explicitly auxiliary, not canonical SRIF state.

### Tests

- [ ] Fix and upgrade the current `tests/matlab/srif_modules/testGivensRotUpdate.m`.
- [ ] Add assertions for:
  - [ ] state solution
  - [ ] SR information matrix/vector
  - [ ] covariance reconstruction
  - [ ] whitening on/off
  - [ ] prior/no-prior modes
  - [ ] non-diagonal cases
- [ ] Add codegen compile tests.
- [ ] Add generated-MEX equivalence tests in the codegen validation flow, not as always-on optional tests in the normal MATLAB suite.

### Exit Criteria

- [ ] Low-level SRIF math is stable, shared-primitive-based, and test-hardened.

## Stage 4: Reorganize Shared Filter Component Ownership

### Implementation

- [ ] Move representation-neutral measurement-editing helpers out of `ekf_modules/components` and into `sharedFiltersModules`.
- [ ] Move observation-model helpers that are not covariance-owned into `sharedFiltersModules/observationModels`.
- [ ] Move representation-neutral sliding-window pose, pose-buffer, and pose-retraction helpers into `sharedFiltersModules`.
- [ ] Keep covariance mutation, EKF update kernels, and covariance-owned sliding-window helpers in `ekf_modules` until SRIF-specific replacements are designed.
- [ ] Do not rename public functions in this stage unless the existing name is actively misleading; this stage is a location/ownership cleanup.
- [ ] Do not promote unresolved leftovers just because they look generic:
  - [ ] `EvaluateDirectionOfMotionModel` overlaps with `ComputeDirOfMotionAndJacobian`; consolidate the contract before moving it.
  - [ ] `ApplyManoeuvreDeltaV` currently mixes mean-state handling with covariance mutation; split it before sharing any representation-neutral part.
  - [ ] `UpdateGlobalQuat` should be reviewed against MathCore quaternion utilities and fixed before moving.
  - [ ] `evalJAC_RayEllipsoidIntersect` is incomplete and must not be promoted until implemented and tested.

### Tests

- [ ] Move or add tests beside the new shared ownership area.
- [ ] Run the EKF observation-update and time-update tests that exercise moved helpers.
- [ ] Run sigma-point tests that already consume shared measurement-editing policy.
- [ ] Run SRIF wrapper tests to ensure the shared path layout remains compatible.
- [ ] Move legacy `.mat` fixture-dependent tests to `tests/matlab/.deprecated/mat_fixture_dependent` before using the active test tree as validation.

### Exit Criteria

- [ ] `ekf_modules/components` no longer owns helpers that are already shared by non-EKF backends or are plainly filter-representation neutral.
- [ ] Active tests used for this stage do not depend on external `.mat` workspace fixtures.

## Stage 5: Rewrite The Legacy `GivensRotEKF` Wrapper

### Implementation

- [ ] Rewrite the legacy `GivensRotEKF` wrapper before backend-selection work depends on it.
- [ ] Replace the ambiguous `ui8FILTER_TYPE` interface with explicit, named wrappers for full covariance, UD covariance, and square-root covariance inputs.
- [ ] Make the square-root covariance input convention explicit in the function name or signature instead of relying on TYPE2 behavior.
- [ ] Preserve current behavior through tests while removing the lower-factor assumption as a hidden long-term contract.
- [ ] Keep the rewritten wrappers thin over the hardened SRIF observation-update primitive.

### Tests

- [ ] Add separate tests for full covariance, UD covariance, and square-root covariance wrappers.
- [ ] Add tests that reject or clearly fail on mismatched square-root covariance orientation.
- [ ] Add parity tests against the direct SRIF wrapper and analytical weighted least-squares references.
- [ ] Add codegen compile coverage for each wrapper that remains public.

### Exit Criteria

- [ ] Legacy TYPE2 behavior is no longer a hidden architecture contract for future SRIF/EKF parity work.

## Stage 6: Add Backend-Selection Plumbing

### Implementation

- [ ] Extend `BuildArchitectureTemplate` and related configs so EKF vs SRIF is chosen at const-config build time.
- [ ] Keep the builder contract shared across backends.
- [ ] Do not add a generic multi-backend mega-pipeline.
- [ ] Add only the minimum backend-selection metadata required for codegen-safe specialization.

### Tests

- [ ] Add builder tests proving EKF and SRIF architecture selection produce coherent structs.
- [ ] Add codegen compile smoke tests for both backend selections.

### Exit Criteria

- [ ] The repo can build configuration data for EKF and SRIF without changing tailoring code.

## Stage 7: Implement `SRIF_SlideWindow_ObsUp`

### Implementation

- [ ] Build the architecture-level SRIF observation update first.
- [ ] Consume the same measurement/tailoring outputs as the EKF observation update.
- [ ] Support:
  - [ ] delayed measurements
  - [ ] multiple measurement types
  - [ ] residual editing/rejection hooks
  - [ ] underweighting hooks
- [ ] Keep SRIF-specific algebra in SRIF-owned modules.

### Tests

- [ ] Add linear-Gaussian EKF-vs-SRIF parity tests for measurement update.
- [ ] Add multi-measurement fusion tests with shared tailoring mocks.
- [ ] Add delayed-measurement tests.
- [ ] Add codegen compile and generated-MEX equivalence tests for the SRIF observation-update entrypoint.

### Exit Criteria

- [ ] SRIF has an architecture-level observation update with EKF-parity inputs and outputs.

## Stage 8: Implement `SRIF_SlideWindow_TimeUp`

### Implementation

- [ ] Implement SRIF-native time update using the shared dynamics/runtime contracts.
- [ ] Use `ComputeTrapzMappedProcessNoiseCov` and the current mapped-noise runtime contract.
- [ ] Keep the propagation in information form.
- [ ] Expose the same orchestration-level products expected by the full step.

### Tests

- [ ] Add linear time-update parity tests against EKF on cases where both should agree.
- [ ] Add process-noise mapping tests in SRIF time-update context.
- [ ] Add codegen compile and generated-MEX equivalence tests for the SRIF time-update entrypoint.

### Exit Criteria

- [ ] SRIF has a native time-update module aligned with the current runtime.

## Stage 9: Audit And Adapt Sliding-Window State Management

### Implementation

- [ ] Audit every EKF sliding-window helper and classify it as:
  - [ ] representation-neutral and reusable
  - [ ] covariance-owned and SRIF-specific
- [ ] Implement SRIF-compatible handling for:
  - [ ] pose augmentation
  - [ ] state buffer updates
  - [ ] state ordering
  - [ ] window marginalization
- [ ] Reuse EKF helpers only where representation ownership is truly neutral.

### Tests

- [ ] Add state-order and pose-buffer tests.
- [ ] Add window augmentation and marginalization tests for SRIF.
- [ ] Add parity tests between EKF and SRIF on simple sliding-window scenarios.

### Exit Criteria

- [ ] SRIF-compatible state/window management exists and is no longer blocked by hidden covariance assumptions.

## Stage 10: Implement `SRIF_SlideWindow_step`

### Implementation

- [ ] Compose the full SRIF step only after Stages 7 to 9 are complete.
- [ ] Match the EKF step sequence:
  - [ ] adaptivity pre-step hook
  - [ ] state management
  - [ ] first time update
  - [ ] observation update
  - [ ] second time update for delayed measurements
  - [ ] manoeuvre handling
  - [ ] output packaging
- [ ] Keep the user-facing tailoring contract unchanged.

### Tests

- [ ] Add step-level smoke tests from shared builders.
- [ ] Add EKF-vs-SRIF parity tests for nominal sliding-window runs.
- [ ] Add delayed-measurement step tests.
- [ ] Add codegen compile and generated-MEX equivalence tests for `SRIF_SlideWindow_step`.

### Exit Criteria

- [ ] SRIF runs end-to-end through the same high-level workflow as EKF.

## Stage 11: Integrate External Adaptivity Layer

### Implementation

- [ ] Integrate SRIF with the standalone external adaptivity layer defined in `doc/developments/adaptivity_algorithms_modernization_plan.md`.
- [ ] Keep all adaptive algorithms, buffers, outage policy state, underweighting policy state, and covariance tuning state outside SRIF estimator kernels.
- [ ] Expose the wrapper-level SRIF outputs required by the external adaptivity layer:
  - [ ] residual/update statistics
  - [ ] innovation covariance or information-form equivalent
  - [ ] gain/update-equivalent statistics when available
  - [ ] timing context
  - [ ] measurement acceptance/rejection masks
- [ ] Apply immediate policy actions at wrapper level.
- [ ] Apply covariance tuning outputs to the next SRIF cycle.

### Tests

- [ ] Add SRIF wrapper contract tests proving it provides the statistics required by the external adaptivity layer.
- [ ] Add SRIF integration tests for immediate adaptive policies.
- [ ] Add SRIF integration tests proving covariance tuning applies to the next cycle.
- [ ] Add codegen coverage for wrapper-level adaptive integration.

### Exit Criteria

- [ ] SRIF consumes the shared external adaptivity layer without owning adaptive algorithms or adaptive buffers inside estimator kernels.

## Stage 12: Add Consider-State Parity

### Implementation

- [ ] Match current EKF consider-state semantics first.
- [ ] Apply consider-state handling coherently in:
  - [ ] SRIF time update
  - [ ] SRIF observation update
  - [ ] sliding-window state management
- [ ] Keep the default semantics conservative and parity-focused.

### Tests

- [ ] Add consider-state propagation tests.
- [ ] Add consider-state observation-update tests.
- [ ] Add outage/adaptivity interaction tests.
- [ ] Add sliding-window consider-state regression tests.

### Exit Criteria

- [ ] SRIF parity includes current EKF consider-state behavior before smoothing begins.

## Stage 13: Implement `SRIS_SlideWindow_Smooth`

### Implementation

- [ ] Implement fixed-lag smoothing over the active sliding window.
- [ ] Reuse SRIF forward products and keep smoothing native in information form.
- [ ] Respect the already-settled contracts for:
  - [ ] state ordering
  - [ ] window pose layout
  - [ ] wrapper-level adaptive outputs
  - [ ] consider-state behavior

### Tests

- [ ] Add fixed-lag linear-Gaussian smoothing reference tests.
- [ ] Add sliding-window end-to-end smoothing tests.
- [ ] Add codegen compile tests where applicable for the smoothing interface.

### Exit Criteria

- [ ] Fixed-lag SRIS works on top of the SRIF sliding-window pipeline.

## Stage 14: Implement `SRIS_FullHistory_Smooth`

### Implementation

- [ ] Add full-history smoothing only after fixed-lag SRIS is stable.
- [ ] Add the minimum history-buffer/storage machinery required.
- [ ] Reuse the same backward-information machinery from fixed-lag smoothing.

### Tests

- [ ] Add full-history consistency tests.
- [ ] Add end-to-end trajectory smoothing tests.
- [ ] Add codegen compile coverage for any public smoothing entrypoint that remains coder-compatible.

### Exit Criteria

- [ ] Full-history SRIS is added without disturbing the fixed-lag design.

## Stage 15: Optional Later Tightening Of Consider-State Semantics

### Implementation

- [ ] Only after parity and smoothing are stable, optionally tighten EKF/SRIF/SRIS consider handling toward a more explicit Schmidt-style contract.
- [ ] Keep this stage out of the critical path for the main SRIF/SRIS migration.

### Tests

- [ ] Add explicit Schmidt-style tests only if this optional stage is started.

### Exit Criteria

- [ ] Optional refinement completed without changing the default parity behavior unexpectedly.

## Global Validation Rules

- [ ] Every stage must add tests before the next stage starts.
- [ ] Every new public SRIF/SRIS entrypoint must have MATLAB tests and codegen compile coverage.
- [ ] Every codegen-capable function touched in a stage must preserve `#codegen` compatibility.
- [ ] No stage is complete if it only works in MATLAB and not in the supported coder path.

## Out Of Scope For The First SRIF/SRIS Rollout

- [ ] Do not redesign the mission-tailoring API.
- [ ] Do not introduce a generic backend-dispatch framework.
- [ ] Do not add compatibility wrappers unless a real downstream caller requires them.
- [ ] Do not fold optional Schmidt-style tightening into the main parity path.
