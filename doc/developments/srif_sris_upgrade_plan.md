# SRIF/SRIS Upgrade Plan

## Summary

- [ ] Keep the mission-tailoring seam unchanged in `matlab/+filter_tailoring`.
- [ ] Bring SRIF to the same architecture level as the EKF sliding-window stack.
- [ ] Preserve codegen compatibility at every stage.
- [ ] Keep SRIF/SRIS native in square-root information form.
- [ ] Avoid helper sprawl and duplicate orchestration.
- [x] Track SR-UKF/SR-USKF modernization separately in `doc/developments/sr_ukf_sruskf_modernization_plan.md` so sigma-point cleanup does not blur SRIF/SRIS ownership.

## Public Contract To Preserve

- [ ] Keep these tailoring entrypoints unchanged:
  - [ ] `BuildArchitectureTemplate`
  - [ ] `BuildInputStructsTemplate`
  - [x] `ComputeMeasPredAndObsJacobian`
  - [ ] `ComputeMeasResiduals`
  - [ ] `ComputeProcessNoiseCov`
- [x] Use `ComputeMeasPredAndObsJacobian` in EKF/SRIF-style updates that need both prediction and Jacobian.
- [x] Keep UKF/SR-UKF updates prediction-only by calling `ComputeMeasPredAndObsJacobian` with prediction outputs only; they do not consume a Jacobian.
- [x] Keep measurement hooks on the EKF-style mutable/constant config interface and use only existing `BuildArchitectureTemplate` state-index fields in `strFilterConstConfig.strStatesIdx`; do not add ad-hoc observed-state fields to `strMeasModelParams`.
- [ ] Keep these shared/public runtime entrypoints outside `+filter_tailoring`:
  - [ ] `ComputeDynFcn`
  - [ ] `ComputeDynMatrix`
  - [ ] `PropagateDyn`
  - [ ] `ManageMeasLatency`

## Public SRIF/SRIS Entry Points To Add

- [ ] `SRIF_SlideWindow_step`
- [ ] `SRIF_SlideWindow_TimeUp`
- [x] `SRIF_SlideWindow_ObsUp`
- [ ] `SRIF_SlideWindow_StateManagementStep`
- [ ] `SRIS_SlideWindow_Smooth`
- [ ] Later: `SRIS_FullHistory_Smooth`

## Const-Config And Interface Changes

- [x] Add a const-config backend selector in `strFilterConstConfig` so EKF and SRIF are chosen at architecture-build time, not by runtime-polymorphic dispatch.
- [x] Keep smoothing selection separate from forward-filter backend selection.
- [ ] Keep adaptivity outside the SRIF/SRIS estimator kernels; use the standalone external adaptivity layer as the source of truth.
- [ ] Add only the SRIF runtime fields that are strictly required for native information-form execution and verification outputs.
- [x] Keep active SR-UKF estimator kernels free of built-in adaptive covariance tuning; wrapper-level adaptivity integration is planned separately.

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

- [x] Repair `ComputeGivensRotMeasUpdate_SRIF` on top of `ApplyGivensRot`.
- [x] Keep `ComputeGivensRotMeasUpdate_SRIF` as an SRIF algorithm wrapper, not the owner of the generic rotation primitive. Rename it appropriately to explain its functionality in the style of the repo. Update test cases as well. Make it codegen friendly and efficient, i.e. "architecture" choices must be coder const,
- [x] Clean the interface, validation, and output contracts.
- [x] Make covariance reconstruction explicitly auxiliary, not canonical SRIF state.

### Tests

- [x] Fix and upgrade the current `tests/matlab/srif_modules/testGivensRotUpdate.m`.
- [ ] Add assertions for:
  - [ ] state solution
  - [ ] SR information matrix/vector
  - [ ] covariance reconstruction
  - [ ] whitening on/off
  - [ ] prior/no-prior modes
  - [ ] non-diagonal cases
- [x] Add codegen compile tests.
- [x] Add generated-MEX equivalence tests in the codegen validation flow, not as always-on optional tests in the normal MATLAB suite.

### Exit Criteria

- [x] Low-level SRIF math is stable, shared-primitive-based, and test-hardened.

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

- [x] Rewrite the legacy `GivensRotEKF` wrapper before backend-selection work depends on it.
- [x] Replace the ambiguous `ui8FILTER_TYPE` interface with explicit, named wrappers for full covariance, UD covariance, and square-root covariance inputs.
- [x] Make the square-root covariance input convention explicit in the function name or signature instead of relying on TYPE2 behavior.
- [x] Preserve current behavior through tests while removing the lower-factor assumption as a hidden long-term contract.
- [x] Keep the rewritten wrappers thin over the hardened SRIF observation-update primitive.

### Tests

- [x] Add separate tests for full covariance, UD covariance, and square-root covariance wrappers.
- [x] Add tests that reject or clearly fail on mismatched square-root covariance orientation.
- [x] Add parity tests against the direct SRIF wrapper and analytical weighted least-squares references.
- [x] Add codegen compile coverage for each wrapper that remains public.

### Exit Criteria

- [x] Legacy TYPE2 behavior is no longer a hidden architecture contract for future SRIF/EKF parity work.

## Stage 6: Add Backend-Selection Plumbing

### Implementation

- [x] Extend `BuildArchitectureTemplate` and related configs so EKF vs SRIF is chosen at const-config build time.
- [x] Keep the builder contract shared across backends.
- [x] Do not add a generic multi-backend mega-pipeline.
- [x] Add only the minimum backend-selection metadata required for codegen-safe specialization.

### Tests

- [x] Add builder tests proving EKF and SRIF architecture selection produce coherent structs.
- [ ] Add codegen compile smoke tests for both backend selections.

### Exit Criteria

- [x] The repo can build configuration data for EKF and SRIF without changing tailoring code.

## Stage 7: Implement `SRIF_SlideWindow_ObsUp`

### Implementation

- [x] Build the architecture-level SRIF observation update first.
- [x] Consume the same measurement/tailoring outputs as the EKF observation update.
- [x] Verify the current implementation convention: `SRIF_SlideWindow_ObsUp` builds
  `y_lin = r_prior + H * x_prior`, `ComputeGivensRotMeasUpdate_SRIF` solves for the
  absolute posterior state, and the caller replaces `dxStatePrior` rather than adding
  a second increment.
- [ ] Migrate the observation-update contract to error-state form before expanding the
  delayed-measurement or multi-measurement surface:
  - [ ] Make the low-level SRIF measurement RHS the prior residual `r_prior`, not the
        pseudo-observation `r_prior + H * x_prior`.
  - [ ] Pass a zero correction prior into the SRIF correction solve while preserving the
        prior square-root information factor.
  - [ ] Return the solved correction/error state from the low-level SRIF kernel.
  - [ ] Let `SRIF_SlideWindow_ObsUp` apply the correction exactly once through the proper
        Euclidean add or manifold/window-pose retraction.
  - [ ] Keep `dAllObservJac` in measurement-model form `dh/dx`; if any future hook returns
        residual Jacobians, convert the sign before calling SRIF.
- [ ] Support:
  - [ ] delayed measurements
  - [ ] multiple measurement types
  - [ ] residual editing/rejection hooks
  - [ ] underweighting hooks
- [x] Keep SRIF-specific algebra in SRIF-owned modules.

### Tests

- [x] Add linear-Gaussian EKF-vs-SRIF parity tests for measurement update.
- [ ] Add an error-state SRIF regression proving the low-level RHS is `r_prior` and the
      caller applies the returned correction exactly once.
- [ ] Add a Euclidean linear-Gaussian equivalence test showing the old absolute-state
      pseudo-observation form and the new error-state residual form produce the same
      posterior state/covariance.
- [ ] Add a Jacobian sign regression that fails if a residual Jacobian is consumed as a
      measurement-model Jacobian.
- [ ] Add multi-measurement fusion tests with shared tailoring mocks.
- [ ] Add delayed-measurement tests.
- [ ] Add codegen compile and generated-MEX equivalence tests for the SRIF observation-update entrypoint.

### Exit Criteria

- [ ] SRIF has an architecture-level observation update with EKF-parity inputs and outputs.

### Current Boundary

- [x] Initial linearized shared-tailoring `SRIF_SlideWindow_ObsUp` exists and is covered by a linear-Gaussian reference test.
- [x] Current pseudo-observation implementation is internally consistent: residuals are
      `z - h(x_prior)`, observation Jacobians are `dh/dx`, the lower square-root
      measurement covariance convention is used for whitening, and no SRIF caller adds
      the returned state as a second increment.
- [ ] Next SRIF implementation priority: replace the pseudo-observation contract with
      the documented residual/error-state contract. The low-level SRIF kernel must
      consume `r_prior` and `H_error`, while `SRIF_SlideWindow_ObsUp` owns residual
      construction, local-coordinate Jacobians, and correction retraction/application.
- [ ] Full delayed-measurement, multi-measurement editing/underweighting, and generated-MEX equivalence coverage remain before this stage is complete.

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

### Design Blocker

- [ ] Do not implement this stage by reconstructing covariance as the canonical algorithm.
- [ ] A native square-root information time update with additive mapped process noise must be designed before implementation proceeds.
- [ ] The no-process-noise deterministic information transform is straightforward, but process-noise injection requires the agreed SRIF augmentation/elimination contract.

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

### Design Blocker

- [ ] EKF window augmentation and marginalization currently mutate covariance-form blocks.
- [ ] SRIF state augmentation and fixed-lag marginalization need an information-form contract before `SRIF_SlideWindow_step` can be claimed architecture-complete.

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

- [x] Record current SR-UKF status separately: `SR_UKF_ObsUp` has square-root Schmidt covariance handling for consider states, and `SR_UKF_TimeUp` preserves consider-state mean/auto-covariance during propagation. This does not complete SRIF consider-state parity.
- [ ] Match current EKF consider-state semantics first.
- [ ] Apply consider-state handling coherently in:
  - [ ] SRIF time update
  - [ ] SRIF observation update
  - [ ] sliding-window state management
- [ ] Upgrade the low-level SRIF measurement-update path to support solved-state selection explicitly:
  - [ ] accept a solved-state mask or solved-state index list
  - [ ] skip Givens rotations for columns that are not solved in the update
  - [ ] keep the wrapper responsible for folding consider-state uncertainty into the effective residual/covariance model before the low-level SRIF kernel is called
  - [ ] do not simply rotate through consider-state columns and zero the correction afterward; that can preserve the wrong algebra for Schmidt-style behavior
- [ ] Keep the default semantics conservative and parity-focused.

### Tests

- [x] SR-UKF covariance-root consider-state observation-update regression exists in `tests/matlab/sigma_points_filters_modules/testSR_UKF_ObsUp.m`.
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

- [ ] Only after parity and smoothing are stable, optionally tighten EKF/SRIF/SRIS consider handling toward a more explicit Schmidt-style contract. SR-UKF covariance-root consider handling is already tracked in `doc/developments/sr_ukf_sruskf_modernization_plan.md`.
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
