# SRIF/SRIS Architecture Instructions

## Purpose

This document defines the hard implementation rules for the SRIF/SRIS upgrade so that the new
stack reaches EKF parity without unnecessary duplication, hidden architecture drift, or codegen
regressions.

## Core Goals

- Keep the same tailoring contract used by EKF.
- Mirror the EKF public architecture for SRIF/SRIS.
- Keep SRIF/SRIS native in square-root information form.
- Reuse only the helpers that are truly representation-neutral.
- Preserve MATLAB Coder compatibility throughout the implementation.

## Tailoring Boundary

- The following functions are the manually tailored mission seam and must remain the same public contract:
  - `BuildArchitectureTemplate`
  - `BuildInputStructsTemplate`
  - `ComputeMeasPredAndObsJacobian`
  - `ComputeMeasResiduals`
  - `ComputeProcessNoiseCov`
- EKF/SRIF-style linearized observation updates must use `ComputeMeasPredAndObsJacobian` when they need
  prediction and Jacobian together. UKF/SR-UKF-style updates call the same hook with prediction outputs
  only and must not be forced to evaluate an observation Jacobian.
- Measurement-tailoring implementations must receive both `strFilterMutabConfig` and
  `strFilterConstConfig`, and must use the EKF-style `strFilterConstConfig.strStatesIdx` contract for
  state indexing. Use only state-index fields already present in `BuildArchitectureTemplate`; do not
  introduce ad-hoc observed-state index fields in `strMeasModelParams` or new shared-template
  `strStatesIdx` fields.
- Do not add an SRIF-specific tailoring package.
- Do not force mission repos to rewrite tailoring code to switch from EKF to SRIF.

## Public Architecture Shape

- SRIF/SRIS must mirror the EKF sliding-window public architecture rather than inventing a new user-facing layout.
- Public SRIF/SRIS entrypoints should be:
  - `SRIF_SlideWindow_step`
  - `SRIF_SlideWindow_TimeUp`
  - `SRIF_SlideWindow_ObsUp`
  - `SRIF_SlideWindow_StateManagementStep`
  - `SRIS_SlideWindow_Smooth`
  - later `SRIS_FullHistory_Smooth`
- Backend selection must be const-config driven.
- Smoothing must remain a separate layer from forward-filter backend selection.

## Shared vs Specific Ownership

- Shared/public runtime entrypoints stay outside `+filter_tailoring`.
- Generic Givens rotation application belongs in MathCore as `ApplyGivensRot`.
- `ComputeGivensRotMeasUpdate_SRIF` is the low-level SRIF measurement-update entrypoint; it is not the owner of the generic primitive.
- `sharedFiltersModules/measurementEditing` owns representation-neutral residual gating/proposal and editing-policy helpers.
- `sharedFiltersModules/observationModels` owns measurement prediction, measurement corrections, measurement covariance models, and observation Jacobians that do not mutate backend covariance/information state.
- `sharedFiltersModules/slidingWindow` owns representation-neutral window pose construction, pose buffers, and pose retraction/update helpers.
- `ekf_modules/components` owns EKF update kernels, covariance mutation, and covariance-owned sliding-window operations.
- Use shared helpers only when the logic is truly representation-neutral or reused materially by more than one architecture.
- If the logic depends on covariance-form ownership, information-form ownership, or smoothing-specific algebra, keep it architecture-specific.

## MathCore Ownership Rule

- The generic Givens primitive must be implemented and owned in MathCore.
- EstimationGears should consume the MathCore primitive rather than becoming the source of truth.
- If multiple MathCore roots or mirrored checkouts exist, confirm the intended owner root before editing the shared primitive.

## Representation Rules

- SRIF and SRIS must carry native square-root information state internally.
- Covariance or square-root covariance reconstruction is auxiliary only:
  - allowed for tests
  - allowed for diagnostics
  - allowed for explicit verification outputs
- Do not build the SRIF architecture on a covariance-backed bridge as the canonical runtime flow.
- Do not implement SRIF time update by reconstructing covariance unless explicitly marked as diagnostic-only or transitional test scaffolding.

## SRIF Observation Residual Contract

- Current implementation status: `SRIF_SlideWindow_ObsUp` still uses a transitional absolute-state
  pseudo-observation, `r_prior + H * x_prior`, and the low-level SRIF update returns an absolute
  posterior state that replaces the prior state directly. This has been verified as internally
  consistent, but it is not the target contract.
- Prioritize migrating this path to the error-state residual contract before adding delayed
  measurements, multi-measurement editing/underweighting, or wider sliding-window behavior.
- The low-level SRIF measurement-update kernel must operate on a linearized error-state residual model, not on absolute measurements:
  - `r_prior = y_meas - h(x_prior)` or the tailored manifold residual equivalent
  - `r_prior ~= H_error * dx + v`
  - `dx` is the local error-state increment
- `ComputeGivensRotMeasUpdate_SRIF` should solve the linear error-state problem and return the estimated correction/error state and post-fit residual products.
- `SRIF_SlideWindow_ObsUp` owns measurement prediction, prior residual construction, valid-measurement masks, measurement editing/underweighting, and manifold retraction/application of the correction.
- Do not pass a pseudo-observation such as `r_prior + H * x_prior` into the low-level SRIF kernel as the long-term contract. That pattern is Euclidean-only scaffolding and hides the required manifold local-coordinate boundary.
- Observation Jacobians consumed by SRIF must be with respect to local error-state coordinates. If a tailoring hook returns global-coordinate Jacobians, the architecture wrapper must explicitly convert them before calling the SRIF kernel.
- On-manifold states, including attitude and sliding-window pose blocks, must be updated through their proper retraction, not by adding the correction directly to the stored global coordinates.

## Helper And Wrapper Rules

- Do not add helpers just to move one short block of code out of sight.
- Do not add wrappers that only rename a single call without clarifying ownership or contract.
- Wrappers are allowed when they:
  - improve interface clarity
  - isolate architecture-specific algebra
  - preserve efficient codegen/runtime behavior
- A shared helper must justify its existence by real reuse, not by style preference alone.

## Naming Rules

- Helper functions never start with `local`.
- Helper function names end with `_`.
- Public entrypoints do not use the trailing underscore suffix.
- Shared math utilities should use precise names that state the real operation and, when useful, the numerical scheme.
- Prefer names that reflect the stable contract, not transient implementation history.

## Anti-Bloat Rules

- Do not create a generic mega-pipeline that dispatches EKF/SRIF/SRIS internals behind a large runtime switch.
- Do not duplicate EKF orchestration just to avoid thinking about helper ownership.
- Do not split one coherent architecture-level step into many tiny private helpers unless that split materially improves reuse or testability.
- Keep builders, state-management utilities, and runtime modules compact and purpose-owned.

## Codegen Rules

- Every new public SRIF/SRIS function must remain `#codegen` compatible unless explicitly documented otherwise.
- Preserve const-config specialization and avoid unsupported runtime-polymorphic behavior.
- Add codegen compile coverage as the feature is introduced, not later.
- If a design looks clean in MATLAB but complicates codegen materially, prefer the simpler coder-safe structure.

## Test Rules

- Every stage must add or update MATLAB tests before the next stage starts.
- Every public SRIF/SRIS entrypoint needs:
  - MATLAB execution tests
  - codegen compile coverage
  - generated-MEX equivalence inside the codegen validation flow when applicable
- Shared primitives moved into MathCore need their own MathCore tests in addition to EstimationGears integration tests.
- EKF-vs-SRIF parity tests are required on simple linear problems before broader scenario tests.
- Tests that depend on external `.mat` workspace fixtures are deprecated until reworked to generated or builder-style test data.

## State-Management Rules

- Audit EKF sliding-window helpers before reusing them.
- Reuse only helpers that are representation-neutral.
- Rewrite or wrap helpers whose logic is covariance-owned.
- Do not discover hidden covariance assumptions late inside `SRIF_SlideWindow_step`; classify them earlier during the state-management audit stage.
- SRIF sliding-window marginalization must be designed in square-root information form before `SRIF_SlideWindow_step` is claimed complete.

## Adaptivity And Consider-State Rules

- The standalone external adaptivity layer is the source of truth for adaptive algorithms, adaptive buffers, outage policy, underweighting policy, and covariance tuning.
- Do not add SRIF-owned adaptive algorithms or an SRIF-specific public adaptivity step.
- Add SRIF adaptivity only after the nominal SRIF step runs end-to-end and the standalone adaptive algorithms are modernized in isolation.
- Keep SR-UKF estimator kernels aligned with the same rule: residual, gain, and innovation statistics may be exposed, but adaptive `R`/`Q` mutation belongs in a wrapper-level `AdaptiveTuning` integration.
- SRIF wrappers may expose the statistics required by the external adaptivity layer and consume its returned policy/tuning outputs.
- Match current EKF consider-state semantics before any smoothing work starts.
- SR-UKF covariance-root observation updates already use square-root Schmidt covariance handling for consider states. Do not copy that implementation blindly into SRIF; information-form consider support needs explicit solved-state selection and wrapper-level handling of consider-state uncertainty.
- Optional Schmidt-style tightening for EKF/SRIF/SRIS remains a later refinement, not part of the first parity rollout.
- Future low-level SRIF consider support should use explicit solved-state selection and skip rotations for unsolved columns. The wrapper must account for consider-state uncertainty in the effective residual/covariance model before calling the low-level kernel.

## Commit Rules

- Commit the current runtime/template/process-noise cleanup separately from SRIF/SRIS work.
- Commit shared MathCore Givens extraction separately from SRIF architecture buildout.
- Keep low-level SRIF hardening separate from architecture-level SRIF modules when possible.
- Avoid mixed commits that combine prerequisite cleanup, shared-math ownership moves, and architecture-level SRIF/SRIS behavior.
