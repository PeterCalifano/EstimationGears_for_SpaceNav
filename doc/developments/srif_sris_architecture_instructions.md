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
  - `ComputeMeasPred`
  - `ComputeObsMatrix`
  - `ComputeMeasResiduals`
  - `ComputeProcessNoiseCov`
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
- `GivensRotSRIF` is an SRIF wrapper/algorithm entrypoint, not the owner of the generic primitive. Rename it to something more representative and respecting the naming style of the repo.
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

## Adaptivity And Consider-State Rules
- The standalone external adaptivity layer is the source of truth for adaptive algorithms, adaptive buffers, outage policy, underweighting policy, and covariance tuning.
- Do not add SRIF-owned adaptive algorithms or an SRIF-specific public adaptivity step.
- Add SRIF adaptivity only after the nominal SRIF step runs end-to-end and the standalone adaptive algorithms are modernized in isolation.
- SRIF wrappers may expose the statistics required by the external adaptivity layer and consume its returned policy/tuning outputs.
- Match current EKF consider-state semantics before any smoothing work starts.
- Optional Schmidt-style tightening is a later refinement, not part of the first parity rollout.

## Commit Rules
- Commit the current runtime/template/process-noise cleanup separately from SRIF/SRIS work.
- Commit shared MathCore Givens extraction separately from SRIF architecture buildout.
- Keep low-level SRIF hardening separate from architecture-level SRIF modules when possible.
- Avoid mixed commits that combine prerequisite cleanup, shared-math ownership moves, and architecture-level SRIF/SRIS behavior.
