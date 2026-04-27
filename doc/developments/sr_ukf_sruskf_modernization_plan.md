# SR-UKF/SR-USKF Modernization Plan

## Goal

- [x] Rename active SR-UKF/SR-USKF entrypoints to the `SR_UKF_*` convention used by the modern filter modules.
- [x] Remove built-in adaptive covariance tuning from estimator kernels.
- [x] Provide current-state SR-UKF time update, observation update, UD-facing observation wrapper, and full-step orchestration.
- [x] Preserve codegen compatibility with one validation script for all public SR-UKF modules.
- [ ] Keep future adaptivity integration at wrapper level through the standalone `AdaptiveTuning` and `AdaptiveLogics` contracts.

## Active Public Entry Points

- [x] `SR_UKF_ObsUp`: square-root covariance observation update. It performs measurement prediction, residual handling, measurement editing, state correction, and square-root posterior covariance factoring. It does not mutate `R` or `Q` through built-in adaptivity.
- [x] `SR_UKF_UDCov_ObsUp`: transitional UD-covariance wrapper over `SR_UKF_ObsUp`. It converts UD to upper square-root covariance and maps the posterior upper square-root factor back to UD with `UDdecompositionFromUpperSqrtCov`, avoiding the dense `S' * S` posterior intermediate.
- [x] `SR_UKF_TimeUp`: current-state square-root sigma-point time update renovated from deprecated SR-USKF logic.
- [x] `SR_UKF_step`: EKF-style orchestration wrapper for time update plus optional observation update.

## Current Boundaries

- [x] Current-state propagation is supported.
- [x] `SR_UKF_ObsUp` uses upper-root Schmidt covariance handling for consider states: update solved states and cross-covariances while preserving consider-state mean and auto-covariance.
- [x] `SR_UKF_TimeUp` preserves consider-state mean and auto-covariance during current-state propagation.
- [x] `SR_UKF_TimeUp` performs one full-interval sigma-point propagation step; EKF-style piecewise propagation flags are ignored by this backend.
- [x] Additive process-noise factoring uses `ComputeFactorProcessNoiseCov`.
- [x] Delayed-measurement backward propagation is explicitly unsupported in `SR_UKF_step`.
- [ ] Native sliding-window SR-UKF state management is not implemented yet.
- [ ] `SR_UKF_UDCov_ObsUp` is not a native UD sigma-point algorithm; it is an interface-compatible wrapper.

## Tests And Codegen

- [x] `tests/matlab/sigma_points_filters_modules/testSR_UKF_ObsUp.m`
- [x] `tests/matlab/sigma_points_filters_modules/testSR_UKF_TimeUp.m`
- [x] `tests/matlab/sigma_points_filters_modules/testSR_UKF_step.m`
- [x] `tests/matlab/sigma_points_filters_modules/testCommonSigmaPointTemplateComparison.m`
- [x] `tests/matlab/ekf_modules/testUDdecompositionFromUpperSqrtCov.m`
- [x] `matlab/sigma_points_filters_modules/codegenSR_UKF_Modules.m` compiles and validates generated MEX outputs for `SR_UKF_ObsUp`, `SR_UKF_UDCov_ObsUp`, `SR_UKF_TimeUp`, and `SR_UKF_step`.
- [x] `tests/matlab/ekf_modules/benchmarkUDdecompositionFromUpperSqrtCov.m` validates the direct UD path and reports dense-path vs direct-path timing.

## Performance Note

- [x] The direct UD recurrence avoids materializing the posterior covariance and is codegen-friendly. The MATLAB benchmark validates equivalence, but MATLAB-level timing can favor the old dense BLAS-backed path for moderate dimensions; do not reintroduce the dense posterior intermediate unless a target-runtime benchmark proves it is required.

## Remaining Work

- [ ] Replace `SR_UKF_UDCov_ObsUp` with a native UD sigma-point update only if strict UD arithmetic becomes a runtime requirement.
- [ ] Add wrapper-level adaptivity adapters after standalone `AdaptiveTuning` and `AdaptiveLogics` contracts are implemented.
- [ ] Design SR-UKF sliding-window state management before adding windowed sigma-point propagation.
- [ ] Add delayed-measurement support only after the sigma-point equivalent of the EKF two-stage/backward-propagation contract is designed.
