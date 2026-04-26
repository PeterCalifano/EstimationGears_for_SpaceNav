# Adaptivity Algorithms Modernization Plan

## Goal

- [x] Modernize the adaptive algorithms first, in isolation from the filter-step implementations.
- [x] Validate adaptive algorithms with algorithmic tests rather than code-shape tests.
- [x] Delay integration into EKF, UKF, UD, and SRIF wrappers until the algorithms themselves are cleaned up and trusted.
- [x] Split the adaptive architecture into an outer `AdaptiveTuning` layer and an inner `AdaptiveLogics` layer.

## Layer Model

- [x] `AdaptiveTuning` is the most external layer and owns adaptation of tuning values such as measurement-noise covariance and process-noise covariance.
- [x] Adaptive `R` and adaptive `Q` are generic tuning categories in `AdaptiveTuning`.
- [x] ASNC is a specific adaptive-`Q` algorithm under `AdaptiveTuning`, not a separate adaptation category.
- [x] `AdaptiveLogics` owns immediate decision logic such as outage handling, underweighting triggers, consider-state triggers, counters, and policy state.
- [x] `AdaptiveLogics` may act in-cycle at wrapper level.
- [x] `AdaptiveTuning` updates tuning values for the next cycle.

## Stage 0: Inventory And Scope Freeze

- [ ] Inventory every adaptive mechanism in the repo: adaptive measurement-noise estimation, adaptive process-noise estimation, ASNC-style process-noise estimation, outage handling, underweighting, consider-state triggering, and adaptive buffer/state update logic.
- [ ] Classify each mechanism as `AdaptiveTuning` or `AdaptiveLogics`.
- [ ] Classify each mechanism as immediate-policy logic or next-cycle tuning logic.
- [ ] Mark each path as active and to be modernized, or legacy/deprecated and review-only.
- [ ] Produce one review table with algorithm, callers, required inputs, internal state, outputs, estimator-family dependencies, and reference basis.

- [ ] Gate: no implementation starts before the inventory table is complete.

## Stage 1: Reference Review And Intended Behavior

- [ ] Review each adaptive algorithm against the relevant references and the repo's current intended use.
- [ ] Define the mathematical purpose, invariants, valid input domain, numerical constraints, expected response direction, and degenerate/failure handling for each algorithm.
- [ ] Decide explicitly what behavior is preserved, corrected, simplified, or removed.

- [ ] Gate: each algorithm has a decision-complete behavioral spec.

## Stage 2: Standalone Interface Design

- [ ] Define standalone contracts for each adaptive algorithm with explicit config parameters, runtime statistics, owned buffers/state, returned updated state, and returned actions or tuning values.
- [x] Keep the split explicit between reusable algorithm core, family-specific adapter, and future wrapper integration.
- [x] Preserve the mixed timing model: immediate policy logic may act in-cycle, while covariance tuning updates apply to the next cycle.
- [ ] Define separate standalone contracts for `AdaptiveTuning` and `AdaptiveLogics`, including their handoff responsibilities.

- [ ] Gate: each adaptive algorithm can be exercised without calling EKF, UKF, UD, or SRIF step implementations directly.

## Stage 3: Refactor And Isolate Algorithms

- [ ] Refactor current implementations into standalone `AdaptiveTuning` and `AdaptiveLogics` modules.
- [ ] Remove hidden coupling to observation-step internals.
- [ ] Make buffer ownership, reset behavior, clipping/bounds, no-op behavior, empty-history handling, and missing-measurement handling explicit.
- [ ] Keep estimator-family specifics in thin adapters only where unavoidable.

- [ ] Gate: the adaptive algorithms run on synthetic inputs alone.

## Stage 4: Deterministic Algorithm Tests

- [ ] Add hand-built tests for nominal and edge cases.
- [ ] Check algorithmic properties, including symmetry, positive semidefiniteness where required, bounds enforcement, counter transitions, reset behavior, no-op behavior, and trigger correctness.
- [ ] Use small interpretable fixtures so failures are easy to diagnose.
- [ ] Test `AdaptiveTuning` value updates separately from `AdaptiveLogics` trigger decisions.

- [ ] Gate: every adaptive algorithm has deterministic nominal and edge-case coverage.

## Stage 5: Randomized Algorithmic Validation

- [ ] Add seeded randomized tests for all adaptive algorithms.
- [ ] Use Monte Carlo and property-style checks over valid parameter ranges and noisy sequences.
- [ ] Validate symmetry and positive semidefiniteness where required, bounded outputs, stable behavior under noisy inputs, expected adaptation direction under systematic residual growth, absence of unphysical drift for nominal-consistent data, and correct outage/underweight/consider behavior across random event sequences.
- [ ] Keep random seeds fixed and acceptance thresholds explicit.
- [ ] Run randomized tests for generic adaptive `R`, generic adaptive `Q`, ASNC as a specific adaptive-`Q` algorithm, and `AdaptiveLogics` event sequences.

- [ ] Gate: randomized tests pass reproducibly with stable thresholds.

## Stage 6: Legacy Cross-Check

- [ ] Compare isolated modernized algorithms against current repo behavior where that behavior is intended to be preserved.
- [ ] Record matched behavior, intentional behavior changes, bug fixes, and legacy artifacts that are not preserved.
- [ ] Avoid preserving incorrect legacy behavior accidentally.

- [ ] Gate: every meaningful delta from legacy behavior is explained.

## Stage 7: Integration Design For Later Landing

- [ ] Define the reintegration contract only after isolated algorithms are stable.
- [ ] Standardize the outputs wrappers must provide to adaptivity: prior/posterior residual data, innovation covariance or factor, gain/update statistics, timing context, and measurement acceptance/rejection masks.
- [ ] Plan wrapper-level integration for full-covariance EKF, UD EKF, SR-UKF/SRUSKF, and SRIF.
- [x] Keep adaptivity external to estimator kernels by design.
- [ ] Define wrapper call order for `AdaptiveLogics` in-cycle actions and `AdaptiveTuning` next-cycle updates.

- [ ] Gate: integration can proceed without redesigning the adaptive algorithms.

## Test Philosophy

- [x] Tests are algorithmic first.
- [x] Functional API tests alone are insufficient.
- [x] Randomized validation is mandatory.
- [x] Acceptance depends on mathematical behavior, numerical robustness, and reproducibility.

## Reference Baseline

- [x] Akhlaghi et al., 2017, adaptive measurement/process-noise covariance adjustment.
- [x] Stacey and D'Amico, 2021, adaptive and dynamically constrained process-noise estimation.
- [x] Repo-local Tapley and Myers-Tapley references for covariance and process-noise structure.

## Assumptions

- [x] The first implementation stage is algorithm-first and intentionally delays estimator integration.
- [x] Existing in-step adaptivity paths are legacy behavior sources to compare against, not the target architecture.
- [x] Same-cycle policy actions belong to `AdaptiveLogics` and stay outside estimator kernels.
- [x] Covariance tuning updates belong to `AdaptiveTuning` and apply to the next cycle.
