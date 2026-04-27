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

- [x] Inventory every adaptive mechanism in the repo: adaptive measurement-noise estimation, adaptive process-noise estimation, ASNC-style process-noise estimation, outage handling, underweighting, consider-state triggering, and adaptive buffer/state update logic.
- [x] Classify each mechanism as `AdaptiveTuning` or `AdaptiveLogics`.
- [x] Classify each mechanism as immediate-policy logic or next-cycle tuning logic.
- [x] Mark each path as active and to be modernized, or legacy/deprecated and review-only.
- [x] Produce one review table with algorithm, callers, required inputs, internal state, outputs, estimator-family dependencies, and reference basis.

- [x] Gate: no implementation starts before the inventory table is complete.

| Mechanism | Current Owner | Current Callers | Layer | Timing | Required Inputs | Owned State | Outputs | Estimator Dependency | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Adaptive measurement-noise covariance | `AdaptMeasCov` | `AdaptRQcovs`; no active SR-UKF kernel caller after `SR_UKF_ObsUp` externalization | `AdaptiveTuning` | next-cycle tuning | base `R`, forgetting factor, post residual, predicted measurement, sigma-point measurement cloud, sigma-point weights, valid-measurement mask | none currently explicit | adapted `R` block | sigma-point statistics today; should accept generic innovation/post-fit statistics | active algorithm, modernize before reintegration |
| Adaptive process-noise covariance | `AdaptProcessCov` | `AdaptRQcovs`; no active SR-UKF kernel caller after `SR_UKF_ObsUp` externalization | `AdaptiveTuning` | next-cycle tuning | base `Q`, forgetting factor, residual, Kalman/update gain | none currently explicit | adapted `Q` block | gain-based filter statistics | active algorithm, modernize before reintegration |
| ASNC process-noise estimation | `AdaptQCovASNC` | no stable architecture-level caller confirmed | `AdaptiveTuning` | next-cycle tuning | timestep, gain buffer, residual-covariance buffer, Q bounds, delta-state buffer, posterior covariance buffer | window buffers supplied externally today | discrete SNC covariance and process-noise covariance | position/velocity state convention | active algorithm candidate, isolate before integration |
| Combined adaptive R/Q wrapper | `AdaptRQcovs` | no active estimator-kernel caller after SR-UKF cleanup | adapter over `AdaptiveTuning` | next-cycle tuning | statistics needed by `AdaptMeasCov` and `AdaptProcessCov`, enable flags | none currently explicit | adapted `R`, adapted `Q` | sigma-point wrapper | transitional legacy adapter, replace after standalone contracts exist |
| Measurement editing/rejection | `EvaluateMeasRejectionProposal`, `ApplyMeasurementEditingPolicy` | EKF and sigma-point observation paths | `AdaptiveLogics` | immediate in-cycle policy | residual, Mahalanobis threshold, counters, max rejection count | editing counters in mutable filter config | accepted/rejected residual policy and updated counters | representation-neutral | active, keep external to kernels |
| Measurement outage consider-state trigger | `EKF_SlideWindow_AdaptivityManagementStep` | EKF sliding-window step | `AdaptiveLogics` | immediate in-cycle policy | measurement-type flags, outage counters, patience, state index config | outage and underweight counters in mutable filter config | consider-state mask, reset bias states, updated counters | EKF state layout today | active behavior source, generalize later |
| Underweighting after outage | `EKF_SlideWindow_AdaptivityManagementStep` | EKF sliding-window step | `AdaptiveLogics` | immediate in-cycle policy | underweight enable flag, counter, max duration, underweight coefficient | underweight counter in mutable filter config | measurement underweight coefficient and counter update | EKF wrapper today | active behavior source, generalize later |

## Stage 1: Reference Review And Intended Behavior

- [x] Review each adaptive algorithm against the relevant references and the repo's current intended use.
- [x] Define the mathematical purpose, invariants, valid input domain, numerical constraints, expected response direction, and degenerate/failure handling for each algorithm.
- [x] Decide explicitly what behavior is preserved, corrected, simplified, or removed.

- [x] Gate: each algorithm has a decision-complete behavioral spec.

| Mechanism | Purpose | Invariants | Valid Domain | Expected Response | Degenerate Handling | Decision |
| --- | --- | --- | --- | --- | --- | --- |
| Adaptive `R` | Tune measurement-noise covariance from post-fit residual statistics | symmetric covariance; updated valid block only; no hidden estimator mutation | nonempty valid-measurement mask, finite residuals, finite weights | residual growth increases relevant `R` entries; nominal residuals decay toward base/statistical value | no valid measurements returns input `R`; invalid covariance must be rejected or projected by the standalone algorithm | preserve purpose, make symmetry/PSD/bounds explicit |
| Adaptive `Q` | Tune process-noise covariance from residual/update statistics | symmetric covariance; bounded to configured subspace | finite residual and update-gain statistics | systematic residual growth increases process-noise estimate in driven states | missing gain/statistics returns input `Q` with no-op status | preserve purpose, separate generic tuning from wrapper statistics |
| ASNC | Estimate SNC process noise from correction/covariance windows | PSD covariance; explicit bounds; fixed state block ownership | filled or partially filled finite buffers with known sample count | larger correction energy increases SNC process noise within bounds | empty window returns prior/default and no-op status | preserve as standalone candidate; do not integrate until deterministic and randomized tests exist |
| Editing/rejection | Reject or override outlier residuals using Mahalanobis policy | counter transitions deterministic; max-rejection override deterministic | finite residual, finite threshold, nonnegative counters | outliers rejected until override threshold; accepted measurements reset counters | invalid residual or disabled editing returns accepted/no-op policy | preserve current behavior, keep representation-neutral |
| Outage/consider/underweight | Trigger conservative consider-state and underweight policies after measurement outage | counters saturate/reset deterministically; state resets are explicit wrapper actions | measurement-type flags and configured state indices present | outage increments counters and eventually raises consider/underweight actions | missing measurement fields should no-op with explicit status in standalone contract | preserve behavior source, move state-layout specifics into adapter |

## Stage 2: Standalone Interface Design

- [x] Define standalone contracts for each adaptive algorithm with explicit config parameters, runtime statistics, owned buffers/state, returned updated state, and returned actions or tuning values.
- [x] Keep the split explicit between reusable algorithm core, family-specific adapter, and future wrapper integration.
- [x] Preserve the mixed timing model: immediate policy logic may act in-cycle, while covariance tuning updates apply to the next cycle.
- [x] Define separate standalone contracts for `AdaptiveTuning` and `AdaptiveLogics`, including their handoff responsibilities.

- [x] Gate: each adaptive algorithm can be exercised without calling EKF, UKF, UD, or SRIF step implementations directly.

### Standalone Contract Sketch

`AdaptiveTuning` algorithms should consume only:

- immutable tuning configuration such as enable flags, forgetting factors, bounds, target blocks, and algorithm IDs
- runtime statistics such as prior/post residuals, predicted measurement statistics, innovation covariance or square-root/information equivalent, update gain or update-equivalent statistics, and timing context
- owned tuning state such as ASNC buffers, previous tuned values, and sample counters

`AdaptiveTuning` algorithms should return:

- updated tuning state
- next-cycle covariance/noise tuning values
- no-op/status flags explaining disabled, missing-data, bounded, or rejected updates

`AdaptiveLogics` algorithms should consume only:

- immutable policy configuration such as thresholds, patience, maximum consecutive rejection/outage counts, underweight coefficients, and state groups
- runtime event statistics such as measurement availability, residual gates, validity masks, and timing context
- owned policy state such as counters and active policy flags

`AdaptiveLogics` algorithms should return:

- updated policy state
- immediate wrapper-level actions such as residual acceptance masks, underweight coefficients, consider-state masks, reset-state masks, and no-op/status flags

Estimator wrappers must provide adapters from EKF, SR-UKF, UD, and SRIF statistics into these standalone inputs. The adaptive algorithms must not call estimator step implementations directly.

## Stage 3: Refactor And Isolate Algorithms

- [x] Remove built-in adaptive `R`/`Q` mutation from the active SR-UKF observation update; `SR_UKF_ObsUp` now exposes residual, gain, error-state, and innovation square-root statistics for an external wrapper-level adaptivity step.
- [ ] Refactor current implementations into standalone `AdaptiveTuning` and `AdaptiveLogics` modules.
- [ ] Remove remaining hidden coupling to observation-step internals.
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
