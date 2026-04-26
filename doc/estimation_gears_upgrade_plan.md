# Runtime Stabilization First, Then Generic EKF Modernization

## Summary

- [x] Close the open runtime defects in EstimationGears before continuing generic EKF work or any MSCKF composition.
- [x] Treat `SR_UKF_Adaptive_ObsUp` as the working example of the new `filter_tailoring.*` observation-hook pattern.
- [x] Keep `adaptiveSRUSKF_ObsUp` as review-only reference and do not import from it.
- [ ] Converge `MSCKF_for_SpaceNav` toward consuming EstimationGears/SimulationGears components instead of duplicating vendored runtime code.

## Confirmed Current Blockers

- [x] Complete [`matlab/sharedFiltersModules/dynamicsModels/ComputeTrapzMappedProcessNoiseCov.m`](../matlab/sharedFiltersModules/dynamicsModels/ComputeTrapzMappedProcessNoiseCov.m), which is currently incomplete.
- [x] Make [`matlab/ekf_modules/full-covariance-sliding/modules/EKF_SlideWindow_FullCov_TimeUp.m`](../matlab/ekf_modules/full-covariance-sliding/modules/EKF_SlideWindow_FullCov_TimeUp.m) coherent with the chosen `PropagateDyn` contract.
- [x] Establish a clean validation path for the repaired runtime before any larger refactor resumes.

## Runtime Ownership Model

- [x] Keep shared/public runtime entrypoints outside `+filter_tailoring`:
- [x] `ComputeDynFcn`
- [x] `ComputeDynMatrix`
- [x] `PropagateDyn`
- [x] `ManageMeasLatency`
- [x] Keep manual tailoring entrypoints inside `+filter_tailoring`:
- [x] `BuildArchitectureTemplate`
- [x] `BuildInputStructsTemplate`
- [x] `ComputeMeasPred`
- [x] `ComputeMeasResiduals`
- [x] `ComputeObsMatrix`
- [x] `ComputeProcessNoiseCov`
- [x] Keep `ComputeDynFcn` and `ComputeDynMatrix` as shared defaults that remain optionally tailorable by mission-specific repos when needed.
- [x] Do not recreate a parallel `+filter_templates_impl` implementation in EstimationGears.

## Stage 1: Runtime Coherence Only

- [x] Finish `ComputeTrapzMappedProcessNoiseCov` with mapped-noise trapezoidal computation and consider-state masking.
- [x] Standardize the `PropagateDyn` output contract and update every EstimationGears caller to that contract.
- [x] Add the missing runtime-side input-noise assembly needed for mapped process noise, either as a shared helper or deterministic builder logic.
- [x] Align `BuildInputStructsTemplate` with the repaired runtime so a template-built filter can execute time update and observation update without hidden field mismatches.
- [x] Keep the current package split intact while repairing the runtime.

## Compatibility Policy

- [ ] Add thin compatibility wrappers only if old plain-function names are actively needed when downstream repos migrate to the updated runtime.
- [ ] Support wrappers only for actively used legacy plain-function names:
- [ ] `computeDynFcn`
- [ ] `computeDynMatrix`
- [ ] `propagateDyn`
- [ ] `manageMeasLatency`
- [ ] `computeMeasPred`
- [ ] `computeMeasResiduals`
- [ ] `computeObsMatrix`
- [ ] `computeProcessNoiseCov`
- [ ] Do not add a package-alias compatibility layer for `filter_templates_impl.*` unless a real caller appears.

## Stage 2: Tests And Native Scenario Builders

- [x] Use the RCS test inventory as a coverage checklist, not as a codebase to import wholesale.
- [x] Build native EstimationGears `matlab.unittest` tests and native scenario builders in this repo.
- [x] Add a direct regression test for `ComputeTrapzMappedProcessNoiseCov`.
- [x] Replace [`tests/matlab/ekf_modules/testEKF_SlideWindow_FullCov_TimeUp.m`](../tests/matlab/ekf_modules/testEKF_SlideWindow_FullCov_TimeUp.m) with a real `matlab.unittest` class.
- [x] Add a builder/runtime smoke test for `BuildArchitectureTemplate`, `BuildInputStructsTemplate`, and one minimal `TimeUp` execution.
- [x] Add a contract test for `PropagateDyn` outputs versus all shared-runtime callers.
- [x] Fix the currently broken measurement-editing helper test.

## Native Scenario-Based Coverage

- [ ] Add minimal linear-Gaussian filter algebra tests.
- [ ] Add nonlinear orbit/FOGM time-update tests.
- [ ] Add simplified lidar + centroid observation-update tests using SimulationGears sublib primitives.
- [ ] Add one end-to-end simulated estimation scenario in this repo.

## After Runtime Stabilization

- [ ] Execute the standalone adaptivity algorithm modernization plan in [`doc/developments/adaptivity_algorithms_modernization_plan.md`](developments/adaptivity_algorithms_modernization_plan.md) before reintegrating adaptive behavior into estimator wrappers.
- [ ] Resume generic observation-update orchestration around `filter_tailoring.*`.
- [ ] Resume generic time-update cleanup on top of the repaired runtime.
- [ ] Keep adaptive algorithms, adaptive buffers, and adaptive policy state outside estimator kernels; wrappers consume the external adaptivity layer after the algorithms are validated in isolation.
- [ ] Avoid new helper sprawl unless justified by reuse or testability.
- [ ] Start MSCKF composition work only after the repaired EKF runtime is stable and validated.

## Analytical Validation

- [ ] Add STM consistency checks.
- [ ] Add mapped process-noise consistency checks.
- [ ] Add FOGM decay consistency checks.
- [ ] Add innovation/NIS sanity checks in simulated runs.
- [ ] Add seeded randomized algorithmic tests for adaptive-noise and adaptive-policy logic as required by the standalone adaptivity plan.

## Assumptions And Defaults

- [x] Leave `SR_UKF_Adaptive_ObsUp` alone except as reference for the new tailoring-hook pattern.
- [x] Make downstream wrapper support conditional on actual consumer migration; do not add a blanket wrapper layer preemptively.
- [x] Keep EstimationGears as the canonical runtime source, with `future-nav`, `nav-backend`, and later MSCKF converging toward consuming it.
