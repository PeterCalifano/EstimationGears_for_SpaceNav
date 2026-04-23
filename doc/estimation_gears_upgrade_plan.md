# Runtime Stabilization First, Then Generic EKF Modernization

## Summary
- [ ] Close the open runtime defects in EstimationGears before continuing generic EKF work or any MSCKF composition.
- [ ] Treat `SR_UKF_Adaptive_ObsUp` as the working example of the new `filter_tailoring.*` observation-hook pattern.
- [ ] Keep `adaptiveSRUSKF_ObsUp` as review-only reference and do not import from it.
- [ ] Converge `MSCKF_for_SpaceNav` toward consuming EstimationGears/SimulationGears components instead of duplicating vendored runtime code.

## Confirmed Current Blockers
- [ ] Complete [`matlab/sharedFiltersModules/dynamicsModels/ComputeLinearizedMappedQcov.m`](../matlab/sharedFiltersModules/dynamicsModels/ComputeLinearizedMappedQcov.m), which is currently incomplete.
- [ ] Make [`matlab/ekf_modules/full-covariance-sliding/modules/EKF_SlideWindow_FullCov_TimeUp.m`](../matlab/ekf_modules/full-covariance-sliding/modules/EKF_SlideWindow_FullCov_TimeUp.m) coherent with the chosen `PropagateDyn` contract.
- [ ] Establish a clean validation path for the repaired runtime before any larger refactor resumes.

## Runtime Ownership Model
- [ ] Keep shared/public runtime entrypoints outside `+filter_tailoring`:
- [ ] `ComputeDynFcn`
- [ ] `ComputeDynMatrix`
- [ ] `PropagateDyn`
- [ ] `ManageMeasLatency`
- [ ] Keep manual tailoring entrypoints inside `+filter_tailoring`:
- [ ] `BuildArchitectureTemplate`
- [ ] `BuildInputStructsTemplate`
- [ ] `ComputeMeasPred`
- [ ] `ComputeMeasResiduals`
- [ ] `ComputeObsMatrix`
- [ ] `ComputeProcessNoiseCov`
- [ ] Keep `ComputeDynFcn` and `ComputeDynMatrix` as shared defaults that remain optionally tailorable by mission-specific repos when needed.
- [ ] Do not recreate a parallel `+filter_templates_impl` implementation in EstimationGears.

## Stage 1: Runtime Coherence Only
- [ ] Finish `ComputeLinearizedMappedQcov` with mapped-noise trapezoidal computation and consider-state masking.
- [ ] Standardize the `PropagateDyn` output contract and update every EstimationGears caller to that contract.
- [ ] Add the missing runtime-side input-noise assembly needed for mapped process noise, either as a shared helper or deterministic builder logic.
- [ ] Align `BuildInputStructsTemplate` with the repaired runtime so a template-built filter can execute time update and observation update without hidden field mismatches.
- [ ] Keep the current package split intact while repairing the runtime.

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
- [ ] Use the RCS test inventory as a coverage checklist, not as a codebase to import wholesale.
- [ ] Build native EstimationGears `matlab.unittest` tests and native scenario builders in this repo.
- [ ] Add a direct regression test for `ComputeLinearizedMappedQcov`.
- [ ] Replace [`tests/matlab/ekf_modules/testEKF_SlideWindow_FullCov_TimeUp.m`](../tests/matlab/ekf_modules/testEKF_SlideWindow_FullCov_TimeUp.m) with a real `matlab.unittest` class.
- [ ] Add a builder/runtime smoke test for `BuildArchitectureTemplate`, `BuildInputStructsTemplate`, and one minimal `TimeUp` execution.
- [ ] Add a contract test for `PropagateDyn` outputs versus all shared-runtime callers.
- [ ] Fix the currently broken measurement-editing helper test.

## Native Scenario-Based Coverage
- [ ] Add minimal linear-Gaussian filter algebra tests.
- [ ] Add nonlinear orbit/FOGM time-update tests.
- [ ] Add simplified lidar + centroid observation-update tests using SimulationGears sublib primitives.
- [ ] Add one end-to-end simulated estimation scenario in this repo.

## After Runtime Stabilization
- [ ] Resume generic observation-update orchestration around `filter_tailoring.*`.
- [ ] Resume generic time-update cleanup on top of the repaired runtime.
- [ ] Avoid new helper sprawl unless justified by reuse or testability.
- [ ] Start MSCKF composition work only after the repaired EKF runtime is stable and validated.

## Analytical Validation
- [ ] Add STM consistency checks.
- [ ] Add mapped process-noise consistency checks.
- [ ] Add FOGM decay consistency checks.
- [ ] Add innovation/NIS sanity checks in simulated runs.

## Assumptions And Defaults
- [ ] Leave `SR_UKF_Adaptive_ObsUp` alone except as reference for the new tailoring-hook pattern.
- [ ] Make downstream wrapper support conditional on actual consumer migration; do not add a blanket wrapper layer preemptively.
- [ ] Keep EstimationGears as the canonical runtime source, with `future-nav`, `nav-backend`, and later MSCKF converging toward consuming it.
