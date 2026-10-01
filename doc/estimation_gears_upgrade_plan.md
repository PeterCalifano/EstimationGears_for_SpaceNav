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
- [x] `PropagateDyn`
- [x] `ManageMeasLatency`
- [x] Treat `ComputeDynFcn` and `ComputeDynMatrix` as consumer-owned plain-function tailoring hooks; EstimationGears tests install temporary fixtures for them.
- [x] Keep manual tailoring entrypoints inside `+filter_tailoring`:
- [x] `BuildArchitectureTemplate`
- [x] `BuildInputStructsTemplate`
- [x] `ComputeMeasPred`
- [x] `ComputeMeasResiduals`
- [x] `ComputeObsMatrix`
- [x] `ComputeProcessNoiseCov`
- [x] Do not ship EstimationGears defaults named `ComputeDynFcn` or `ComputeDynMatrix`, because those names shadow mission-specific tailoring in consumer repos.
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

## Optional SRP LUT consolidation - 30 September 2026

Owner: `feature/upgrade-filter-smoother-interface`. Reviewer: Codex gpt-6.
Keep this batch independent of the unfinished filter/smoother modernization.
Stage the SRP adapters, their existing orbit dispatch, code generator, focused
tests and consolidated documentation. Leave dependency revisions unchanged.

### Ownership correction requested after staged review

Route force composition through SimulationGears `EvalRHS_InertialDynOrbit`.
Remove LUT force evaluation and residual injection from both EstimationGears
orbit wrappers. Keep filter state mapping and consider-mode extraction here;
keep physical force, unit conversion and bias partials in SimulationGears.

- [x] Inspect callers, current indexes and the shared RHS ownership boundary
- [x] Move the generic LUT/bias kernel into SimulationGears and add optional
  SRP inputs to its orbital RHS while preserving existing positional calls
- [x] Share filter input preparation; retain a thin filter adapter for direct
  force/Jacobian code generation and remove force composition from wrappers
- [x] Validate active/consider/absent bias, units, eclipse, force reporting,
  unchanged cannonball behavior and generated interfaces
- [x] Review and restage only the EstimationGears correction; preserve the
  SimulationGears panel index and leave its required SRP provider edits unstaged

Keep one implementation of bias physics. Pass resolved numeric spacecraft
data and an immutable table across the library boundary; do not make
SimulationGears depend on filter configuration or EstimationGears callbacks.

### Stage 1: Establish scope

- [x] Inspect the branch, complete index, working changes and dependency checkouts
- [x] Preserve unrelated plans, `CLAUDE.md`, nested checkout changes and the
  existing four-file SimulationGears panel-visibility index
- [x] Identify the SimulationGears LUT provider and committed MathCore hash
  utility required by the focused validation

### Stage 2: Review source and documentation

- [x] Review all twelve SRP MATLAB files for documentation, formatting,
  descriptive names, imperative comments and unnecessary complexity
- [x] Consolidate the standalone SRP note into `doc/main_page.md`
- [x] Preserve cannonball selection, bias units and consider-state semantics;
  retain one-output callbacks and generated Jacobian defaults

### Stage 3: Validate and stage

- [x] Run focused source contracts, shared-evaluation profiling and Code
  Analyzer on every candidate; verify the committed cannonball behavior
- [x] Verify generated RHS/Jacobian interfaces and existing orbit dispatch
- [x] Restage the fourteen-path allowlist and only this section of the mixed plan
- [x] Review the complete corrected index and verify protected files/indexes are unchanged
- [x] Prepare the corrected staged batch for the user's review and commit

### Stage 4: Integrate dependencies separately

- [ ] Consolidate the SimulationGears LUT provider
- [ ] Update dependency revisions after the relevant provider commits
- [ ] Verify ordinary recorded-dependency setup and downstream consumers

Validation currently resolves the sibling SimulationGears consolidation
checkout explicitly. Its recorded MathCore dependency predates the shared
hash utility, so resolve the committed canonical utility explicitly for the
focused harness. This does not qualify ordinary submodule-only setup. Keep
that integration gate separate from the source and code-generation checks.

#### Initial consolidation validation (before the ownership correction)

Pass twelve unit/mode/bias cases, sixteen pole cases and the 60-second RK4/STM
contract. Compare all three cannonball entry points with committed `59a3217`
in 24 cases across units, consider states, pressure and eclipse. Profile two
LUT evaluations for separate calls and one for the combined interface.
Pass sixteen fresh MEX/C++ builds and 48 runtime parity queries, including
the existing orbit dispatch, absent-bias layout and generated output prefixes.
Code Analyzer reports zero findings in all eleven MATLAB files. No mission
simulation or flight-hardware timing qualification ran.

The fixed-ephemeris adapter also aligns its function declaration with the
existing filename and bounds interpolation storage by coefficient capacity.
Pass one additional fixed-capacity MEX build with three runtime degrees
(2, 3 and 4), with exact source/MEX agreement. The initial probe requested
degree zero, below the pre-existing Chebyshev evaluator's minimum of two;
correct the probe and retain that input constraint without changing algorithms.

The first review run stopped on three Code Analyzer requests for an explicit
`Input` attribute in test-helper argument blocks. Apply those attributes and
rerun the complete gate successfully; change no numerical behavior. Retain
both logs and the tested source snapshots in the external evidence folder:
`/tmp/estimation-srp-consolidation-review-20260930-8w856r3z`.

#### Corrected ownership validation

Replace the filter-owned force kernel with `EvalFilterSRPLutWithBias`, a
filter input adapter to the SimulationGears physical kernel. Share input
preparation with both orbit wrappers through `BuildFilterSrpLutInputs`.
The orbital wrappers forward numerical inputs without evaluating LUT force
or adding SRP to residual acceleration. Keep physical bias mathematics and
unit conversion in the provider; keep state indices and consider flags here.

Pass 30 standalone orbital composition/reporting cases, twelve filter cases,
sixteen pole cases, the RK4/STM contract and sixteen adjacent orbital/gravity
tests. Compare the three cannonball entry points in 24 cases with committed
EstimationGears and SimulationGears sources. Profile two table evaluations
for separate calls, one for combined force/Jacobian and one through the
shared orbital RHS. Pass sixteen fresh filter MEX/C++ builds and 48 runtime
queries, plus two standalone SimulationGears MEX builds and eight LUT queries
with EstimationGears absent from its path. Preserve the default generated
cannonball diagnostic layout. Code Analyzer reports zero findings in all
twelve filter and three provider MATLAB files.
Pass one additional fixed-capacity MEX at runtime polynomial degrees 2, 3 and
4 with exact source/generated agreement, for nineteen fresh builds in total.

Correct two MATLAB Coder issues introduced by input preparation: assign empty
and selected struct schemas in separate compile-time branches, and call the
existing constant-unit pressure helper with a literal selector in each runtime
unit branch. Preserve force mathematics and retain both failed logs. Keep
the provider's four-file panel index unchanged and its required SRP source
unstaged. Current evidence:
`/tmp/srp-orbital-rhs-ownership-20260930-x2x77wan`.

## RHS and Jacobian naming consolidation - 30 September 2026

Standardize shared MATLAB providers and their active consumers on `EvalRHS_…`
and `EvalJac_…`. Preserve model equations, argument contracts and numerical
results. Rename source files, declarations, handles, code-generation entry
points, test spies and current documentation together. Do not add aliases.

### Stage 1 - Preserve and map

- [x] Record provider and consumer branches, source snapshots and Git indexes.
- [x] Identify 26 EstimationGears and three SimulationGears provider renames.
- [x] Identify active callers, including the Bennu/Apophis validation worktree.
- [x] Exclude independent RCS-1 filters, archived code, frozen outputs and
      unrelated worktrees with their own provider revisions.

### Stage 2 - Rename

- [x] Rename provider files and public functions without changing equations.
- [x] Update active callers, local Jacobian wrappers, test spies and MEX builders.
- [x] Update current naming references and verify no old symbols remain in scope.
- [x] Compare each edited source with its snapshot after reversing the renames.

### Stage 3 - Validate

- [x] Run affected existing MATLAB harnesses and generated-interface checks.
- [x] Verify consumer resolution against the reviewed provider checkouts.
- [x] Record unrelated failures separately; do not repair them in this batch.
- [x] Skip MATLAB execution for consumers without existing MATLAB tests.

### Stage 4 - Review and stage

- [x] Review complete diffs, signatures, documentation, comments and readability.
- [x] Stage only the naming batch and existing EstimationGears SRP corrections.
- [x] Preserve protected SimulationGears and COSMICA indexes and all gitlinks.
- [x] Report remaining provider-integration requirements before committing.

Keep dependency revisions fixed. Updated consumer sources require the renamed
provider sources when integrated; do not claim that old recorded dependencies
supply the new names. No simulation campaign or numerical redesign is required.
Evidence: `/tmp/matlab-eval-case-standardization-20260930-v1il4f89`.

Reviewed by Codex (GPT-6). Preserve existing historical changelog attribution.

#### Validation and integration limits

Resolve all 29 renamed provider symbols. Pass 58 existing provider unit tests,
the existing filter and orbital LUT harnesses, 13 backend unit tests, two
FUTURE dependency tests, the guidance source/delegation/integration harnesses,
and COSMICA truth/approach/interface harnesses. Pass three max-fidelity MEX
contract tests, the complete existing filter code-generation harness (including
48 runtime queries), and guidance MEX/runtime-degree checks. Check 105 edited
MATLAB files with Code Analyzer; introduce no findings relative to the saved
source. Confirm that every source edit is an identifier substitution.

Export exact working sources to an external integration snapshot for these
checks. Compose consumers with the reviewed EstimationGears, SimulationGears
and MathCore sources there; preserve workspace dependency revisions. The first
consumer attempt failed on missing snapshot kernel paths and provider path
ordering. Correct only the temporary snapshot and test runner, then pass all
checks. Retain the initial logs. Keep the wrapper teardown warning separate
from the passing assertions. No unresolved failure caused by the rename.

Skip MATLAB execution for the GTSAM triangulation and MSCKF callers, which
have no matching MATLAB harnesses. Do not hand-edit tracked FUTURE generated
C++ outputs (`cxx/` and `cxx_armv8/`); their old symbols form self-contained
generated implementations. Regenerate deployment artifacts from the updated
MATLAB sources during the consumer release. Skip the FUTURE time-update
script: its unconditional early return prevents exercising its assertions.

Integrate renamed providers and consumers together. Old nested dependency
revisions do not expose the new names, so normal setup against those revisions
remains an integration gate. Rebuild affected generated artifacts after that
update. Do not add casing aliases or move gitlinks in this consolidation pass.

#### Final staged review

- [x] Complete the existing parent-GUI suites with `runtests`: seven unit tests
      pass. Do not count the earlier suite-factory calls as test execution.
- [x] Pass the SSTO discovery scaffold in the external snapshot; keep its
      pre-existing untracked package files outside the index.
- [x] Pass 83 unit tests and eleven existing assertion/code-generation
      harnesses in total; introduce no new tests or configuration-value checks.
- [x] Review every staged blob against the prior index; preserve unrelated
      source hunks, documentation edits and all gitlinks.
- [x] Remove trailing spaces from two renamed declaration/help lines; preserve
      every equation and the remaining legacy formatting.
- [x] Pass whitespace checks for all twelve inspected indexes.

Stage 51 EstimationGears paths, retaining its prior SRP batch. Stage naming
hunks in nav-backend (7), its embedded consumer checkout (5), GUI-System (3),
gui-trajectory-generation (10), FUTURE (8), GTSAM SpaceNav (2) and MSCKF (1).
Leave three pre-existing untracked EstimationGears plans unstaged. Preserve
the protected SimulationGears and both COSMICA indexes; leave their naming
edits unstaged. Keep the untracked SSTO package outside staging.

The validation worktree advanced externally from `81ad4fd0` to `0132cf2`
during this pass. Preserve that commit and its index content; make no commits
here. Verify that all owned source edits remain mechanical substitutions.

The consumer-test MATLAB process remained in wrapper teardown for more than
fourteen minutes after saving passing results. Terminate only that owned
process after checking its command and completed result files. Treat this
exit behavior separately from assertion results; change no wrapper code.
Keep the initial snapshot-setup errors and the suite-discovery-only attempt
in the evidence. The remaining integration gate belongs to provider commits,
dependency updates and generated deployment artifacts, outside this pass.

## SRP Jacobian composition correction - 30 September 2026

Keep model selection and matrix assembly in `EvalJac_InertialPosVelDyn`,
matching the orbital RHS composition. Keep `EvalJac_SRPwithBias` limited to
cannonball partials and use `EvalJac_SRPLutWithBias` for LUT partial mapping.
Preserve shared LUT force/derivative evaluation and both existing signatures.

### Stage 1 - Inspect and preserve

- [x] Inspect the staged batch and both provider boundaries.
- [x] Trace optional bias: exclude absent/disabled states, retain sensitivity
      for configured active or consider states, and preserve nominal modes.
- [x] Save source/index snapshots; preserve unrelated staged batches.

### Stage 2 - Integrate selection

- [x] Remove LUT dispatch and its early return from the cannonball Jacobian.
- [x] Select both SRP models in the orbital Jacobian and share column assembly.
- [x] Clarify optional bias documentation without adding selectors or physics.
- [x] Update existing harnesses for component and composition ownership;
      cover absent and zero bias indexes without mode-vector access.

### Stage 3 - Validate and stage

- [x] Run focused source/cannonball regressions and existing MEX/C++ harnesses.
- [x] Verify considered-state sensitivity, absent-bias shape and LUT reuse.
- [x] Complete documentation, formatting and full staged-diff review.
- [x] Restage only this correction in EstimationGears and preserve other indexes.

Do not repair unrelated algorithm defects or change physical provider kernels.
Reviewer: Codex GPT-6. Evidence: `/tmp/srp-jacobian-model-composition-20260930-9baib3ho`.

Validation: pass thirteen existing cannonball unit tests and four associated
MEX builds. Pass the LUT source harness across twelve unit/mode cases,
sixteen pole cases, sixty-second propagation and the added absent/disabled
bias contracts. Pass seventeen fresh filter MEX/C++ builds and 48 runtime
parity queries, including both units, existing orbital dispatch and a disabled
bias index without a consider-mode vector. Rerun the source harness after
replacing its literal bias column with the fixture's configured index.

Review the sectioned public headers and imperative block comments. Remove
uncommitted changelog entries for the superseded component dispatcher.
Retain two analyzer suppressions only for unused inputs required by existing
callback/component signatures; report no analyzer findings in the five
changed MATLAB files. Preserve the cannonball numerical body and shared
LUT force/partial computation. Do not change provider source or legacy
cannonball consider-state rules.

Stage only the seven-file EstimationGears correction, retaining its 51-path
batch and leaving the unrelated plan prefix unstaged. Preserve the four-file
SimulationGears panel index and nineteen-file COSMICA config index exactly.
Keep the COSMICA plan update unstaged. Validate against exact reviewed source
snapshots; integration with recorded dependency revisions remains separate.
Create no commits or simulation runs.

### SRP acronym spelling - 1 October 2026

- [x] Complete the user replacement with `EvalJac_SRPLutWithBias`,
      `EvalRHS_SRPLutWithBias` and the filter adapter `EvalFilterSRPLutWithBias`.
- [x] Rename files and update active callers, builders, examples and docs;
      preserve equations, signatures, configuration keys and other SRP names.
- [x] Pass existing source/code-generation harnesses and refresh only the
      affected EstimationGears index; preserve the provider and COSMICA indexes.

Reviewer: Codex GPT-6.

Pass both existing source harnesses: twelve filter unit/mode cases, sixteen
pole cases, optional bias/force reuse checks and thirty orbital composition
cases. Pass seventeen fresh MEX/C++ builds and 48 runtime parity queries.
Report no Code Analyzer findings in nine touched MATLAB files. Confirm all
source changes are identifier substitutions plus changelog updates, with
filenames matching their primary functions. Validate composed provider sources;
keep recorded-dependency integration pending.

Refresh eight affected EstimationGears entries, including two file renames;
retain its 51-path batch and leave the unrelated plan prefix unstaged. Keep
SimulationGears changes behind its existing four-file panel index. Preserve
both COSMICA indexes. Create no commits or simulation runs.

Resolve a staging retry by adding `--add` for the renamed index paths.
Preserve source and all unrelated index entries during the retry.

### Filter-module code-generation builder rename - 1 October 2026

- [x] Rename `CodegenSrpLutFilter` to `CodegenSrpLutFilterModules`, including
      its file, declaration, examples, callers and error identifiers.
- [x] Validate the existing generated-interface harness and Code Analyzer.
- [x] Review documentation, imperative comments, formatting and complete index.
- [x] Stage and commit only the builder rename, linked test and documentation;
      preserve unrelated changes and the SimulationGears/COSMICA indexes.

Reviewer: Codex GPT-6. Preserve numerical kernels, argument contracts and
generated model names. Use the existing SRP filter code-generation harness;
introduce no configuration-value tests or simulation runs. The user explicitly
authorized a commit for this correction. Keep dependency integration separate
from validation against the reviewed provider sources.

Pass seventeen fresh MEX/C++ builds and 48 runtime parity queries through
`testSrpLutFilterCodegen`, including both units, optional bias and existing
orbital dispatch. Report zero analyzer findings in the renamed builder and
its caller harness. Confirm that source changes are identifier substitutions,
one declaration-alignment correction and a changelog entry. Preserve generated
model names and output contracts. Keep external validation artifacts under
`/tmp/srp-filter-builder-rename-20261001-sgs4src5`.


### Transverse SRP table correction - 1 October 2026

- [x] Move only `bIncludeTransverseSrp` to constant filter configuration; keep
      consider-state flags and numerical inputs mutable.
- [x] Pass inclusion separately through the existing filter/shared orbital
      interfaces; consume SimulationGears nodal transverse samples.
- [x] Validate both scalar/transverse source and generated targets; review and
      stage scoped source, test and documentation updates.

Reviewer: Codex GPT-6. Pass twelve source unit/mode/bias cases, sixteen pole
cases, the existing RK/STM checks and 34 fresh filter/dispatch MEX/C++ builds.
Preserve optional-bias layouts, units, output pruning and runtime consider flags.
Keep scalar generated targets free of vector storage. Report zero analyzer
findings in the changed MATLAB files. Record complete commands, refinement
metrics, resolved Coder discrepancies and storage checks in the SimulationGears
worktree consolidation plan, section "Correct transverse SRP tables"; retain
evidence under `/tmp/srp-transverse-correction-20261001-0160gjc4`.

Stage only this correction and preserve unrelated plans and existing index
entries. Keep recorded-dependency integration separate from composed-source
validation. Run no mission simulation and create no commit.


### Selected SRP acceleration diagnostics - 1 October 2026

- [x] Rename the LEO acceleration record to `dAccSRP`; preserve its force law.
- [x] Pass all eight existing LEO/third-body unit tests and the independent SRP
      diagnostic check; introduce no MATLAB analyzer findings.
- [x] Review and stage only the source/header change and this progress note.

Reviewer: Codex GPT-6. Record cross-repo consumers and generated validation in
the SimulationGears worktree consolidation plan's "Unify selected SRP
acceleration diagnostics" section. Keep evidence under
`/tmp/srp-diagnostics-unification-20261001-kgidqwta`; preserve unrelated work
and dependency pointers. Create no commit.

### Review EstimationGears after SimulationGears consolidation - 1 October 2026

- [x] Compare all twelve staged files with the qualified SRP batch and verify
      the SRP provider against SimulationGears commit `361076a`.
- [x] Review the complete index for documentation, imperative comments,
      formatting, readability and unnecessary duplication. Align continuation
      arguments and specify transverse inclusion in the documented example.
- [x] Rerun the filter source harness against the committed provider: pass twelve
      unit/mode/bias cases, sixteen pole cases and RK/STM checks. Pass all eight
      LEO/third-body tests and the independent selected-SRP diagnostic check.
- [x] Report zero Code Analyzer findings in ten staged MATLAB files and clean
      cached whitespace. Reuse the 34 generated filter/dispatch checks after
      verifying that only comments and whitespace changed executable sources.
- [ ] Complete user review and commit the twelve-file SRP batch.

Reviewer: Codex GPT-6. Keep evidence under
`/tmp/estimationgears-srp-consolidation-20261001-bt4i9_qn`. Use the reviewed
SimulationGears worktree for these checks; keep dependency-pointer integration
separate. Preserve unrelated plans, sibling indexes and dependency pointers.
Run no mission simulation and create no commit.
