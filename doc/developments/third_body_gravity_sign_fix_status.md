# Third-Body Gravity Sign Fix Status

**Status:** EstimationGears-owned MATLAB correction implemented and focused-verified in the working tree; the
SimulationGears dependency update, commits, releases, and downstream migrations remain pending.

This status file records only the independently stageable B1 third-body slice. The broader
`third_body_gravity_sign_bug_report_and_fix_plan.md` also covers process-noise PSD semantics and remains a separate
multi-stage planning artifact.

## EstimationGears-owned correction

- [x] Correct the generic-body and Sun differential-gravity signs in `evalRHS_DynLEO`.
- [x] Correct `evalJAC_3rdBodyGrav` to evaluate
  `mu*(3*u*u'/|u|^5-I/|u|^3)` for `u = r_body-r_sc`.
- [x] Document the main-body-relative inertial-frame convention in the RHS and Jacobian headers.
- [x] Add independent axial and arbitrary three-dimensional RHS oracles.
- [x] Exercise the separate Sun path and configured Sun-plus-body accumulation.
- [x] Add direct axial and arbitrary three-dimensional analytic Jacobian oracles.
- [x] Retain an independent central finite-difference Jacobian check as secondary consistency coverage.

## Verification evidence

- [x] The focused `testThirdBodyGravityConsistency` suite passes 7/7 tests using the standalone EstimationGears
  RHS and Jacobian paths resolved by `SetupPaths_EstimationGears`.
- [ ] Run the wider affected EstimationGears regression suites when the mixed worktree is separated into reviewable
  branches; unrelated adaptive-tuning, STM, PSD, SRIF, and dependency changes are already present in this checkout.

## Dependency and consumer boundary

- [ ] Update `lib/SimulationGears_for_SpaceNav` only after the standalone source-owner correction has a reviewed
  commit; do not stage the currently dirty nested checkout as a replacement for that dependency revision.
- [ ] Publish compatible SimulationGears and EstimationGears revisions under separate release authorization.
- [ ] Keep FUTURE-nav source and generated C++ unchanged until those reviewed upstream revisions exist.
