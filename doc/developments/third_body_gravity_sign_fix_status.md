# Third-Body Gravity Sign Fix and Release Evidence

**Status:** The EstimationGears B1 correction is implemented, independently verified, generated for both FUTX C++
targets, and integrated into `future-onboard-sw`. Version selection, tagging, and release publication remain
user-owned and pending.

This report covers only B1: the sign of main-body-relative differential third-body acceleration and its
spacecraft-position Jacobian. B6 process-noise PSD semantics, simulator changes, and general SimulationGears work
remain separate.

## Finding and corrected model

Let `r` be spacecraft position and `s` the perturbing-body position. Both are measured from the main body in one
inertial frame. The physical acceleration of the spacecraft relative to the main body is

```text
a_3b(r,s) = mu * ((s-r)/|s-r|^3 - s/|s|^3).
```

The first term is the direct attraction of the spacecraft by the third body. The second removes the acceleration
of the main-body-centred frame origin. Before the correction, the implementation evaluated the exact negative:

```text
a_legacy = mu * ((r-s)/|r-s|^3 + s/|s|^3) = -a_3b.
```

With `u=s-r`, the derivative with respect to spacecraft position is

```text
da_3b/dr = mu * (3*u*u'/|u|^5 - I/|u|^3).
```

The indirect term has no spacecraft-position derivative. The legacy Jacobian was the derivative of the legacy
RHS, so ordinary finite differences of production RHS against production Jacobian could not reveal the shared
physical sign error.

## External cross-check

The corrected convention agrees with independent published implementations and references:

- JPL's modified-equinoctial-elements note defines `d=r-s` and gives the secondary-body term as
  `-mu*(d/|d|^3+s/|s|^3)`, algebraically identical to the corrected implementation.
- Orekit's `ThirdBodyAttraction` constructs `satToBody = centralToBody - spacecraftPosition` and returns positive
  direct attraction along `satToBody` plus negative indirect acceleration along `centralToBody`.
- JPL Publication 78-40 describes N-body perturbation as the direct spacecraft acceleration combined with the
  opposite indirect acceleration of the central body.
- The same direct-minus-indirect construction is the standard model described in Vallado, *Fundamentals of
  Astrodynamics and Applications*, 4th edition, and Battin, *An Introduction to the Mathematics and Methods of
  Astrodynamics*, revised edition.

Primary online references:

- [JPL modified equinoctial elements](https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf)
- [Orekit ThirdBodyAttraction source](https://www.orekit.org/static/xref/org/orekit/forces/gravity/ThirdBodyAttraction.html)
- [JPL Publication 78-40](https://ntrs.nasa.gov/api/citations/19780019213/downloads/19780019213.pdf)

## Source ownership and revision chain

| Boundary | Revision | Evidence |
| --- | --- | --- |
| EstimationGears source correction | `db9ee029217bb691eb75716ed159567425cf008d` | RHS, Jacobian, seven independent tests, and initial status report |
| EstimationGears revision consumed by FUTX | `8110959088e53a58c15806b9eca5648f02c20465` | Contains the B1 parent plus repository-property and nested-gitlink updates |
| future-nav dependency/validation | `1cb1ce2` and preceding `50a8d68` | Frozen 7200 s filter case, physical oracle, active-path proof, updated gitlink |
| future-nav generated C++ | `aa18419dbea19e48257a789a42a48e4155a2b81a` | Four generated host/ARM files |
| future-onboard-sw integration | `97c8c0d317803eebda5db6812414c69d8a5a7a97` | The same four generated files transferred by the standard export script |
| future-onboard-sw release naming | `b88b026cc90dcf204b9ca015f2153271f54b6009` | Restores final generation timestamp and removes the legacy `hil` name token |

The active FUTX path is:

```text
ComputeDynFcn -> filterDynLEO -> evalRHS_DynLEO
ComputeDynMatrix -> evalJAC_DynLEO -> evalJAC_3rdBodyGrav
```

`testThirdBodyGravityActivePath` proves both runtime resolution and static dependency closure select the embedded
EstimationGears checkout. It also proves `ComputeDynFcn` does not depend on SimulationGears'
`evalRHS_InertialDynOrbit`.

## RHS trace

`evalRHS_DynLEO` stores `dPos3rdBodiesToSC = r-s`. It now applies a leading minus to both the stored direct vector
and the main-to-body vector:

```text
-mu * ((r-s)/|r-s|^3 + s/|s|^3).
```

The separate Sun block builds the same two stored-vector terms and multiplies their sum by `-mu_sun`. Generated
`ComputeDynFcn.cpp` preserves these signs in both x86_64 and ARMv8 outputs. The generated host and ARM files are
not expected to be byte-identical to each other because MATLAB Coder emits target-specific outputs; each is instead
checked against its corresponding onboard copy.

## Jacobian and covariance trace

`evalJAC_3rdBodyGrav` stores `dBodyPosToSC=r-s` and accumulates

```text
mu * (-I/|r-s|^3 + 3*(r-s)*(r-s)'/|r-s|^5).
```

This is identical to the `u=s-r` expression because the outer product is invariant under `u=-(r-s)`. The block is
inserted into velocity rows and position columns by `ComputeDynMatrix`. During the full-covariance time update:

```text
ComputeDynMatrix -> getDiscreteTimeSTM -> dDeltaFlowSTM
                 -> P_next = dDeltaFlowSTM * P * dDeltaFlowSTM' + Q
```

The accumulated flow STM is updated with the same step matrix. Therefore, the sign correction affects both the
mean trajectory and covariance mapping, as intended.

## Units and configuration review

- EstimationGears' generic formula is valid for any internally consistent length scale. Its test documentation now
  states that contract without incorrectly pinning the generic fixtures to metres.
- FUTX explicitly sets `bUseKilometersScale=true`; state position is in kilometres, velocity in kilometres per
  second, GM in `km^3/s^2`, acceleration in `km/s^2`, and the acceleration-position Jacobian in `1/s^2`.
- The frozen physical parameters are Earth `398600.44180000003`, Sun `132712440041.27942`, and Moon
  `4904.8695000000007`, all in `km^3/s^2`. No gravitational parameter is scaled for the test.
- The 7200 s fixture retains Earth point-mass plus Sun/Moon differential gravity. J2, drag, SRP, residual
  acceleration, measurement updates, process noise, and non-orbit covariance blocks are disabled to isolate B1.
- B6 process-noise migration is absent. The validation explicitly observes zero integrated process noise.

## Verification evidence

Fresh Stage 5 verification on 07-Aug-2026 produced:

| Layer | Result |
| --- | --- |
| EstimationGears physical/Jacobian suite | 7/7 passed |
| future-nav active dependency-path suite | 2/2 passed |
| future-nav 7200 s filter-performance suite | 1/1 passed |
| future-onboard-sw native CTest | 55/55 passed |
| future-onboard-sw ARMv8 CTest under QEMU | 55/55 passed |
| ARM artifact inspection | `acquisitionWindowTool`, `moonIP`, and `navFilter` are ELF64 AArch64 |

The frozen legacy baseline and corrected candidate used fixture SHA-256
`a6b0833b35475686b00f116352cfb27fa5d56175afe29af135bf94faa0dfe364`:

| Metric | Legacy sign | Corrected sign | Acceptance |
| --- | ---: | ---: | ---: |
| Final position error | `0.00832826883754 km` | `2.46288260342e-10 km` | corrected <= 50% of legacy |
| Final velocity error | `8.93765413186e-06 km/s` | `1.66567238054e-13 km/s` | corrected <= 50% of legacy |
| Orbit covariance trace | `32860.7566606` | `32860.7310551` | relative change <= 5% |
| Covariance-trace relative change | - | `7.79210281967e-07` | pass |
| Covariance symmetry error | - | `1.07416936496e-09` | limit `3.27224197236e-06` |
| Minimum covariance eigenvalue | - | `2.7504690684e-10` | lower limit `-2.45807650256e-08` |
| Integrated process-noise norm | `0` | `0` | required `0` |

The position- and velocity-error ratios are approximately `2.96e-08` and `1.86e-08`, respectively. The corrected
production propagation is therefore materially closer to the independent oracle, not merely internally
RHS/Jacobian-consistent.

## Generated-code provenance and transfer evidence

- Both targets were generated with MATLAB Coder 24.2 / MATLAB R2024b on 07-Aug-2026.
- Host generation timestamp: `17:24:49`; ARMv8 generation timestamp: `17:25:10`.
- The corresponding future-nav and future-onboard-sw files are byte-identical:

| Generated file | SHA-256 |
| --- | --- |
| x86_64 `ComputeDynFcn.cpp` | `8d49767f110272b377a94b3b9678a25a5baad6cdfbd8c43b0d0080781dfdabcf` |
| x86_64 `ComputeDynMatrix.cpp` | `bde132493fcc64dfaa80d8648d3ee8e0c1afa580fb78092d1f77a5ffd0e9f84d` |
| ARMv8 `ComputeDynFcn.cpp` | `522f322a93a1cf8cba8b23a7f1723175b860ecb8c5db83d5aea7c3b51c756461` |
| ARMv8 `ComputeDynMatrix.cpp` | `018503e3c3e7866f541a8700d8aec7314b48c2408938bd209389303300de1595` |

The transfer used future-nav's standard `export_codegen_src.sh`. No generated C++ was hand-edited.

## Scope audit and residual risk

- The EstimationGears B1 commit changes only the two production functions, the independent test, and this report.
- The future-nav generated commit changes only four generated filter files. Its semantic dynamics changes are the
  RHS and Jacobian signs; remaining diff is regenerated comments, source-line metadata, and algebraically identical
  zero-product spelling from the current dependency set.
- The future-onboard-sw integration commit changes only the corresponding four generated files. Simulator,
  SimulationGears implementation, Moon-IP, B6 process-noise, and unrelated filter algorithms do not enter that
  commit.
- The following `b88b026` release-process commit is algorithm-neutral. It changes the source-delivery filename to
  `polimi_apps_delivery_<release-id>_l4t-r36.4.4_aarch64_<YYYYMMDD_HHMM>.tar.gz`, updates CI consumers and release
  documentation, and passes all 21 packaging-contract tests.
- Revision `8110959` does contain a nested SimulationGears gitlink advance and `.vscode` changes after the isolated
  B1 parent. The FUTX runtime and generated-code closure do not consume that nested implementation, so this does not
  alter the onboard algorithm. It remains a repository-provenance caveat if EstimationGears itself is released as a
  source package.
- Per explicit project decision, future-onboard-sw does not duplicate the MATLAB B1 physical oracle. Its 55-test
  native and ARM suites prove integration/build/runtime regression only. Physical correctness remains owned by the
  independent EstimationGears and future-nav tests.
- The 7200 s filter case is a realistic FUTX Earth/Sun/Moon propagation using tracked ephemerides, but it is not a
  full closed-loop measurement-update campaign. It directly exposes the sign defect by more than 8 m over two hours
  and is sufficient for the isolated algorithm gate; later mission-level HIL acceptance remains valuable release
  evidence.

## Release handoff

The B1 source, generated-code, transfer, native, and ARM gates are complete. The remaining actions are deliberately
outside this report's authority:

- choose and record the release version;
- merge through the project workflow;
- create the user-selected tag;
- run the official tag pipeline and inspect all release assets;
- perform Jetson/HIL acceptance and publish the release.

No agent-owned tag, release version, or publication action is authorized.
