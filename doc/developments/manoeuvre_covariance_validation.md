# Manoeuvre covariance models and validation

`ComputeManoeuvreInputNoise` now implements the four-source Gates model and has separate
validation oracles for the nonlinear polar-angle and linear Gaussian models. Existing calls
retain their model and zero fixed-error defaults. The navigation backend still selects `MAG_DIR_THR`.

## Model contracts

| Model | Error convention | Validation |
|---|---|---|
| `MAG_DIR_THR` | Gaussian fractional magnitude and polar angle, uniform azimuth; thrust along TH +X | Nonlinear samples, covariance and mean |
| `HERA_GNC` | Existing GMV approximation; thrust along TH +X | Numerical reference, limiting cases, frame and attitude checks |
| `MAG_DIR_DIRECT` | Gaussian proportional magnitude and independent per-axis small-angle errors | Linear Gaussian samples with an arbitrary command |
| `GATES` | DIRECT proportional errors plus fixed axial and per-transverse-axis errors | Four independent Gaussian error sources, references and limiting cases |

The first covariance output is in W, as defined by `dDCM_WfromSC`. Only this output includes the
optional spacecraft attitude contribution. The second output is execution covariance in TH.
Attitude covariance is interpreted as right-local uncertainty in SC, in radians squared.

HERA's existing magnitude-only variance is half that of DIRECT and THR. The reference tests
preserve this coefficient; they do not establish its calibration or derive it from a physical
sampler. The prior shared Monte Carlo test did not resolve these differences between models.

## Gates implementation

For a nonzero command `v`, write `u = v / norm(v)`. Conditional on that command, the execution
covariance is

```text
Q = (sigmaFixedMagnitude^2 + norm(v)^2 * sigmaMagnitudeFraction^2) * u*u'
  + (sigmaFixedPointing^2 + norm(v)^2 * sigmaAngle^2) * (I - u*u').
```

The four independent sources are shutoff, resolution, pointing, and autopilot errors, following
sections II-III, pp. 2-3 of [Gates, JPL TR 32-504 (1963)](https://ntrs.nasa.gov/citations/19640003365).
This is a linear execution-error model with zero mean, conditional on a deterministic command.
It does not integrate over uncertainty in the command itself.

- `dSigmaFixedMagnitudeDV` and `dSigmaFixedPointingDV` are optional ninth and tenth arguments,
  in the same velocity units as the command; both default to zero
- Fixed pointing sigma is per transverse axis, not a polar-angle sigma
- With zero fixed errors, GATES reproduces MAG_DIR_DIRECT covariance
- A zero command with nonzero fixed errors is rejected because its direction is undefined
- Other models reject nonzero fixed-error arguments rather than ignoring them
- `bUseAveragePerturbDeltaV` leaves the Gates command unchanged; its additive errors have zero mean

Example from the repository root, after adding the runtime MATLAB paths:

```matlab
[dCov_W, dCov_TH, dMean_W] = ComputeManoeuvreInputNoise([2; 0; 0], ...
    0.1, 0.2, eye(3), eye(3), zeros(3), EnumManCovModel.GATES, true, 0.3, 0.4);
% dCov_W = dCov_TH = diag([0.13, 0.32, 0.32]); dMean_W = [2; 0; 0].
```

The nonlinear polar-angle sampler follows the random-vector construction in
[Laurens et al. (2021), sections 2 and 4](https://conference.sdo.esoc.esa.int/proceedings/sdc8/paper/121/SDC8-paper121.pdf).
Uniform azimuth on `[0, 2*pi)` and a symmetric Gaussian polar angle give the same distribution as
the paper's half-circle azimuth convention. This sampler is not used for DIRECT or GATES.

## Test changes and evidence

- [x] Replace the old direction-as-rotation-vector sampler with direct polar-angle samples
- [x] Replace the `1e-4` absolute covariance tolerance with checks on normalized covariance entries
- [x] Cover burns scaled by 0.01, 1, and 10; noncommuting rotations; covariance cross terms; and means
- [x] Check correlated attitude uncertainty with rotation finite differences and nonlinear samples
- [x] Cover all models' zero-noise, zero-burn, magnitude-only, frame and scaling behavior
- [x] Cover Gates fixed terms, reduction to DIRECT, unchanged mean, input errors and eight-argument calls
- [x] Reject eight deliberately introduced defects in disposable copies
- [x] Upgrade `EnumManCovModel` to uint8 with codes 0-3 in existing declaration order
- [x] Generate and execute Gates and legacy eight-argument THR MEX functions; match MATLAB outputs
- [x] Run Code Analyzer on the function, enum, and test file: zero findings
- [ ] Validate HERA's physical error convention against its source and operational calibration
- [ ] Validate actual campaign command/thruster alignment and combined uncertainty calibration

The sampled covariance is normalized by the analytical Cholesky factor and compared with identity.
The maximum entry error must be below 0.02, with 300,000-400,000 fixed-seed samples. There is no
absolute covariance floor that can admit zero covariance. The rank-two attitude-only case is also
checked in its observable tangent plane. Local random streams preserve the caller's RNG state.

Validation on 2026-09-08:

- 18 manoeuvre tests plus two neighboring process-noise tests passed, MATLAB exit 0
- 15 isolated backend Delta-V/manoeuvre compatibility tests passed, MATLAB exit 0, with the
  existing wrapper teardown warning
- Eight deliberate defects were rejected: zero covariance, missing frame rotation, omitted
  attitude noise, half DIRECT direction variance, wrong HERA magnitude coefficient, omitted
  Gates fixed magnitude error, Gates mean shrinking, and the former direction sampler
- MATLAB Coder generated both MEX variants; their covariance outputs matched MATLAB within
  `1e-14` absolute Frobenius error, with matching command outputs

Earlier failed attempts are retained: the original GATES rejection had no error identifier;
Coder first rejected the plain enum, then the redundant string-array validator. The enum now
inherits uint8, and the typed enum argument supplies validation without that redundant validator.
The Gates tests were run against the unimplemented interface before adding the implementation.

Scripts, logs, source snapshots for deliberate defects, and generated MEX files are under
`/tmp/manoeuvre_cov_validation_20260908/`. No closed-loop campaign was run. Existing serialized
MATLAB objects were not used to validate the enum representation change.
