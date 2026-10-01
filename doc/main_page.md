# EstimationGears_for_SpaceNav {#mainpage}

EstimationGears provides MATLAB estimation algorithms and a C++20/CUDA library scaffold for spacecraft-navigation applications.

## Entry points

- [Repository overview](../README.md)
- [C++ and CUDA builds](cpp_cuda_build.md)
- [Build script](build_script_doc.md)
- [Python and MATLAB wrappers](wrappers.md)
- [Testing and CI](testing_and_ci.md)
- [Versioning](versioning.md)
- [Logging](logging.md)
- [Documentation workflow](documentation_workflow.md)

## Installed CMake API

```cmake
find_package(EstimationGears_for_SpaceNav REQUIRED CONFIG)
target_link_libraries(
    application
    PRIVATE EstimationGears_for_SpaceNav::EstimationGears_for_SpaceNav)
```

The native build is independent of MATLAB-only source/runtime checkouts under `lib/`. CUDA, Python wrapping, and MATLAB wrapping are optional.

## Optional SRP response-table model

Select the response table with constant `strFilterConstConfig.bUseSrpLut=true`
and supply the validated numeric payload in `strResponseLut`. The existing
`EvalFilterDynOrbit` and `EvalFilterDynOrbit_FixedEph` prepare numerical inputs
with `BuildFilterSrpLutInputs` and forward them to SimulationGears
`EvalRHS_InertialDynOrbit` for force selection and composition. Use
`EvalFilterSRPLutWithBias` for direct force/partial requests and
`EvalJac_SRPLutWithBias` for filter Jacobian mapping. The existing
`EvalJac_InertialPosVelDyn` selects the SRP component and assembles its state
columns for either model. Keep `EvalJac_SRPwithBias` specific to cannonball SRP.
An absent or false selector retains the existing
cannonball model. Set runtime `strFilterMutabConfig.bIncludeTransverseSrp`
to include the table's transverse response.

SimulationGears owns table preparation, interpolation and physical SRP
evaluation, unit conversion, additive-bias physics and orbital force composition.
EstimationGears owns state-index/consider-mode input preparation and Jacobian
mapping. Gravity and residual acceleration continue through their existing
dynamics paths. See the SimulationGears
[response-table contract](../lib/SimulationGears_for_SpaceNav/doc/main_page.md#spacecraft-srp-response-tables)
for payload fields, interpolation and derivative conventions.

### Runtime inputs and units

Supply resolved Sun-first `strDynParams.dBodyEphemerides`, the eclipse flag,
reference pressure at 1 AU in `strSRPdata.dP_SRP0`, spacecraft mass in
`strSCdata.dSCmass` and explicit `strSrpPointing` fields:

| Field | Convention and units |
| --- | --- |
| `dDcm` | Body-to-inertial proper rotation, 3 by 3 |
| `dAttitudePositionPartials` | Three `dR/dr_j` slices, 3 by 3 by 3, per length unit |
| `dIllumination` | External illumination fraction in [0, 1] |
| `dIlluminationGradient` | Position gradient, 1 by 3, per length unit |

Use consistent metre or kilometre filter units (LU). Pressure is in
kg/(LU s²); table geometry and area remain in metres and square metres.
The shared physical kernel converts once to the SI evaluator and returns acceleration
in LU/s² and its position partial in 1/s². Hold the supplied attitude fixed
unless its position partial is provided. Velocity partials are zero.

Retain the optional `strStatesIdx.ui8CoeffSRPidx` as an additive acceleration
bias in LU/s². Its nominal contribution is illumination times the bias times
the unit Sun-to-spacecraft direction. Include Sun-line rotation and the
illumination gradient in position partials. Eclipse, zero pressure or
unavailable Sun geometry suppress nominal force and all partials.

For a consider bias, exclude its stored value from nominal acceleration and
position partials while retaining its sensitivity for uncertainty propagation.
That sensitivity describes the uncertain parameter; it is not a derivative
with respect to the ignored nominal state value.

### Shared force and Jacobian evaluation

```matlab
strConstant.bUseSrpLut = true;
strConstant.strResponseLut = strPreparedLut;
% Request mapped partials and force together for identical inputs.
[dSrpJacobian, dAcceleration] = EvalJac_SRPLutWithBias(dxState, strDynParams, ...
                                                      strMutable, strConstant);
```

The Jacobian has six orbital rows and three position columns, plus a bias
column when mapped. Request the optional acceleration output to obtain both
results from one LUT evaluation. Separate RHS and Jacobian calls evaluate
the table twice. Reuse results only for identical state, dynamics inputs and
configuration. Keep distinct integration stages separate.

### Code generation and validation

Use `CodegenSrpLutFilter` to generate a fixed-size MEX or C++ library with the
table and state mapping embedded in constant configuration. Constant inputs
are removed from MEX signatures. Numerical kernels disable variable sizing
and dynamic allocation; MATLAB gateways still allocate output storage.

| Entry point | Default MEX outputs | Default C++ outputs | Explicit output prefixes |
| --- | --- | --- | --- |
| `EvalFilterSRPLutWithBias` | Acceleration | All four | Acceleration, position partial, bias sensitivity, derivative regularity |
| `EvalJac_SRPLutWithBias` | Mapped Jacobian | Mapped Jacobian | Mapped Jacobian, acceleration |

Set `ui8OutputCount` to the number of leading outputs needed. Request
`uint8(2)` explicitly for the combined Jacobian/acceleration interface.
Keep `bForceOnly=true` as the one-output RHS shorthand. Bias columns are
compiled out when omitted from the state mapping. Transverse selection and
consider mode remain runtime choices.

Run `testSrpLutFilterModels` for units, independent force/partial references,
bias modes, inactive radiation, poles and the existing RK4/STM path. Run
`testSrpLutFilterCodegen(charNewOutputRoot)` for MEX parity, C++ libraries
and compilation of existing orbit dispatch. Add the standalone SimulationGears
source after EstimationGears setup for sibling-repository development; keep
recorded dependency updates separate.

The LUT is piecewise smooth. Knots use one-sided derivatives; adjusted pole
lookups return a finite Jacobian with a false regularity flag. Keep the
physical Sun direction and additive bias unchanged. The filter consumes that
finite Jacobian; the regularity output reports the lookup's smoothness.
These checks establish numerical and interface behavior. Flight-hardware
timing and mission navigation performance require separate validation.

### MATLAB dynamics and Jacobian names

Use `EvalRHS_…` for right-hand-side kernels and `EvalJac_…` for Jacobians.
Match each primary function name to its source filename. Update callers,
function handles and generated entry points with the provider; rebuild MEX
artifacts compiled from these entry points. The previous case spellings have
no forwarding aliases.
