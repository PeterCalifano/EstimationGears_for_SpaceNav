# Filter Tailoring Package
This package contains the manually tailored entrypoints used to adapt EstimationGears and
SimulationGears filter implementations to a specific mission or measurement stack.

The functions expected to be edited on a case-by-case basis live here:
`BuildArchitectureTemplate`, `BuildInputStructsTemplate`, `ComputeMeasPredAndObsJacobian`,
`ComputeMeasResiduals`, and `ComputeProcessNoiseCov`.

EKF/SRIF-style updates that need both measurement prediction and observation Jacobian should use
`ComputeMeasPredAndObsJacobian` with three outputs so both products come from the same measurement-model
branch. UKF/SR-UKF updates should call the same hook with prediction outputs only because they do not
consume a Jacobian.

Measurement hooks receive both `strFilterMutabConfig` and `strFilterConstConfig`. They must use the
EKF-style `strFilterConstConfig.strStatesIdx` state-index contract, with only fields already present in
`BuildArchitectureTemplate`; do not add ad-hoc observed-state index fields to `strMeasModelParams` or
new state-index fields to `strStatesIdx` in the shared template layer.

Generic shared functions do not live here. `ComputeDynFcn`, `ComputeDynMatrix`, `PropagateDyn`, and
`ManageMeasLatency` are shared/public EstimationGears entrypoints and stay outside this package.

This separation keeps the manually tailored layer explicit without changing the filter mechanization code,
while remaining compatible with MATLAB Coder and Simulink.
