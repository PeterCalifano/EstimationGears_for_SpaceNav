# Filter Tailoring Package
This package contains the manually tailored entrypoints used to adapt EstimationGears and
SimulationGears filter implementations to a specific mission or measurement stack.

The functions expected to be edited on a case-by-case basis live here:
`BuildArchitectureTemplate`, `BuildInputStructsTemplate`, `ComputeMeasPred`, `ComputeObsMatrix`,
`ComputeMeasResiduals`, and `ComputeProcessNoiseCov`.

Generic shared functions do not live here. `ComputeDynFcn`, `ComputeDynMatrix`, `PropagateDyn`, and
`ManageMeasLatency` are shared/public EstimationGears entrypoints and stay outside this package.

This separation keeps the manually tailored layer explicit without changing the filter mechanization code,
while remaining compatible with MATLAB Coder and Simulink.
