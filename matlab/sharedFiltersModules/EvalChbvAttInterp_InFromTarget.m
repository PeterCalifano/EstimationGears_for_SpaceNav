function dTargetDCM_INfromTF = EvalChbvAttInterp_InFromTarget(dTimestamp, strAttData) %#codegen
%% SIGNATURE
% dTargetDCM_INfromTF = EvalChbvAttInterp_InFromTarget(dTimestamp, strAttData)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate nominal target attitude from scalar-first passive quaternion ephemerides.
% Return the TF-to-IN rotation without applying estimated target-attitude bias.
% MATLAB calls validate the data layout and timestamp. Generated code uses a fixed
% struct schema; the MathCore interpolator retains its MEX timestamp assertion.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dTimestamp    Evaluation epoch in the same time units as the ephemeris bounds.
% strAttData    ui32PolyDeg; stacked dChbvPolycoeffs (degree+1 entries per component);
%               dsignSwitchIntervals; dTimeLowBound; dTimeUpBound.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dTargetDCM_INfromTF    Rotation mapping target-fixed vectors into inertial axes.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Move the backend evaluator into EstimationGears.
% 09-09-2026  Pietro Califano, Codex gpt-6    Keep dynamic schema validation on the MATLAB path.
% 09-09-2026  Pietro Califano, Codex gpt-6    Derive polynomial capacity from fixed coefficient storage.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalAttQuatChbvPolyWithCoeffs, Quat2DCM (MathCore).
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dTimestamp (1,1) double
    strAttData (1,1) struct
end
arguments (Output)
    dTargetDCM_INfromTF (3,3) double
end

% Generated code has a fixed input schema and must not compile dynamic field checks.
if coder.target('MATLAB')
    ValidateAttData_(strAttData, dTimestamp);
end

% Four equally sized coefficient blocks determine the compiled polynomial degree.
% Keep the metadata check explicit so a mismatched degree cannot be ignored.
ui32PolyDegree = coder.const(uint32(numel(strAttData.dChbvPolycoeffs) / 4 - 1));
assert(strAttData.ui32PolyDeg == ui32PolyDegree, ...
    'Attitude polynomial degree does not match coefficient storage.');

dQuaternion = evalAttQuatChbvPolyWithCoeffs(ui32PolyDegree, uint32(4), ...
    dTimestamp, strAttData.dChbvPolycoeffs, strAttData.dsignSwitchIntervals, ...
    strAttData.dTimeLowBound, strAttData.dTimeUpBound);

dTargetDCM_INfromTF = Quat2DCM(dQuaternion, false);
end

%% Local helper
function ValidateAttData_(strAttData, dTimestamp)
% Check the struct before evaluating its component-wise coefficient blocks.
cellRequiredFields = {'ui32PolyDeg', 'dChbvPolycoeffs', 'dsignSwitchIntervals', ...
    'dTimeLowBound', 'dTimeUpBound'};

bPresentFields = isfield(strAttData, cellRequiredFields);
if ~all(bPresentFields)
    cellMissingFields = cellRequiredFields(~bPresentFields);
    error('EvalTarget:MissingField', 'strAttData.%s is required.', cellMissingFields{1});
end

validateattributes(strAttData.ui32PolyDeg, {'numeric','uint32'}, ...
    {'scalar','integer','>=',0}, mfilename, 'strAttData.ui32PolyDeg');

dCoeffsPerComponent = double(strAttData.ui32PolyDeg) + 1;

validateattributes(strAttData.dChbvPolycoeffs, {'numeric'}, ...
    {'real','2d','nrows',4*dCoeffsPerComponent,'ncols',1}, ...
    mfilename, 'strAttData.dChbvPolycoeffs');

% An out-of-domain request must fail rather than silently clamp the epoch.
validateattributes(strAttData.dTimeLowBound, {'numeric'}, {'scalar'}, ...
    mfilename, 'strAttData.dTimeLowBound');
validateattributes(strAttData.dTimeUpBound, {'numeric'}, ...
    {'scalar','>=',strAttData.dTimeLowBound}, mfilename, 'strAttData.dTimeUpBound');

if dTimestamp < strAttData.dTimeLowBound || dTimestamp > strAttData.dTimeUpBound
    error('EvalTarget:TimeOutOfRange', ...
        'dTimestamp (%.3f) must be within [%g, %g].', ...
        dTimestamp, strAttData.dTimeLowBound, strAttData.dTimeUpBound);
end
end
