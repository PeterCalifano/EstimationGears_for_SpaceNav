function strScenario = BuildFullCovObservationTestProblem()
%% SIGNATURE
% strScenario = BuildFullCovObservationTestProblem()
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Build an ellipsoidal LiDAR update with nonzero target bias and one camera clone.
% The prior includes current / clone correlations. Tests can extend its allocated window
% or activate centroid and direction measurements while retaining the same geometry.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% strScenario    State, covariance, timestamps, configuration and sensor inputs.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Share the full-covariance observation fixture.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% filter_tailoring.BuildArchitectureTemplate, filter_tailoring.BuildInputStructsTemplate,
% RotationVectorToDCM, DCM2quat.
% -------------------------------------------------------------------------------------------------------------
arguments (Output)
    strScenario (1, 1) struct
end

strConstant = filter_tailoring.BuildArchitectureTemplate('bWriteBusDefs', false);
[strMutable, strDynamics, strModel, strMeasurements] = ...
    filter_tailoring.BuildInputStructsTemplate(strConstant);
strConstant.bEstimateGravParam = false;
strMutable.ui16WindowStateCounter = uint16(1);
strMutable.bEnableEditing = false;
strMutable.bConsiderStatesMode(:) = false;
strMutable.bConsiderStatesMode(15:16) = true;
strMutable.dMeasUnderweightCoeff = 0;
strMutable.dRangeLidarSigma = 0.2;
strMutable.ui8LidarShapeModelMode = uint8(2);
strMutable.dEllipsoidInvDiagShapeCoeffs = [1/9;1/4;1/1.44];
strMutable.bEnableLidarFallbackPrediction = false;

% Use nontrivial spacecraft and target attitudes and a displaced ray origin.
dBeam_IN = [-1;0.2;-0.1];
dBeam_IN = dBeam_IN / norm(dBeam_IN);
dSpacecraftRotation = RotationVectorToDCM([0.2;-0.3;0.1]);
strModel.dDCM_SCBiFromIN(:, :, :) = repmat(dSpacecraftRotation, ...
    1, 1, size(strModel.dDCM_SCBiFromIN, 3));
strMutable.dLidarBeamDirection_SCB = dSpacecraftRotation * dBeam_IN;
dNominalRotation = RotationVectorToDCM([0.1;0.05;-0.15]);
dQuaternion = DCM2quat(dNominalRotation', false);
strDynamics.strMainData.strAttData.dChbvPolycoeffs(:) = 0;
strDynamics.strMainData.strAttData.dChbvPolycoeffs(1:3:12) = dQuaternion;

ui32StateSize = uint32(strConstant.ui16StateSize);
dxState = zeros(ui32StateSize + 7, 1);
dxState(1:3) = [4;-1;0.5];
dxState(7:9) = [0.12;-0.08;0.04];
dxState(14) = 0.03;
dxState(ui32StateSize + uint32(1:7)) = [5;-1;0.5;1;0;0;0];
ui32CovSize = ui32StateSize + 6;
dFactor = diag(linspace(0.2, 0.5, ui32CovSize)) + ...
    0.002 * reshape(sin(1:double(ui32CovSize)^2), ui32CovSize, ui32CovSize);
strMeasurements.bMeasTypeFlags(:) = false;
strMeasurements.bMeasTypeFlags(3) = true;
strMeasurements.dRangeLidarCentroid(1) = 1.1;

strScenario = struct('dxState', dxState, 'dCovariance', dFactor * dFactor', ...
    'dTimestamps', [0;-0.2], 'strConstant', strConstant, 'strMutable', strMutable, ...
    'strDynamics', strDynamics, 'strModel', strModel, 'strMeasurements', strMeasurements, ...
    'dNominalRotation', dNominalRotation, 'dBeam_IN', dBeam_IN);
end
