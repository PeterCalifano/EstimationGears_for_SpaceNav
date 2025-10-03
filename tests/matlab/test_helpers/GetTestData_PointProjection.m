function [dKcam, dxStatePost, dDCM_CiFromIN, strMeasModelParams, ...
            strFilterConstConfig, strFilterMutabConfig] = GetTestData_PointProjection(dTargetPosition_IN, ...
                                                                                      dDCM_CamFromSCB)
arguments
    dTargetPosition_IN (3,1) double = zeros(3,1);
    dDCM_CamFromSCB (3,3) double    = eye(3);
end


% Mock camera intrinsics (fx, fy, cx, cy)
dKcam = [5300, 0, 1024; ...
    0, 5300, 768; ...
    0,   0,   1];

% Mock attitude/position states
strFilterConstConfig = struct();
strFilterConstConfig.ui16StateSize = 15;
strFilterConstConfig.strStatesIdx.ui8posVelIdx = 1:6;

dxStatePost = zeros(strFilterConstConfig.ui16StateSize, 1);
dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)) = [75; -60; 1200];

% Dummy Sun position
dSunPosition_IN = zeros(3,1);
dSunPosition_IN(2) = 1E7;
dSunPosition_IN(3) = 1E6 * dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(3));

objAttPointGenerator = CAttitudePointingGenerator(dxStatePost(strFilterConstConfig.strStatesIdx.ui8posVelIdx(1:3)), ...
                                                                dTargetPosition_IN, ...
                                                                dSunPosition_IN);

% Mock direction cosine matrices
[objAttPointGenerator, dDCM_INfromCi] = objAttPointGenerator.pointToTarget("enumConstraintType", "YorthogonalSun"); %#ok<ASGLU>
dDCM_CiFromIN = transpose(dDCM_INfromCi);

strMeasModelParams = struct();
strMeasModelParams.dDCM_SCBiFromIN =  transpose(dDCM_CamFromSCB) * dDCM_CiFromIN; % Generate attitude of spacecraft from camera one

% Mock mutable configuration
strFilterMutabConfig = struct();
strFilterMutabConfig.dDCM_CamFromSCB = dDCM_CamFromSCB;

end
