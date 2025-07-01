%% testTriangulateFeaturesFromMotionFcn.m
% Function-based tests using the custom fixture
function tests = testTriangulateFeaturesFromMotionFcn
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % arguments
    %     testCase matlab.unittest.TestCase
    % end
    % Apply fixture with default parameters
    fixture = testCase.applyFixture(CTestFixture_TriangulateFeaturesFromMotion(185.0, 0.05, uint32(1)));
    testCase.TestData.fix = fixture;
end

function teardownOnce(~)
    % nothing to clean up
end

function testSingleFeature(testCase)
    % arguments
    %     testCase matlab.unittest.TestCase
    % end
    fixture = testCase.TestData.fix;
    testCase.verifyNotEmpty(fixture.dLandmarkPosGT_TB);
    % Prepare inputs
    ui32Nfeatures = size(fixture.dLandmarkPosGT_TB,2);
    ui32Nposes    = numel(fixture.dTimegridID);
    idAnchor      = 1;
    % Inverse-depth guess in anchor frame (use GT + noise)
    dFeatInvGuess = transformEPtoIDP(fixture.dRefAttitudes_TBfromCAM(:,:,idAnchor)' * (fixture.dLandmarkPosGT_TB(:,1) - fixture.dPosRefStates_TB(idAnchor,1:3)') + 100*randn(3,1));
    % Pack normalized measurements
    dyPix = reshape(fixture.dMeasKeypoints_uv(:,1:ui32Nposes),2,[]);
    dyNorm = transformPixelsToNormCoords(dyPix, fixture.dKcam, ui32Nposes, ui32Nposes+1);
    dyNorm = reshape(dyNorm(:,1:ui32Nposes),[],1);
    % Compute relative poses
    [dRelPos, dDCM_CiFromCk, dPosCk, dDCM_NavFromCk] = ComputeCamRelPoses(...
        fixture.dRefAttitudes_TBfromCAM, fixture.dPosRefStates_TB', idAnchor, ui32Nposes, ui32Nposes+1);
    % Call triangulation
    [dFeatNav,~,~,~] = TriangulateFeaturesFromMotion(dyNorm, fixture.dMeasCovSigma, dDCM_CiFromCk, dRelPos, dDCM_NavFromCk, dPosCk, dFeatInvGuess, fixture.dPrincipalPoint_UV, ui32Nposes, uint32(1), uint32(5), 1e-3, 1e-12, ui32Nposes+1, ui32Nfeatures+1);
    % Verify against ground truth
    dTol = double(max(eps('single'), 4*fixture.dMeasCovSigma_GT));
    testCase.verifyEqual(dFeatNav(:,1), fixture.dLandmarkPosGT_TB(:,1), 'AbsTol', dTol);
end
