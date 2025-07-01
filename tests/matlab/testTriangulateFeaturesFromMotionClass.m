%% Class-based tests: testTriangulateFeaturesFromMotionClass.m
classdef testTriangulateFeaturesFromMotionClass < matlab.unittest.TestCase
    properties (TestParameter)
        % You can parameterize class here if desired
    end
    properties (Access=private)
        fixture CTestFixture_TriangulateFeaturesFromMotion
    end
    methods (TestClassSetup)
        function attachFixture(testCase)
            testCase.fixture = testCase.applyFixture(CTestFixture_TriangulateFeaturesFromMotion(185.0, 0.05, uint32(1)));
        end
    end

    methods (Test)
        function testSingleFeatureClass(testCase)
            
            fix = testCase.fixture;
            testCase.verifyNotEmpty(fix.dLandmarkPosGT_TB);
            ui32Nfeatures = size(fix.dLandmarkPosGT_TB,2);
            ui32Nposes    = numel(fix.dTimegridID);
            idAnchor      = 1;

            dFeatInvGuess = transformEPtoIDP(fix.dRefAttitudes_TBfromCAM(:,:,idAnchor)' * ...
                (fix.dLandmarkPosGT_TB(:,1) - fix.dPosRefStates_TB(idAnchor,1:3)') + 100 * randn(3,1));

            dyPix  = reshape(fix.dMeasKeypoints_uv(:,1:ui32Nposes),2,[]);
            dyNorm = transformPixelsToNormCoords(dyPix, fix.dKcam, ui32Nposes, ui32Nposes+1);
            dyNorm = reshape(dyNorm(:,1:ui32Nposes),[],1);

            [dRelPos, dDCM_CiFromCk, dPosCk, dDCM_NavFromCk] = ComputeCamRelPoses( fix.dRefAttitudes_TBfromCAM, fix.dPosRefStates_TB', idAnchor, ui32Nposes, ui32Nposes+1);
            [dFeatNav,~,~,~] = Triangulate3dPointInvDepth_GaussNewton(dyNorm, ...
                                                        fix.dMeasCovSigma, ...
                                                        dDCM_CiFromCk, ...
                                                        dRelPos, ...
                                                        dDCM_NavFromCk, ...
                                                        dPosCk, ...
                                                        dFeatInvGuess, ...
                                                        fix.dPrincipalPoint_UV, ...
                                                        ui32Nposes, ...
                                                        uint32(1), ...
                                                        uint32(5), ...
                                                        1e-3, ...
                                                        1e-12, ...
                                                        ui32Nposes+1, ...
                                                        ui32Nfeatures+1);

            dTol = double(max(eps('single'), 4 * fix.dMeasCovSigma_GT));
            testCase.verifyEqual(dFeatNav(:,1), fix.dLandmarkPosGT_TB(:,1), 'AbsTol', dTol);
        end
    end
end
