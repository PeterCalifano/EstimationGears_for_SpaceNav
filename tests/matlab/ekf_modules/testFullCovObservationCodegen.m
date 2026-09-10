function tests = testFullCovObservationCodegen
%% SIGNATURE
% tests = testFullCovObservationCodegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Compile the full observation update with dynamic allocation and variable sizing disabled.
% Cross-compile standalone C for ARM64 and inspect the linked allocation contract.
% Compare MATLAB and MEX outputs while measurement flags, window count, editing,
% underweighting and consider mode change at runtime. Exercise two compiled capacities.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None. MATLAB Coder and a configured C compiler are required for this acceptance test.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Verify fixed-size generated observation updates.
% 09-09-2026  Pietro Califano, Codex gpt-6    Exercise zero and nonzero camera lever arms.
% 10-09-2026  Pietro Califano, Codex gpt-6    Cross-compile both constant window-frame modes.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% MATLAB Coder, BuildFullCovObservationTestProblem, EKF_SlideWindow_FullCov_ObsUp.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOriginalPath = path;
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','..'));
SetupPaths_EstimationGears;
addpath(fullfile(fileparts(mfilename('fullpath')),'..','test_helpers'));
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testFixedStorageAcrossRuntimeModes(testCase)
assumeTrue(testCase, ~isempty(which('codegen')),'MATLAB Coder is required.');
charBuildRoot = tempname;
mkdir(charBuildRoot);
testCase.addTeardown(@() RemoveBuild_(charBuildRoot));
addpath(charBuildRoot);
objConfig = coder.config('mex');
objConfig.EnableDynamicMemoryAllocation = false;
objConfig.EnableVariableSizing = false;

for bCorrelated = [false, true]
    for ui16Capacity = uint16([3, 13])
        strScenario = BuildFullWindowScenario_(ui16Capacity);
        strScenario.strConstant.ui8RelDirDesign = uint8(bCorrelated);
        strScenario.strConstant.bUseMeasNoiseCrossCov = bCorrelated;
        cellInputs = BuildInputs_(strScenario);
        cellCodegenArgs = cellInputs;
        cellCodegenArgs{8} = coder.Constant(strScenario.strConstant);
        clear FullCovObsFixed_test_mex
        codegen('-config', objConfig,'EKF_SlideWindow_FullCov_ObsUp','-args', cellCodegenArgs, ...
            '-d', fullfile(charBuildRoot,'build'), ...
            '-o', fullfile(charBuildRoot,'FullCovObsFixed_test_mex'));

        % One binary must accept each runtime combination without changing array shapes.
        ui16Counts = uint16([0, 1, ceil(double(ui16Capacity)/2), double(ui16Capacity)]);
        for ui16WindowCount = ui16Counts
            for ui8Mask = uint8(0):uint8(7)
                if ui16WindowCount == 0 && bitget(ui8Mask, 1)
                    continue
                end
                for ui8Policy = uint8(0):uint8(3)
                    strScenario.strMutable.ui16WindowStateCounter = ui16WindowCount;
                    strScenario.strMutable.bEnableEditing = logical(bitget(ui8Policy, 1));
                    bConsider = logical(bitget(ui8Policy, 2));
                    strScenario.strMutable.bConsiderStatesMode(7:9) = bConsider;
                    strScenario.strMutable.dMeasUnderweightCoeff = 0.4*double(bConsider);
                    strScenario.strMutable.dCameraPosition_SCB = [.7;-.3;.2]*double(bConsider);
                    strScenario.strMutable.ui8CenMeasCovModel = uint8(bConsider);
                    strScenario.strMutable.i8CentroidingAlgorithmMode = uint8(bConsider);
                    strScenario.strMeasurements.bMeasTypeFlags = logical(bitget(ui8Mask, uint8(1:3)))';
                    cellInputs = BuildInputs_(strScenario);
                    cellMatlab = cell(1, 10);
                    cellMex = cell(1, 10);
                    [cellMatlab{:}] = EKF_SlideWindow_FullCov_ObsUp(cellInputs{:});
                    [cellMex{:}] = FullCovObsFixed_test_mex(cellInputs{:});
                    for ui32Output = [1, 2, 3, 6, 7, 8, 9, 10]
                        verifyEqual(testCase, cellMex{ui32Output}, cellMatlab{ui32Output}, ...
                            'AbsTol', 2e-9,'RelTol', 2e-12);
                    end
                    verifyEqual(testCase, cellMex{4}, cellMatlab{4});
                    verifyEqual(testCase, cellMex{5}, cellMatlab{5});
                end
            end
        end

        % A received LiDAR sample may fail prediction while the other sensor blocks remain usable.
        strScenario.strMutable.ui16WindowStateCounter = uint16(1);
        strScenario.strMutable.bEnableEditing = false;
        strScenario.strMutable.dLidarBeamDirection_SCB = -strScenario.strMutable.dLidarBeamDirection_SCB;
        for bFallback = [false, true]
            strScenario.strMutable.bEnableLidarFallbackPrediction = bFallback;
            for ui8Mask = uint8([4, 5, 6, 7])
                strScenario.strMeasurements.bMeasTypeFlags = logical(bitget(ui8Mask, uint8(1:3)))';
                VerifyNativeOutputs_(testCase, BuildInputs_(strScenario));
            end
        end
    end
end
end

function testCentroidWithoutBiasUsesStaticStorage(testCase)
assumeTrue(testCase, ~isempty(which('codegen')),'MATLAB Coder is required.');
charBuildRoot = tempname;
mkdir(charBuildRoot);
testCase.addTeardown(@() RemoveBuild_(charBuildRoot));
addpath(charBuildRoot);
objConfig = coder.config('mex');
objConfig.EnableDynamicMemoryAllocation = false;
objConfig.EnableVariableSizing = false;

strScenario = BuildFullWindowScenario_(uint16(3));
strScenario.strMutable.ui16WindowStateCounter = uint16(1);
strScenario.strMutable.i8CentroidingAlgorithmMode = uint8(0);
strScenario.strMeasurements.bMeasTypeFlags = logical([0;1;0]);
strScenario.strConstant.strStatesIdx.ui8CenMeasBiasIdx = uint8([]);
cellCodegenArgs = BuildInputs_(strScenario);
cellCodegenArgs{8} = coder.Constant(strScenario.strConstant);
clear FullCovObsFixed_test_mex
codegen('-config', objConfig,'EKF_SlideWindow_FullCov_ObsUp','-args', cellCodegenArgs, ...
    '-d', fullfile(charBuildRoot,'build'), ...
    '-o', fullfile(charBuildRoot,'FullCovObsFixed_test_mex'));
VerifyNativeOutputs_(testCase, BuildInputs_(strScenario));
end

function VerifyNativeOutputs_(testCase, cellInputs)
cellMatlab = cell(1, 10);
cellMex = cell(1, 10);
[cellMatlab{:}] = EKF_SlideWindow_FullCov_ObsUp(cellInputs{:});
[cellMex{:}] = FullCovObsFixed_test_mex(cellInputs{:});
for ui32Output = [1, 2, 3, 6, 7, 8, 9, 10]
    verifyEqual(testCase, cellMex{ui32Output}, cellMatlab{ui32Output},'AbsTol', 2e-9,'RelTol', 2e-12);
end
verifyEqual(testCase, cellMex{4}, cellMatlab{4});
verifyEqual(testCase, cellMex{5}, cellMatlab{5});
end

function testArm64StandaloneBuild(testCase)
assumeTrue(testCase, ~isempty(which('codegen')),'MATLAB Coder is required.');
[i32CompilerStatus, ~] = system('command -v aarch64-linux-gnu-gcc');
assumeEqual(testCase, i32CompilerStatus, 0,'An AArch64 GCC cross-compiler is required.');
for enumFrame = [EnumWindowRefFrame.INERTIAL, EnumWindowRefFrame.TARGET_FIXED]
    for bCorrelated = [false, true]
        charBuildRoot = tempname;
        mkdir(charBuildRoot);
        testCase.addTeardown(@() rmdir(charBuildRoot,'s'));
        strScenario = BuildFullWindowScenario_(uint16(13));
        strScenario.strConstant.enumWindowRefFrame = enumFrame;
        strScenario.strConstant.ui8RelDirDesign = uint8(bCorrelated);
        strScenario.strConstant.bUseMeasNoiseCrossCov = bCorrelated;
        cellCodegenArgs = BuildInputs_(strScenario);
        cellCodegenArgs{8} = coder.Constant(strScenario.strConstant);
        objConfig = coder.config('lib');
        objConfig.GenCodeOnly = true;
        objConfig.EnableDynamicMemoryAllocation = false;
        objConfig.EnableVariableSizing = false;
        objConfig.HardwareImplementation.ProdHWDeviceType = 'ARM Compatible->ARM Cortex-A (64-bit)';
        objConfig.HardwareImplementation.TargetHWDeviceType = 'ARM Compatible->ARM Cortex-A (64-bit)';
        codegen('-config', objConfig,'EKF_SlideWindow_FullCov_ObsUp','-args', cellCodegenArgs, ...
            '-d', charBuildRoot,'-o','FullCovObsFixed_arm64');

        % Compile every generated algorithm source with the target compiler. The MEX
        % gateway is absent from this standalone build; its MATLAB allocations are irrelevant here.
        strSources = dir(fullfile(charBuildRoot,'*.c'));
        verifyNotEmpty(testCase, strSources);
        for ui32Source = 1:numel(strSources)
            charSourcePath = fullfile(charBuildRoot, strSources(ui32Source).name);
            [~, charSourceName] = fileparts(charSourcePath);
            charObjectPath = fullfile(charBuildRoot,[charSourceName,'.o']);
            charCommand = sprintf(['aarch64-linux-gnu-gcc -std=c99 -O2 -fPIC ', ...
                '-Werror=implicit-function-declaration -I"%s" -I"%s/extern/include" ', ...
                '-c "%s" -o "%s"'], charBuildRoot, matlabroot, charSourcePath, charObjectPath);
            [i32Status, charOutput] = system(charCommand);
            assertEqual(testCase, i32Status, 0, charOutput);
        end

        % Require a complete target link and verify the artifact architecture. Execution
        % on ARM hardware or an emulator is a separate acceptance step.
        charLibraryPath = fullfile(charBuildRoot,'libFullCovObservation_arm64.so');
        [i32Status, charOutput] = system(sprintf( ...
            'aarch64-linux-gnu-gcc -shared -Wl,--no-undefined -o "%s" "%s"/*.o -lm', ...
            charLibraryPath, charBuildRoot));
        assertEqual(testCase, i32Status, 0, charOutput);
        [i32Status, charOutput] = system(sprintf('aarch64-linux-gnu-readelf -h "%s"', charLibraryPath));
        assertEqual(testCase, i32Status, 0, charOutput);
        verifyTrue(testCase, contains(charOutput,'AArch64'));
        [i32Status, charSymbols] = system(sprintf('aarch64-linux-gnu-nm -u "%s"', charLibraryPath));
        assertEqual(testCase, i32Status, 0, charSymbols);
        verifyEmpty(testCase, regexp(charSymbols,'\<(malloc|calloc|realloc|free)\>','once'));
    end
end
end

function cellInputs = BuildInputs_(strScenario)
cellInputs = {strScenario.dxState, strScenario.dCovariance, strScenario.dTimestamps, ...
    strScenario.strMeasurements, strScenario.strDynamics, strScenario.strModel, ...
    strScenario.strMutable, strScenario.strConstant};
end

function RemoveBuild_(charBuildRoot)
clear FullCovObsFixed_test_mex
rmpath(charBuildRoot);
rmdir(charBuildRoot,'s');
end

function strScenario = BuildFullWindowScenario_(ui16WindowCapacity)
strScenario = BuildFullCovObservationTestProblem();
% Exercise nonzero backward N; augmentation must leave this interval noise in the joint prior.
strScenario.strModel.dIntegrProcessNoiseCovQ = 1e-5*eye(17);
strScenario.strModel.dFlowSTM(1:3, 4:6) = eye(3);
strConst = strScenario.strConstant;
strConst.ui16NumWindowPoses = ui16WindowCapacity;
strConst.ui32WindowMaxSize = uint32(ui16WindowCapacity)*uint32(strConst.ui16WindowPoseSize);
strConst.ui32FullStateSize = uint32(strConst.ui16StateSize)+strConst.ui32WindowMaxSize;
strConst.ui32FullCovSize = uint32(strConst.ui16StateSize) + ...
    uint32(ui16WindowCapacity)*uint32(strConst.ui16WindowStateCovSize);
strScenario.strConstant = strConst;
ui32StateSize = uint32(strConst.ui16StateSize);
strScenario.dxState = zeros(strConst.ui32FullStateSize, 1);
strScenario.dxState(1:3) = [4;-1;0.5];
strScenario.dxState(7:9) = [0.01;-0.02;0.03];
strScenario.dxState(14) = 0.03;
strScenario.dCovariance = 0.01*eye(strConst.ui32FullCovSize);
strScenario.dTimestamps = linspace(0,-0.5, double(strConst.ui16NumWindowPoses)+1)';
strScenario.strMutable.ui16WindowStateCounter = strConst.ui16NumWindowPoses;
strScenario.strMutable.dKcam = [500, 0, 512;0, 500, 512;0, 0, 1];
strScenario.strMutable.i8CentroidingAlgorithmMode = uint8(0);
strScenario.strMutable.dReferenceMetricRadius = 2;
strScenario.strDynamics.strMainData.dRefRadius = 2;
strScenario.strMutable.dCenMeasApparentSizeLawCoeff = .05;
strScenario.strMutable.dMeanInstFOVinRadPx = 1e-3;
strScenario.strMutable.dDirOfMotionMeasCov = 1e-3*eye(3);
strScenario.strMutable.ui32MaxMeasEditingOccurrence = uint32(3);
strScenario.strMutable.dMahaDist2MeasThr = 9;
strScenario.strDynamics.strBody3rdData(1).strOrbitData.dChbvPolycoeffs(:) = 0;
strScenario.strDynamics.strBody3rdData(1).strOrbitData.dChbvPolycoeffs(1:3:9) = [20;5;3];
dCameraFromIN = strScenario.strMutable.dDCM_CamFromSCB*strScenario.strModel.dDCM_SCBiFromIN(:, :, 1);
dCurrentRotation = ComputeTargetAttitudeBias(strScenario.dxState(7:9))*strScenario.dNominalRotation;
dPrevPosition = strScenario.dxState(1:3)-[0.2;0.1;0.1];
dDirection = dCameraFromIN*(strScenario.dxState(1:3)-dPrevPosition);
dDirection = dDirection/norm(dDirection);
strScenario.strMeasurements.dDirectionOfMotion_CurrentCamFromPrevCam_Cam = dDirection;
strScenario.strMeasurements.dRangeLidarCentroid(2:3) = pinholeProjectHP( ...
    strScenario.strMutable.dKcam, dCameraFromIN, strScenario.dxState(1:3), zeros(3, 1));
for ui32Pose = uint32(1):uint32(strConst.ui16NumWindowPoses)
    ui32Start = ui32StateSize+(ui32Pose-1)*uint32(7);
    strScenario.dxState(ui32Start+uint32(1:3)) = ...
        strScenario.dxState(1:3)-double(ui32Pose)*[0.2;0.1;0.1];
    strScenario.dxState(ui32Start+uint32(4:7)) = DCM2quat(dCurrentRotation*dCameraFromIN', false);
end
end
