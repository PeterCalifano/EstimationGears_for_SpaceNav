function tests = testWindowRefFrame
%% SIGNATURE
% tests = testWindowRefFrame
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Verify constant ownership of clone coordinates and static MEX specializations for both frames.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB tests.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 10-09-2026  Pietro Califano, Codex gpt-6    Validate constant window-frame selection.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% SetupPaths_EstimationGears, ComputeWindowPose, MATLAB Coder.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOriginalPath = path;
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','..'));
SetupPaths_EstimationGears;
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testBuilderOwnsFrame(testCase)
strConstant = filter_tailoring.BuildArchitectureTemplate('bWriteBusDefs',false);
strMutable = filter_tailoring.BuildInputStructsTemplate(strConstant);
verifyEqual(testCase,strConstant.enumWindowRefFrame,EnumWindowRefFrame.INERTIAL);
verifyTrue(testCase,isa(strConstant.enumWindowRefFrame,'uint8'));
verifyFalse(testCase,isfield(strMutable,'charWindowRefFrame'));
verifyFalse(testCase,isfield(strMutable,'enumWindowRefFrame'));
strTarget = filter_tailoring.BuildArchitectureTemplate('bWriteBusDefs',false, ...
    'enumWindowRefFrame',EnumWindowRefFrame.TARGET_FIXED);
verifyEqual(testCase,strTarget.enumWindowRefFrame,EnumWindowRefFrame.TARGET_FIXED);
end

function testStaticFrameSpecializations(testCase)
assumeFalse(testCase,isempty(which('codegen')),'MATLAB Coder is required.');
charBuildDir = tempname;
mkdir(charBuildDir);
charOriginalDir = pwd;
addpath(charBuildDir);
cd(charBuildDir);
objCleanup = onCleanup(@() Cleanup_(charBuildDir,charOriginalDir));
objConfig = coder.config('mex');
objConfig.EnableVariableSizing = false;
objConfig.EnableDynamicMemoryAllocation = false;
dTargetRotation = [0,-1,0;1,0,0;0,0,1];

% Each frame has its own compiled constant; positions and attitudes remain runtime inputs.
for enumFrame = [EnumWindowRefFrame.INERTIAL,EnumWindowRefFrame.TARGET_FIXED]
    strConstant.enumWindowRefFrame = enumFrame;
    cellInputs = {[3;-2;1],[1;0;0;0],[sqrt(.5);0;0;-sqrt(.5)], ...
        [1;0;0;0],strConstant,zeros(3,1)};
    cellCodegenInputs = cellInputs;
    cellCodegenInputs{5} = coder.Constant(strConstant);
    clear WindowFramePose_mex
    codegen('-config',objConfig,'ComputeWindowPose','-args',cellCodegenInputs, ...
        '-o',fullfile(charBuildDir,'WindowFramePose_mex'),'-d',fullfile(charBuildDir,'build'));
    for dScale = [1,2]
        cellInputs{1} = dScale*[3;-2;1];
        dExpectedPosition = cellInputs{1};
        if enumFrame == EnumWindowRefFrame.TARGET_FIXED
            dExpectedPosition = dTargetRotation*dExpectedPosition;
        end
        [dPosition,dQuaternion] = WindowFramePose_mex(cellInputs{:});
        verifyEqual(testCase,dPosition,dExpectedPosition,'AbsTol',2e-14);
        verifyEqual(testCase,Quat2DCM(dQuaternion,false),dTargetRotation,'AbsTol',2e-14);
    end
end
end

function Cleanup_(charBuildDir,charOriginalDir)
clear WindowFramePose_mex
cd(charOriginalDir);
rmpath(charBuildDir);
rmdir(charBuildDir,'s');
end
