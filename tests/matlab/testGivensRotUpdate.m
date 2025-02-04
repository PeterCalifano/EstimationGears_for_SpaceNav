close all
clear
clc

% Add paths for tests
addpath(genpath("../../matlab/"))
addpath(genpath("../../lib/UnitTesting4MATLAB/utils"))

% TEST SETUP
% Reference: example from Tapley chapter 5.6.4
% Inputs
bNPRIOR_INFO = false;
dYobs = [-1.1; 1.2; 1.8];
dHobsMatrix = [1, -2;
               2, -1;
               1, 1];

dxPrior = [2; 2];
dPcov = diag([100, 100]);

dSroot = chol(dPcov, 'lower');
dSRInfoMatPrior = dSroot^-1;
dMeasCovSR = diag([1, 1, 1]);
bRUN_WHITENING = false;
bENABLE_EDITING = false;

% Correct estimate after ith residual (FIXTURES)
dxState_r0 = [2.17964071856287; 1.64071856287425]; 
dxState_r1 = [1.1750640102856; 1.14176767288272]; 
dxState_r2 = [1.00335913215659; 0.970062794753707];

% Missing covariance fixture?
% TODO

dDefaultTol = 1e-9;

%% test_GivensRotSRIF_withMexEquivalence

tic
% Algorithm test
[dxPost, dSRInfoMatPost, dInfoVecPost, dErrorVec] = GivensRotSRIF(dxPrior, ...
    dSRInfoMatPrior, ...
    dYobs, ...
    dHobsMatrix, ...
    bNPRIOR_INFO, ...
    bRUN_WHITENING,...
    dMeasCovSR);
toc

assertDifference(dxPost, dxState_r2, dDefaultTol);

if exist('GivensRotSRIF_MEX', 'file')
    % MEX EQUIVALENCE TEST
    tic
    [dxPost_MEX, dSRInfoMatPost_MEX, dInfoVecPost_MEX, dErrorVec_MEX] = GivensRotSRIF_MEX(dxPrior, ...
        dSRInfoMatPrior, ...
        dYobs, ...
        dHobsMatrix, ...
        bNPRIOR_INFO, ...
        bRUN_WHITENING,...
        dMeasCovSR);
    toc

    % Test assertions
    assertDifference(dxPost         , dxPost_MEX);
    assertDifference(dSRInfoMatPost , dSRInfoMatPost_MEX);
    assertDifference(dInfoVecPost   , dInfoVecPost_MEX  );
    assertDifference(dErrorVec      , dErrorVec_MEX     );

else
    warning('MEx equivalence test skipped due to missing function.')
end


%% test_GivensRotEKFvsSRIF_TYPE0
% TODO
ui8FILTER_TYPE = uint8(0);
% Test to do for complete coverage:
% 0) Full covariance passing Pcov
% 1) UD passing UD decomposition
% 2) Square root covariance passing Sroot
% All must return the same result

tic
if ui8FILTER_TYPE == 0
    dInputTest = dPcov;
    dDxPrior = zeros(2);

elseif ui8FILTER_TYPE == 2
    dInputTest = dSroot;
    dDxPrior = zeros(2);

elseif ui8FILTER_TYPE == 1
    [dInputTest, dDxPrior] = UDdecomposition(dPcov);
end


% Call GivensRotEKF
[dxPost, dPxPost, dDxPost] = GivensRotEKF(dxPrior, ...
    dInputTest, ...
    dYobs, ...
    dHobsMatrix, ...
    bNPRIOR_INFO, ...
    bRUN_WHITENING, ...
    dMeasCovSR, ...
    ui8FILTER_TYPE, ...
    dDxPrior);

toc

% Algorithm test
[dxPost_SRIF, dSRInfoMatPost_SRIF, dInfoVecPost_SRIF, dErrorVec_SRIF] = GivensRotSRIF(dxPrior, ...
    dSRInfoMatPrior, ...
    dYobs, ...
    dHobsMatrix, ...
    bNPRIOR_INFO, ...
    bRUN_WHITENING,...
    dMeasCovSR);

assertDifference(dxPost, dxState_r2, dDefaultTol);

% Test assertions
assertDifference(dxPost, dxState_r2, dDefaultTol);
assertDifference(dxPost, dxPost_SRIF, dDefaultTol);

% TODO add assertion on covariance

dSpost = dSRInfoMatPost_SRIF^-1;
dPxPost_ref = dSpost * dSpost';

% Discrepancy between GivensRotSRIF and GivensRotEKF should be a numerical zero
if ui8FILTER_TYPE == 2
    dPxPost* dPxPost'- dPxPost_ref %#ok<*NOPTS>
elseif ui8FILTER_TYPE == 0
    dPxPost - dPxPost_ref
elseif ui8FILTER_TYPE == 1
    dPxPost * dDxPost * dPxPost' - dPxPost_ref
end

%% test_GivensRotEKF_ModifiedAgeeTurner_AgeeTurnerRank1Update
ui8FILTER_TYPE = uint8(1);

[dInputTest, dDxPrior] = UDdecomposition(dPcov);
dUprior = dInputTest;

% GivensEKF Algorithm test (reference)
tic

[dxPost, dPxPost, dDxPost] = GivensRotEKF(dxPrior, ...
    dInputTest, ...
    dYobs, ...
    dHobsMatrix, ...
    bNPRIOR_INFO, ...
    bRUN_WHITENING, ...
    dMeasCovSR, ...
    ui8FILTER_TYPE, ...
    dDxPrior);
toc

% Algorithm test: Modified Agee-Turner and (basic, No consider states) UD EKF update module
dzPrior = dxPrior;
dSRmeasCov = chol(dMeasCovSR);

dyRes = dYobs - dHobsMatrix * dzPrior; % Evaluate residuals

% Kalman gain check
Kcheck = dPcov* dHobsMatrix' / (dHobsMatrix* dPcov *dHobsMatrix' + dMeasCovSR);

[dzPost, dUpost, dK, dDpost] = UDobsUpdate_ModAgeeTurner(dzPrior, ...
    dUprior, ...
    dDxPrior, ...
    dSRmeasCov, ...
    dHobsMatrix, ...
    dyRes, ...
    bENABLE_EDITING)

dzPrior + Kcheck * dyRes

dErrState = dxPost - dzPost %% DEVNOTE: state update is incorrect. Error may be in the Kalman gain,
% but does not influence the update of the UD factors, which are correct
dErrCov = dPxPost* dDxPost * dPxPost' - dUpost * dDpost * dUpost' % PASSED
% TODO: add test assertions

% Algorithm test Agee-Turner
% Execute Rank 1 update of UD factors through  Agee-Turner
% algorithm and mean estimate update
% bValidResiduals = true(Nres, 1);
%
% for idRes = 1:Nres
%
%     if bValidResiduals(idRes) == true
%
%         % Update UD factors of covariance and compute Kalman gain column
%         dHrow = dHobs(idRes, :);
%
%         [dUpost, dDpost, dK(:, idRes)] = UDRank1Up_ModAgeeTurner(dUpost, dDpost, measCovRdiag(idRes), dHrow);
%
%
%         dK(:, idRes)
%     end
%
%     % Update mean state estimate using accepted residuals
%     dzPost = dzPost + dK(:, idRes) * dyRes(idRes);
% end
