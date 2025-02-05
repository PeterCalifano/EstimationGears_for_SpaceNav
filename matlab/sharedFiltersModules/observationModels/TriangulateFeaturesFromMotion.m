function [dFeatPosVec_NavFrame, dFeatPosVec_Ck, dRelPos_CkFromCi_Ci, dDCM_CiFromCk, bConvergenceFlag] = ...
                            TriangulateFeaturesFromMotion(dyMeasVec, ...
                            dMeasCovSigma, ...
                                                          dDCM_CiFromCk, ...
                                                          dRelPos_CkFromCi_Ci, ...
                                                          dDCM_NavFrameFromCk, ...
                                                          dPositionCk_NavFrame, ...
                                                          dFeatInverseDepthGuess_Ck, ...
                                                          dPrincipalPoint_UV, ...
                                                          ui32NumOfPoses,...
                                                          ui32NumOfFeatures, ...
                                                          ui32MaxIter,...
                                                          dDeltaResNormRelTol, ...
                                                          ui32MaxNumOfPoses, ...
                                                          ui32MaxNumOfFeatures) %#codegen
arguments
    % DEVNOTE: comment validation args for speed
    dyMeasVec                   (:,:)   double  {ismatrix, isnumeric}   % Actual size: (2*ui32NumOfPoses, MaxNumberOfFeatures)
    dMeasCovSigma               (1,1)   double  {isscalar, isnumeric}
    dDCM_CiFromCk               (3,3,:) double  {ismatrix, isnumeric}
    dRelPos_CkFromCi_Ci         (3,:)   double  {ismatrix, isnumeric} 
    dDCM_NavFrameFromCk         (3,3)   double  {ismatrix, isnumeric}   % Attitude of anchor pose
    dPositionCk_NavFrame        (3,1)   double  {ismatrix, isnumeric}   % Position of anchor pose
    dFeatInverseDepthGuess_Ck   (3,:)   double  {ismatrix, isnumeric}   % Actual size: (3, MaxNumberOfFeatures)
    dPrincipalPoint_UV          (2,1)   double  {ismatrix, isnumeric}
    ui32NumOfPoses              (1,1)   uint32  {isscalar, isnumeric}
    ui32NumOfFeatures           (1,1)   uint32  {isscalar, isnumeric}
    ui32MaxIter                 (1,1)   uint32  {isscalar, isnumeric} = uint32(5);
    dDeltaResNormRelTol         (1,1)   double  {isscalar, isnumeric} = 1E-6;
    ui32MaxNumOfPoses           (1,1)   uint32  {isscalar, isnumeric} = ui32NumOfPoses; % TBC if not needed
    ui32MaxNumOfFeatures        (1,1)   uint32  {isscalar, isnumeric} = ui32NumOfFeatures;
end
%% PROTOTYPE
% 
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% REFERENCE:
% 1) High-precision, consistent EKF-based visual-inertial odometry, Li, Mourikis, 2023
% 2) Vision-Aided Inertial Navigation for Spacecraft Entry, Descent, and Landing, Mourikis, 2009
% 3) Statistical Orbit determination, Chapter 4, Tapley 2004
% 4) Optimal State estimation: Kalman, H Infinity, and Nonlinear Approaches, Dan Simon, 2006
% 5) Hartley, R. and Zisserman, A., 2003. Multiple view geometry in computer vision. Cambridge university press.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description

% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-12-2023    Pietro Califano     First prototype coded.
% 29-12-2023    Pietro Califano     Major reworking. Validation needed.
% 03-01-2024    Pietro Califano     Update of function to new coding rules; upgrade to support estimation 
%                                   of multiple features together.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% pinholeProjectIDP()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Replace WLS with Givens Rotations for improved efficiency TBC
% 2) Modify code for static memory allocation (WIP)
% 3) MAJOR TODO: modify code to reduce memory footprint in constructing the problem when features are not
% coupled through noise (see DEVNOTE below)
% -------------------------------------------------------------------------------------------------------------
%% Function code

% DEVNOTE: since we are solving for the features given the cameras, there is not share of information and block problems can be solved
% independently, unless noise couples them (e.g. through the cameras?). Do this upgrade for memory efficiency!

ui32SOLVED_FOR_SIZE = uint32(3 * ui32NumOfFeatures); % Depends on number of features. 1 feature <-> 3 params.
ui32SOLVED_FOR_MAX_SIZE = uint32(3 * ui32MaxNumOfFeatures); % Depends on number of features. 1 feature <-> 3 params.

% Input validation checks
% assert(ui32NumOfPoses <= size(dDCM_NavFrameFromC, 3));
% assert()

% Covariance weight
dMeasInformWeight = 1./dMeasCovSigma; % ACHTUNG, here we assume scalar to reduce (a lot) computations

%% Gauss-Newton iteration loop
ui8IterCounter          = uint8(0);
ui8NonConvergingCounter = uint8(0);
bConvergenceFlag        = false;

% Initialize inverse depth parameterization at estimation time tk
% NOTE: do division only ONCE (like in transformEPtoIDP!
dFeatInverseDepth_Ck = reshape(dFeatInverseDepthGuess_Ck, 3*size(dFeatInverseDepthGuess_Ck, 2), 1);
% TODO: line above MUST be modified for static size allocation

% Compute squared relative change threshold
dDeltaResNormRelTol2 = dDeltaResNormRelTol*dDeltaResNormRelTol; 

% Initialize relative residual norm delta
dDeltaResRelNorm2 = 1E2; 

while dDeltaResRelNorm2 > dDeltaResNormRelTol2
   
    % Initialize LAMBDA Info matrix and dProjectResidualVec (residuals projection)
    dLAMBDA                 = zeros(ui32SOLVED_FOR_MAX_SIZE, ui32SOLVED_FOR_MAX_SIZE); % TBC TO UPDATE
    dLeastSquaresRHS        = zeros(ui32SOLVED_FOR_MAX_SIZE, 1); % TO UPDATE
    dFeatInvDepUpdateDelta  = zeros(ui32SOLVED_FOR_MAX_SIZE, 1); % TO UPDATE

    % Initialize residual vector norm squared
    dResVecNorm2 = 0.0;

    % Compute Observation matrix and residual vector (linerized problem)
    dFeatHobs = coder.nullcopy(zeros(2, 3));

    for ui32IdPose = 1:ui32NumOfPoses % Loop through Ci 

        % Reset Measurement extraction index
        ui32IdMeas = uint32(1); % DEVNOTE: can be easily made uint16 (up to 65535 if memory bound)
        ui32IdFeat = uint32(1);

        dDCM_CkFromCi           = transpose( dDCM_CiFromCk(:,:, ui32IdPose) ); % TBD remove transpose!
        dTmpRelPos_CkFromCi_Ci  = dRelPos_CkFromCi_Ci(:, ui32IdPose);

        % Compute dh/dInvDep jacobian for camera Ci
        % dDCM required in Jacobian --> dDeltaDCM_CkFromCi % DEVNOTE review this in detail, my view of the
        % problem seems different from papers?
        dJacNormCoordsWrtIDP = [transpose(dDCM_CkFromCi(1,:)), transpose(dDCM_CkFromCi(2,:)), dTmpRelPos_CkFromCi_Ci];

        for ui32FeatID = 1:ui32NumOfFeatures % Loop through features
    
            % ACHTUNG: TO REWORK to handle multiple features
            ui32TmpFeatIDs = ui32IdFeat : ui32IdFeat + uint32(2);

            % Compute ith predicted measurement [2x1]
            % [dFeatPosPred_Ci, dFeatProjPred_UVi] =  pinholeProjectIDP(dDCM_CkFromCid, ...
            %                                                          dRelPos_CkFromCid_Cid, ...
            %                                                          dFeatInverseDepth_Ck(ui32TmpFeatIDs, 1), ...
            %                                                          dPrincipalPoint_UV);

            [dPredictFeatPos_Ci, dPredictFeatProj_Ci] = normalizedProjectIDP(dDCM_CkFromCi, ...
                                                                     dTmpRelPos_CkFromCi_Ci, ...
                                                                     dFeatInverseDepth_Ck(ui32TmpFeatIDs, 1) );

            % Compute relative attitude as DCM
            % DEVNOTE: original implementation uses normalized coordinates, but working in pixel coordinats
            % is likely better to reduce computations (avoid conversion to normalized coordinates, i.e. lots
            % of divisions).

            % Compute ith Observation Matrix [2x3] per feature
            % Hobs = dzdh * dh/dInvDep; (Ref.[1] notation, normalized coordinates)
            % dDCM required in Jacobian --> dDeltaDCM_CiFromCk(:, ui32IdPose);

            dFeatHobs(1,:) = [ 1/dPredictFeatPos_Ci(3), 0, - dPredictFeatPos_Ci(1) / dPredictFeatPos_Ci(3) ^ 2];
            
            % TODO: verify that dDCM_CkFromCid is the correct matrix here
            dFeatHobs(2,:) = [0, 1/dPredictFeatPos_Ci(3), -dPredictFeatPos_Ci(2)/dPredictFeatPos_Ci(3)^2];
                                              
            % Chain jacobian of jth feature for ith camera (Ci)
            dFeatHobs = dFeatHobs * dJacNormCoordsWrtIDP;

            % Accumulate LAMBDA Information matrix and Nmatrix residuals projection vector
            % NOTE: NO measurement noise covariance
            dFeatureProjResVec = dyMeasVec( ui32IdMeas:ui32IdMeas + uint32(1) , 1) - dPredictFeatProj_Ci;
            
            % Accumulate LAMBDA matrix (A^T * W * A)
            dLAMBDA(ui32TmpFeatIDs, ui32TmpFeatIDs) = dLAMBDA(ui32TmpFeatIDs, ui32TmpFeatIDs) + dMeasInformWeight * transpose(dFeatHobs) * dFeatHobs; 
            % Actual matrix is band diagonal here --> features are independent of each other
            
            % Accumulate LS RHS (A^T * W * b)
            dLeastSquaresRHS(ui32TmpFeatIDs) = dLeastSquaresRHS(ui32TmpFeatIDs) + ...
                    transpose(dFeatHobs) * dMeasInformWeight * dFeatureProjResVec;

            % Accumulate residual vector norm squared
            dResVecNorm2 = dResVecNorm2 + dFeatureProjResVec' * dFeatureProjResVec;

            % Store residuals (mainly for debug: remove for mem. optim.)
            % dPixMeasResVec( ui32IdMeas:ui32IdMeas + uint32(2) , 1) = dFeatureProjResVec;

            % Move to pointers next entry
            ui32IdMeas = ui32IdMeas + uint32(2);
            ui32IdFeat = ui32IdFeat + uint32(3);

        end % END Loop through features

    end % END Loop through poses (i.e., time instants in the window)
    
    % FIND LOCAL LINEAR Least Squares solution (all features together)
    dFeatInvDepUpdateDelta(1:ui32SOLVED_FOR_SIZE) = dLAMBDA(1:ui32SOLVED_FOR_SIZE, 1:ui32SOLVED_FOR_SIZE)\dLeastSquaresRHS(1:ui32SOLVED_FOR_SIZE); 
    
    assert( all(size(dFeatInvDepUpdateDelta(1:ui32SOLVED_FOR_SIZE)) == size(dFeatInverseDepth_Ck), 'all') );

    % Apply correction to reference state for next iteration
    dFeatInverseDepth_Ck = dFeatInverseDepth_Ck + dFeatInvDepUpdateDelta(1:ui32SOLVED_FOR_SIZE);


    if ui8IterCounter > 0
        % Compute relative square residual change
        dDeltaResRelNorm2 = abs(dResVecNorm2 - dResVecNorm2_prev) ./ dResVecNorm2_prev;
        if dResVecNorm2 > dResVecNorm2_prev
            warning('Detected increase in residual norm. Solver is not converging!')
            ui8NonConvergingCounter = ui8NonConvergingCounter + uint8(1);
        end
    end

    % SOLVER LOOP COUNTER
    if ui8IterCounter >= ui32MaxIter
        disp('SOLVER STOP: MAX ITER REACHED.')
        break;
    end

    % Store previous residual vector norm squared
    dResVecNorm2_prev = dResVecNorm2;

    ui8IterCounter = ui8IterCounter + uint8(1);
end


% Output definition
dFeatPosVec_NavFrame = zeros(3, ui32MaxNumOfFeatures);

% Check convergence and return output
if ui8NonConvergingCounter == ui8IterCounter - uint8(1)
    bConvergenceFlag = false;
    return % RETURN NON CONVERGED
    
elseif ui8NonConvergingCounter < ui8IterCounter - uint8(1) && dResVecNorm2 <= dResVecNorm2_prev
    
    bConvergenceFlag = true;

    % Compute output in Camera frame at tk time
    dFeatPosVec_Ck = transformIDPtoEP(dFeatInverseDepth_Ck, ui32NumOfFeatures, ui32MaxNumOfFeatures);

    % Compute landmarks positions in Navigation frame at tk time (TODO validated)
    for ui32FeatID = 1:ui32NumOfFeatures
        dFeatPosVec_NavFrame(1:3, ui32FeatID) =  ( dDCM_NavFrameFromCk * dFeatPosVec_Ck(1:3, ui32FeatID) ) ...
            +  dPositionCk_NavFrame;
    end

    return % RETURN CONVERGED
% else TODO: verify that there is not need for this else (make sure no undefined behaviour!) The two if
% conditions above should cover any feasible case.
end



end
