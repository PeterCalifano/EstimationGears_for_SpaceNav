function [xhatPost, Spost, Ppost, Pprior, WmWc] = SRUKF_ObsUpdate(xhat, Scov, i_dyObs, Q, R, params) %# codegen
%% PROTOTYPE
% [xhatPost, Spost, Ppost] = SR_UKF_kernel(xhat, Scov, ymeas, Q, R)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% % DEV IN PROGRESS: NOT YET WORKING CORRECTLY. PAUSED: SR-USKF is already
%                    more general.
%
% ACHTUNG: ScaledUTsquareRoot function requires tuning of the parameters of
% the UT propagation to properly work. Too small or too large values lead
% to incorrect estimation of the propagated mean state.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 06-07-2023    Pietro Califano     Template for SR-UKF coded. Applied to 
%                                   SGN Assignment exe 3.
% 05-08-2023    Pietro Califano     Verified. Process noise may be critical
%                                   for the stability of the filter due to
%                                   negative Cholupdate (mean Sigma point).
%                                   The test bench (SGN exe 2.3) is still
%                                   not handled by the filter.
% 17-09-2023    Pietro Califano     Template reworking and split in Time
%                                   and Observation update.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Time Update
% Input Scov here must be: UPPER triangular

i_bComputeWeights = 1;
i_dWmWc = 0;
i_dPertubStep = 0;

[xhatProp, Sprop, dCsi, WmWc] = ScaledUTsquareRoot(xhat, Scov, Q, i_bComputeWeights, i_dWmWc, i_dPertubStep, params);

% FOR DEBUG
dt = 1;
nonlinfcn = @(x) CW_analytical(dt, params.nOrb, x);
k = 0;
alpha = 1e-3;
beta = 2;
Pcov = Scov' * Scov;
[xhatPriorTest, PpriorTest] = ScaledUT(nonlinfcn, xhat, Pcov, k, alpha, beta);

Sprior = chol(PpriorTest);

% diffNormMean = norm(xhatPrior - xhatPriorTest);

% NOTE: Sprior seems completely wrong!
% diffNormCov = norm(Sprior - SpriorTest);

% assert(diffNormMean < 1e-8, 'Mean vector greatly differs between UT and SR version.')
% assert(diffNormCov < 1e-8, 'Square Root covariance greatly differs between UT and SR version.')


% Variables allocation
Spost = Sprior;

if not(isempty(i_dyObs))
    %% Measurement Update
    % Compute Predicted measurements
    % ycsi = MeasModel(csi, params);
    yCsi = zeros(length(i_dyObs), size(dCsi, 1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Retrieve info from params
    corners_cell = params.corners_cell;
    vertices_b = params.vertices_b;
    q_state = params.q_state;
    idt = params.idt;
    Nviscorners = params.Nviscorners;
    u0 = params.u0;
    v0 = params.v0;
    pixel_d = params.pixel_d;
    baseline = params.baseline;
    focal_l = params.focal_l;
    elapsed_time = params.elapsed_time;
    nOrb = params.nOrb;

    % Application specific measurement model
    hmeas_fcn = @(Xcam) [u0 - pixel_d*focal_l* Xcam(2)./Xcam(3);
        v0 + pixel_d*focal_l* Xcam(1)./Xcam(3);
        baseline*pixel_d*focal_l ./ Xcam(3)];

    idfirst = 1;
    idlast = 3;

    for idp = 1:Nviscorners

        for csi_id = 1:size(dCsi, 2)
            % Extract vertex id
            vertex_id = corners_cell{idt, 2}(idp);

            % Get vector identifiyng the vertex in the target body frame
            vertex_body = vertices_b(vertex_id, :);

            % Conversion from target body frame to Chaser Cam frame
            vertex_cam = TargetBRFtoCamRF(vertex_body, dCsi(:, csi_id), elapsed_time, q_state(idt, :), nOrb);
            % Generate expected measurements for csi sigma point
            yCsi(idfirst:idlast, csi_id) = hmeas_fcn(vertex_cam);

        end
        idfirst = idlast + 1;
        idlast = idlast + 3;

    end

    % Do not change code below (filter formulation is the same in all cases)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yPredMean = sum(WmWc(:, 1)' .*yCsi, 2);

    % Compute Measurement SR covariance
    % RsquareRoot = chol(R, 'upper'); % Compute SR of the Process Noise matrix
    RsquareRoot = sqrtm(R); 


    [~, Syk] = qr( [sqrt(WmWc(2:end, 2)) .* (yCsi(:, 2:end) - yPredMean),  RsquareRoot]', 'econ');
    % Note: The Wc needs to be "embedded" in X, due to how the MATLAB
    % cholupdate works. Therefore, if Wc is negative --> apply sqrt(abs()),
    % then use "downdate" instead of "update".
    Syk = cholupdate(Syk, sqrt(abs(WmWc(1, 2))).*(yCsi(:, 1) - yPredMean), '-'); % cholupdate

    % Compute Cross Covariance
    Pxyk = WmWc(:, 2) .* (dCsi - xhatPrior) * (yCsi - yPredMean)'; % to check carefully

    % ACHTUNG: Syk is upper diagonal as it exits from cholupdate. Reference paper
    % assumes it as lower diagonal.

    % Compute Kalman Gain
%     Kk = (Pxyk/Syk)/Syk'; % Use back-substitution (Syk upper triangular here)
    Kk = Pxyk/(Syk'*Syk);

    % Update Mean
    xhatPost = xhatPrior + Kk*(i_dyObs - yPredMean);

    % Update SR covariance (sequentially)
    Uk = Kk * Syk';

    % Warning: delicate step next. The filter may fail due to loose of the
    % positive-definiteness of the covariance.
    for i = 1:size(Uk, 2)
        Spost = cholupdate(Spost, Uk(:, i), '-');
    end

else
    % Time update only
    xhatPost = xhatPrior;

end

% Optionally compute Full Covariance Matrix
if nargout > 2
    Ppost = Spost' * Spost;
    if nargout > 3
        Pprior = Sprior' * Sprior;
    end
end



end

