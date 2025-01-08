close all
clear
% clc

%% Test cases definition
testcase = 16;

% TODO: 
% Exponential fitting (Non-linear LS)
% Isserlis formula for High-Order moments

% Test 1: WLS function:                                                                         PASSED
% Test 2: LOESS:                                                                                PASSED
% Test 3: Unscented transform (Scaled):                                                         PASSED
% Test 4: EKF/UKF modules for UDEKF-SRUKF hybrid:                                               TODO
% Test 5: GivensRotEKF to validate Mod. Agee-Turner. UD update module using Mod. Agee-Turner.   PASSED                                       
% Test 6: Chol, cholupdate, cholRank1Update, Agee-Turner Rank-1 Update: 
%         Custom cholRank1Update: PASSED; Agee-Turner Rank-1:                                   PASSED
% Test 7: CVX L1, L2 regression filters:                                                        PASSED
% Test 8: Test and validation of POLRLS:                                                        FAILING
% Test 9: GivensRotSRIF validation + MEX equivalence test                                       PASSED
% Test 10: GivensRotEKF validation                                                              PASSED
% Test 11: Feature init. problem using Newton-Raphson triangulation:                            FAILING
% Test 12: Test AdaptQCovASNC runtime (verification):                                           TODO
% Test 13: EM Algorithm for GMM fitting and testing of Covariance generator                     TODO                            
% Test 14: Test Weighted Modified Gram-Schmidt, UD Covariance Time Update                       TODO
% Test 15: Schur Marginalization                                                                PASSED
% Test 16: Covariance recovery algorithm by (Kaess, Dellaert 2009)                              IN PROGRESS, FAILING
% Test 17: Mutual Information and Information gain computation algorithm (Kaess, 2009)          TODO


switch testcase

    case 1
        %% TEST 1: LS algorithms validation
        I = [0.2, 0.3, 0.4, 0.5, 0.6];
        V = [1.23, 1.38, 2.06, 2.47, 3.17];

        % Test LS and WLS
        % LS % Validated
        H_LS = I';
        y = V';
        W_LS = eye(length(I));

        [xhat_LS, P_LS, res_LS, RMS_LS] = WLS(y, H_LS, W_LS);

        % WLS % Validated
        H_LS = I';
        y = V';
        R_WLS = diag(var(V));

        [xhat_WLS, P_WLS, res_WLS, RMS_WLS] = WLS(y, H_LS, R_WLS);

        % Test WLS by accumulation
        y = V;
        H_LS = reshape(I, 1, 1, length(y));
        R_WLS = ones(1, 1, length(y)) * var(V);
        [xhat_WLSaccum, P_WLSaccum, res_WLSaccum, RMS_WLSaccum] = WLSaccum(y, H_LS, R_WLS);

        
        % Test RWLS
        % NOT COMPLETED
        %         H = I;
        %         y = V;
        %         xhat0 = V(1)/I(1);
        %         W = var(V);

        %         [xhat, P, xhist] = RecursiveWLS(y, H, W, xhat0);

    case 2
        %% TEST 2 LOESS validation
        xmesh = linspace(-10, 10, 200)';
        Nsamples = length(xmesh);

        th1 = 0.6;
        th3 = 0.4;
        qintercept = -6.4;

        % Build custom testing function
        sumSines = zeros(length(xmesh), 1);
        f = [1, 5, 16, 29, 99, 589];
        for idf = 1:length(f)
            sumSines = sumSines + sin(f(idf) * xmesh);
        end

        trueCurve = th3 * xmesh.^3 + th1*xmesh + qintercept;
        trueCurve = trueCurve + mean(abs(trueCurve)) * sumSines;

        ytestSignal = trueCurve + mean(abs(trueCurve))*randn(1, Nsamples)'.^3;

        testPlot = figure;
        hold on;
        plot(xmesh, trueCurve, 'k-', 'Linewidth', 1.05, 'DisplayName', 'True')
        plot(xmesh, ytestSignal, 'b.', 'MarkerSize', 10, 'DisplayName', 'Sample points')
        xlabel('Indep. variable X')
        ylabel('Sample points Y')
        legend()
        DefaultPlotOpts();

         
        % LOESS options 
        f = 0.3;
        d = 1; % 0: constant, 1: linear, 2: quadratic
        t = 3;

        % Test LOESS
        [theta, Yhat, FirstDer] = LOESS(ytestSignal, xmesh, f, d, t);

        % Apply global LS to fit a polynomial
        degree = 3;
        H = zeros(length(xmesh), degree + 1);

        for d = 0:degree
            H(:, d+1) = xmesh.^d;
        end

        [coeffLS, P_LS] = WLS(ytestSignal, H);
        
        % Evaluate polynomial on xmesh
        YhatLS = zeros(length(xmesh), 1);
        for j = 1:degree + 1
            YhatLS = YhatLS + coeffLS(j) .* H(:, j);
        end

        figure(testPlot);
        hold on
        plot(xmesh, Yhat, 'r.-', 'Linewidth', 1.05, 'DisplayName', 'Smoothed curve by LOESS')
        plot(xmesh, YhatLS, 'g.-', 'Linewidth', 1.05, 'DisplayName', ['Order ', num2str(degree),' curve by LS'])

    case 3
        %% TEST 3: Unscented Transform validation
        xhat = [15, 1, 3];
        Rot = Rot1(0.05);
        Pcov = Rot' * diag([100, 10, 50]) * Rot;

        nlinfun = @(x) [x(1).^3 + x(3) - x(2); x(3)^2; x(2)^3 - x(1)^2];

        %         N = 1e5;
        %         MCsamples = mvnrnd(xhat, Pcov, N);
        %         MCsamplesPost = zeros(N, 3);
        %         for i = 1:N
        %             MCsamplesPost(i, :) = nlinfun(MCsamples(i, :));
        %         end
        %
        %         xhatPost = mean(MCsamplesPost, 1);
        %         PcovPost = 1/N * (MCsamplesPost - xhatPost)' * (MCsamplesPost - xhatPost);
        Nsize = 10000;
        avgtime = zeros(Nsize, 2);

        k = 0;
        alpha = 1e-3;
        beta = 2;

        for i = 1:Nsize

            tic
            [xhat_UT_chol, P_UT_chol] = SigmaPointTransform(nlinfun, xhat, Pcov, k, alpha, beta);
            avgtime(i, 1) = toc;

            tic
            [xhat_UT, P_UT] = SigmaPointTransform_sqrtm(nlinfun, xhat, Pcov, k, alpha, beta);
            avgtime(i, 2) = toc;

            assert(sum(abs(xhat_UT - xhat_UT_chol) > 1e-6) == 0, 'Results are not equivalent');
        end

        timePerRun_chol = mean(avgtime(:, 1));
        timePerRun_sqrtm = mean(avgtime(:, 2));

        fprintf('\nAvg time per run. Sqrtm: %1.8f vs chol: %1.8f\n', timePerRun_sqrtm, timePerRun_chol);
        %         fprintf('Mean MC: %6.3g, Mean UT: %6.3g \nP_MC: %6.3g, P_UT: %6.3g\n',...
        %             xhatPost, xhat_UT, PcovPost, P_UT)


    case 4
        %% TEST 4: EKF/UKF modules for UDEKF-SRUKF hybrid
        % Example from Cholesky decomposition Wikipedia page
        P = blkdiag([4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98], ...
                    [4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98], ...
                    [4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98], ...
                    [4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98], ...
                    [4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98], ...
                    [4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98], ...
                    [4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98], ...
                    [4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98], ...
                    [4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98], ...
                    [4, 12, -16;
                    12, 37, -43;
                    -16, -43, 98]);

        Ntests = 1E4;
        zeros(1, Ntests);

        for idN = 1:Ntests

            tic
            [U, D] = UDdecomposition(P); % Custom function
            runtime1(idN) = toc;

        end

        errCheck = max(abs(P - U*D*U'), [], 'all');

        runtimePerCycle1 = mean(runtime1)
        stdRuntimePerCycle1 = std(runtime1);
        medianRuntimePerCycle1 = median(runtime1)

        runtime2 = zeros(1, Ntests);

        for idN = 1:Ntests
            tic
            % Test computation of Cholesky factor S
            S = U*sqrt(D); % PASSED

            % err = max(abs(U*D*U' - S*S'), [], 'all');

            % Retrieve UD from S
            Cvec = diag(S);
            Nx = uint16(length(Cvec));
            Urecomp = zeros(Nx);

            for idc = 1:Nx
                Urecomp(1:idc, idc) = S(1:idc, idc)/Cvec(idc);
            end

            Drecomp = diag(Cvec.^2);

            % Check matrix again
            % Precomp = Urecomp * Drecomp * Urecomp';
            % errP = max(abs(Precomp - P), [], 'all'); % PASSED
            runtime2(idN) = toc;

        end

        runtimePerCycle2 = mean(runtime2)
        stdRuntimePerCycle2 = std(runtime2);
        medianRuntimePerCycle2 = median(runtime2)



        
    case 5
        %% TEST 5: GivensRotEKF to validate Modified Agee-Turner, Agee-Turner Rank-1 Updates, EKF templates: TODO
        % Example from Tapley chapter 5.6.4 (using GivensRotSRIF)

        % Test inputs
        i_bNO_PRIOR_INFO = false;
        i_dYobs = [-1.1; 1.2; 1.8]; 
        % Re-check example by Tapley. Not sure Givens algorithm takes the observation, but rather the
        % residual. If so, rename variable. --> No change needed, it is correct.
        i_dHobsMatrix = [1, -2;
                         2, -1;
                         1, 1];

        i_dxPrior = [2; 2];
        Pcov = diag([100, 100]);
    
        i_ui8FILTER_TYPE = 1;
       
        [inputTest, i_dDxPrior] = UDdecomposition(Pcov);
        i_dUprior = inputTest;

        % GivensEKF Algorithm test (reference)
        tic
        i_bRUN_WHITENING = false;
        i_dMeasCovSR = diag([1, 1, 1]);

        [o_dxPost, o_dPxPost, o_dDxPost] = GivensRotEKF(i_dxPrior, ...
            inputTest, ...
            i_dYobs, ...
            i_dHobsMatrix, ...
            i_bNO_PRIOR_INFO, ...
            i_bRUN_WHITENING, ...
            i_dMeasCovSR, ...
            i_ui8FILTER_TYPE, ...
            i_dDxPrior);
        toc

        % Algorithm test: Modified Agee-Turner and (basic, No consider states) UD EKF update module
        i_dzPrior = i_dxPrior;
        i_dSRmeasCov = chol(i_dMeasCovSR);

        i_dyRes = i_dYobs - i_dHobsMatrix * i_dzPrior; % Evaluate residuals

        i_bENABLE_EDITING = false;

        % Correct estimate after ith residual:
        x1 = [2.17964071856287; 1.64071856287425]; % OK, correct
        x2 = [1.1750640102856; 1.14176767288272]; % ERROR
        x3 = [1.00335913215659; 0.970062794753707];

        % Kalman gain check
        Kcheck = Pcov* i_dHobsMatrix' / (i_dHobsMatrix* Pcov *i_dHobsMatrix' + i_dMeasCovSR);

        [o_dzPost, o_dUpost, o_dK, o_dDpost] = UDobsUpdate_ModAgeeTurner(i_dzPrior, ...
            i_dUprior, ...
            i_dDxPrior, ...
            i_dSRmeasCov, ...
            i_dHobsMatrix, ...
            i_dyRes, ...
            i_bENABLE_EDITING)

        i_dzPrior + Kcheck * i_dyRes

        errState = o_dxPost - o_dzPost %% DEVNOTE: state update is incorrect. Error may be in the Kalman gain, 
        % but does not influence the update of the UD factors, which are correct
        errCov = o_dPxPost* o_dDxPost * o_dPxPost' - o_dUpost * o_dDpost * o_dUpost' % PASSED

        % Algorithm test Agee-Turner
        % Execute Rank 1 update of UD factors through  Agee-Turner
        % algorithm and mean estimate update
        % o_bValidResiduals = true(Nres, 1);
        % 
        % for idRes = 1:Nres
        % 
        %     if o_bValidResiduals(idRes) == true
        % 
        %         % Update UD factors of covariance and compute Kalman gain column
        %         i_dHrow = i_dHobs(idRes, :);
        % 
        %         [o_dUpost, o_dDpost, o_dK(:, idRes)] = UDRank1Up_ModAgeeTurner(o_dUpost, o_dDpost, measCovRdiag(idRes), i_dHrow);
        % 
        % 
        %         o_dK(:, idRes)
        %     end
        % 
        %     % Update mean state estimate using accepted residuals
        %     o_dzPost = o_dzPost + o_dK(:, idRes) * i_dyRes(idRes);
        % end



    case 6
        %% TEST 6: Chol Update and Agee-Turner Rank-1 Update test

        P = Rot1(0.1)*Rot3(0.5)*diag([2, 3, 4])*Rot3(0.5)'*Rot1(0.1)';
        R = chol(P);
        
        i_dUpdateCoeff = 0.5;
        i_dUpdateVec = [-3; 1; -0.4];
        
        % Compute referene result
        Pcheck = P + i_dUpdateCoeff * (i_dUpdateVec*i_dUpdateVec');
        
        % Check cholupdate
        Rcheck = chol(Pcheck);

        Ntests = 1e8;
        timevec1 = zeros(1, Ntests);
        timevec2 = zeros(1, Ntests);

        for idt = 1:Ntests
            tic
            RcholUp = cholupdate(R, sqrt(i_dUpdateCoeff)*i_dUpdateVec, "+");
            timevec1(idt) = toc;
        end

        avgtimeCustom1 = mean(timevec1)

        resRank1Chol = Pcheck - RcholUp' * RcholUp

        % Testing of custom cholRank1Update
        for idt = 1:Ntests
            tic
            Lnew = cholRank1Update(R', i_dUpdateVec, i_dUpdateCoeff, 1) % VALIDATED
            timevec2(idt) = toc;
        end

        avgtimeCustom2 = mean(timevec2)

        resCholRank1Custom = Pcheck - Lnew*Lnew'

        % Algorithm test
        [i_dU, i_dD] = UDdecomposition(P);
        [o_dUnew, o_dDnew] = UDrank1Up_AgeeTurner(i_dU, i_dD, i_dUpdateCoeff, i_dUpdateVec);
        
        % Check result
        resAgeeTurner = Pcheck - o_dUnew * o_dDnew * o_dUnew'

    case 7
        %% TEST 7: CVX L1, L2 regression filters
        % Generate a random signal with WN
        sigma = 5;
        slope = -0.8;
        Nsamples = 48;
        q = 13000;

        Xset = transpose(linspace(0, Nsamples-1, Nsamples));

        y = slope*Xset + q + sigma * randn(Nsamples, 1);

        % Reference solution
        ytruth = slope*Xset + q;

        % Pick random points and make them outliers randomly
        ids = randi(length(y), round(Nsamples/15), 1);
        y(ids) = y(ids) + mean(y);

        % Test L1, L2 regressions and smoothdata
        Norder = 1;
        H = zeros(length(Xset), Norder+1);

        H(:, 1) = 1;
        H(:, 2) = Xset;

        tic
        cL2_MATLAb = H\y;
        toc

        tic
        cL2_cvx = Regress1DLpNorm(y, Xset, 2);
        toc
        
        tic
        cL1_cvx = Regress1DLpNorm(y, Xset, 1);
        toc

        yL2_MATLAB = cL2_MATLAb(2) * Xset + cL2_MATLAb(1);
        yL2cvx = cL2_cvx(2) * Xset + cL2_cvx(1);
        yL1cvx = cL1_cvx(2) * Xset + cL1_cvx(1);


        figure;
        plot(Xset, ytruth, 'k-', 'LineWidth', 1.05, 'DisplayName', 'Truth')
        hold on;
        scatter(Xset, y, 8, "black", '*', 'DisplayName', 'Samples')
        plot(Xset, yL2_MATLAB, 'r-', 'LineWidth', 1.05, 'DisplayName', 'L2 MATLAB')
        plot(Xset, yL2cvx, 'b-', 'LineWidth', 1.05, 'DisplayName', 'L2 CVX')
        plot(Xset, yL1cvx, 'g-', 'LineWidth', 1.05, 'DisplayName', 'L1 CVX')

        title('Test L2 MATLAB vs L2 and L1 CVX')
        legend()
        DefaultPlotOpts();
        axis padded
        xlabel('X')   
        ylabel('Y')

    case 8
        %% TEST 8: Test and validation of POLRLS
        % Setup experiment as in paper:
        % Tracking time varying parameters with local regression
        % ARX model 
        
        % True parameters
        a = 0.7;
        b = @(i) 5 + 4*sin(2*pi/1000 * i);


        dt = 1;
        Nt = 1000;
        yi = zeros(Nt, 1);
        yprev = 0;

        % POLRLS options and initialization
        aOrder = 0;
        bOrder = 1;
        regressOrder = 1;

        alpha = 1e5;
        paramsDescriptorsOrders = [aOrder; bOrder];
        Ncoeffs = regressOrder+1 + sum(paramsDescriptorsOrders);
        P0 = alpha * eye(Ncoeffs);
        Mk = zeros(Ncoeffs, Ncoeffs, Nt);
        Mk(:, :, 1) = P0^-1;
        paramsHatk(:, 1) = [0; 0];

        PhiHatkPrev = [0; 0; 0];

        % Simulation loop
        for i = 1:Nt-1
            % Simulate noise
            ei = randn(1);
            zi(i) = randn(1); %#ok<*SAGROW> 

            % Generate simulated response (ARX model)
            yi(i) = a*yprev + b(i) * zi(i) + ei;
            yprev = yi(i);

            [paramsHatk(:, i+1), PhiHatk, Mk(:,:,i+1)] = POLRLSstep(yi(i), zi(i),...
                i, Mk(:, :, i), paramsHatk(:, i), PhiHatkPrev, regressOrder,...
                paramsDescriptorsOrders, 0.9999);

            PhiHatkPrev = PhiHatk;
        end


        figure;
        plot(1:Nt, yi, 'Color', 'r', 'LineWidth', 1.05)
        xlabel('Time sample [-]')
        ylabel('Response value')
        DefaultPlotOpts();
        title('Response signal')

        figure;
        subplot(2, 1, 1)
        plot(1:Nt, paramsHatk(1, :), '--', 'Color', 'r', 'LineWidth', 1.05, 'DisplayName', 'Estimated')
        hold on;
        plot(1:Nt, a*ones(Nt, 1), 'Color', 'b', 'LineWidth', 1.05, 'DisplayName', 'True value');
        xlabel('Time sample [-]')
        ylabel('Parameter a [-]')
        title('Parameter a constant')
        DefaultPlotOpts();
        axis padded

        subplot(2, 1, 2)
        plot(1:Nt, paramsHatk(2, :), '--', 'Color', 'r', 'LineWidth', 1.05, 'DisplayName', 'Estimated')
        hold on;
        plot(1:Nt, b(1:Nt), 'Color', 'b', 'LineWidth', 1.05, 'DisplayName', 'True value');
        xlabel('Time sample [-]')
        ylabel('Parameter b [-]')
        title('Parameter b time-varying')
        DefaultPlotOpts();
        axis padded
    case 9
        %% TEST 9: GivensRotSRIF
        % Example from Tapley chapter 5.6.4

        % Inputs
        i_bNO_PRIOR_INFO = false;
        i_dYobs = [-1.1; 1.2; 1.8];
        i_dHobsMatrix = [1, -2;
                         2, -1;
                         1, 1];

        i_dxPrior = [2; 2];
        Pcov = diag([100, 100]);

        Sroot = chol(Pcov, 'lower');      
        i_dSRInfoMatPrior = Sroot^-1;

        i_bRUN_WHITENING = false;
        i_dMeasCovSR = diag([1, 1, 1]);

        tic
        % Algorithm test
        [o_dxPost, o_dSRInfoMatPost, o_dInfoVecPost, o_dErrorVec] = GivensRotSRIF(i_dxPrior, ...
            i_dSRInfoMatPrior, ...
            i_dYobs, ...
            i_dHobsMatrix, ...
            i_bNO_PRIOR_INFO, ...
            i_bRUN_WHITENING,...
            i_dMeasCovSR);

        toc

        % MEX EQUIVALENCE TEST: PASSED
        tic
        [o_dxPost_MEX, o_dSRInfoMatPost_MEX, o_dInfoVecPost_MEX, o_dErrorVec_MEX] = GivensRotSRIF_MEX(i_dxPrior, ...
            i_dSRInfoMatPrior, ...
            i_dYobs, ...
            i_dHobsMatrix, ...
            i_bNO_PRIOR_INFO, ...
            i_bRUN_WHITENING,...
            i_dMeasCovSR);
        toc

        Spost = o_dSRInfoMatPost^-1;
        o_dPxPost = Spost * Spost';

        o_dxPost         - o_dxPost_MEX
        o_dSRInfoMatPost - o_dSRInfoMatPost_MEX
        o_dInfoVecPost   - o_dInfoVecPost_MEX
        o_dErrorVec      - o_dErrorVec_MEX



    case 10
        %% TEST 10: GivensRotEKF
        % Example from Tapley chapter 5.6.4 (using GivensRotSRIF)

        % Inputs
        i_bNO_PRIOR_INFO = false;
        i_dYobs = [-1.1; 1.2; 1.8];
        i_dHobsMatrix = [1, -2;
                         2, -1;
                         1, 1];

        i_dxPrior = [2; 2];
        Pcov = diag([100, 100]);

        Sroot = chol(Pcov, 'lower'); 
        
        i_ui8FILTER_TYPE = 2;
        % Test to do for complete coverage:
        % 0) Full covariance passing Pcov
        % 1) UD passing UD decomposition
        % 2) Square root covariance passing Sroot
        % All must return the same result


        tic
        if i_ui8FILTER_TYPE == 0
            inputTest = Pcov;
            i_dDxPrior = zeros(2);
        elseif i_ui8FILTER_TYPE == 2
            inputTest = Sroot;
            i_dDxPrior = zeros(2);
        elseif i_ui8FILTER_TYPE == 1
            [inputTest, i_dDxPrior] = UDdecomposition(Pcov);
        end

        i_bRUN_WHITENING = false;
        i_dMeasCovSR = diag([1, 1, 1]);

        % Algorithm test
        [o_dxPost, o_dPxPost, o_dDxPost] = GivensRotEKF(i_dxPrior, ...
            inputTest, ...
            i_dYobs, ...
            i_dHobsMatrix, ...
            i_bNO_PRIOR_INFO, ...
            i_bRUN_WHITENING, ...
            i_dMeasCovSR, ...
            i_ui8FILTER_TYPE, ...
            i_dDxPrior);

        toc

        i_dSRInfoMatPrior = Sroot^-1;

        i_bRUN_WHITENING = false;
        i_dMeasCovSR = diag([1, 1, 1]);

        % Algorithm test
        [o_dxPost_ref, o_dSRInfoMatPost_ref, o_dInfoVecPost_ref, o_dErrorVec_ref] = GivensRotSRIF(i_dxPrior, ...
            i_dSRInfoMatPrior, ...
            i_dYobs, ...
            i_dHobsMatrix, ...
            i_bNO_PRIOR_INFO, ...
            i_bRUN_WHITENING,...
            i_dMeasCovSR);

       
        o_dxPost - o_dxPost_ref
        Spost = o_dSRInfoMatPost_ref^-1;
        o_dPxPost_ref = Spost * Spost';

        % Discrepancy between GivensRotSRIF and GivensRotEKF should be a numerical zero
        if i_ui8FILTER_TYPE == 2
            o_dPxPost* o_dPxPost'- o_dPxPost_ref %#ok<*NOPTS> 
        elseif i_ui8FILTER_TYPE == 0
            o_dPxPost - o_dPxPost_ref
        elseif i_ui8FILTER_TYPE == 1
            o_dPxPost * o_dDxPost * o_dPxPost' - o_dPxPost_ref
        end

    case 11
        % TODO: check implementation of cost function for lsqnonlin
        rng("default");

        bADD_ATT_NOISE = 1;
        bADD_PIX_NOISE = 1;

        %% TEST 11: Feature LS-based initialization using Gauss-Newton iterations
        % Truth position 
        cspice_kclear()

        % Lunar Reconnaissance Orbiter kernel
        if strcmpi(computer("arch"), "glnxa64")
            LRO_KERNELS_PATH = fullfile("/mnt", "c", "Users", "pietr", "OneDrive - Politecnico di Milano",...
                "PoliMi - LM","MATLABwideCodes","SPICE_KERNELS","LRO_kernels");

        else
            LRO_KERNELS_PATH = 'C:\\Users\\pietr\\OneDrive - Politecnico di Milano\\PoliMi - LM\\MATLABwideCodes\\SPICE_KERNELS\\LRO_kernels';

        end
        cspice_furnsh( char(fullfile(LRO_KERNELS_PATH, 'lrorg_2023166_2023258_v01.bsp')) );
        cspice_furnsh( char(fullfile(LRO_KERNELS_PATH, 'moon_pa_de440_200625.bpc')) );
        cspice_furnsh( char(fullfile(LRO_KERNELS_PATH, 'moon_de440_200625.tf')) );

        cspice_furnsh(fullfile(which('SPICE_MK_DEFAULT.tm'))) % THIS ONLY WORKS FOR WINDOWS DUE TO MK FORMAT
        cspice_ktotal('all')

        R_Moon = mean(cspice_bodvrd('MOON', 'RADII', 3)) % [km]
        mu_Moon = cspice_bodvrd('MOON', 'GM', 1) % [km^2/s^3]

        
        % LRO WAC parameters
        Npixels = 1008;
        focalLength = 6; % [mm]
        pixSize = 9e-3; % [mm]
        UVopticalCentre = [Npixels/2; Npixels/2];
        skewFactor = 0;

        % Field of View (full): 90°

        % Camera parameters (intrinsic and pose)
        i_dKcam = [focalLength/pixSize, skewFactor, UVopticalCentre(1);
                   0, focalLength/pixSize, UVopticalCentre(2);
                   0,     0,  1]; 

%         i_dKcam = [10560, 0, 510;
%             0, 10560, 510;
%             0,     0,  1]; % Hera NAVCAM

%         SCaltitude = 400; % [km]
%         i_dRSC_IN = (R_Moon + SCaltitude) * posFeatureDir_TB; % [km] Above feature
        
        % Generate random velocity vector assuming random plane circular orbit

%         Vnorm = sqrt(mu_Moon/norm(i_dRSC_IN + SCaltitude)); 
%         randOrthogonal = randn(3, 1);
%         Vdir = randOrthogonal - dot(posFeatureDir_TB, randOrthogonal) * posFeatureDir_TB;
%         Vdir = Vdir/norm(Vdir); % Direction orthogonal to position vector
%         i_dVSC_TB = Vnorm * Vdir;

        % Simulate centroiding measurements of features from their true
        % position on the body.
        i_dqCAMwrtSCB = [0; 0; 0; 1];
        i_dSigmaAPE = deg2rad(0.1); % [rad]
        sigmaIP = 2; % [pix]

        MAX_NUM_IMAGES = 1000;
        IMG_ACQ_DELTA_TIME = 1; % [s]

        tspan = 0:IMG_ACQ_DELTA_TIME:( (MAX_NUM_IMAGES-1)*IMG_ACQ_DELTA_TIME); % [s]

        % Generate Ephemerides of SC and Moon attitude
        dateString = 'Aug 30 9:45:10 PDT 2023';
        et0 = cspice_str2et(dateString);
        timeETgrid = tspan + et0;

        i_dXSC_IN_t0 = cspice_spkezr('LUNAR RECONNAISSANCE ORBITER', timeETgrid(1), 'J2000', 'NONE', 'MOON')';
      
        i_dXSC_IN_t0 = [2 * i_dXSC_IN_t0(1:3), sqrt(mu_Moon/norm(i_dXSC_IN_t0(1:3))) * i_dXSC_IN_t0(4:6)/norm(i_dXSC_IN_t0(4:6))]';

        dDCM_Moon_fromTBtoIN = cspice_pxform('MOON_PA_DE440', 'J2000', timeETgrid);

        [~, i_dXSC_IN] = ode113(@(t, x) RHS_2BP(t, x, mu_Moon), tspan, i_dXSC_IN_t0, odeset('RelTol', 1e-12, 'AbsTol', 1e-12));
       

        % Feature measurements simulation
        % Construct feature position unit vector from position at 1/3 time
        % 
        tFeatureLocSim = round(length(i_dXSC_IN)/2);

        posFeatureDir_TB = dDCM_Moon_fromTBtoIN(:, :, tFeatureLocSim)' *...
        (i_dXSC_IN(tFeatureLocSim, 1:3)'/norm(i_dXSC_IN(tFeatureLocSim, 1:3)) ); 

        % Translate to Moon surface assumed spherical in TB (fixed)
        posFeature_TB = R_Moon * posFeatureDir_TB;
  
        % Initialize data arrays
        UVpixCoordsMeas = zeros(length(tspan), 2);
        o_dqCAMwrtIN_hat = zeros(length(tspan), 4);
        o_dqCAMwrtIN_ref = zeros(length(tspan), 4);
        dDCM_fromINtoCAM = zeros(3, 3, length(tspan));
        i_drFeature_IN = zeros(length(tspan), 3);

        Nimages = 0;
        idt = 0;
        FEATURE_IN_VIEW = true;

        % Introduce "is in view" flag to stop the 
        while FEATURE_IN_VIEW && Nimages < MAX_NUM_IMAGES 
            
            idt = idt + 1;

            % Use the feature position as if it was the centre of the body
            % to simulate pointing to the it
            i_drFeature_IN(idt, :) = dDCM_Moon_fromTBtoIN(:, :, idt) * posFeature_TB;

            % Simulate pointing to feature and attitude error
            i_drTargetPoint_IN = [0; 0; 0]; % Nadir pointing since origin is Moon centre

            [o_dqSCBwrtIN_ref, o_dqSCBwrtIN_hat, o_dqCAMwrtIN_ref(idt, :), o_dqCAMwrtIN_hat(idt, :)] =...
                      simulateTBPointing(i_dXSC_IN(idt, :)', ...
                      i_drTargetPoint_IN, ...
                      i_dqCAMwrtSCB, ...
                      bADD_ATT_NOISE*i_dSigmaAPE);

            % Project feature position from Inertial frame to image plane
            % knowing the true camera pose and calibration matrix

            %% DEVNOTE: error here: (u, v) coordinates cannot result outside 
            % the image plane if the camera is pointing to the feature!
            [trueFeatureProj, dDCM_fromINtoCAM(:, :, idt)] = pinholeProject(i_dKcam, o_dqCAMwrtIN_ref(idt, :)',...
                i_dXSC_IN(idt, 1:3)', i_drFeature_IN(idt, :)');

            if or(trueFeatureProj(1) > 1020, trueFeatureProj(1) < 0) || ...
                    or(trueFeatureProj(2) > 1020, trueFeatureProj(2) < 0)

                FEATURE_IN_VIEW = false;
                break;
            end

            UVpixCoordsMeas(idt, :) = trueFeatureProj + bADD_PIX_NOISE * (sigmaIP * rand(2, 1));

            Nimages = Nimages + 1;
        end


        %% DEBUG PLOTS
        
        i_dXSC_TB = zeros(size(i_dXSC_IN(:, 1:3) ));
        featurePhaseAngle = zeros(length(tspan), 1);

        for idt = 1:length(tspan)
            % Compute trajectory in TB
            i_drSC_TB(idt, :) = dDCM_Moon_fromTBtoIN(:, :, idt)' * i_dXSC_IN(idt, 1:3)';

            % Compute angle between camera boresight and feature position wrt SCB
            featureRelPos = i_drFeature_IN(idt, 1:3)' - i_dXSC_IN(idt, 1:3)';

            boresightCamAxis = dDCM_fromINtoCAM(3, 1:3, idt)';
            featurePhaseAngle(idt) = acosd(dot(featureRelPos./norm(featureRelPos), boresightCamAxis));
        end

        opts.Units = 'km';
        planet3D('Moon', opts);
        hold on;
        % Plot trajectory
        plot3(i_drSC_TB(1, 1), i_drSC_TB(1, 2), i_drSC_TB(1, 3), 'ko', 'MarkerSize', 5, 'LineWidth', 3);
        plot3(i_drSC_TB(:, 1), i_drSC_TB(:, 2), i_drSC_TB(:, 3), 'b-', 'LineWidth', 1.1);
        % Plot feature in target body
        plot3(posFeature_TB(1), posFeature_TB(2), posFeature_TB(3), 'rx', 'MarkerSize', 5, 'LineWidth', 1.05)

        DefaultPlotOpts();
        axis equal;
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('Visualization in Moon Centred Moon Fixed frame')
        % Visualize spacecraft motion
        %         plotAttitudeQuat(o_dqCAMwrtIN_ref', ...
        %             i_dXSC_IN(:, 1:3)', ...
        %             true, ...
        %             true, ...
        %             true, ...
        %             0.001);


        figure;
        cmap = jet(size(UVpixCoordsMeas, 1));
        scatter(UVpixCoordsMeas(:, 1), UVpixCoordsMeas(:, 2), 1, cmap)
        DefaultPlotOpts();
        axis equal;
        xlabel('X [pix]')
        ylabel('Y [pix]');     
        ax = gca;
        ax.YDir = 'reverse';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get ID of last image in which feature was visible
        lastImgID = Nimages;

        % MEASUREMENT FUNCTION 
        Yobs = reshape(UVpixCoordsMeas(1:lastImgID, :)', [], 1); % Reshape to column vector [x1; y1; ...; xN; yN]
        % Compute estimated 3D position of feature using lsqnonlin
        lsqnonlinOpts = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxIterations', 2000, ...
            'FunctionTolerance', 1e-12, 'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', 2e5);
        
        % Initial guess: [point projection at t0]

        dist2FeatureGuess = norm(i_dXSC_IN(1, 1:3)) - R_Moon;
        dir2FeatureGuess = Pixel2LoS_NoDistorsion(UVpixCoordsMeas(1, :)', i_dKcam);

%         x0FeaturePosEP_CAMk = dist2FeatureGuess * dir2FeatureGuess; % Initial guess of position of feature in CAM frame at t0
   
        % Set "truth" as initial guess
        x0FeaturePosEP_CAMk = transpose(Quat2DCM(o_dqCAMwrtIN_ref(idt, :), true)') * ...
                              (dDCM_Moon_fromTBtoIN(:, :, 1) * posFeature_TB - i_dXSC_IN(1, 1:3)');

        x0FeaturePosIDP_CAMk = transformEPtoIDP(x0FeaturePosEP_CAMk);

        Nt = size(Yobs, 1)/2
        % Relative displacements wrt initial position and attitude
        dr_CkwrtCi_Ci = zeros(3, Nt);
        DeltaQuat_CkwrtCi = zeros(4, Nt);

        for idt = 1:Nt
            % Compute relative displacement wrt initial position in each Ci frame
            dr_CkwrtCi_Ci(:, idt) = Quat2DCM(o_dqCAMwrtIN_hat(idt, :), true)' * (i_dXSC_IN(1, 1:3) - i_dXSC_IN(idt, 1:3))';

            % Compute relative attitudes Ck wrt Ci
            % q_CkwrtCi = q_INwrtCi qCross q_CkwrtIN
            DeltaQuat_CkwrtCi(:, idt) = qCross( qInvert(o_dqCAMwrtIN_hat(idt, :)', true), o_dqCAMwrtIN_hat(1, :)');
        end

        tic
        % Solve for optimal inverse depth parameters (for camera at initial
        % position t0)
        [invDepthParamsEst_LSQNONLIN, resnorm, residual, exitflag, output, lambda, jacobian] =...
            lsqnonlin(@(invDepthParamsEst) findFeatureCostFcn(invDepthParamsEst, Yobs,...
            Nt, DeltaQuat_CkwrtCi, dr_CkwrtCi_Ci), x0FeaturePosIDP_CAMk, [], [], lsqnonlinOpts);

        lsqnonlinTOC = toc

        % Convert inverse depth parameters to global frame position at t0
        featurePos_CAMk = transformEPtoIDP(invDepthParamsEst_LSQNONLIN);
        estFeaturePos_IN = Quat2DCM(o_dqCAMwrtIN_hat(1, :), true) * featurePos_CAMk + i_dXSC_IN(1, 1:3)';
        
        %% MODULE UNDER TESTING
        i_bIS_JPL_CONV = true;
        i_ui16EstTimeID = 1;
        i_dDeltaResNormRelTol = 1e-3;
        i_dCamCentreCoord = [Npixels/2; Npixels/2];
        i_ui8MaxIter = 10;

        i_dFtPosVecGuess_CAMk = x0FeaturePosEP_CAMk;

%         i_dyMeas = reshape(Yobs, 2, []);
        i_dyMeas = Yobs; 
   
        tic

        [o_dFeaturePosVec_CAM, o_dFeaturePosVec_IN, o_dMeasVec] = locateFeature_NewtonRaphson(i_dyMeas, ...
            o_dqCAMwrtIN_hat(1:Nt, :), ...
            i_dXSC_IN(1:Nt, 1:3)', ...
            i_ui16EstTimeID, ...
            i_dFtPosVecGuess_CAMk, ...
            i_dCamCentreCoord, ...
            i_bIS_JPL_CONV, ...
            i_dDeltaResNormRelTol, ...
            i_ui8MaxIter);

        locateFeatureTOC = toc

        diffwrtLSQNONLIN = o_dFeaturePosVec_IN - estFeaturePos_IN

        % Position error of estimated feature
        errFeaturePos_TB = posFeature_TB - dDCM_Moon_fromTBtoIN(:, :, 1)' * estFeaturePos_IN

    case 12
        %% TEST 12: Test AdaptQCovASNC runtime (verification)
        
        % Load data from Hera EXP test bench for testing
        load("testdataASNC.mat");

        i_dProcessNoiseCovSNC = zeros(3); % Process noise covariance matrix (acceleration)
        i_dTimeStep = 1; % [s]
        windowSize = 100;
        i_dQlowerBound = 0;
        i_dQupperBound = 10;
        i_dPxxPostBuffer = zeros(6, 6, 1); % Not needed for now
        idCount = 1;
        runtime = 0;

        for idt = tspan(30000:40000)


            i_dDeltaStatesBuffer = zeros(6, windowSize);
            i_dKalmanGainBuffer = zeros(6, 3, windowSize);
            i_dPyyResCovBuffer = zeros(3, 3, windowSize);
            idAlloc = 1;

            % Compute state corrections buffer
            dt = 1;
            dtMax = dt;
            for idtBuffer = idt:-1:idt-windowSize+1

                if  any(priorRes_ts(idtBuffer, :), 'all')
                     i_dDeltaStatesBuffer(:, idAlloc)    = KalmanGain_ts(:, :, idtBuffer)*priorRes_ts(idtBuffer, :)';
                     i_dKalmanGainBuffer(:, :, idAlloc) = KalmanGain_ts(:, :, idtBuffer);
                     i_dPyyResCovBuffer(:, :, idAlloc)  = resSyy_ts(:, :, idtBuffer) * resSyy_ts(:, :, idtBuffer)';

                    idAlloc = idAlloc + 1;
                    if dtMax < dt
                        dtMax = dt;
                    end
                    dt = 1;
                else
                    dt = dt + 1;
                end
            end

            tic
            if any(i_dKalmanGainBuffer, 'all') 

                i_dTimeStep = dtMax;
                [o_dDiscreteQposVelSNC, o_dProcessNoiseCovSNC] = adaptQCovASNC(i_dProcessNoiseCovSNC, ...
                    i_dTimeStep, ...
                    i_dKalmanGainBuffer(:, :, 1:idAlloc-1), ...
                    i_dPyyResCovBuffer(:, :, 1:idAlloc-1), ...
                    i_dQlowerBound, ...
                    i_dQupperBound, ...
                    i_dDeltaStatesBuffer(:, 1:idAlloc-1), ...
                    i_dPxxPostBuffer);

            else
                o_dProcessNoiseCovSNC = i_dProcessNoiseCovSNC;
            end
            runtime = runtime + toc;
            idCount = idCount + 1;
        end

        runtime = runtime/idCount;
        
    case 13
        %% TEST 13: Expectation Maximization for GMM fitting and testing of Covariance generator
        % Generate 2 Gaussian distributions realizations
        Npdfs = 2;
        dimSpace = 2;

        MeanVec = zeros(dimSpace, Npdfs);
        CovMatArray = zeros(dimSpace, dimSpace, Npdfs);

        % Gaussian 1
        MeanVec(:, 1) = [4, -5];
        CovMatArray(:, :, 1) = eye(2);
        % Gaussian 2
        MeanVec(:, 2) = [-14, 6.7];
        CovMatArray(:, :, 2) = eye(2);

        % Randomly generate a correlated covariance
        % CHECK IF I ALREADY CODED A FUNCTION

    case 14
        %% TEST 14: Test Weighted Modified Gram Schmidt and UD Time Update
        % Test case used: 
        
        A = randn(50);

        % Decompose matrix in UD factors
        [i_dW, i_dD] = UDdecomposition(A);

        % LEVEL-1 TEST
        [o_dU, o_dV] = orthogonalizeUD_WMGS(i_dW, i_dD); % Seems ok
    
        all((eye(50) - o_dV' * o_dV == 0), 'all')
     
        % LEVEL-2 TEST
        % [o_dUprior, o_dDprior] = UDCov_TimeUp(i_dSTM, i_dUpost, i_dDpost, i_dProcessCov)


        return
        % LEVEL-3 TEST
        % [outputArg1, outputArg2] = EKF_UDcov_TimeUp(inputArg1,inputArg2)
        
        % LEVEL-4 TEST
        % UDCov_MaxEffTimeUp

    case 15
        %% TEST 15: Schur Marginalization
        i_dCovMatrix = getRandomCov([10.2, 15, 3.4, 3.5, 6.7, 9.25, 0.5, 1.2]);
        
        i_ui16FirstMargStateIdx = 3;
        [o_dMarginalCov] = ShurMarginalization(i_dCovMatrix, i_ui16FirstMargStateIdx);

    case 16
        %% TEST 16: Covariance recovery algorithm by (Kaess, Dellaert 2009)
        % Example from Tapley chapter 5.6.4
        diagSigma = [2.1, 5, 4.3, 8.9, 3.8];

        covMatrix = getRandomCov(diagSigma);
        InfoMatrix = eye(length(diagSigma))/covMatrix;

        i_ui16ExtrIndices = [1, 2];

        marginalCov = ShurMarginalization(covMatrix, 4); % TO FIX! --> the second argument is non-sense like this

        o_dCovSubMatrix = RecoverCovFromSRInfoMat(InfoMatrix, i_ui16ExtrIndices) % Not working 
end

%% LOCAL FUNCTION
% Temporary for testing: Inverse Depth Model
function UVpixPrediction = inverseDepthPinhole(DeltaQuat_CkwrtCi, dr_CkwrtCi_Ci, invDepthParam)

alpha = invDepthParam(1);
beta = invDepthParam(2);
rho = invDepthParam(3);

% Inverse depth model
Hvector = Quat2DCM(DeltaQuat_CkwrtCi, true) * [alpha; beta; 1] + rho * dr_CkwrtCi_Ci;
UVpixPrediction = 1/Hvector(3) * [Hvector(1); Hvector(2)] + 1008/2;

end

% Cost function for Gauss-Newton algorithm
function resY = findFeatureCostFcn(invDepthParam, Yobs, Nt, DeltaQuat_CkwrtCi, dr_CkwrtCi_Ci)

% DeltaQuat_CkwrtCi, dr_CkwrtCi_Ci

% Predict feature location in image plane
IDpos = 1;
UVpixPrediction = nan(2*Nt, 1);

for idt = 1:Nt
    UVpixPrediction(IDpos:IDpos+1) = inverseDepthPinhole(DeltaQuat_CkwrtCi(:, idt), dr_CkwrtCi_Ci(:, idt), invDepthParam);
    IDpos = IDpos + 2;
end

if any(isnan(UVpixPrediction))
    error('nan detected')
end

% Compute residuals
resY = Yobs - UVpixPrediction;

end





