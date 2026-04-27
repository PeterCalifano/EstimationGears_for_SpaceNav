classdef CDensityFcnPropagator 
    %% DESCRIPTION
    % Class to collect and implement methods to perform uncertainty propagation, with specific focus to
    % the first two moments (mean and covariance). Current implementation provides sigma points transform
    % and Monte Carlo propagation methods.
    % -------------------------------------------------------------------------------------------------------------
    %% CHANGELOG
    % 18-08-2025    Pietro Califano     First implementation prototype.
    % -------------------------------------------------------------------------------------------------------------
    %% METHODS
    % Method1: Description
    % -------------------------------------------------------------------------------------------------------------
    %% PROPERTIES
    % Property1: Description, dtype, nominal size
    % -------------------------------------------------------------------------------------------------------------
    %% DEPENDENCIES
    % [-]
    % -------------------------------------------------------------------------------------------------------------
 
    properties (Access = public)
        
    end
    
    methods (Access = public)
        % CONSTRUCTOR
        function self = CDensityFcnPropagator()
            
        end

        % GETTERS


        % SETTERS
        
   
    end

    methods(Static, Access = private)
            function [dxCurrentCov] = PropagateLinCov_internal_(dxCurrentCov, dFlowSTM)%#codegen
            arguments
                dxCurrentCov    (:,:) double {ismatrix, mustBeReal}
                dFlowSTM        (:,:) double {ismatrix, mustBeReal}
            end
            
            assert( all(size(dxCurrentCov) == size(dFlowSTM), 'all'), 'ERROR: STM and covariance does not match in size');
            dxCurrentCov = dFlowSTM * dxCurrentCov * transpose(dxCurrentCov);
        end
    end
    
    methods (Static, Access = public)

        %%% Linear Covariance methods
        function [dxMeanOut, dxCovarianceOut] = PropagateFiniteDiffLinCovPropagateDyn(dxMeanIn, ...
                                                                                      dxCovarianceIn, ...
                                                                                      strDynParams, ...
                                                                                      strFilterMutabConfig, ...
                                                                                      strFilterConstConfig, ...
                                                                                      dTimestampStart, ...
                                                                                      dTimestampEnd, ...
                                                                                      dDeviationSize, ...
                                                                                      dIntegrTimestep)
            arguments
                dxMeanIn          (:,1) double {isvector, mustBeReal}
                dxCovarianceIn    (:,:) double {ismatrix, mustBeReal}
                strDynParams         (1,1) {isstruct}
                strFilterMutabConfig (1,1) {isstruct}
                strFilterConstConfig (1,1) {isstruct}
                dTimestampStart   (1,1) double {isscalar, mustBeReal} = -1.0
                dTimestampEnd     (1,1) double {isscalar, mustBeReal, mustBeGreaterThan(dTimestampEnd, dTimestampStart)} = 0.0
                dDeviationSize    (1,1) double {mustBeGreaterThan(dDeviationSize, 0.0)} = 1e-5;
                dIntegrTimestep   (1,1) double {isscalar, mustBeReal} = 1.0
            end

            % Define function handle for PropagateDyn function
            fcnDynHandle = @(dxSamplePoint, dTimestart, dDeltaTime, dIntegrTimestep) PropagateDyn(dxSamplePoint, ...
                                                                                        dTimestart, ...
                                                                                        dDeltaTime, ...
                                                                                        dIntegrTimestep, ...
                                                                                        strDynParams, ...
                                                                                        strFilterMutabConfig, ...
                                                                                        strFilterConstConfig);

            % Call handle version
            [dxMeanOut, dxCovarianceOut] = CDensityFcnPropagator.PropagateFiniteDiffLinCovHandle(fcnDynHandle, ...
                                                                            dxMeanIn, ...
                                                                            dxCovarianceIn, ...
                                                                            dTimestampStart, ...
                                                                            dTimestampEnd, ...
                                                                            dDeviationSize, ...
                                                                            dIntegrTimestep);
        end
        
        function [dxMeanOut, dxCovarianceOut] = PropagateFiniteDiffLinCovHandle(fcnHandle, ...
                                                                                dxMeanIn, ...
                                                                                dxCovarianceIn, ...
                                                                                dTimestampStart, ...
                                                                                dTimestampEnd, ...
                                                                                dDeviationSize, ...
                                                                                varargin)%#codegen
            arguments (Input)
                fcnHandle         (1,1) {mustBeA(fcnHandle, "function_handle")}
                dxMeanIn          (:,1) double {isvector, mustBeReal}
                dxCovarianceIn    (:,:) double {ismatrix, mustBeReal}
                dTimestampStart   (1,1) double {isscalar, mustBeReal} = -1.0
                dTimestampEnd     (1,1) double {isscalar, mustBeReal, mustBeGreaterThan(dTimestampEnd, dTimestampStart)} = 0.0
                dDeviationSize    (1,1) double {mustBeGreaterThan(dDeviationSize, 0.0)} = 1e-5;
            end
            arguments (Input, Repeating)
                varargin
            end
            arguments (Output)
                dxMeanOut          (:,1) double {ismatrix, mustBeReal}
                dxCovarianceOut    (:,:) double {ismatrix, mustBeReal}
            end

            % Initialize output variables
            dxMeanOut       = zeros(size(dxMeanIn));
            dxCovarianceOut = zeros(size(dxCovarianceIn));

            dIntegrTimestep = 0.0; % Default value to allow parfor usage
            if ~isempty(varargin)
                dIntegrTimestep = varargin{1};
            end

            dFiniteDiffSTM = zeros(length(dxMeanOut), length(dxMeanOut));
            d2ndOrdFiniteDiffDenom = 1/(2*dDeviationSize);

            % Propagate mean state
            if dTimestampStart > 0 - eps
                % Input point + timestamps
                dxMeanOut(:)  = fcnHandle(dxMeanIn, dTimestampStart, ...
                    dTimestampEnd - dTimestampStart, ...
                    dIntegrTimestep);

            else
                % Input point only signature
                dxMeanOut(:)  = fcnHandle(dxMeanIn);

            end

            % Loop over input states to evaluate function
            parfor ii = 1:length(dxMeanIn)

                dEpsVec = zeros(size(dxMeanIn));
                dEpsVec(ii) = dDeviationSize;

                % Compute plus-minus function vector values
                if dTimestampStart > 0 - eps
                    % Input point + timestamps
                    dfPlus  = fcnHandle(dxMeanIn + dEpsVec, dTimestampStart, ...
                                        dTimestampEnd - dTimestampStart, ...
                                        dIntegrTimestep);

                    dfMinus = fcnHandle(dxMeanIn - dEpsVec, dTimestampStart, ...
                                        dTimestampEnd - dTimestampStart, ...
                                        dIntegrTimestep);
                else
                    % Input point only signature
                    dfPlus  = fcnHandle(dxMeanIn + dEpsVec);
                    dfMinus = fcnHandle(dxMeanIn - dEpsVec);
                end


                if isvector(dfPlus)
                    dfPlus = dfPlus(:);
                end

                if isvector(dfMinus)
                    dfMinus = dfMinus(:);
                end

                dDiffVec = reshape(dfPlus - dfMinus, [], 1);
                dFiniteDiffSTM(:, ii) = d2ndOrdFiniteDiffDenom * (dDiffVec);
            end

            % Propagate covariance
            dxCovarianceOut(:,:) = CDensityFcnPropagator.PropagateLinCov_internal_(dxCovarianceIn, dFiniteDiffSTM);

        end


        %%% Sigma Points method
        function [dxMeanOut, dxCovarianceOut] = PropagateSigmaPointTransformPropagateDyn(dxMeanIn, ...
                                                                                        dxCovarianceIn, ...
                                                                                        strDynParams, ...
                                                                                        strFilterMutabConfig, ...
                                                                                        strFilterConstConfig, ...
                                                                                        dTimestampStart, ...
                                                                                        dTimestampEnd, ...
                                                                                        dIntegrTimestep, ...
                                                                                        enumVariantName, ...
                                                                                        dAlphaCoeff, ...
                                                                                        dBetaCoeff, ...
                                                                                        dKappaCoeff)%#codegen
            arguments (Input)
                dxMeanIn          (:,1) double {isvector, mustBeReal}
                dxCovarianceIn    (:,:) double {ismatrix, mustBeReal}
                strDynParams         (1,1) {isstruct}
                strFilterMutabConfig (1,1) {isstruct}
                strFilterConstConfig (1,1) {isstruct}
                dTimestampStart   (1,1) double {isscalar, mustBeReal} = 0.0
                dTimestampEnd     (1,1) double {isscalar, mustBeReal} = 1.0
                dIntegrTimestep   (1,1) double {isscalar, mustBeReal} = 1.0
                enumVariantName    (1,:) {mustBeMember(enumVariantName, ["scaled", "original"]), mustBeA(enumVariantName, ["char", "string"])} = "original"
                dAlphaCoeff        (1,1) double {mustBeGreaterThan(dAlphaCoeff, 0.0)} = 2e-3;
                dBetaCoeff         (1,1) double {mustBeGreaterThan(dBetaCoeff , 0.0)} = 2.0
                dKappaCoeff        (1,1) double {mustBeGreaterThanOrEqual(dKappaCoeff , 0.0)} = 0.0
            end
            arguments (Output)
                dxMeanOut          (:,1) double {ismatrix, mustBeReal}
                dxCovarianceOut    (:,:) double {ismatrix, mustBeReal}
            end

            % Define function handle for PropagateDyn function
            fcnDynHandle = @(dxSigmaPoint, dTimestart, dDeltaTime, dIntegrTimestep) PropagateDyn(dxSigmaPoint, ...
                                                                                        dTimestart, ...
                                                                                        dDeltaTime, ...
                                                                                        dIntegrTimestep, ...
                                                                                        strDynParams, ...
                                                                                        strFilterMutabConfig, ...
                                                                                        strFilterConstConfig);

            % Call handle version
            [dxMeanOut, dxCovarianceOut] = CDensityFcnPropagator.PropagateSigmaPointTransformHandle(fcnDynHandle, ...
                                                                                               dxMeanIn, ...
                                                                                               dxCovarianceIn, ...
                                                                                               dTimestampStart, ...
                                                                                               dTimestampEnd, ...
                                                                                               enumVariantName, ...
                                                                                               dAlphaCoeff, ...
                                                                                               dBetaCoeff, ...
                                                                                               dKappaCoeff, ...
                                                                                               dIntegrTimestep);%#codegen

        end

        function [dxMeanOut, dxCovarianceOut] = PropagateSigmaPointTransformHandle(fcnHandle, ...
                                                                                   dxMeanIn, ...
                                                                                   dxCovarianceIn, ...
                                                                                   dTimestampStart, ...
                                                                                   dTimestampEnd, ...
                                                                                   enumVariantName, ...
                                                                                   dAlphaCoeff, ...
                                                                                   dBetaCoeff, ...
                                                                                   dKappaCoeff, ...
                                                                                   varargin)%#codegen
            arguments (Input)
                fcnHandle         (1,1) {mustBeA(fcnHandle, "function_handle")}
                dxMeanIn          (:,1) double {isvector, mustBeReal}
                dxCovarianceIn    (:,:) double {ismatrix, mustBeReal}
                dTimestampStart   (1,1) double {isscalar, mustBeReal} = -1.0
                dTimestampEnd     (1,1) double {isscalar, mustBeReal, mustBeGreaterThan(dTimestampEnd, dTimestampStart)} = 0.0
                enumVariantName    (1,:) {mustBeMember(enumVariantName, ["scaled", "original"]), mustBeA(enumVariantName, ["char", "string"])} = "original"
                dAlphaCoeff        (1,1) double {mustBeGreaterThan(dAlphaCoeff, 0.0)} = 2e-3;
                dBetaCoeff         (1,1) double {mustBeGreaterThan(dBetaCoeff , 0.0)} = 2.0
                dKappaCoeff        (1,1) double {mustBeGreaterThanOrEqual(dKappaCoeff , 0.0)} = 0.0
            end
            arguments (Input, Repeating)
                varargin
            end
            arguments (Output)
                dxMeanOut          (:,1) double {ismatrix, mustBeReal}
                dxCovarianceOut    (:,:) double {ismatrix, mustBeReal}
            end

            % Initialize output variables
            dxMeanOut       = zeros(size(dxMeanIn));
            dxCovarianceOut = zeros(size(dxCovarianceIn));

            ui32DomainSize = uint32(size(dxMeanIn,1));

            dIntegrTimestep = 0.0; % Default value to allow parfor usage
            if ~isempty(varargin)
                dIntegrTimestep = varargin{1};
            end

            % Generate sigma points from input moments
            [dWeightsMean, dWeightsCov, dPerturbScale] = CDensityFcnPropagator.GenerateUnscentedWeights(ui32DomainSize, ...
                                                                                                    enumVariantName, ...
                                                                                                    dAlphaCoeff, ...
                                                                                                    dBetaCoeff, ...
                                                                                                    dKappaCoeff);

            ui32NumOfSigmaPoints = uint32(length(dWeightsMean));

           [dSigmaPointsSet] = CDensityFcnPropagator.GenerateSigmaPointsSet(dxMeanIn, ...
                                                     dxCovarianceIn, ...
                                                     ui32NumOfSigmaPoints, ...
                                                     dPerturbScale);

            % Run propagation through function handle
            parfor idCsi = 1:ui32NumOfSigmaPoints
                % Template:   dxSigmaPointsPrior(:, idCsi) = F(dxSigmaPoints(:, idCsi));

                if dTimestampStart > 0 - eps
                    % Sigma points + timestamps
                    dSigmaPointsSet(:, idCsi) = fcnHandle(dSigmaPointsSet(:, idCsi), ...
                                                        dTimestampStart, ...
                                                        dTimestampEnd - dTimestampStart, ...
                                                        dIntegrTimestep);

                else
                    % Sigma points only signature
                    dSigmaPointsSet(:, idCsi) = fcnHandle(dSigmaPointsSet(:, idCsi));

                end

            end
        

            % Compute sample mean and covariance
            [dxMeanOut(:), dxCovarianceOut(:,:)] = SumSigmaPointsToMoments(dSigmaPointsSet, dWeightsMean, dWeightsCov);

        end

        function [dWeightsMean, dWeightsCov, dPerturbScale] = GenerateUnscentedWeights(ui32DomainSize, ...
                                                                        enumVariantName, ...
                                                                        dAlphaCoeff, ...
                                                                        dBetaCoeff, ...
                                                                        dKappaCoeff)%#codegen
            arguments
                ui32DomainSize  (1,1) uint32 {mustBeGreaterThan(ui32DomainSize, 0)}
                enumVariantName (1,:) {mustBeMember(enumVariantName, ["scaled", "original"]), mustBeA(enumVariantName, ["char", "string"])}
                dAlphaCoeff     (1,1) double {mustBeGreaterThan(dAlphaCoeff, 0.0)} = 2e-3;
                dBetaCoeff      (1,1) double {mustBeGreaterThan(dBetaCoeff , 0.0)} = 2.0
                dKappaCoeff     (1,1) double {mustBeGreaterThanOrEqual(dKappaCoeff , 0.0)} = 0.0
            end

            % Define output
            NumSigmaPoints    = 2*double(ui32DomainSize) + 1;
            dWeightsMean = zeros(NumSigmaPoints, 1);
            dWeightsCov = zeros(NumSigmaPoints, 1);

            switch enumVariantName
                case "scaled"
                    % Scaled version Unscented Transform weights
                    % Reference: TODO


                    dLambda = dAlphaCoeff^2*(NumSigmaPoints + k) ;
                    dPerturbScale = sqrt(dLambda); % Square root covariance perturbation scale factor

                    dWeightsMean(1) = 1 - NumSigmaPoints / (dAlphaCoeff^2 * (NumSigmaPoints + dKappaCoeff));
                    dWeightsCov(1) = (2 - dAlphaCoeff^2 + beta) - NumSigmaPoints / (dAlphaCoeff^2 * (NumSigmaPoints + dKappaCoeff));
                    
                    dWeightsCov(2:end) = ( 1/(2* dAlphaCoeff^2 *(NumSigmaPoints + dKappaCoeff)) ) .* ones(1, 2*NumSigmaPoints); % * ones(1, 2*ui32Ncsi);
                    dWeightsMean(2:end) = dWeightsCov(2:end);


                case "original"
                    % Original Unscented Transform weights
                    % Reference: TODO

                    dDomainSize = double(ui32DomainSize);

                    dKappaCoeff      = 3 - dDomainSize;
                    dLambda = dAlphaCoeff^2 * (dDomainSize + dKappaCoeff) - dDomainSize;
                    dPerturbScale = sqrt(dDomainSize + dLambda); % Square root covariance perturbation scale factor

                    dWeightsMean(1) = dLambda./(dDomainSize + dLambda);
                    dWeightsCov(1) = dWeightsMean(1) + (1.0 - dAlphaCoeff^2 + dBetaCoeff);

                    dWeightsCov(2:end) = 1.0/(2.0* (dDomainSize + dLambda));
                    dWeightsMean(2:end) = dWeightsCov(2:end);

                otherwise
                    error('Invalid transform selection.')
            end
        end

        function [dSigmaPointsSet] = GenerateSigmaPointsSet(dxMean, ...
                                                            dxCovariance, ...
                                                            ui32NumOfSigmaPoints, ...
                                                            dPerturbScale)%#codegen
            arguments
                dxMean                  (:,1) double {mustBeReal, mustBeNumeric}
                dxCovariance            (:,:) double {mustBeReal, mustBeNumeric}
                ui32NumOfSigmaPoints    (1,1) uint32 {mustBeGreaterThan(ui32NumOfSigmaPoints, 0)}
                dPerturbScale           (1,1) double {mustBeReal, mustBeGreaterThan(dPerturbScale, 0.0)}
            end
            %% PROTOTYPE
            % [dSigmaPointsSet] = GenerateSigmaPointsSet(dxMean, ...
            %                                           dxCovariance, ...
            %                                           ui32NumOfSigmaPoints, ...
            %                                           dPerturbScale)%#codegen
            % -------------------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % 
            % -------------------------------------------------------------------------------------------------------------
            %% INPUT
            % dxMean                  (:,1) double {mustBeReal, mustBeNumeric}
            % dxCovariance            (:,:) double {mustBeReal, mustBeNumeric}
            % ui32NumOfSigmaPoints    (1,1) uint32 {mustBeGreaterThan(ui32NumOfSigmaPoints, 0)}
            % dPerturbScale           (1,1) double {mustBeReal, mustBeGreaterThan(dPerturbScale, 0.0)}
            % -------------------------------------------------------------------------------------------------------------
            %% OUTPUT
            % dSigmaPointsSet
            % -------------------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 17-08-2025    Pietro Califano     Renewed implementation from deprecated code.
            % 27-04-2026    Pietro Califano     Avoid resizing the covariance input in-place so the
            %                                   generator remains MATLAB Coder compatible.
            % -------------------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % [-]
            % -------------------------------------------------------------------------------------------------------------

            %% Function code
            coder.inline('always');
            dSqrtCovariance = zeros(size(dxCovariance));

            if istriu(dxCovariance)
                % Already factorized as upper triangular, transpose
                dSqrtCovariance(:,:) = transpose(dxCovariance);
                
            elseif not(istril(dxCovariance))
                % Full covariance matrix, factorize cholesky
                dSqrtCovariance(:,:) = chol(dxCovariance, 'lower');
            else
                dSqrtCovariance(:,:) = dxCovariance;
            end

            % Initialize output
            ui32StateSize = uint32(size(dxMean,1));            
            dSigmaPointsSet = zeros(ui32StateSize, ui32NumOfSigmaPoints);

            % Compute perturbations vectors
            dDeltaX = dPerturbScale .* dSqrtCovariance;

            % Assign first Sigma point (Mean state)
            dSigmaPointsSet(:, 1) = dxMean(1:ui32StateSize);

            % Perturb state vector to get Sigma Points
            for idState = 1:double(ui32StateSize)
                dSigmaPointsSet(:, 1 + idState) = dxMean(1:ui32StateSize) + dDeltaX(:, idState);
                dSigmaPointsSet(:, 1 + double(ui32StateSize) + idState) = dxMean(1:ui32StateSize) - dDeltaX(:, idState);
            end


        end
    
        function [dxMean, dxCovariance] = SumSigmaPointsToMoments(dSigmaPoints, ...
                                                                dMeanWeights, ...
                                                                dCovWeights)%#codegen
            arguments (Input)
                dSigmaPoints (:,:) double {ismatrix, mustBeReal}
                dMeanWeights (1,:) double {ismatrix, mustBeReal}
                dCovWeights  (1,:) double {ismatrix, mustBeReal}
            end
            arguments (Output)
                dxMean          (:,1) double {ismatrix, mustBeReal}
                dxCovariance    (:,:) double {ismatrix, mustBeReal}
            end
            %% PROTOTYPE
            % [dxMean, dxCovariance] = SumSigmaPointsToMoments(dSigmaPoints, dMeanWeights, dCovWeights)%#codegen
            % -------------------------------------------------------------------------------------------------------------
            %% DESCRIPTION
            % Function computing the weighted sample mean and covariance of a set of Sigma Points with generic
            % weights (including Cubature and UT variants).
            % -------------------------------------------------------------------------------------------------------------
            %% INPUT
            % dSigmaPoints (:,:) double {ismatrix, mustBeReal}
            % dMeanWeights (1,:) double {ismatrix, mustBeReal}
            % dCovWeights  (1,:) double {ismatrix, mustBeReal}
            % -------------------------------------------------------------------------------------------------------------
            %% OUTPUT
            % dxMean          (:,1) double {ismatrix, mustBeReal}
            % dxCovariance    (:,:) double {ismatrix, mustBeReal}
            % -------------------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 17-08-2025    Pietro Califano     Renewed implementation from deprecated code.
            % -------------------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % [-]
            % -------------------------------------------------------------------------------------------------------------


            %% Input handling
            if coder.target('MATLAB') || coder.target('MEX')
                % TODO: asserts validity of sizes
            end

            ui32StateSize = uint32(size(dSigmaPoints, 1));
            dxMean        = coder.nullcopy(zeros(ui32StateSize, 1));
            dxCovariance  = coder.nullcopy(zeros(ui32StateSize, ui32StateSize));

            %% Mean computation
            % Compute Weighted mean of Sigma Points
            dxMean(:) = sum(dMeanWeights .* dSigmaPoints, 2);

            %% Covariance computation
            dxDevFromMean = dSigmaPoints - dxMean;

            dxCovariance(:,:) = dCovWeights(1) .* (dxDevFromMean(:,1)) * transpose(dxDevFromMean(:,1)) + ...
                dCovWeights(2) .* (dxDevFromMean(:,2:end)) * transpose(dxDevFromMean(:,2:end));

            % Ensure symmetry by averaging
            dxCovariance(:,:) = 0.5 * (dxCovariance + dxCovariance');

        end
    
    
        %%% ------------------------------------------------------------------------------------------------- %%%
        %%% Monte Carlo methods
        function [dxMeanOut, dxCovarianceOut] = PropagateMonteCarloTransformPropagateDyn(dxMeanIn, ...
                                                                                        dxCovarianceIn, ...
                                                                                        strDynParams, ...
                                                                                        strFilterMutabConfig, ...
                                                                                        strFilterConstConfig, ...
                                                                                        dTimestampStart, ...
                                                                                        dTimestampEnd, ...
                                                                                        dIntegrTimestep, ...
                                                                                        ui32NumOfSamples, ...
                                                                                        bStartParpool)%#codegen
            arguments (Input)
                dxMeanIn          (:,1) double {isvector, mustBeReal}
                dxCovarianceIn    (:,:) double {ismatrix, mustBeReal}
                strDynParams         (1,1) {isstruct}
                strFilterMutabConfig (1,1) {isstruct}
                strFilterConstConfig (1,1) {isstruct}
                dTimestampStart   (1,1) double {isscalar, mustBeReal} = 0.0
                dTimestampEnd     (1,1) double {isscalar, mustBeReal} = 1.0
                dIntegrTimestep   (1,1) double {isscalar, mustBeReal} = 1.0
                ui32NumOfSamples  (1,1) uint32 {isscalar, mustBeReal} = 100
                bStartParpool     (1,1) logical {isscalar} = false
            end
            arguments (Output)
                dxMeanOut          (:,1) double {ismatrix, mustBeReal}
                dxCovarianceOut    (:,:) double {ismatrix, mustBeReal}
            end

            % Define function handle for PropagateDyn function
            fcnDynHandle = @(dxSamplePoint, dTimestart, dDeltaTime, dIntegrTimestep) PropagateDyn(dxSamplePoint, ...
                                                                                        dTimestart, ...
                                                                                        dDeltaTime, ...
                                                                                        dIntegrTimestep, ...
                                                                                        strDynParams, ...
                                                                                        strFilterMutabConfig, ...
                                                                                        strFilterConstConfig);

            % Call handle version
            [dxMeanOut, dxCovarianceOut] = CDensityFcnPropagator.PropagateMonteCarloTransformHandle(fcnDynHandle, ...
                                                                                               dxMeanIn, ...
                                                                                               dxCovarianceIn, ...
                                                                                               dTimestampStart, ...
                                                                                               dTimestampEnd, ...
                                                                                               ui32NumOfSamples, ...
                                                                                               bStartParpool, ...
                                                                                               dIntegrTimestep);%#codegen

        end

        function [dxMeanOut, dxCovarianceOut] = PropagateMonteCarloTransformHandle(fcnHandle, ...
                                                                                    dxMeanIn, ...
                                                                                    dxCovarianceIn, ...
                                                                                    dTimestampStart, ...
                                                                                    dTimestampEnd, ...
                                                                                    ui32NumOfSamples, ...
                                                                                    bStartParpool, ...
                                                                                    varargin)%#codegen
            %PROPAGATEMONTECARLOTRANSFORMHANDLE Monte Carlo propagation of state moments through a nonlinear map
            % -------------------------------------------------------------------------------------------------
            % DESCRIPTION
            %   Draws ui32NumOfSamples samples from a Gaussian N(dxMeanIn, dxCovarianceIn), propagates each
            %   sample through the provided function handle fcnHandle, and returns the empirical mean and
            %   covariance of the propagated set.
            %
            %   The function mirrors the design and calling style of the Sigma-Points counterpart, including:
            %   - Optional time-augmented signature for fcnHandle
            %   - Parallel propagation over samples (parfor-safe)
            %   - Arguments validation and codegen readiness
            %
            % IMPROVABLE SAMPLING OPTIONS (not implemented here; consider for accuracy/variance reduction)
            %   1) Antithetic Variates: for each standard normal z, also use -z before correlation injection.
            %      Preserves mean exactly and reduces variance of estimators.
            %   2) Latin Hypercube / Stratified Sampling (LHS): stratify each dimension in [0,1] (via inverse
            %      CDF to normal) to produce better space-filling than pure i.i.d. normals.
            %   3) Quasi-Monte Carlo (Sobol/Halton): generate low-discrepancy points u in [0,1]^n and map
            %      via the Gaussian inverse CDF; then correlate with Cholesky. Faster convergence in practice.
            %   4) Moment-Matching / Spherical-Radial MC: enforce zero sample mean and unit covariance on
            %      the standard normal set before correlation, then inject correlation to better match
            %      target moments with finite N (useful for small/medium sample sizes).
            %   5) Importance Sampling: when fcnHandle induces strong nonlinearity in specific regions, bias
            %      sampling toward such regions and reweight samples accordingly.
            %
            %   All the above can be adapted to respect your naming and codegen constraints (e.g., prebuilt
            %   Sobol sequences mapped through icdf('Normal',...); or deterministic antithetic pairing).
            %
            % CHANGELOG
            %   19-08-2025    Pietro Califano, GPT-5 Thinking     First implementation.
            %
            % DEPENDENCIES
            %   [-]
            %
            % FUTURE UPGRADES
            %   - Add antithetic/Stratified/Quasi-MC paths (compile-time switches friendly to codegen)
            %   - Add moment-matching normalization of standard normals before correlation injection
            %   - Optional return of full propagated cloud
            % -------------------------------------------------------------------------------------------------

            arguments (Input)
                fcnHandle         (1,1) {mustBeA(fcnHandle, "function_handle")}
                dxMeanIn          (:,1) double {isvector, mustBeReal}
                dxCovarianceIn    (:,:) double {ismatrix, mustBeReal}
                dTimestampStart   (1,1) double {isscalar, mustBeReal} = -1.0
                dTimestampEnd     (1,1) double {isscalar, mustBeReal, mustBeGreaterThan(dTimestampEnd, dTimestampStart)} = 0.0
                ui32NumOfSamples  (1,1) {mustBeInteger, mustBePositive} = uint32(512)
                bStartParpool     (1,1) logical {isscalar} = false
            end
            arguments (Input, Repeating)
                varargin
            end
            arguments (Output)
                dxMeanOut          (:,1) double {ismatrix, mustBeReal}
                dxCovarianceOut    (:,:) double {ismatrix, mustBeReal}
            end

            % Initialize outputs
            dxMeanOut       = zeros(size(dxMeanIn));
            dxCovarianceOut = zeros(size(dxCovarianceIn));
            ui32DomainSize  = uint32(size(dxMeanIn,1));

            % Optional integration timestep 
            dIntegrTimestep = 0.0; % Default value to allow parfor usage
            if ~isempty(varargin)
                dIntegrTimestep = varargin{1};
            end

            % -------------------------------------------------------------------------
            % Generate Monte Carlo samples ~ N(dxMeanIn, dxCovarianceIn)
            % -------------------------------------------------------------------------
            % Symmetrize covariance and build a numerically robust Cholesky factor.
            dSymCov_tmp = 0.5 * ( dxCovarianceIn + transpose(dxCovarianceIn) );
            dTraceAbs   = sum(abs(diag(dSymCov_tmp)));
            dJitter     = max(1e-12, 1e-12 * dTraceAbs);

            % Try Cholesky; if it fails (semi-definite), fall back to eig-clip.
            % try
            dLchol = chol(dSymCov_tmp + dJitter * eye(double(ui32DomainSize)), 'lower');
            % catch %#ok<CTCH>
            %     % Eigenvalue clipping to nearest SPD approximation
            %     [dVEig_tmp, dDEig_tmp] = eig(dSymCov_tmp);
            %     dLam = max(diag(dDEig_tmp), 0.0);
            %     dDEig_clip = diag(dLam + dJitter);
            %     dSPD_tmp   = dVEig_tmp * dDEig_clip * dVEig_tmp.';
            %     dSPD_tmp   = 0.5 * (dSPD_tmp + dSPD_tmp.'); % enforce symmetry
            %     dLchol     = chol(dSPD_tmp, 'lower');
            % end

            % Standard normal draws (population size = ui32NumOfSamples)
            % Note: For reproducibility in MATLAB, set rng(...) BEFORE calling this function.
            dStdNormalSamples = randn(double(ui32DomainSize), double(ui32NumOfSamples));

            % Inject correlation and mean shift
            dxSamplesSet = dLchol * dStdNormalSamples;
            dxSamplesSet = dxSamplesSet + dxMeanIn; 

            % -------------------------------------------------------------------------
            % Propagate samples through the user function handle
            % -------------------------------------------------------------------------
            if bStartParpool
                objCurrentPool = gcp('nocreate');
                if isempty(objCurrentPool)
                    parpool('Threads');
                end
            end

            tic
            parfor idSample = 1:double(ui32NumOfSamples)
                if dTimestampStart > 0 - eps
                    % Time-augmented signature
                    dxSamplesSet(:, idSample) = fcnHandle(dxSamplesSet(:, idSample), ...
                                                        dTimestampStart, ...
                                                        dTimestampEnd - dTimestampStart, ...
                                                        dIntegrTimestep);
                else
                    % State-only signature
                    dxSamplesSet(:, idSample) = fcnHandle(dxSamplesSet(:, idSample));
                end
            end

            dElapsedTimeMC = toc;
            fprintf('\nMC propagation of %d samples completed in %2.4g [s].\n', ui32NumOfSamples, dElapsedTimeMC);

            % -------------------------------------------------------------------------
            % Empirical moments (population versions: divide by N)
            % -------------------------------------------------------------------------
            [dxMeanOut(:), dxCovarianceOut(:,:)] = CDensityFcnPropagator.SumSamplesToMoments(dxSamplesSet);

        end

        function [dMeanOut, dCovOut, dCentredDataMatrix] = SumSamplesToMoments(dxSamplesSet)
            arguments
                dxSamplesSet (:,:) double {mustBeReal, ismatrix}
            end
            %SUMSAMPLESTOMOMENTS Compute empirical mean and covariance from sample matrix

            dMeanOut           = mean(dxSamplesSet, 2);
            dCentredDataMatrix = dxSamplesSet - dMeanOut;

            ui32SampleSize = size(dxSamplesSet, 2);

            % Population covariance (consistent with sigma-point weighting style)
            dCovOut = (dCentredDataMatrix * transpose(dCentredDataMatrix)) / ui32SampleSize;
        end

    end

end
