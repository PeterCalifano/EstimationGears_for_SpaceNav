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

    methods (Static, Access = public)

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
                dxMean                  (:,1) double {isvector, mustBeReal}
                dxCovariance            (:,:) double {ismatrix, mustBeReal}
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
            % dxMean                  (:,1) double {isvector, mustBeReal}
            % dxCovariance            (:,:) double {ismatrix, mustBeReal}
            % ui32NumOfSigmaPoints    (1,1) uint32 {mustBeGreaterThan(ui32NumOfSigmaPoints, 0)}
            % dPerturbScale           (1,1) double {mustBeReal, mustBeGreaterThan(dPerturbScale, 0.0)}
            % -------------------------------------------------------------------------------------------------------------
            %% OUTPUT
            % dSigmaPointsSet
            % -------------------------------------------------------------------------------------------------------------
            %% CHANGELOG
            % 17-08-2025    Pietro Califano     Renewed implementation from deprecated code.
            % -------------------------------------------------------------------------------------------------------------
            %% DEPENDENCIES
            % [-]
            % -------------------------------------------------------------------------------------------------------------

            if istriu(dxCovariance)
                % Already factorized as upper triangular, transpose
                dxCovariance = transpose(dxCovariance);
            elseif not(istril(dxCovariance))
                % Full covariance matrix, factorize cholesky
                dxCovariance = chol(dxCovariance, 'lower');
            end

            % Initialize output
            ui32StateSize = uint32(size(dxMean,1));            
            dSigmaPointsSet = zeros(ui32StateSize, ui32NumOfSigmaPoints);

            % Compute perturbations vectors
            dDeltaX = dPerturbScale .* dxCovariance;

            % Assign first Sigma point (Mean state)
            dSigmaPointsSet(:, 1) = dxMean(1:ui32StateSize);

            % Perturb state vector to get Sigma Points
            % "Right-side" Sigma Points
            dSigmaPointsSet(:, 2:ui32StateSize+1) = dxMean(1:ui32StateSize) + dDeltaX;

            % "Left-side" Sigma Points
            dSigmaPointsSet(:, ui32StateSize + 2:ui32NumOfSigmaPoints) = dxMean(1:ui32StateSize) - dDeltaX;


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
    end

end

