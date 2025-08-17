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
                                                                                        dTimestampEnd)%#codegen
            arguments (Input)
                dxMeanIn          (:,1) double {isvector, mustBeReal}
                dxCovarianceIn    (:,:) double {ismatrix, mustBeReal}
                strDynParams         (1,1) {isstruct}
                strFilterMutabConfig (1,1) {isstruct}
                strFilterConstConfig (1,1) {isstruct}
                dTimestampStart   (1,1) double {isscalar, mustBeReal} = 0.0
                dTimestampEnd     (1,1) double {isscalar, mustBeReal} = 1.0
                % TODO add params
            end
            arguments (Output)
                dxMeanOut          (:,1) double {ismatrix, mustBeReal}
                dxCovarianceOut    (:,:) double {ismatrix, mustBeReal}
            end

            % Define function handle for PropagateDyn function
            fcnDynHandle = @(dxSigmaPoint, dTimestart, dDeltaTime, dIntegrStep) PropagateDyn(dxSigmaPoint, ...
                                                                                        dTimestart, ...
                                                                                        dDeltaTime, ...
                                                                                        dIntegrStep, ...
                                                                                        strDynParams, ...
                                                                                        strFilterMutabConfig, ...
                                                                                        strFilterConstConfig);

            % Call handle version
            [dxMeanOut, dxCovarianceOut] = CDensityFcnPropagator.PropagateSigmaPointTransformHandle(fcnDynHandle, ...
                                                                                               dxMeanIn, ...
                                                                                               dxCovarianceIn, ...
                                                                                               dTimestampStart, ...
                                                                                               dTimestampEnd);%#codegen

        end

        function [dxMeanOut, dxCovarianceOut] = PropagateSigmaPointTransformHandle(fcnHandle, ...
                                                                                   dxMeanIn, ...
                                                                                   dxCovarianceIn, ...
                                                                                   dTimestampStart, ...
                                                                                   dTimestampEnd)%#codegen
            arguments (Input)
                fcnHandle         (1,1) {mustBeA(fcnHandle, "function_handle")}
                dxMeanIn          (:,1) double {isvector, mustBeReal}
                dxCovarianceIn    (:,:) double {ismatrix, mustBeReal}
                dTimestampStart   (1,1) double {isscalar, mustBeReal} = 0.0
                dTimestampEnd     (1,1) double {isscalar, mustBeReal} = 1.0
                % TODO add other settings
            end
            arguments (Output)
                dxMeanOut          (:,1) double {ismatrix, mustBeReal}
                dxCovarianceOut    (:,:) double {ismatrix, mustBeReal}
            end

            % Initialize output variables
            dxMeanOut       = zeros(size(dxMeanIn));
            dxCovarianceOut = zeros(size(dxCovarianceIn));

            ui32DomainSize = uint32(size(dxMeanIn,1));

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

                dSigmaPointsSet(:, idCsi) = fcnHandle(dSigmaPointsSet(:, idCsi), ...
                                                      dTimestampStart, ...
                                                      dTimestampEnd - dTimestampStart, ...
                                                      dIntegrTimestep);

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
            ui32NumSigmaPoints    = 2*ui32DomainSize + 1;
            dWeightsMean = zeros(ui32NumSigmaPoints, 1);
            dWeightsCov = zeros(ui32NumSigmaPoints, 1);

            switch enumVariantName
                case "scaled"
                    % Scaled version Unscented Transform weights
                    % Reference: TODO


                    dLambda = dAlphaCoeff^2*(ui32NumSigmaPoints + k) ;
                    dPerturbScale = sqrt(dLambda); % Square root covariance perturbation scale factor

                    dWeightsMean(1) = 1 - ui32NumSigmaPoints / (dAlphaCoeff^2 * (ui32NumSigmaPoints + dKappaCoeff));
                    dWeightsCov(1) = (2 - dAlphaCoeff^2 + beta) - ui32NumSigmaPoints / (dAlphaCoeff^2 * (ui32NumSigmaPoints + dKappaCoeff));
                    
                    dWeightsCov(2:end) = ( 1/(2* dAlphaCoeff^2 *(ui32NumSigmaPoints + dKappaCoeff)) ) .* ones(1, 2*ui32NumSigmaPoints); % * ones(1, 2*ui32Ncsi);
                    dWeightsMean(2:end) = dWeightsCov(2:end);


                case "original"
                    % Original Unscented Transform weights
                    % Reference: TODO

                    dKappaCoeff      = 3 - ui32DomainSize;
                    dLambda = dAlphaCoeff^2 * (ui32DomainSize + dKappaCoeff) - ui32DomainSize;
                    dPerturbScale = sqrt(ui32DomainSize + dLambda); % Square root covariance perturbation scale factor

                    dWeightsMean(1) = dLambda./(ui32DomainSize + dLambda);
                    dWeightsCov(1) = dWeightsMean(1) + (1.0 - dAlphaCoeff^2 + dBetaCoeff);

                    dWeightsCov(2:end) = 1.0/(2.0* (ui32DomainSize + dLambda));
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
            elseif not(istril)
                % Full covariance matrix, factorize cholesky
                dxCovariance = chol(dxCovariance, 'lower');
            end

            % Initialize output
            ui16StateSize = uint16(size(dxMean,1));            
            dSigmaPointsSet = zeros(ui16StateSize, ui32NumOfSigmaPoints);

            % Compute perturbations vectors
            dDeltaX = dPerturbScale .* dxCovariance;

            % Assign first Sigma point (Mean state)
            dSigmaPointsSet(:, 1) = dxMean(1:ui16StateSize);

            % Perturb state vector to get Sigma Points
            % "Right-side" Sigma Points
            dSigmaPointsSet(:, 2:ui16StateSize+1) = dxMean(1:ui16StateSize) + dDeltaX;

            % "Left-side" Sigma Points
            dSigmaPointsSet(:, ui16StateSize + 2:ui32NumOfSigmaPoints) = dxMean(1:ui16StateSize) - dDeltaX;


        end
    end

end

