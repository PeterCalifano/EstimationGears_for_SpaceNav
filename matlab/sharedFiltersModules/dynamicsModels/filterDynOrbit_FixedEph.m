function dDrvDt = filterDynOrbit_FixedEph(dStateTimetag, ...
                                dxState, ...
                                strDynParams, ...
                                strFilterConstConfig) %#codegen
arguments
    dStateTimetag         (1, 1) double
    dxState               (:, 1) double
    strDynParams          {isstruct}
    strFilterConstConfig  {isstruct}
end
%% PROTOTYPE
% dDrvDt = filterDynOrbit_FixedEph(dStateTimetag, ...
%                                 dxState, ...
%                                 strDynParams, ...
%                                 strFilterConstConfig) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCurrentTime
% dxState
% strDynParams
% strStatesIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dxdt
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025    Pietro Califano     Implement version taking from legact filterDynOrbit and
%                                   for compatibility with evalRHS_InertialDynOrbit
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalRHS_InertialDynOrbit()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% TODO make more general purpose
% TODO convert to use static sized arrays
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Get configuration flags to configure dynamics 
% TODO if any
dDrvDt = zeros(6, 1);

% Get ephemerides from dynamic params struct
dBodyEphemerides = strDynParams.dBodyEphemerides;
dDCMmainAtt_INfromTF = zeros(3,3);

if isfield(strDynParams, "dDCMmainAtt_INfromTF")
    dDCMmainAtt_INfromTF(:,:) = strDynParams.dDCMmainAtt_INfromTF;
end

% Compute SRP coefficient 
dBiasCoeffSRP = 0.0;
if isfield(strFilterConstConfig.strStatesIdx, "ui8CoeffSRPidx")
    dBiasCoeffSRP(:) = dxState( strFilterConstConfig.strStatesIdx.ui8CoeffSRPidx);
end

dCoeffSRP = (strDynParams.strSRPdata.dP_SRP * strDynParams.strSCdata.dReflCoeff * ...
             strDynParams.strSCdata.dA_SRP)/strDynParams.strSCdata.dSCmass; % Move to compute outside, since this

dCoeffSRP = dCoeffSRP + dBiasCoeffSRP;

% Get residual acceleration if any
dResidualAccel = zeros(3,1);

if isfield(strFilterConstConfig.strStatesIdx, "ui8ResidualAccelIdx")
    dResidualAccel(:) = dxState( strFilterConstConfig.strStatesIdx.ui8ResidualAccelIdx);
end

%% Evaluate RHS
dDrvDt(strFilterConstConfig.ui8posVelIdx) = evalRHS_InertialDynOrbit(dxState, ...
                                                                     dDCMmainAtt_INfromTF, ...
                                                                     strDynParams.strMainData.dGM, ...
                                                                     strDynParams.strMainData.dRefRadius, ...
                                                                     dCoeffSRP, ...
                                                                     d3rdBodiesGM, ...
                                                                     dBodyEphemerides, ...
                                                                     strDynParams.strMainData.dSHcoeff, ...
                                                                     strDynParams.strMainData.ui16MaxSHdegree, ...
                                                                     ui16StatesIdx, ...
                                                                     dResidualAccel);


end
