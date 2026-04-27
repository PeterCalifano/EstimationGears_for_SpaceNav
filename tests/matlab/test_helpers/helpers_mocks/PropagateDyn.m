function [dxOut, dTimeOut, strDynParams] = PropagateDyn(dxSigmaPoint, ...
                                                        dTimestart, ...
                                                        dDeltaTime, ...
                                                        dIntegrStep, ...
                                                        strDynParams, ...
                                                        strFilterMutabConfig, ...
                                                        strFilterConstConfig)
arguments
    dxSigmaPoint (:,1) double {mustBeReal}
    dTimestart   (1,1) double {mustBeReal}
    dDeltaTime   (1,1) double {mustBeReal}
    dIntegrStep  (1,1) double {mustBeReal}
    strDynParams (1,1) struct
    strFilterMutabConfig (1,1) struct
    strFilterConstConfig (1,1) struct
end

% Identity-like propagation with a deterministic drift; keeps shapes consistent for smoke testing.
if isempty(dIntegrStep) || isempty(strFilterMutabConfig) || isempty(strFilterConstConfig)
end

dxOut = dxSigmaPoint + dDeltaTime;
dTimeOut = dTimestart + dDeltaTime;
end
