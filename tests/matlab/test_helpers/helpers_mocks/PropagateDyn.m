%%% Local helper with the exact signature expected by PropagateSigmaPointTransformPropagateDyn()
function dxOut = PropagateDyn(dxSigmaPoint, ...
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
% Identity-like propagation with tiny drift to avoid degenerate paths
% (keeps shapes consistent for smoke testing)

dxOut = dxSigmaPoint + 1.0*dDeltaTime;
end
