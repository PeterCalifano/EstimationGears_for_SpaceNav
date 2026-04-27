function [dUFactor, dDFactor] = UDdecompositionFromUpperSqrtCov(dSqrtCovUpper) %#codegen
%% PROTOTYPE
% [dUFactor, dDFactor] = UDdecompositionFromUpperSqrtCov(dSqrtCovUpper)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing the U-D decomposition of the covariance represented by an upper-triangular
% square-root factor S, without explicitly forming P = S' * S.
%
% The output follows the repository convention:
%
%   P = U * D * U'
%
% where U is unit upper triangular and D is diagonal.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dSqrtCovUpper:       [N, N]       Upper-triangular covariance square-root factor, P = S' * S
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dUFactor:            [N, N]       Unit upper triangular U factor of covariance P
% dDFactor:            [N, N]       Diagonal D factor of covariance P
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-04-2026    Pietro Califano     Add direct U-D decomposition from an upper square-root covariance factor.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code

ui16StateSize = uint16(size(dSqrtCovUpper, 1));

if coder.target('MATLAB') || coder.target('MEX')
    assert(ui16StateSize == size(dSqrtCovUpper, 2), ...
        'ERROR: upper square-root covariance factor must be square.');
end

dDFactor = zeros(double(ui16StateSize), double(ui16StateSize));
dUFactor = zeros(double(ui16StateSize), double(ui16StateSize));
dDDiagVec = zeros(double(ui16StateSize), 1);

% Match the current SR-UKF wrapper regularization policy for pathological
% non-positive covariance diagonals without materializing the full covariance.
dDiagJitter = 0.0;
if coder.target('MATLAB') || coder.target('MEX')
    dMinCovDiag = 0.0;

    for ui16IdDiag = uint16(1):ui16StateSize
        dCovDiag = 0.0;
        for ui16IdDot = uint16(1):ui16IdDiag
            dCovDiag = dCovDiag + dSqrtCovUpper(ui16IdDot, ui16IdDiag) * dSqrtCovUpper(ui16IdDot, ui16IdDiag);
        end

        if ui16IdDiag == 1 || dCovDiag < dMinCovDiag
            dMinCovDiag = dCovDiag;
        end
    end

    if dMinCovDiag <= 0.0
        dDiagJitter = abs(dMinCovDiag) + 1.0e-12;
    end
end

for ui16IdCol = ui16StateSize:-1:uint16(1)
    for ui16IdRow = ui16IdCol:-1:uint16(1)
        dCovEntry = 0.0;

        % Because S is upper triangular and ui16IdRow <= ui16IdCol:
        % P(row,col) = S(1:row,row)' * S(1:row,col).
        for ui16IdDot = uint16(1):ui16IdRow
            dCovEntry = dCovEntry + dSqrtCovUpper(ui16IdDot, ui16IdRow) * dSqrtCovUpper(ui16IdDot, ui16IdCol);
        end

        if ui16IdRow == ui16IdCol
            dCovEntry = dCovEntry + dDiagJitter;
        end

        dTmp = dCovEntry;
        for ui16IdAuxCol = ui16IdCol + uint16(1):ui16StateSize
            dTmp = dTmp - dUFactor(ui16IdRow, ui16IdAuxCol) * ...
                dDDiagVec(ui16IdAuxCol) * dUFactor(ui16IdCol, ui16IdAuxCol);
        end

        if ui16IdRow == ui16IdCol
            dDDiagVec(ui16IdCol) = dTmp;
            dUFactor(ui16IdCol, ui16IdCol) = 1.0;
        else
            dUFactor(ui16IdRow, ui16IdCol) = dTmp / dDDiagVec(ui16IdCol);
        end
    end
end

for ui16IdDiag = uint16(1):ui16StateSize
    dDFactor(ui16IdDiag, ui16IdDiag) = dDDiagVec(ui16IdDiag);
end

end
