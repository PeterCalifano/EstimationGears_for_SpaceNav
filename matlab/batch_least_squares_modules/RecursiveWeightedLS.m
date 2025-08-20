function [xhat, P, xhist, Phist] = RecursiveWLS(y, H, W, xhat0) %# codegen
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21-06-202     Pietro Califano     Code structure coded. To be completed.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades

% -------------------------------------------------------------------------------------------------------------
%% Function code


%% Input handling

% Get size of H
[nrowH, ncolH, n3rdH] = size(H);
% Get size of W
[nrowW, ncolW, n3rdW] = size(W);
% Get size of y
[nrowy, ncoly, n3rdy] = size(y);

%% Data array checks and formatting
assert(n3rdy == 1, 'Error: data array y with > 2 dimensions NOT supported.')

% Format y to column vector
if ~iscolumn(y)
    % Measurement vector dimension assumed as number of rows
    measDim = Nrowy;
    y = reshape(y, Nrowy*Ncoly, []);
end

%% Model H and Weight W checks and formatting


assert(nrowH == nrowW, 'Error: size mismatch between H and W.')


% If W has 3rd dimension --> change over time

% If x0 guess not provided, use first N meas in y to solve LS and get x0, P0
if nargin <= 3
    % x0 not provide: apply WLS to get first guess
    idsFirst = 5;
    [xhat0, P] = WLS(y(:, 1:idsFirst), H, W);

else
    % For Minimum Variance Estimator, use first R as Covariance
    P = W(:, :, 1);
end

% Initialize estimate xhat 
xkprev = xhat0;
Pkprev = P;

if nargout > 2
    idsFirst = 2;
    xhist = zeros(length(xkprev), Nsamples);
    xhist(1) = xhat0;

    if nargout > 3
        Phist = zeros(length(xkprev), lenght(xkprev), Nsamples);
        Phist(:, :, 1) = W(:, :, 1);
    end
end



% Main iteration loop
for ids = idsFirst:Nsamples

    % Apply RWLSkernel for each step
    [xhat, P] = RWLSkernel(y(:, ids), H, xkprev, Pkprev, Wk);

    % Save current variable for next iteration
    Pkprev = P;
    xkprev = xhat;

    if nargout > 2
        xhist(:, ids) = xhat;
        if nargout > 3
            Phist(:, :, ids) = P;
        end
    end

end


    function [xkhat, Pk] = RWLSkernel(yk, Hk, xkprev, Pkprev, Wk)


        % Model: yk = Hk*xk + ek;
        % Compute optimal gain
        Kk = Pkprev * Hk' / (Hk' * Pkprev * Hk + Wk);

        % Update estimate using Innovation and Kk gain
        innk = yk - Hk * xkprev;
        xkhat = xkprev + Kk*innk;

        % Update covariance (Joseph Formula)
        UpdateTerm = eye(length(xkhat)) - Kk * Hk;
        Pk = UpdateTerm * Pkprev * UpdateTerm' + Kk * Rk * Kk';

    end


end


