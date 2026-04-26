clear
close all
clc

%% CODEGEN SCRIPT for GivensRotSRIF() function
% Created by PeterC 25-04-2024
% 24-04-2026    Pietro Califano     Replaced the missing makeCodegen helper with a direct MATLAB
%                                   Coder invocation.
targetFcnName = 'GivensRotSRIF';
targetMexName = 'GivensRotSRIF_MEX';

nMaxStates = 100;

dxPrior         = coder.typeof(0, [nMaxStates, 1], [1,0]);
dSRInfoMatPrior = coder.typeof(0, [nMaxStates, nMaxStates], [1,1]);
dYobs           = coder.typeof(0, [Inf, 1], [1,0]);
dHobsMatrix     = coder.typeof(0, [Inf, nMaxStates], [1,1]);
bNPRIOR_INFO    = coder.typeof(false, [1,1]);
bRUN_WHITENING  = coder.typeof(false, [1,1]);
dMeasCovSR      = coder.typeof(0, [Inf, Inf], [1,1]);

args_cell = {dxPrior, ...
             dSRInfoMatPrior, ...
             dYobs, ...
             dHobsMatrix, ...
             bNPRIOR_INFO, ...
             bRUN_WHITENING, ...
             dMeasCovSR};

cfg = coder.config('mex');
cfg.GenerateReport = false;
cfg.RowMajor = false;

codegen('-config', cfg, targetFcnName, '-args', args_cell, '-o', targetMexName);
