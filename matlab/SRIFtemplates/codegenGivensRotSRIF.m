clear
close all
clc

%% CODEGEN SCRIPT for GivensRotSRIF() function
% Created by PeterC 25-04-2024
targetFcnName = "GivensRotSRIF.m";


% Input types and sizes definition

% [o_dxPost, o_dSRInfoMatPost, o_dInfoVecPost, o_dErrorVec, o_dSqrtPxPost, o_dJcost] =
% GivensRotSRIF(i_dxPrior, ...
%     i_dSRInfoMatPrior, ...
%     i_dYobs, ...
%     i_dHobsMatrix, ...
%     i_bNO_PRIOR_INFO,...
%     i_bRUN_WHITENING,...
%     i_dMeasCovSR)

nMaxStates = 100;

i_dxPrior         = coder.typeof(0, [nMaxStates, 1], [1,0]); 
i_dSRInfoMatPrior = coder.typeof(0, [nMaxStates, nMaxStates], [1,1]); 
i_dYobs           = coder.typeof(0, [Inf, 1], [1,0]); 
i_dHobsMatrix     = coder.typeof(0, [Inf, nMaxStates], [1,1]); 
i_bNO_PRIOR_INFO  = coder.typeof(false, [1,1]); 
i_bRUN_WHITENING  = coder.typeof(false, [1,1]); 
i_dMeasCovSR      = coder.typeof(0, [Inf, Inf], [1,1]);  


args_cell{1} = i_dxPrior        ;
args_cell{2} = i_dSRInfoMatPrior;
args_cell{3} = i_dYobs          ;
args_cell{4} = i_dHobsMatrix    ;
args_cell{5} = i_bNO_PRIOR_INFO ;
args_cell{6} = i_bRUN_WHITENING ;
args_cell{7} = i_dMeasCovSR     ;


% Call automatic maker for codegen
makeCodegen(targetFcnName, args_cell);
