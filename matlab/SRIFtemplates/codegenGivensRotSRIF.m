clear
close all
clc

%% CODEGEN SCRIPT for GivensRotSRIF() function
% Created by PeterC 25-04-2024
targetFcnName = "GivensRotSRIF.m";

% Input types and sizes definition

% [dxPost, dSRInfoMatPost, dInfoVecPost, dErrorVec, dSqrtPxPost, dJcost] = GivensRotSRIF(dxPrior, ...
%                                                                                        dSRInfoMatPrior, ...
%                                                                                        dYobs, ...
%                                                                                        dHobsMatrix, ...
%                                                                                        bNPRIOR_INFO,...
%                                                                                        bRUN_WHITENING,...
%                                                                                        dMeasCovSR)

nMaxStates = 100;

dxPrior         = coder.typeof(0, [nMaxStates, 1], [1,0]); 
dSRInfoMatPrior = coder.typeof(0, [nMaxStates, nMaxStates], [1,1]); 
dYobs           = coder.typeof(0, [Inf, 1], [1,0]); 
dHobsMatrix     = coder.typeof(0, [Inf, nMaxStates], [1,1]); 
bNPRIOR_INFO  = coder.typeof(false, [1,1]); 
bRUN_WHITENING  = coder.typeof(false, [1,1]); 
dMeasCovSR      = coder.typeof(0, [Inf, Inf], [1,1]);  


args_cell{1} = dxPrior        ;
args_cell{2} = dSRInfoMatPrior;
args_cell{3} = dYobs          ;
args_cell{4} = dHobsMatrix    ;
args_cell{5} = bNPRIOR_INFO ;
args_cell{6} = bRUN_WHITENING ;
args_cell{7} = dMeasCovSR     ;


% Call automatic maker for codegen
makeCodegen(targetFcnName, args_cell);
