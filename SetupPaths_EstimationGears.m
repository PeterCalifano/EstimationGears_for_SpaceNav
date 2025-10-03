% Reset and add path for this repository only
charCallDir = cd(fileparts(mfilename('fullpath')));

restoredefaultpath;

addpath(genpath('matlab'))
addpath(genpath('simulink'))
addpath(genpath('lib'))
addpath(genpath('tests'))
addpath(genpath('.'))

cd(charCallDir);
