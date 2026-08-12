% Reset and add runtime paths for this repository only.
charRepoRoot = fileparts(mfilename('fullpath'));
charCallDir = cd(charRepoRoot);

restoredefaultpath;

addpath(charRepoRoot)
addpath(genpath(fullfile(charRepoRoot, 'matlab')))
addpath(genpath(fullfile(charRepoRoot, 'simulink')))
addpath(genpath(fullfile(charRepoRoot, 'lib')))

cd(charCallDir);
