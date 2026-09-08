function initSSLIP()
% Initialize SSLIP after starting the desired MTEX version.
% Add the repository root to the MATLAB path, then call initSSLIP.
% This does not change the working folder or save the MATLAB path.
%
% Adapted from Philipp (PhilKro), PR #2, commit d3ee1a5:
% https://github.com/Tijmenvermeij/SSLIP/pull/2

baseDir = fileparts(mfilename('fullpath'));
% Add only the library folders; examples add their own utilities.
addpath(fullfile(baseDir,'src'), ...
    fullfile(baseDir,'src','plotting'),fullfile(baseDir,'src','utils'));
fprintf('SSLIP library initialized from: %s\n',baseDir);
end
