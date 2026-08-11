function initSSLIP()
% initSSLIP Initializes the SSLIP library by adding all source folders to the MATLAB path.
%
% Usage:
%   Call this function at the top of your analysis scripts instead of manually adding paths.
%   Example: initSSLIP();

    % Get the directory where this initialization script is located
    baseDir = fileparts(mfilename('fullpath'));
    
    % Add the root 'src' directory and all of its subdirectories recursively
    srcDir = fullfile(baseDir, 'src');
    addpath(genpath(srcDir));
    
    fprintf('SSLIP Library initialized successfully from: %s\n', baseDir);
end
