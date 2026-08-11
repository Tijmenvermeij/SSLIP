function writeSSLIP_Report(filepath, reportData)
% writeSSLIP_Report Writes a CSV report of the slip activities and bundling logic.
%
% reportData is expected to be a struct array containing:
%   BundleID (double or string)
%   SystemID (double)
%   Activity (double)
%   SchmidFactor (double)
%   H11, H12, H21, H22 (double)
%   IsRepresentative (logical)
%   n_x, n_y, n_z, b_x, b_y, b_z (double)

if nargin < 2
    error('writeSSLIP_Report requires a filepath and reportData struct array.');
end

if isempty(reportData)
    warning('No report data to write.');
    return;
end

% Convert struct array to table
try
    T = struct2table(reportData);
catch
    warning('Failed to convert reportData to table. Make sure it is a uniform struct array.');
    return;
end

% Write table to CSV
try
    writetable(T, filepath);
    fprintf('Successfully exported slip activity report to:\n  %s\n', filepath);
catch ME
    warning('Failed to export slip activity report to %s.\nError: %s', filepath, ME.message);
end

end
