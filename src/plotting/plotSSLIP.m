function plotSSLIP(ebsdID, sSLocal, opt)
% plotSSLIP High-level orchestrator for plotting SSLIP results.
% It now delegates plotting to specific sub-functions based on 'opt' flags.
% Note: To plot context boundaries (e.g. grain boundaries) over the maps,
% ensure opt.plotBoundariesOverDisps = 1 and provide opt.grains. This will
% automatically trigger plotSSLIP_Context.m inside the sub-functions.

if nargin < 3
    opt = struct;
end

% Set default orchestrator flags
if ~isfield(opt, 'plotSlipActivities'), opt.plotSlipActivities = 1; end
if ~isfield(opt, 'plotResidual'), opt.plotResidual = 1; end
if ~isfield(opt, 'plotRotation'), opt.plotRotation = 1; end
if ~isfield(opt, 'plotDominant'), opt.plotDominant = 0; end
if ~isfield(opt, 'plotDefGrad'), opt.plotDefGrad = 1; end

% 1. Deformation Fields
if opt.plotDefGrad
    try
        plotSSLIP_DeformationFields(ebsdID, opt);
    catch ME
        warning(ME.identifier, 'Failed to plot Deformation Fields: %s:\n%s\n', ME.message, getReport(ME));
    end
end

% 2. Slip Activities
if opt.plotSlipActivities && ~opt.plotDominant
    try
        plotSSLIP_SlipActivities(ebsdID, sSLocal, opt);
    catch ME
        warning(ME.identifier, 'Failed to plot Slip Activities: %s:\n%s\n', ME.message, getReport(ME));
    end
end

% 3. Dominant Activities (Top N + Leftover)
if opt.plotDominant
    try
        plotSSLIP_DominantActivities(ebsdID, sSLocal, opt);
    catch ME
        warning(ME.identifier, 'Failed to plot Dominant Activities: %s:\n%s\n', ME.message, getReport(ME));
    end
end

% 4. Residual Error
if opt.plotResidual
    try
        plotSSLIP_Residual(ebsdID, opt);
    catch ME
        warning(ME.identifier, 'Failed to plot Residual: %s:\n%s\n', ME.message, getReport(ME));
    end
end

% 5. Rotation Field
if opt.plotRotation && isfield(ebsdID.prop, 'rotationIDcor')
    try
        plotSSLIP_Rotation(ebsdID, opt);
    catch ME
        warning(ME.identifier, 'Failed to plot Rotation: %s:\n%s\n', ME.message, getReport(ME));
    end
end

end
