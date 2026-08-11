function plotSSLIP_Context(opt)
% plotSSLIP_Context Plots the background grain context if opt.grains is provided.
% This should be called immediately before plotting the main EBSD data field
% on a given axis, so that the context lies behind the data.

if ~isfield(opt, 'plotGrainContext')
    opt.plotGrainContext = false;
end

if opt.plotGrainContext && isfield(opt, 'grains') && ~isempty(opt.grains)
    % Default styling for the context grains
    if ~isfield(opt, 'grainFaceColor'), opt.grainFaceColor = [0.8 0.8 0.8]; end
    if ~isfield(opt, 'grainLineColor'), opt.grainLineColor = 'w'; end
    if ~isfield(opt, 'grainLineWidth'), opt.grainLineWidth = 1; end

    % Plot background
    plot(opt.grains, 'facecolor', opt.grainFaceColor); 
    hold on;
    plot(opt.grains.boundary, 'linecolor', opt.grainLineColor, 'linewidth', opt.grainLineWidth);
    hold on;
end
end
