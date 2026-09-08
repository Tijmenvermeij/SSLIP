function rotationFigure = plotSSLIP_Rotation(ebsdID, opt)
% Plot stored rotation in degrees and return the figure handle.
% rotationIDcor remains in radians. Missing rotation data creates no figure.
% Adapted from Philipp (PhilKro), PR #2, original commit d3ee1a5:
% https://github.com/Tijmenvermeij/SSLIP/pull/2
% Retains the current examples' colors, limits, labels, and saved filename.

if nargin < 2, opt = struct; end
if ~isfield(opt,'saveFig'), opt.saveFig = 0; end
rotationFigure = [];
if ~isfield(ebsdID.prop,'rotationIDcor'), return; end

figure;
rotationDegrees = ebsdID.prop.rotationIDcor / degree;
plot(ebsdID,rotationDegrees,'micronbar','off');
title('Inferred rotation correction');
mtexColorMap(blue2redColorMap);
rotationLimit = max(abs(rotationDegrees),[],'omitnan');
if isfinite(rotationLimit) && rotationLimit > 0
    clim([-rotationLimit rotationLimit]);
end
mtexColorbar('title','Rotation [deg]');
rotationFigure = gcf;
if opt.saveFig
    saveFigure([opt.plotname '_rotation.png']);
end
end
