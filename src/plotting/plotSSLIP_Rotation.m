function plotSSLIP_Rotation(ebsdID, opt)
% plotSSLIP_Rotation Plots the solved pseudo-slip system activities for rigid body rotation.
% NOTE: Supports context boundaries via plotSSLIP_Context (requires opt.plotGrainContext=1 and opt.grains)

if nargin < 2
    opt = struct;
end

if ~isfield(ebsdID.prop, 'rotationIDcor')
    % No rotation field exists
    return;
end

if ~isfield(opt, 'clim_percentile'), opt.clim_percentile = 99; end
if ~isfield(opt, 'saveFig'), opt.saveFig = 0; end

figure;
fRot = newMtexFigure;
plotSSLIP_Context(opt);
plot(ebsdID, ebsdID.prop.rotationIDcor, 'micronbar', 'off'); title('Rotation');

validRot = ebsdID.prop.rotationIDcor(~isinf(ebsdID.prop.rotationIDcor) & ~isnan(ebsdID.prop.rotationIDcor));

% symmetric limits for rotation
maxRot = prctile(abs(validRot(:)), opt.clim_percentile);
if isempty(maxRot) || maxRot <= 0 || isnan(maxRot), maxRot = 0.01; end

caxis([-maxRot maxRot])
mtexColorMap(jet(512))
mtexColorbar('title', 'Rotation')

if isfield(opt, 'sizeAdjust')
    fRot.figSizeFactor = opt.sizeAdjust;
    fRot.drawNow;
end

if opt.saveFig && isfield(opt, 'plotname') && ~isempty(opt.plotname)
    saveFigure([opt.plotname, '_Rotation.png'])
    if isfield(opt, 'saveExt')
        saveFigure([opt.plotname, '_Rotation', opt.saveExt])
    end
end
end
