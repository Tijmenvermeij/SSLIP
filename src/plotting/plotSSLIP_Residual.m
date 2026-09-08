function plotSSLIP_Residual(ebsdID, opt)
% Plot the stored residual field with the existing SSLIP color scales.
% When residualScaleSame is supplied, opt.caxisMinMax carries the activity
% limits; plotSSLIP passes these explicitly, retaining its existing behavior.
% Adapted from Philipp (PhilKro), PR #2, original commit d3ee1a5:
% https://github.com/Tijmenvermeij/SSLIP/pull/2

if nargin < 2, opt = struct; end
if ~isfield(opt,'cmap'), opt.cmap = viridis(256); end
if ~isfield(opt,'logscale'), opt.logscale = 0; end
if ~isfield(opt,'saveFig'), opt.saveFig = 0; end
residualEeff = ebsdID.prop.residualEeff;

figure;
meanResidual = mean(residualEeff(:),'omitnan');
plot(ebsdID, residualEeff,'micronbar','off' ); title(['residual Eeff,mean=',num2str(meanResidual)]);
if isfield(opt,'residualScaleSame')
    caxis(opt.caxisMinMax)
    if opt.logscale
        set(gca,'colorscale','log')
    end
else
    % Show the full finite residual range, including constant residuals.
    maxResidual = max(residualEeff(isfinite(residualEeff)));
    if isempty(maxResidual) || maxResidual <= 0
        maxResidual = 1; % Non-degenerate limits for zero or missing data.
    end
    caxis([0 maxResidual])
end
mtexColorMap(opt.cmap)
mtexColorbar
if opt.saveFig
    saveFigure([opt.plotname, '_Residual.png'])
    if isfield(opt,'saveExt')
        saveFigure([opt.plotname, '_Residual',opt.saveExt])
    end
end
end
