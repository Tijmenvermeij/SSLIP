function plotSSLIP_Residual(ebsdID, opt)
% plotSSLIP_Residual Plots the residual error (Eeff) map.
% NOTE: Supports context boundaries via plotSSLIP_Context (requires opt.plotGrainContext=1 and opt.grains)

if nargin < 2
    opt = struct;
end

if ~isfield(opt, 'cmap'), opt.cmap = parula(256); try opt.cmap = viridis(256); catch, end; end
if ~isfield(opt, 'saveFig'), opt.saveFig = 0; end

residualEeff = ebsdID.prop.residualEeff;

figure;
meanResidual = mean(residualEeff(:),'omitnan');
stdResidual = std(residualEeff(:),'omitnan');

plotSSLIP_Context(opt);
plot(ebsdID, residualEeff,'micronbar','off' ); 
title(['residual Eeff, mean=',num2str(meanResidual)]);

if isfield(opt, 'caxisMinMax') && isfield(opt, 'residualScaleSame') && opt.residualScaleSame
    caxis(opt.caxisMinMax)
    if isfield(opt, 'logscale') && opt.logscale
        set(gca,'colorscale','log')
    end
else
    if stdResidual <= 0 || isnan(stdResidual), stdResidual = 0.05; end
    caxis([0 2*stdResidual])
end

mtexColorMap(opt.cmap)
mtexColorbar

if opt.saveFig && isfield(opt, 'plotname') && ~isempty(opt.plotname)
    saveFigure([opt.plotname, '_Residual.png'])
    if isfield(opt,'saveExt')
        saveFigure([opt.plotname, '_Residual',opt.saveExt])
    end
end
end
