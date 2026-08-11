function plotSSLIP_DeformationFields(ebsdID, opt)
% plotSSLIP_DeformationFields Plots the DIC deformation gradients (Hxx, Hxy, Hyx, Hyy) or displacements (U, V).
% NOTE: Supports context boundaries via plotSSLIP_Context (requires opt.plotGrainContext=1 and opt.grains)  used internally by the SSLIP solver.

if nargin < 2
    opt = struct;
end

if ~isfield(opt, 'cmap'), opt.cmap = parula(256); try opt.cmap = viridis(256); catch, end; end
if ~isfield(opt, 'clim_percentile'), opt.clim_percentile = 99; end
if ~isfield(opt, 'saveFig'), opt.saveFig = 0; end

figure;
f1 = newMtexFigure('layout',[3 2]);

% U_x
plotSSLIP_Context(opt);
plot(ebsdID, ebsdID.prop.U, 'micronbar', 'off'); title('U_x');

% E_eff
nextAxis;
plotSSLIP_Context(opt);
plot(ebsdID, ebsdID.prop.Eeff, 'micronbar', 'off'); title('E_{eff}');

Eeff = ebsdID.prop.Eeff;
validData = Eeff(~isinf(Eeff) & ~isnan(Eeff));
cmin = 0;
cmax = prctile(validData, opt.clim_percentile);
if isempty(cmax) || isnan(cmax) || cmax <= cmin, cmax = cmin + 0.05; end
caxis([cmin cmax]);

% H components
nextAxis; plotSSLIP_Context(opt); plot(ebsdID, ebsdID.prop.Hxx, 'micronbar', 'off'); title('H_{11}'); caxis(robustCaxis(ebsdID.prop.Hxx, opt.clim_percentile));
nextAxis; plotSSLIP_Context(opt); plot(ebsdID, ebsdID.prop.Hxy, 'micronbar', 'off'); title('H_{12}'); caxis(robustCaxis(ebsdID.prop.Hxy, opt.clim_percentile));
nextAxis; plotSSLIP_Context(opt); plot(ebsdID, ebsdID.prop.Hyx, 'micronbar', 'off'); title('H_{21}'); caxis(robustCaxis(ebsdID.prop.Hyx, opt.clim_percentile));
nextAxis; plotSSLIP_Context(opt); plot(ebsdID, ebsdID.prop.Hyy, 'micronbar', 'off'); title('H_{22}'); caxis(robustCaxis(ebsdID.prop.Hyy, opt.clim_percentile));

mtexColorbar;

f1.children(1).Colormap = jet(512);
f1.children(2).Colormap = opt.cmap;
for i = 3:6, f1.children(i).Colormap = jet(512); end

if isfield(opt, 'sizeAdjust')
    f1.figSizeFactor = opt.sizeAdjust;
    f1.innerPlotSpacing = f1.innerPlotSpacing * opt.sizeAdjust;
    f1.drawNow;
end

if opt.saveFig && isfield(opt, 'plotname') && ~isempty(opt.plotname)
    saveFigure([opt.plotname, '_FilteredGradients.png']);
    if isfield(opt, 'saveExt')
        saveFigure([opt.plotname, '_FilteredGradients', opt.saveExt]);
    end
end

end

function clim = robustCaxis(X, p_val)
    valid = X(~isinf(X) & ~isnan(X));
    lim = prctile(abs(valid(:)), p_val);
    if isempty(lim) || lim <= 0 || isnan(lim), lim = 0.01; end
    clim = [-lim lim];
end
