function plotSSLIP_DeformationFields(ebsdID, opt)
% Plot processed displacements, effective strain, and displacement gradients.
% Use the options returned by SSLIP to reproduce its scales and saved names.
% Adapted from Philipp (PhilKro), PR #2, original commit d3ee1a5:
% https://github.com/Tijmenvermeij/SSLIP/pull/2
% Retains the displacement-input layout, limits, and export conventions.
% Gradient-only input omits the unavailable displacement panel. Constant zero
% or missing fields use non-degenerate limits without modifying stored data.

if nargin < 2, opt = struct; end
if ~isfield(opt,'cmap'), opt.cmap = viridis(256); end
if ~isfield(opt,'logscale'), opt.logscale = 0; end
if ~isfield(opt,'saveFig'), opt.saveFig = 0; end

Eeff = ebsdID.prop.Eeff;
Hxx = ebsdID.prop.Hxx;
Hxy = ebsdID.prop.Hxy;
Hyx = ebsdID.prop.Hyx;
Hyy = ebsdID.prop.Hyy;

figure;
f1=newMtexFigure('layout',[3 2]);

hasDisplacements = isfield(ebsdID.prop,'U');
if hasDisplacements
    plot(ebsdID,ebsdID.prop.U,'micronbar','off'); title('U_x')
    nextAxis
end

plot(ebsdID,Eeff,'micronbar','off'); title('E_{eff}');

if opt.logscale
    set(gca,'colorscale','log')
    cmin = opt.logmin;
else
    cmin = 0;
end

if isfield(opt,'maxE')
    cmax = opt.maxE;
else
    cmax = max(Eeff(:));
end
if isempty(cmax) || ~isfinite(cmax) || cmax <= cmin
    cmax = cmin + 1;
end
caxis([cmin cmax])

components = {Hxx,Hxy,Hyx,Hyy};
titles = {'H_{11}','H_{12}','H_{21}','H_{22}'};
for k = 1:4
    nextAxis
    values = components{k};
    plot(ebsdID,values,'micronbar','off'); title(titles{k});
    limit = max(abs(values(isfinite(values))));
    if isempty(limit) || limit <= 0, limit = 1; end
    caxis([-limit limit]);
end
mtexColorbar

strainAxis = 1 + hasDisplacements;
if hasDisplacements, f1.children(1).Colormap = jet(512); end
f1.children(strainAxis).Colormap = opt.cmap;
for k = strainAxis+(1:4)
    f1.children(k).Colormap = jet(512);
    if isfield(opt,'DefGradLim'), f1.children(k).CLim = opt.DefGradLim; end
end

if isfield(opt,'fontSize')
    set(findall(gcf,'-property','FontSize'),'FontSize',opt.fontSize)
end

if isfield(opt,'sizeAdjust')
    f1.figSizeFactor = opt.sizeAdjust;
    f1.innerPlotSpacing = f1.innerPlotSpacing * opt.sizeAdjust;

    f1.drawNow;
end

if opt.saveFig
    plotName = [opt.casename '_' opt.comment '_'];
    saveFigure(['SSLIP_CGR_' num2str(opt.coarsegrain),'_Filt_' num2str(opt.filterSize), '_' plotName '_gradients.png'])
    if isfield(opt,'saveExt')
        saveFigure(['ssAnalysis_CoarseGr_' num2str(opt.coarsegrain),'_Filt_' num2str(opt.filterSize), '_' plotName '_disp_grad_tensor',opt.saveExt])
    end
end
end
