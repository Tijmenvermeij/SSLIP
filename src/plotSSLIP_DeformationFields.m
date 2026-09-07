function plotSSLIP_DeformationFields(ebsdID, opt)
% Plot processed displacements, effective strain, and displacement gradients.
% Use the options returned by SSLIP to reproduce its scales and saved names.
% Adapted from Philipp (PhilKro), PR #2, original commit d3ee1a5:
% https://github.com/Tijmenvermeij/SSLIP/pull/2
% Retains the existing SSLIP figure layout, limits, and export conventions.

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

plot(ebsdID,ebsdID.prop.U,'micronbar','off'); title('U_x')

nextAxis

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
caxis([cmin cmax])

nextAxis
plot(ebsdID,Hxx,'micronbar','off'); title('H_{11}'); caxis([-1*max([max(Hxx(:)),abs(min(Hxx(:)))]) max([max(Hxx(:)),abs(min(Hxx(:)))])]);
nextAxis
plot(ebsdID,Hxy,'micronbar','off'); title('H_{12}'); caxis([-1*max([max(Hxy(:)),abs(min(Hxy(:)))]) max([max(Hxy(:)),abs(min(Hxy(:)))])]);
nextAxis
plot(ebsdID,Hyx,'micronbar','off'); title('H_{21}'); caxis([-1*max([max(Hyx(:)),abs(min(Hyx(:)))]) max([max(Hyx(:)),abs(min(Hyx(:)))])]);
nextAxis
plot(ebsdID,Hyy,'micronbar','off'); title('H_{22}'); caxis([-1*max([max(Hyy(:)),abs(min(Hyy(:)))]) max([max(Hyy(:)),abs(min(Hyy(:)))])]);
mtexColorbar

f1.children(1).Colormap = jet(512);
f1.children(2).Colormap = opt.cmap;
f1.children(3).Colormap = jet(512);
f1.children(4).Colormap = jet(512);
f1.children(5).Colormap = jet(512);
f1.children(6).Colormap = jet(512);

if isfield(opt,'DefGradLim')
    f1.children(3).CLim = opt.DefGradLim ;
    f1.children(4).CLim = opt.DefGradLim ;
    f1.children(5).CLim = opt.DefGradLim ;
    f1.children(6).CLim = opt.DefGradLim ;
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
