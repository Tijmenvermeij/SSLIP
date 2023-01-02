function plotSSLIP(slipIDcor,residualEeff,ebsdID,sSLocal,opt)

% define stress tensor, only used to report SF
if ~isfield(opt,'stress')
    warning('be careful, stress is assumed to be in x direction (used to calc SF for plotting), define opt.stress = ... if necessary')
    opt.stress = stressTensor.uniaxial(xvector);
end

% determine if activity fields need to be plotted in individual (single)
% plots, instead of combinining them
if ~isfield(opt,'plotSingle')
    opt.plotSingle = 0;
end

% plot slip trace and direction in activity maps
if ~isfield(opt,'plotTraces')
    opt.plotTraces = 1;
end

% choose to plot empty fields to use same structure (immediately assign num
% of systems)
if ~isfield(opt,'plotEmptyFields')
    opt.plotEmptyFields = 0;
end


% max fields plotted per figure (for cases of more than 24 slip systems, to
% avoid that figures become too large
maxFields = 24;



% prepare for plotting

NoSs = opt.NoSs;

clear caxis
if isempty(opt.layout)
    figure;
    f1=newMtexFigure;
else
    figure;
    f1=newMtexFigure('layout',opt.layout);
end

% initialize this variable, which keeps track of the color limits
caxisMinMax = [0 0];


% calc SF, used in the title of each slip activity field
SF = abs(SchmidFactor(sSLocal(NoSs),opt.stress));
SF = round(SF*100)/100;

% if logaritmic scale needs to be used for plotting, and positive constraint was not used, plot the absolute
% value of slip activities!
if opt.logscale && ~opt.posConstr
    slipIDcor = abs(slipIDcor);
    warning('Since logaritmic plotting is required (opt.logscale==1), while pos. constr. was not used (opt.posConstr==0), absolute values of slip activities are plotted')
end

%% start loop for plotting. Loop over every slip system.
%%%
if ~opt.plotSingle % in this case, all activities are plotted in one figure (using Mtex "nextAxis" command)
    
    for i = 1:min([length(NoSs) , maxFields]) % only plot maximum "maxFields" fields per figure
        nextAxis
        % plot activity field
        ha(i) = plot(ebsdID, slipIDcor(i,:),'micronbar','off' ); title(['#',num2str(NoSs(i)),'-SF=',num2str(SF(i))]);
        % if required, plot slip system trace and slip system direction
        % (in-plane) over the map
        if isfield(opt,'plotTraces')
            if opt.plotTraces
                ebsdTrace = ebsdID(round(size(ebsdID,1)/2),round(size(ebsdID,2)/2));
                hold on
                quiver(ebsdTrace,0.2*(ebsdID.xmax-ebsdID.xmin) * sSLocal(NoSs(i)).trace,'color','r');
                quiver(ebsdTrace,0.2*(ebsdID.xmax-ebsdID.xmin) * sSLocal(NoSs(i)).b.normalize,'color','r');
                hold off
            end
        end
        
%         % extract the colorscale of this plot
%         caxis_i = caxis;
%         caxis_i(caxis_i==0 | caxis_i==1) = NaN;
% 
%         % remember max colorscale for later usage
%         if caxis_i(2) > caxisMinMax(2)
%             caxisMinMax(2) = caxis_i(2);
%         end
%         if caxis_i(1) < caxisMinMax(1)
%             caxisMinMax(1) = caxis_i(1);
%         end

    end
    
    % set all colorbars to same scale, depending on defined plotting options
    if opt.logscale
        caxisMinMax(1) = opt.logmin;
    elseif isfield(opt,'logmin')
        caxisMinMax(1) = opt.logmin;
    elseif isfield(opt,'maxE')
        caxisMinMax(1) = -opt.maxE;
    else
        caxisMinMax(1) = min(slipIDcor(~isinf(slipIDcor)),[],'all','omitnan');
    end

    if isfield(opt,'maxE')
        caxisMinMax(2) = opt.maxE;
    else
        caxisMinMax(2) = max(slipIDcor(~isinf(slipIDcor)),[],'all','omitnan');
    end

    % loop over all axes to adjust color scale and if needed logaritmic
    % scale
    for i=1:min([length(NoSs) , maxFields])

        ha(i).Parent.CLim = caxisMinMax;
        if opt.logscale
            ha(i).Parent.ColorScale = 'log';
        end
    end

    % use colorbar
    mtexColorMap(opt.cmap)
    if opt.logscale && ~opt.posConstr
        mtexColorbar('title','|\gamma|')
    else
        mtexColorbar('title','\gamma')
    end
    
    
    % if necessary, adjust size of figure
    if isfield(opt,'sizeAdjust')
        f1.figSizeFactor = opt.sizeAdjust;
        f1.innerPlotSpacing = f1.innerPlotSpacing * opt.sizeAdjust;
        
        f1.drawNow;
    end
    
    % if necessary, save figure
    if opt.saveFig
        saveFigure([opt.plotname, '.png'])
        if isfield(opt,'saveExt')
            saveFigure([opt.plotname, opt.saveExt])
        end
    end

    %%% second plot, if more than 24 systems. Same thing as before
    if length(NoSs) > maxFields
        % plot results
        clear caxis
        if isempty(opt.layout)
            figure;
            f1=newMtexFigure;
        else
            f1=newMtexFigure('layout',opt.layout);
        end
        caxisMinMax = [0 0];

        for i = maxFields+1 : length(NoSs) % plot the rest of the fields
            nextAxis
            ha(i-maxFields) = plot(ebsdID, slipIDcor(i,:),'micronbar','off' ); title(['#',num2str(NoSs(i)),'-SF=',num2str(SF(i))]);
            if isfield(opt,'plotTraces')
                if opt.plotTraces
                    ebsdTrace = ebsdID(round(size(ebsdID,1)/2),round(size(ebsdID,2)/2));
                    hold on
                    quiver(ebsdTrace,0.2*(ebsdID.xmax-ebsdID.xmin) * sSLocal(NoSs(i)).trace,'color','r');
                    quiver(ebsdTrace,0.2*(ebsdID.xmax-ebsdID.xmin) * sSLocal(NoSs(i)).b.normalize,'color','r');
                    hold off
                end
            end

            caxis_i = caxis;
            caxis_i(caxis_i==0 | caxis_i==1) = NaN;

            % remember max colorscale for later adjusting
            if caxis_i(2) > caxisMinMax(2)
                caxisMinMax(2) = caxis_i(2);
            end
            if caxis_i(1) < caxisMinMax(1)
                caxisMinMax(1) = caxis_i(1);
            end

        end
        % set all colorbars to same scale, depending on defined plotting options
        if opt.logscale
            caxisMinMax(1) = opt.logmin;
        elseif isfield(opt,'logmin')
            caxisMinMax(1) = opt.logmin;
        elseif isfield(opt,'maxE')
            caxisMinMax(1) = -opt.maxE;
        end

        if isfield(opt,'maxE')
            caxisMinMax(2) = opt.maxE;
        end


        for i=1:length(NoSs)-maxFields

            ha(i).Parent.CLim = caxisMinMax;
            if opt.logscale
                ha(i).Parent.ColorScale = 'log';
            end
        end

        mtexColorMap(opt.cmap)
        
        if opt.logscale && ~opt.posConstr
            mtexColorbar('title','|\gamma|')
        else
            mtexColorbar('title','\gamma')
        end
        
        if isfield(opt,'sizeAdjust')
            f1.figSizeFactor = opt.sizeAdjust;
            f1.innerPlotSpacing = f1.innerPlotSpacing * opt.sizeAdjust;

            f1.drawNow;
        end
        if opt.saveFig
            saveFigure([opt.plotname, '24plus.png'])
            if isfield(opt,'saveExt')
                saveFigure([opt.plotname, '24plus',opt.saveExt])
            end
        end
    end
else
    %%% if opt.plotSingle == 1, plot all fields in single plots, and store
    %%% in a folder
    close all
    
    if opt.saveFig
        if ~exist(opt.plotname,'dir')
            mkdir(opt.plotname)
        end
    end
    
    % plot single figures
    for i = 1:min(length(NoSs)) 
        figure;
        f1=newMtexFigure;
        ha = plot(ebsdID, slipIDcor(i,:),'micronbar','off' ); title(['#',num2str(NoSs(i)),'-SF=',num2str(SF(i))]);
        if isfield(opt,'plotTraces')
            if opt.plotTraces
                ebsdTrace = ebsdID(round(size(ebsdID,1)/2),round(size(ebsdID,2)/2));
                hold on
                quiver(ebsdTrace,0.2*(ebsdID.xmax-ebsdID.xmin) * sSLocal(NoSs(i)).trace,'color','r','linewidth',5);
                quiver(ebsdTrace,0.2*(ebsdID.xmax-ebsdID.xmin) * sSLocal(NoSs(i)).b.normalize,'color','r','linewidth',5);
                hold off
            end
        end
        
        if opt.logscale
            set(gca,'colorscale','log')
            cmin = opt.logmin;
        elseif isfield(opt,'logmin')
            cmin = opt.logmin;
        elseif isfield(opt,'maxE')
            cmin = -opt.maxE;
        end
        
        if isfield(opt,'maxE')
            cmax = opt.maxE;
        else
            cmax = max(slipIDcor,[],'all');
        end
        
        caxis([cmin cmax])
        mtexColorMap(opt.cmap)
        
        if opt.logscale && ~opt.posConstr
            mtexColorbar('title','|\gamma|')
        else
            mtexColorbar('title','\gamma')
        end

        set(gca,'fontSize',50)
        f1.drawNow;
        
        
        if isfield(opt,'sizeAdjust')
            f1.figSizeFactor = opt.sizeAdjust;
            f1.drawNow;
        end
        

        if opt.saveFig
            saveFigure([opt.plotname,'/sS',num2str(NoSs(i)) '.png'])
            if isfield(opt,'saveExt')
                saveFigure([opt.plotname,'/sS',num2str(NoSs(i)),opt.saveExt])
            end
            
        end

    end

    
end








% if required, plot residual
if opt.plotResidual && (opt.IDMethod == 1 || opt.IDMethod == 2)
    figure;
    meanResidual = mean(residualEeff(:),'omitnan');
    stdResidual = std(residualEeff(:),'omitnan');
    plot(ebsdID, residualEeff,'micronbar','off' ); title(['residual Eeff,mean=',num2str(meanResidual)]);
    if isfield(opt,'residualScaleSame')
        caxis(caxisMinMax)
        if opt.logscale
            set(gca,'colorscale','log')
        end
    else
        caxis([0 2*stdResidual])
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

end

