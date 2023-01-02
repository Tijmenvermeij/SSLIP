function [ebsdID,opt] = SSLIP(ebsd,U,V,sSLocal,opt)
%% Function to Perform SSLIP (Slip System based Identification of Local Plasticity)
% For a list of slip systems (of a single crystal), with
% input a displacement field, compute slip system activity fields.
%
% Syntax
%   [PLOTEBSD,opt] = SSLIP(ebsd,U,V,sSLocal,opt)
% Input
%   ebsd        - Mtex ebsd variable, used predominantly for the position
%               grid
%   U           - X-component of the displacement field (in um), same size as EBSD
%   V           - Y-component of the displacement field (in um), same size as EBSD
%   sSLocal     - list of slipSystems (MTex slipSystem objects), used for
%               identification, should already be rotated into local
%               crystal orientation
%   opt         - struct with options, see the defaults below. 

% Output
%   ebsdID      - Updated Mtex ebsd variable, with slip system activity
%   fields in "prop" field
%   opt         - struct with options, updated with defaults where applicable. 

% This function contains the SSLIP method as proposed in the paper 
% "T. Vermeij et al., Automated identification of slip system activity
% fields from digital image correlation data, Acta Mater. 243, 2022"
% DOI: https://doi.org/10.1016/j.actamat.2022.118502
% Please consider citing this paper when you use this code.
%
%%%
% Author: T. Vermeij
% // Eindhoven University of Technology, Hoefnagels Group
% Date: 30-11-2022
% the latest version of this code can be found on
% www.github.com/TijmenVermeij/SSLIP
%
% MTEX is required to use this code


%% Set default options, if needed

%%%
% set SSLIP method
% 1: constrained and minimized slip ID (As used predominantly in the SSLIP paper. 
% Contraint: || H^exp - H^their || < H_thresh. Minimzation of sum of absolute value of slip activities ) 

% 2: constrained slip ID (only solve problem based on the || H^exp - H^their || < H_thresh  constraint, no minimization)


% 3: single slip system ID (check for each pixel if a SINGLE system fits the activity. Used to "initialize" the SSLIP id om Figure 11 of the paper)
% Recommended as trial for uncertain/complex situations.

if ~isfield(opt,'IDMethod')
    opt.IDMethod = 1;
end
%%%

% threshold for residual
if ~isfield(opt,'threshResidual')
    opt.threshResidual = 0.01;
end

% set minimum effective strain for which SSLIP needs to be performed at a datapoint (i.e. skip
% pixels with low strain). For improved speed
if ~isfield(opt,'minEeff')
    opt.minEeff = 0.01;
end

% gaussian blurring filter size, applied to displacement field before computing gradients and
% performing SSLIP, to reduce noise (but also reduced spatial resolution)
%
% Defined in datapoints
% use 0 for no filtering
if ~isfield(opt,'filterSize')
    opt.filterSize = 1;
end

% coarse graining setting, 1 = no coarse graining, 2 = 2x2 pixels 
% coarsegrained into 1 pixel, ...
% (for improved speed)
if ~isfield(opt,'coarsegrain')
    opt.coarsegrain = 1;
end

% slip "numbers" systems to be used for SSLIP (in order of the sSLocal variable),
% default is all of them. 
if ~isfield(opt,'NoSs')
    opt.NoSs = 1:length(sSLocal);
end

% Positive constraint: constrain the slip amplitudes to be positive. This
% only works well if the slip systems are "configured" to have a positive
% amplitude under a certain load (which is normally assured in the main
% script, assuming e.g. uniaxial tension).
% How to "reconfigure" the slip system under complex loads is T.B.D.
if ~isfield(opt,'posConstr')
    opt.posConstr = 0;
end



%%%%
%%%% plotting options
%%%%


% whether or not to plot the def grad tensor and eq strain field, before
% performing slip ID
if ~isfield(opt,'plotDefGrad')
    opt.plotDefGrad = 0;
end

% colormap for strain and activity plots
if ~isfield(opt,'cmap')
    opt.cmap = viridis(256);
end

% layout for plotting multiple activity fields, e.g. [4 3] 
% means 4 rows and 3 columns
if ~isfield(opt,'layout')
    opt.layout = [];
end
% 
% % max strain/activity to plot
% if ~isfield(opt,'maxE')
%     opt.maxE = 0.1;
% end

% extra comments, maybe for plotting
if ~isfield(opt,'comment')
    opt.comment = '';
end

% casename, maybe for plotting
if ~isfield(opt,'casename')
    opt.casename = 'Nameless';
end

% use logscale for plotting?
if ~isfield(opt,'logscale')
    opt.logscale = 0;
end

% % min value for log plotting?
% if ~isfield(opt,'logmin')
%     opt.logmin = 0.01;
% end

if ~isfield(opt,'plotResidual')
    opt.plotResidual = 1;
end


%% some checks
% check if ebsd data is same size as U and V
if size(ebsd) ~= size(U) | size(ebsd) ~= size(V)
    error('ebsd data is not same as U and/or V')
end

% transpose sSLocal if needed 
if length(sSLocal) > 1
    if size(sSLocal,2) ~= 1
        sSLocal = transpose(sSLocal);
    end
end

%% prepare data for SSLIP
% include U and V, as a vector field, in ebsd and make sure data is gridded
ebsd.prop.U = vector3d(U,V,zeros(size(U)));
ebsd = ebsd.gridify;

% define plotting name
plotName = [opt.casename '_' opt.comment '_' ];

% extract position grid from ebsd variable
X = ebsd.x;
Y = ebsd.y;

% apply filtering on displacement field
if opt.filterSize ~= 0
    % filtersize:
    options.filt_std = opt.filterSize;
    
    % replace 0 values in disp field by NaNs (some DIC codes export 0
    % instead of NaN on non-correlated points)
    ebsd.prop.U.x(ebsd.prop.U.x == 0) = NaN;
    ebsd.prop.U.y(ebsd.prop.U.y == 0) = NaN;

    % filter displacements
    data = filterDisplacements(ebsd.prop.U.x,ebsd.prop.U.y,options);
else
    data.U = ebsd.prop.U.x;
    data.V = ebsd.prop.U.y;
end

% coarse graining to increase speed
crs = coarsegrainDisp(data.U,data.V,X(1,:),Y(:,1)',opt.coarsegrain);

%create dummy EBDS for plotting
ebsdID = dummyEBSDSimple(ebsd.orientations(1),crs.X,crs.Y);

% store displacement field after coarsegraining
data.U = crs.f;
data.V = crs.g;

% calculate numerical gradients (displacement gradient tensor components)
[Hxx, Hxy] = gradient(data.U,crs.pixelsize(1),crs.pixelsize(2));
[Hyx, Hyy] = gradient(data.V,crs.pixelsize(1),crs.pixelsize(2));


% calc effective shear strain (for plotting purposes)
Eeff = calcEffectiveE(Hxx,Hxy,Hyx,Hyy);

% store data in the coarsegrained ebsd variable
ebsdID.prop.Eeff = Eeff;

% store other fields in PLOTEBSD
ebsdID.prop.U = data.U;
ebsdID.prop.V = data.V;

ebsdID.prop.Hxx = Hxx;
ebsdID.prop.Hxy = Hxy;
ebsdID.prop.Hyx = Hyx;
ebsdID.prop.Hyy = Hyy;

% plot some fields (just for visualization and to check filtering)
if opt.plotDefGrad
    figure;
    f1=newMtexFigure('layout',[3 2]);
    
    plot(ebsdID,data.U,'micronbar','off'); title('U_x')
    
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
        saveFigure(['SSLIP_CGR_' num2str(opt.coarsegrain),'_Filt_' num2str(opt.filterSize), '_' plotName '_gradients.png'])
        if isfield(opt,'saveExt')
            saveFigure(['ssAnalysis_CoarseGr_' num2str(opt.coarsegrain),'_Filt_' num2str(opt.filterSize), '_' plotName '_disp_grad_tensor',opt.saveExt])
        end
    end
end


%% perform SSLIP analysis
fprintf(['Now running SSLIP, using method ',num2str(opt.IDMethod),'\n'])

% take the required slip systems
NoSs = opt.NoSs; 
sSAnalysis = sSLocal(NoSs);

% define plotting name
opt.plotname = ['SSLIP_CGR_' num2str(opt.coarsegrain),'_Filt_' num2str(opt.filterSize), '_' plotName];

% select and perform SSLIP method
if opt.IDMethod == 1 % combined & minimized slip ID
    [slipIDcor,residualEeff] = SSLIPConeprogConstrMinAbs(sSAnalysis,Hxx,Hxy,Hyx,Hyy,opt);

    opt.plotname = [opt.plotname,'_constr_min'];

elseif opt.IDMethod == 2 % constrained slip ID
    [slipIDcor,residualEeff] = SSLIPConstr(sSAnalysis,Hxx,Hxy,Hyx,Hyy,opt);
        opt.plotname = [opt.plotname,'_constr'];

elseif opt.IDMethod == 3 % single slip ID
    opt.plotname = [opt.plotname,'_singleSlip'];

    % initialize matrices
    slipIDcor = zeros(length(NoSs),length(ebsdID));
    residualEeff = zeros(length(NoSs),length(ebsdID));

    % don't use constraints at single slip ID (since it has no benefit
    % usually, and is much faster without it)
    if opt.posConstr == 1
        warning('option "posConstr" changed to 0, since it has no added value for single slip ID, and is much slower')
        opt.posConstr = 0;
    end

    % loop over single slip systemsm, and perform SSLIP (without
    % minimization)
    for j = 1:length(NoSs)
        [slipIDcor(j,:),residualEeff(j,:)] = SSLIPConstr(sSAnalysis(j),Hxx,Hxy,Hyx,Hyy,opt);
    end
   
    % use residual to "filter" the fields, only leaving activities with low
    % residual
    if ~isfield(opt,'threshResidualFraction')
        % "clean" slipID by using theshold on residual, based on single systems
        goodData = residualEeff < opt.threshResidual;
    else
        % "clean" slipID by using theshold on residual, based on single
        % systems, based on factor of total effective shear (not used
        % often)
        threshResidual = Eeff(:) * opt.threshResidualFraction;
        goodData = residualEeff < threshResidual';
    end

    %%% not used a lot:
    % this option leaves only 1 slip system active at each datapoint, which has the lowest residual.
    % Although this might make sense, noise can be detrimental here
    if isfield(opt,'singleSlipPerPixel')
        if opt.singleSlipPerPixel
            [minThres,minThreshInd] = min(residualEeff,[],1);
            slipIDcor = zeros(size(slipIDcor));
            for j=1:length(ebsdID)
                slipIDcor(minThreshInd(j),j) = slipIDcor(minThreshInd(j),j);
            end
        end
    end
    
    
    % make all slip activities, with residual abovet threshold, zero.
    % (might be better to make it NaN, but not nice for plotting)
    slipIDcor(~goodData) = 0;
    
else
    error('IDMethod unknown (should be 1, 2, 3)')
end


% plot SSLIP results
if opt.plotSSLIP
    plotSSLIP(slipIDcor,residualEeff,ebsdID,sSLocal,opt)
end

ebsdID.prop.residualEeff = residualEeff;
ebsdID.prop.slipIDcor = slipIDcor;

% save(IDoptions.plotname,'slipIDcor','residualEeff','PLOTEBSD','IDoptions');


end




