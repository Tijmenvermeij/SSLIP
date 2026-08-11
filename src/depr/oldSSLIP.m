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

if isfield(opt,'posConstr')
    if opt.posConstr == 1
        % sSLocal.Schmid  
        warning('Make sure that the slip systems are "configured" to have a positive SF amplitude in the assumed loading. This has to be done in the main script.')
        
    end
end
% Positive constraint: constrain the slip amplitudes to be positive. This
% only works well if the slip systems are "configured" to have a positive
% amplitude under a certain load (which is normally assured in the main
% script, assuming e.g. uniaxial tension).
% How to "reconfigure" the slip system under complex loads is T.B.D.
if ~isfield(opt,'posConstr')
    opt.posConstr = 0;
end

% enable rotation as pseudo-slip system
if ~isfield(opt,'enableRotation')
    opt.enableRotation = 0;
end



% Use pre-computed gradient fields directly from opt if provided
if ~isfield(opt, 'plotResidual')
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
if isfield(opt, 'Hxx') && isfield(opt, 'Hxy') && isfield(opt, 'Hyx') && isfield(opt, 'Hyy')
    fprintf('  [SSLIP Bypass] Utilizing pre-computed displacement gradient fields from opt structure.\n');
    Hxx = opt.Hxx; Hxy = opt.Hxy;
    Hyx = opt.Hyx; Hyy = opt.Hyy;
    
    % If raw grid gradients were passed, apply identical Gaussian filtering and superpixel coarse-graining
    if isequal(size(Hxx), size(ebsd.x))
        if opt.filterSize ~= 0
            filt_opts.filt_std = opt.filterSize;
            fH1 = filterDisplacements(Hxx, Hxy, filt_opts);
            fH2 = filterDisplacements(Hyx, Hyy, filt_opts);
            Hxx = fH1.U; Hxy = fH1.V;
            Hyx = fH2.U; Hyy = fH2.V;
        end
        if opt.coarsegrain > 1
            crsH1 = coarsegrainDisp(Hxx, Hxy, X(1,:), Y(:,1)', opt.coarsegrain);
            crsH2 = coarsegrainDisp(Hyx, Hyy, X(1,:), Y(:,1)', opt.coarsegrain);
            Hxx = crsH1.f; Hxy = crsH1.g;
            Hyx = crsH2.f; Hyy = crsH2.g;
        end
    end
else
    error('Please use pre-computed displacement gradient fields from opt structure.')
    [Hxx, Hxy] = gradient(data.U,crs.pixelsize(1),crs.pixelsize(2));
    [Hyx, Hyy] = gradient(data.V,crs.pixelsize(1),crs.pixelsize(2));
end


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

% perform SSLIP analysis
fprintf(['Now running SSLIP, using method ',num2str(opt.IDMethod),'\n'])

% take the required slip systems
NoSs = opt.NoSs; 
sSAnalysis = sSLocal(NoSs);

% define plotting name only if not already specified by user
custom_plotname = isfield(opt, 'plotname') && ~isempty(opt.plotname);
if ~custom_plotname
    opt.plotname = ['SSLIP_CGR_' num2str(opt.coarsegrain),'_Filt_' num2str(opt.filterSize), '_' plotName];
end

% select and perform SSLIP method
if opt.IDMethod == 1 % combined & minimized slip ID
    [slipIDcor,residualEeff] = SSLIPConeprogConstrMinAbs(sSAnalysis,Hxx,Hxy,Hyx,Hyy,opt);

    if ~custom_plotname, opt.plotname = [opt.plotname,'_constr_min']; end

elseif opt.IDMethod == 2 % constrained slip ID
    [slipIDcor,residualEeff] = SSLIPConstr(sSAnalysis,Hxx,Hxy,Hyx,Hyy,opt);
    if ~custom_plotname, opt.plotname = [opt.plotname,'_constr']; end

elseif opt.IDMethod == 3 % single slip ID
    if ~custom_plotname, opt.plotname = [opt.plotname,'_singleSlip']; end

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
        threshResidual(threshResidual<opt.threshResidual) = opt.threshResidual;
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
if isfield(opt, 'plotSSLIP') && opt.plotSSLIP
    try
        plotSSLIP(ebsdID,sSLocal,opt);
    catch ME
        % Display a warning message containing the error identifier and message
        warning('SSLIP Plotting Error (%s): %s', ME.identifier, ME.message);
    end
end

ebsdID.prop.residualEeff = residualEeff;

if isfield(opt, 'enableRotation') && opt.enableRotation
    % Explicitly extract the rotation pseudo-slip system and store separately
    ebsdID.prop.rotationIDcor = slipIDcor(end, :);
    slipIDcor = slipIDcor(1:end-1, :);
end
ebsdID.prop.slipIDcor = slipIDcor;

end

function clim = robustCaxis(X, p_val)

    valid = X(~isinf(X) & ~isnan(X));
    valid = valid(:);

    % Robust upper limit based on absolute value percentile
    lim = prctile(abs(valid), p_val);
    if lim <= 0 || isnan(lim), lim = 0.01; end

    clim = [-lim lim];

end


