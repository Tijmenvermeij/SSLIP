function [ebsdID,opt] = SSLIP(ebsd,U,V,sSLocal,opt)
%% Function to Perform SSLIP (Slip System based Identification of Local Plasticity)
% For a list of slip systems (of a single crystal), with
% input a displacement field, compute slip system activity fields.
%
% Syntax
%   [PLOTEBSD,opt] = SSLIP(ebsd,U,V,sSLocal,opt)
%   [PLOTEBSD,opt] = SSLIP(ebsd,DeformationData,sSLocal,opt)
%   DeformationData contains either U,V or Hxx,Hxy,Hyx,Hyy.
%   Supplied gradients are ready for fitting, aligned with ebsd; use
%   filterSize = 0 (the gradient-input default) and coarsegrain = 1.
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
%   fields in "prop" field. For method 1, prop.solverExitFlag contains each
%   coneprog exit condition (NaN for pixels not solved).
%   With opt.enableRotation, prop.rotationIDcor contains signed small-angle
%   rotation in radians (one value per pixel); prop.slipIDcor retains only
%   physical slip systems, in opt.NoSs order.
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

% Accept Philipp's structured input alongside the original five-input call.
% Outputs and the flat options structure retain their existing conventions.
if nargin == 4 && isstruct(U)
    opt = sSLocal;
    sSLocal = V;
    DeformationData = U;
elseif nargin == 5
    DeformationData = struct('U',U,'V',V);
else
    error('SSLIP:InvalidInputs', ...
        'Use SSLIP(ebsd,U,V,sSLocal,opt) or SSLIP(ebsd,DeformationData,sSLocal,opt).');
end
isGradientInput = any(isfield(DeformationData,{'Hxx','Hxy','Hyx','Hyy'}));

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

% Absolute residual tolerance. If opt.threshResidualFraction is supplied,
% methods 1 and 3 use max(threshResidual, Eeff * threshResidualFraction).
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
    if isGradientInput
        opt.filterSize = 0;
    else
        opt.filterSize = 1;
    end
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

% Optional small-angle rotation correction, adapted from Philipp (PhilKro),
% PR #2: https://github.com/Tijmenvermeij/SSLIP/pull/2 (commit d3ee1a5).
% Only the combined coneprog solver implements this additional basis.
% minEeff still applies: use minEeff = 0 to include pure-rotation pixels.
if ~isfield(opt,'enableRotation')
    opt.enableRotation = 0;
end
if opt.enableRotation && opt.IDMethod ~= 1
    error('SSLIP:RotationRequiresMethod1', ...
        'enableRotation is supported only with IDMethod = 1.');
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
% transpose sSLocal if needed 
if length(sSLocal) > 1
    if size(sSLocal,2) ~= 1
        sSLocal = transpose(sSLocal);
    end
end

%% prepare data for SSLIP
% Preprocessing separation adapted from Philipp (PhilKro), PR #2.
[data,ebsdID] = preprocessSSLIP(ebsd,DeformationData,opt);
Hxx = data.Hxx;
Hxy = data.Hxy;
Hyx = data.Hyx;
Hyy = data.Hyy;

% define plotting name
plotName = [opt.casename '_' opt.comment '_' ];

% plot some fields (just for visualization and to check filtering)
if opt.plotDefGrad
    plotSSLIP_DeformationFields(ebsdID,opt);
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
    [slipIDcor,residualEeff,solverExitFlag] = SSLIPConeprogConstrMinAbs(sSAnalysis,Hxx,Hxy,Hyx,Hyy,opt);
    ebsdID.prop.solverExitFlag = solverExitFlag;
    if opt.enableRotation
        ebsdID.prop.rotationIDcor = slipIDcor(end,:)';
        slipIDcor = slipIDcor(1:end-1,:);
    end

    opt.plotname = [opt.plotname,'_constr_min'];

elseif opt.IDMethod == 2 % constrained slip ID
    [slipIDcor,residualEeff] = SSLIPConstr(sSAnalysis,Hxx,Hxy,Hyx,Hyy,opt);
        opt.plotname = [opt.plotname,'_constr'];

elseif opt.IDMethod == 3 % single slip ID
    opt.plotname = [opt.plotname,'_singleSlip'];
    % Return updated options as well: method 3 switches posConstr off.
    [slipIDcor,residualEeff,opt] = solveSSLIP_SingleSlip(sSAnalysis,Hxx,Hxy,Hyx,Hyy,opt);
    
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
