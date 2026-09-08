function [ebsdID,opt,sSLocal] = SSLIP(ebsd,U,V,sSLocal,opt)
%% Function to Perform SSLIP (Slip System based Identification of Local Plasticity)
% For a list of slip systems (of a single crystal), with
% input a displacement field, compute slip system activity fields.
%
% Syntax
%   [PLOTEBSD,opt,sSLocal] = SSLIP(ebsd,U,V,sSLocal,opt)
%   [PLOTEBSD,opt,sSLocal] = SSLIP(ebsd,DeformationData,sSLocal,opt)
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
%   sSLocal     - full slip-system list after any stress alignment. Save this
%                 with the fit; activity rows correspond to sSLocal(opt.NoSs).

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

% Positive constraint: constrain slip amplitudes to be nonnegative.
% For methods 1 and 2, opt.stress aligns each slip direction with the resolved
% shear stress. Supply one stressTensor in specimen coordinates or a vector3d
% tension direction. Without stress, the caller must orient the systems.
% Method 3 always uses signed activities and does not align the systems.
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

%% orient slip directions for positive-constrained fits
% Adapted from Philipp (PhilKro), PR #2, commit 9c393a3:
% https://github.com/Tijmenvermeij/SSLIP/pull/2
% Return the full aligned list so NoSs still refers to the original labels.
% Method 3 switches posConstr off, so it must keep the supplied directions.
if opt.posConstr && any(opt.IDMethod == [1 2])
    if isfield(opt,'stress') && ~isempty(opt.stress)
        if isa(opt.stress,'stressTensor')
            sig = opt.stress;
        elseif isa(opt.stress,'vector3d')
            direction = opt.stress;
            if ~isscalar(direction) || any(~isfinite([direction.x direction.y direction.z])) || norm(direction) == 0
                error('SSLIP:InvalidStressValue','Supply one finite, nonzero tension direction.');
            end
            sig = stressTensor.uniaxial(direction);
        else
            error('SSLIP:InvalidStressType','opt.stress must be a stressTensor or vector3d.');
        end
        if ~isscalar(sig) || any(~isfinite(sig.matrix),'all')
            error('SSLIP:InvalidStressValue','Supply one finite stress tensor in specimen coordinates.');
        end
        % Use the physical tensor contraction, without Schmid-factor or CRSS
        % normalization. Zero resolved shear leaves the direction unchanged.
        resolvedShear = sig : transpose(sSLocal.deformationTensor);
        directionSigns = reshape(sign(resolvedShear),size(sSLocal));
        directionSigns(directionSigns == 0) = 1;
        sSLocal.b = directionSigns .* sSLocal.b;
    else
        warning('SSLIP:MissingStress', ...
            'posConstr is enabled without opt.stress; using the supplied slip directions.');
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
