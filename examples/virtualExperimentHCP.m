%% Script to generate a (HCP) virtual experiment with some slips steps, including identification using the SSLIP method
% This script replicates the HCP virtual experiment and its identification,
% as shown in Figures 2,3,4 in the paper "T. Vermeij et al., Automated identification of slip system activity fields from digital image correlation data, Acta Mater. 243, 2022
% doi: https://doi.org/10.1016/j.actamat.2022.118502
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


%% Initiallize Mtex

clear
close all

addpath(fullfile(pwd,'..','src'));

setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoplane');

%% define a crystal orientation and slip systems

% HCP crystal symmetry
CS = crystalSymmetry('6/mmm',[0.266 0.266 0.495]);

nameOri = 'ori1';

% define orientation
ori = orientation(rotation.byAxisAngle(zvector,45*degree)*rotation.byAxisAngle(xvector,45*degree)*rotation.byAxisAngle(zvector,15*degree),CS);

% define the slip systems for each slip family
sSBas = symmetrise(slipSystem.basal(CS),'antipodal');
sSPris = symmetrise(slipSystem.prismaticA(CS),'antipodal');
sSPyr = symmetrise(slipSystem.pyramidalA(CS),'antipodal');
sSPyrca = symmetrise(slipSystem.pyramidalCA(CS),'antipodal');
sSPyr2ca = symmetrise(slipSystem.pyramidal2CA(CS),'antipodal');

% put the (here-considered) slip systems into one list
sS = [sSBas; sSPris; sSPyrca; sSPyr2ca];

% get the rotated slip systems
sSLocal = ori * sS;

% assuming uniaxial tension in x-direction, reconfigure slip systems (by changing signs of b and n) to be
% "positive". i.e., slip amplitudes will be positive under uniaxial
% tension and Schmid factors will be positive
loadDir = xvector;
sSLocal.b = sign(cos(angle(loadDir,sSLocal.b))) .* sSLocal.b;
sSLocal.n = sign(sSLocal.SchmidFactor(loadDir)) .* sSLocal.n;

%% create (dummy) EBSD and DIC maps
% some dimensions, in um
sizeArea = [7 7];
psizeEBSD = 0.03;

% create EBSD variable
[prop.x,prop.y] =  meshgrid(0:psizeEBSD:sizeArea(1),0:psizeEBSD:sizeArea(2));
ebsd = EBSD(repmat(ori,[numel(prop.x),1]),ones(size(prop.x)),CS,prop);
ebsd = ebsd.gridify;

%% Define the virtual slip system activities: system numbers, amplitudes, positions
% in this case, there will be 4 slip "activities" (traces), therefore the
% variables have length 4

% slip system numbers
systems = [3 3 14 24]';

% name of the virtual exp
nameSlip = 'test1';

% slip system amplitudes, in um
amplitudes =[0.02 0.04 0.03 0.02]';

% locations of the "center" of the activities, first column x-position,
% second column y-position, in um
locations = [4.25, sizeArea(2)/2
             2.75, sizeArea(2)/2
             3.25, sizeArea(2)/2
             3.75, sizeArea(2)/2];

% define how much displacement field noise to add (this is actually done
% differently in the paper, that's why it's 0 here)
noise_level = 0;

%% Generate the slip steps

% make sure locations are on full pixel positions
locations = round(locations / psizeEBSD) * psizeEBSD;

slipSteps = struct;

for i=1:length(systems)
    slipSteps(i).System = systems(i); % slip system number
    slipSteps(i).traceLocation = locations(i,:); % trace location
    slipSteps(i).Magnitude = amplitudes(i); % total slip magnitude
    slipSteps(i).sSLocal = sSLocal(slipSteps(i).System);

    
    [slipSteps(i).U] = generateSlipStepField(ebsd,slipSteps(i).sSLocal,slipSteps(i).traceLocation,slipSteps(i).Magnitude);
end


% combine slipsteps and add data to dic

% if desired, add a rigid body motion to the disp field (does not influence
% the results)
% RBM = vector3d(0.8,0.2,0); % add RBM to all disp fields
RBM = vector3d(0,0,0); % add RBM to all disp fields
RBM = repmat(RBM,size(slipSteps(1).U));

% add noise to RBM
noise_x = RBM.x + randn(size(RBM.x)) * noise_level;
noise_y = RBM.y + randn(size(RBM.y)) * noise_level;

RBM.x = noise_x;
RBM.y = noise_y;


% combine RBM and slip steps into single displacement field
ebsd.prop.U = RBM;
for i=1:length(slipSteps)
    ebsd.prop.U = ebsd.prop.U + slipSteps(i).U;
end

% as a check, plot x component of displacement field
figure;
plot(ebsd,ebsd.prop.U.x)
mtexColorbar('title','X-Disp [um]')


%% slip ID
clear IDoptions
close all


% Batches for identification (name of batch and slip system numbers to use
% for ID). Can be used to easily try some things. Not currently capable of
% running multiple batches in a loop.

% idName = name used to same figures (if IDoptions.saveFig == 1)
% idName = name used to same figures
idName{1} = 'HCP_24'; NoSs_batch{1} = 1:24;
idName{2} = 'HCP_3_14_24'; NoSs_batch{2} = [3 14 24]; % the systems which are actually active...

% choose which "case/batch" to run
Ss_batch = [1];

%%%%%%% INPUT OPTIONS %%%%%%%
%%% for more information, see the SSLIP.m file, in which most options are
%%% explained a bit better

%%%
% set SSLIP method
% 1: constrained and minimized slip ID (As used predominantly in the SSLIP paper. 
% Contraint: || H^exp - H^their || < H_thresh. Minimzation of sum of absolute value of slip activities ) 

% 2: constrained slip ID (only solve problem based on the || H^exp - H^their || < H_thresh  constraint, no minimization)


% 3: single slip system ID (check for each pixel if a SINGLE system fits the activity. Used to "initialize" the SSLIP id om Figure 11 of the paper)
% Recommended as trial for uncertain/complex situations.
IDoptions.IDMethod = 1;
%%%

% filtering and coarse graining options (performed on displacement fields)
IDoptions.filterSize = 4;
IDoptions.coarsegrain = 2;

% define threshold under which the "residual" should be to stop solving
% (note that this is an "absolute" value, so in case of low strains (e.g.
% 0.02), the default value (0.01) is rather high.
IDoptions.threshResidual = 0.01;

% set minimum strain to perform SSLIP at a datapoint (i.e. skip
% pixels with low strain)
IDoptions.minEeff = 0.02;

% Positive constraint: choose to constrain the slip amplitudes to be positive or not. This
% only works well if the slip systems are "configured" to have a positive
% amplitude under a certain load (which is normally assured in the main
% script, assuming e.g. uniaxial tension).
% How to "reconfigure" the slip system under complex loads is T.B.D.
IDoptions.posConstr = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% PLOTTING OPTIONS %%%%%%%

% plot SSLIP results?
IDoptions.plotSSLIP = 1;

% max strain to plot, do not specify if no max is required (only required
% in some cases where something weird happens)
% IDoptions.maxE = 0.2;

% layout of plotting (leave empty if unknown) [rows columns]
IDoptions.layout = [4 6];

% plot slip traces and slip directions over the activity fields
IDoptions.plotTraces = 1;

% give this a name for saving
IDoptions.casename = ['HCP_Virt_',nameOri,'_',nameSlip,'_',idName{Ss_batch}];

% specify wether or not to plot defGrad fields and Eff shear strain fields,
% before performing SSLIP
IDoptions.plotDefGrad = 1;

% save the figures?
IDoptions.saveFig = 0;

% extra saved figure with different extension? (only used when saveFig==1)
IDoptions.saveExt = '.eps';

% change "size" of figure, e.g. to enlarge font etc (only used for
% individual plots)
IDoptions.sizeAdjust = 0.5;

% plot residual field after ID
IDoptions.plotResidual = 1;

% use logscale for plotting
IDoptions.logscale = 0;
% IDoptions.logmin = 0.01;

% assign colormap
IDoptions.cmap = viridis(512);

% put in a comment, will be used in filenames when saving figures
IDoptions.comment = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% END OF INPUT %%%%%%%%%




%%% set some variables

%%% choose which analysis to do, based on batch number (could be used to try
% a lot of things in a loop)
IDoptions.NoSs = NoSs_batch{Ss_batch};

% define the displacement field which should be used
U = ebsd.prop.U.x;
V = ebsd.prop.U.y;

%%% perform SSLIP analysis
[ebsdID,optOut] = SSLIP(ebsd,U,V,sSLocal,IDoptions);

%%% save the results in a matfile
save(optOut.plotname,'ebsdID','sSLocal','optOut');




% plot one activity field (system 3), just to demonstrate how data is
% structured
figure;
plot(ebsdID,ebsdID.prop.slipIDcor(find(optOut.NoSs==3),:))
mtexColorbar




% %%% potentially, for replotting:
% plotSSLIP(ebsdID.prop.slipIDcor,ebsdID.prop.residualEeff,ebsdID,sSLocal,optOut)










