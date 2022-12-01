%% Script to perform slip system identification, using the SSLIP method, on SEM-DIC/EBSD data of a Ni based super alloy.
% This script replicates the identification,
% as shown in Figures 5&6 in the paper "T. Vermeij et al., Automated identification of slip system activity fields from digital image correlation data, Acta Mater. 243, 2022
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
%
% The data was retrieved from Zenodo (https://doi.org/10.5281/zenodo.3691903), after which we aligned it using our Nanomechanical testing framework (https://doi.org/10.1007/s11340-022-00884-0) 
% The Ni-based superalloy RR1000 was deformed under uniaxial tension to a global strain of ∼ 0.02, data was acquired by Harte et al.


%% Initiallize Mtex etc

clear
close all
addpath('./src')

% load aligned data
load('./data/NiSuperAloy_Aligned.mat');

% set MTex preferences
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoplane');

%% plot some things to show the data

% plot effective strain map with grain boundary overlay to check alignment
figure;
plot(newEBSD,newEBSD.prop.Eeff)
colormap viridis
caxis([0.01 0.1])
set(gca,'colorscale','log')

hold on
plot(grains.boundary,'linewidth',2,'linecolor','r')
hold off

% plot grains with grain ID numbers
figure;
plot(grains,grains.meanOrientation)
hold on
text(grains,grains.id)
hold off


%% choose a grain and define slip systems
grainId = 141; %grain used in Figure 5 in SSLIP paper

% define FCC slip systems
sS = slipSystem.fcc(CS{2});
sS = sS.symmetrise('antipodal');
% reorder systems for later plotting of 3 rows and 4 cols, with each column
% having the same slip plane
sS = sS([1 4 7 10 2 5 8 11 3 6 9 12]);

% get local slip systems
sSLocal = grains(grainId).meanOrientation * sS;

% extract ebsd data (which included the disp data) for the choosen grain
ebsd = newEBSD(grains(grainId));
ebsd = ebsd.gridify;

%% slip ID
clear IDoptions
close all


% Batches for identification (name of batch and slip system numbers to use
% for ID). Can be used to easily try some things. Not currently capable of
% running multiple batches in a loop.

% idName = name used to same figures (if IDoptions.saveFig == 1)
% idName = name used to same figures
idName{1} = 'FCC_all'; NoSs_batch{1} = 1:12;
idName{2} = 'FCC_1_8'; NoSs_batch{2} = [1 8]; % the systems which are actually active...

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
IDoptions.filterSize = 1;
IDoptions.coarsegrain = 3;

% define threshold under which the "residual" should be to stop solving
% (note that this is an "absolute" value, so in case of low strains (e.g.
% 0.02), the default value (0.01) is rather high.
IDoptions.threshResidual = 0.005;

% set minimum strain to perform SSLIP at a datapoint (i.e. skip
% pixels with low strain)
IDoptions.minEeff = 0.01;

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
IDoptions.layout = [3 4];

% plot slip traces and slip directions over the activity fields
IDoptions.plotTraces = 1;

% give this a name for saving
IDoptions.casename = ['Ni_grain',grainId,'_'];

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
IDoptions.logscale = 1;
IDoptions.logmin = 0.01;

% assign colormap
IDoptions.cmap = viridis(512);

% put in a comment, will be used in filenames when saving figures
IDoptions.comment = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% END OF INPUT %%%%%%%%%




%%% set some variables

%%% choose which analysis to do, based on batch number (could be used to try
% a lot of things in a loop)
IDoptions.casename = idName{Ss_batch}; 
IDoptions.NoSs = NoSs_batch{Ss_batch};

% define the displacement field which should be used
U = ebsd.prop.U;
V = ebsd.prop.V;

%%% perform SSLIP analysis
[ebsdID,optOut] = SSLIP(ebsd,U,V,sSLocal,IDoptions);

%%% save the results in a matfile
save(optOut.plotname,'ebsdID','sSLocal','optOut');




% plot one activity field (system 1), just to demonstrate how data is
% structured
figure;
plot(ebsdID,ebsdID.prop.slipIDcor(find(optOut.NoSs==1),:))
mtexColorbar




% %%% potentially, for replotting:
% plotSSLIP(ebsdID.prop.slipIDcor,ebsdID.prop.residualEeff,ebsdID,sSLocal,optOut)










