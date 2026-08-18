function [slipIDcor,residualEeff] = solveSSLIP_Combined(sS,Hxx,Hxy,Hyx,Hyy,cfg_solver)
    
%%% Function to perform slip id analysis using def grad tensor component
%%% solving. Uses constraints and minimized the sum of absolute values of
%%% slip amplitudes

% This function contains the SSLIP method as proposed in the paper 
% "T. Vermeij et al., Automated identification of slip system activity
% fields from digital image correlation data, Acta Mater. 243, 2022"
% DOI: https://doi.org/10.1016/j.actamat.2022.118502
% Please consider citing this paper when you use this code.

% Date: 30-11-2022
% the latest version of this code can be found on
% www.github.com/TijmenVermeij/SSLIP

%%% Tijmen Vermeij / TUe / t.vermeij@tue.nl


%%
if nargin < 6
    cfg_solver = struct;
end

%% set some options if not given



%% start setting up the ID
% get theoretical disp grad tensor components from the slip systems:
Hslip = sS.deformationTensor.matrix;

if isfield(cfg_solver, 'enableRotation') && cfg_solver.enableRotation
    Hslip(:,:,end+1) = [0 -1 0
                1 0 0
                0 0 0];
    if isfield(cfg_solver, 'posConstr') && cfg_solver.posConstr
        Hslip(:,:,end+1) = [0 1 0
                    -1 0 0
                    0 0 0];
    end
end

% extract the theoretical 2D components
Hslip11 = reshape( Hslip(1,1,:) , 1, []);
Hslip12 = reshape( Hslip(1,2,:) , 1, []);
Hslip21 = reshape( Hslip(2,1,:) , 1, []);
Hslip22 = reshape( Hslip(2,2,:) , 1, []);

% put the theoretical components in matrix form for A*x = B type solving  
% (rows: different components, columns: different slip systems)
A = [Hslip11
    Hslip12
    Hslip21
    Hslip22];

% if necessary, normalize inplane
if cfg_solver.normalizeInplane
    A = A./sqrt(sum(A.^2,1));
end

% put the experimental disp grad components in a 3D-matrix (2x2xn), with n
% number of datapoints
HExp = zeros(2,2,length(Hxx(:)));
HExp(1,1,:) = Hxx(:);
HExp(1,2,:) = Hxy(:);
HExp(2,1,:) = Hyx(:);
HExp(2,2,:) = Hyy(:);

% set options for coneprog
coneprogoptions = optimoptions('coneprog','Display','none');

% initialize some variables

gamma = zeros(size(Hslip,3),size(HExp,3)); 
res = zeros(size(HExp,3),1);
fobj = zeros(size(HExp,3),1);
flag = zeros(size(HExp,3),1);
N = size(Hslip,3); % no of considered slip systems

% calc eff strains, which will be used to determine on which points to
% perform ID
Eeff = calcEffectiveE(Hxx(:),Hxy(:),Hyx(:),Hyy(:));

% calc number of points on which ID will be performed
numAnalysis = sum(Eeff>cfg_solver.minEeff);

% create PARFOR waitmessage to monitor progress
WaitMessage = parfor_wait(numAnalysis,'ReportInterval',ceil(numAnalysis/20));

%%% loop over all points
parfor i=1:length(Hxx(:))

    % define threshold based on input
    if isfield(cfg_solver,'threshResidualFraction')
        threshResidual = Eeff(i) * cfg_solver.threshResidualFraction;
        threshResidual(threshResidual<cfg_solver.threshResidual) = cfg_solver.threshResidual;
    else
        threshResidual = cfg_solver.threshResidual;
    end

    % do some checks to see if ID needs to be performed
    
    % skip NaNs
    if any(isnan(HExp(:,:,i)),'all')
        gamma(:,i) = NaN;
        res(i) = NaN;
    % skip points with low Eeff, assign 0 activity and residual
    elseif Eeff(i) < cfg_solver.minEeff
        gamma(:,i) = 0;
        res(i) = 0;
    % perform ID
    else
        % extract exp H, for point i
        HExpi = HExp(:,:,i);
        HExpi = HExpi';
        HExpi = HExpi(:);
        
        % minimization, constraining activities to be positive
        if cfg_solver.posConstr
            % define the lower bound of the activity (i.e. the postive
            % constraint)
            gamma0 = zeros(N,1);
            
            % define the constraints: || H^exp - H^their || < H_thresh
            % (see documentation of coneprog for clarifications)
            socContraints = secondordercone(A,HExpi,gamma0,-1*threshResidual);
            
            % run coneprog (see documentation of coneprog for clarifications)
            [x,f,flag] = coneprog(ones(size(gamma0)),socContraints,[],[],[],[],gamma0,[],coneprogoptions);
        else
            % minimization, not constraining activities to be positive (by
            % minimizing sum of absolute values, in an indirect way, by
            % splitting the slip systems in their positive and negative
            % "parts"
            
            % each (positive or negative) slip system will still have a
            % positive constraint:
            gamma0 = zeros(2*N,1);
            
            % extend the A matrix with the negative slip systems
            A2 = [A -1*A];
            
            % define the constraints: || H^exp - H^their || < H_thresh
            % (see documentation of coneprog for clarifications)
            socContraints = secondordercone(A2,HExpi,gamma0,-1*threshResidual);
            
            % run coneprog (see documentation of coneprog for clarifications)
            [x,f,flag] = coneprog(ones(size(gamma0)),socContraints,[],[],[],[],gamma0,[],coneprogoptions);
            
            % recombine the slip activities, by substracting the "negative
            % slip system amplitudes", from the amplitudes of their
            % positive counterpart
            if flag == 1
                x = x(1:N) - x(N+1:end);
            end
        end
        
        % check the solution. if no good flag, put in a NaN
        if flag == 1
            gamma(:,i) = x;
            fobj(i) = f;
        else
            gamma(:,i) = NaN;
            fobj(i) = NaN;
        end
        
        % calculate residual 
        res(i) = norm(A*gamma(:,i)-HExpi);

        WaitMessage.Send;
    end
    
end
WaitMessage.Destroy
end

slipIDcor = gamma;
if isfield(cfg_solver, 'enableRotation') && cfg_solver.enableRotation
    if isfield(cfg_solver, 'posConstr') && cfg_solver.posConstr
        slipID_rot = slipIDcor(end-1,:) - slipIDcor(end,:);
        slipIDcor = [slipIDcor(1:end-2,:); slipID_rot];
    end
end

residualEeff = res;





end



