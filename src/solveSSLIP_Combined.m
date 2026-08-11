function [slipIDcor,residualEeff] = solveSSLIP_Combined(sS,Hxx,Hxy,Hyx,Hyy,options)
    
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
    options = struct;
end

%% set some options if not given



%% start setting up the ID
% get theoretical disp grad tensor components from the slip systems:
Hslip = sS.deformationTensor.matrix;

if isfield(options, 'enableRotation') && options.enableRotation
    Hslip(:,:,end+1) = [0 -1 0
                1 0 0
                0 0 0];
    if isfield(options, 'posConstr') && options.posConstr
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
if options.normalizeInplane
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
numAnalysis = sum(Eeff>options.minEeff);

if isfield(options, 'ADMM') && options.ADMM
    % =========================================================================
    % VECTORIZED ADMM SOLVER (Modern Implementation)
    % =========================================================================
    % PROS:
    % - Extremely fast. Converts 100,000+ independent optimization problems into 
    %   a single batched matrix operation.
    % - Avoids MATLAB's heavy 'coneprog' overhead.
    % - Exact same mathematical objective (Basis Pursuit Denoising / L1-Min).
    % CONS / LIMITATIONS:
    % - It is an iterative numerical approximation (1st-order method), so results 
    %   may differ from 'coneprog' (Interior Point method) by very small margins.
    % - Requires tuning of hyperparameters (rho1, rho2, iterations) for specific 
    %   noise profiles, though defaults work well for general DIC.
    % =========================================================================
    
    fprintf('  [SSLIP] Solving %d points using Vectorized ADMM (%d iters)...\n', numAnalysis, options.ADMM_iters);
    
    % 1. Prepare global thresholds
    if isfield(options,'threshResidualFraction')
        threshResidual = Eeff(:) * options.threshResidualFraction;
        threshResidual(threshResidual < options.threshResidual) = options.threshResidual;
    else
        threshResidual = ones(length(Eeff(:)), 1) * options.threshResidual;
    end
    
    % 2. Identify Valid Pixels
    % Reshape HExp to [4 x N_pixels]
    B_all = reshape(HExp, 4, []);
    valid_idx = find(Eeff(:) >= options.minEeff & ~any(isnan(B_all), 1)');
    num_valid = length(valid_idx);
    
    if num_valid > 0
        B_valid = B_all(:, valid_idx);
        eps_valid = threshResidual(valid_idx)'; % [1 x num_valid]
        
        % Hyperparameters
        % Activities are typically O(0.01). kappa = 1/rho1 is the soft-threshold cutoff.
        % Therefore, rho1 must be very large (e.g., 5000) so kappa = 0.0002.
        rho1 = 5000;
        rho2 = rho1 * 5;
        
        % 3. Precompute Matrix Inverse for exact X-update
        % Because A is constant across all pixels, we only invert this ONCE!
        M = inv(rho1 * eye(N) + rho2 * (A' * A));
        At = A';
        
        % =========================================================================
        % DEBUG DATA DUMP
        % =========================================================================
        if isfield(options, 'ADMM') && options.ADMM == 1
            save('C:\Users\phkr\MTEX_workEnv\debug_ADMM_matrices.mat', 'A', 'HExp', 'Eeff', 'options');
            fprintf('  [DEBUG] Saved A, HExp, Eeff, options to debug_ADMM_matrices.mat\n');
        end
        
        % Initialize outputsADMM variables
        x  = zeros(N, num_valid);
        z1 = zeros(N, num_valid);
        u1 = zeros(N, num_valid);
        z2 = zeros(4, num_valid);
        u2 = zeros(4, num_valid);
        
        kappa = 1 / rho1;
        
        % 4. ADMM Iteration Loop
        for k = 1:options.ADMM_iters
            % X-Update (Matrix Multiplications)
            rhs = rho1 * (z1 - u1) + rho2 * At * (z2 - u2);
            x = M * rhs;
            
            % Z1-Update (Soft Thresholding / Positivity)
            v1 = x + u1;
            if options.posConstr
                z1 = max(v1 - kappa, 0);
            else
                z1 = sign(v1) .* max(abs(v1) - kappa, 0);
            end
            
            % Z2-Update (Projection onto L2 Residual Ball)
            v2 = A * x + u2;
            diff_v = v2 - B_valid;
            norms = sqrt(sum(diff_v.^2, 1));
            
            proj_mask = norms > eps_valid;
            z2 = v2;
            if any(proj_mask)
                scale = eps_valid(proj_mask) ./ norms(proj_mask);
                z2(:, proj_mask) = B_valid(:, proj_mask) + diff_v(:, proj_mask) .* scale;
            end
            
            % Dual Updates
            u1 = u1 + x - z1;
            u2 = u2 + A * x - z2;
        end
        
        % Extract solved activities and calculate exact residuals
        gamma(:, valid_idx) = z1;
        res(valid_idx) = sqrt(sum((A * z1 - B_valid).^2, 1));
    end
    
    % Clean up NaNs
    nan_mask = any(isnan(B_all), 1);
    gamma(:, nan_mask) = NaN;
    res(nan_mask) = NaN;

else
    % =========================================================================
    % 🐌 LEGACY LOOP-BASED SOLVER (parfor + coneprog)
    % =========================================================================
    % create PARFOR waitmessage to monitor progress
    WaitMessage = parfor_wait(numAnalysis,'ReportInterval',ceil(numAnalysis/20));
    
    %%% loop over all points
    parfor i=1:length(Hxx(:))
    
        % define threshold based on input
        if isfield(options,'threshResidualFraction')
            threshResidual = Eeff(i) * options.threshResidualFraction;
            threshResidual(threshResidual<options.threshResidual) = options.threshResidual;
        else
            threshResidual = options.threshResidual;
        end
    
        % do some checks to see if ID needs to be performed
        
        % skip NaNs
        if any(isnan(HExp(:,:,i)),'all')
            gamma(:,i) = NaN;
            res(i) = NaN;
        % skip points with low Eeff, assign 0 activity and residual
        elseif Eeff(i) < options.minEeff
            gamma(:,i) = 0;
            res(i) = 0;
        % perform ID
        else
            % extract exp H, for point i
            HExpi = HExp(:,:,i);
            HExpi = HExpi';
            HExpi = HExpi(:);
            
            % minimization, constraining activities to be positive
            if options.posConstr
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
if isfield(options, 'enableRotation') && options.enableRotation
    if isfield(options, 'posConstr') && options.posConstr
        slipID_rot = slipIDcor(end-1,:) - slipIDcor(end,:);
        slipIDcor = [slipIDcor(1:end-2,:); slipID_rot];
    end
end

residualEeff = res;





end



