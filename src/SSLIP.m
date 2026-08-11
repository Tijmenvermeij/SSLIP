function [ebsdID, prepData, optOut] = SSLIP(ebsd, DeformationData, sSLocal, cfg)
% SSLIP orchestrator for the SSLIP optimization process.
%
% Inputs:
%   ebsd            - MTEX EBSD object (used for grid / spatial context)
%   DeformationData - struct containing either (U, V) or (Hxx, Hxy, Hyx, Hyy)
%   sSLocal         - List of slip systems rotated to local crystal orientation
%   cfg             - Nested configuration struct with the following optional fields:
%                     cfg.preprocess.filterSize  (default: 1) - Gaussian filter std dev
%                     cfg.preprocess.coarsegrain (default: 1) - Factor for coarse-graining grid
%                     cfg.solver.IDMethod        (default: 1) - 1: Combined, 2: Constrained, 3: Single Slip
%                     cfg.solver.minEeff         (default: 0.01) - Minimum effective strain to process a pixel
%                     cfg.solver.threshResidual  (default: 0.01) - Maximum allowed residual error
%                     cfg.solver.threshResidualFraction - (Optional) Dynamic residual threshold factor based on Eeff
%                     cfg.solver.posConstr       (default: 0) - If 1, constrains slip activities to positive values
%                     cfg.solver.enableRotation  (default: 0) - If 1, includes pseudo-slip system for rotation
%                     cfg.solver.normalizeInplane (default: 0) - If 1, normalizes theoretical displacement gradient
%                     cfg.solver.singleSlipPerPixel (default: 0) - (Method 3 only) Forces only 1 active system per pixel
%                     cfg.solver.NoSs            (default: all) - Indices of slip systems to evaluate
%
% Outputs:
%   ebsdID          - The original MTEX EBSD object, enriched with new properties
%                     in 'ebsdID.prop' containing the solver results:
%                     .slipIDcor    : The solved slip activities (array)
%                     .residualEeff : The residual error map of the fit
%                     .rotationIDcor: (Optional) Solved rotation pseudo-slip activities
%   prepData        - Struct containing the final preprocessed displacement/gradient 
%                     fields used by the solver (e.g. .Hxx, .Hxy, .Hyx, .Hyy)
%   optOut          - The exact 'cfg.solver' options that were actually executed, 
%                     useful to pass down to downstream plotting functions.

    if nargin < 4, cfg = struct; end
    if ~isfield(cfg, 'preprocess'), cfg.preprocess = struct; end
    if ~isfield(cfg, 'solver'), cfg.solver = struct; end
    
    % --- PREPROCESS DEFAULTS ---
    if ~isfield(cfg.preprocess, 'filterSize'), cfg.preprocess.filterSize = 1; end                   % Gaussian blur size (pixels) applied to displacement field to reduce noise. Used in preprocessSSLIP.m
    if ~isfield(cfg.preprocess, 'coarsegrain'), cfg.preprocess.coarsegrain = 1; end                 % Coarse graining factor (1=none, 2=2x2) for performance. Used in preprocessSSLIP.m

    % --- SOLVER DEFAULTS ---
    if ~isfield(cfg.solver, 'IDMethod'), cfg.solver.IDMethod = 1; end                               % 1=constrained L1-minimized, 2=constrained only, 3=single slip. Used to route to solveSSLIP_*.m
    if ~isfield(cfg.solver, 'minEeff'), cfg.solver.minEeff = 0.01; end                              % Minimum effective strain required to evaluate a pixel (skips low strain). Used in solveSSLIP_*.m
    if ~isfield(cfg.solver, 'threshResidual'), cfg.solver.threshResidual = 0.01; end                % L2-norm threshold for the optimization residual (||H_exp - H_theor|| < H_thresh). Used in solveSSLIP_*.m
    if ~isfield(cfg.solver, 'threshResidualFraction'), cfg.solver.threshResidualFraction = 0; end   % (Optional) Dynamic residual threshold (Eeff * fraction), bounded below by threshResidual. Deactivate by omitting or setting to 0.
    if ~isfield(cfg.solver, 'posConstr'), cfg.solver.posConstr = 0; end                             % 1=constrain slip amplitudes to be positive. Used in solveSSLIP_*.m
    if ~isfield(cfg.solver, 'enableRotation'), cfg.solver.enableRotation = 0; end                   % 1=add rigid body rotation as a pseudo-slip system. Used in solveSSLIP_*.m
    if ~isfield(cfg.solver, 'normalizeInplane'), cfg.solver.normalizeInplane = 0; end               % 1=normalize the 2D projected theoretical slip deformation gradients. Used in solveSSLIP_*.m
    if ~isfield(cfg.solver, 'singleSlipPerPixel'), cfg.solver.singleSlipPerPixel = 0; end           % 1=only allow a single active slip system per pixel (forces ultimate sparsity). Used in solveSSLIP_Combined.m
    % Ensure sSLocal is a column vector
    if length(sSLocal) > 1 && size(sSLocal, 2) ~= 1
        sSLocal = transpose(sSLocal);
    end
    
    if ~isfield(cfg.solver, 'NoSs'), cfg.solver.NoSs = 1:length(sSLocal); end               % Indices of slip systems to include in the solver. Used in runSSLIP.m
    
    % Apply posConstr warning
    if cfg.solver.posConstr == 1
        warning('Make sure that the slip systems are configured to have a positive SF amplitude in the assumed loading. This can be done by running: sf_vals = sSLocal.SchmidFactor(stressState); sf_signs = sign(sf_vals); sf_signs(sf_signs == 0) = 1; sSLocal.b = sf_signs .* sSLocal.b;');
    end
    
    % Subselect slip systems based on NoSs
    sSLocal = sSLocal(cfg.solver.NoSs);

    % --- 1. Preprocessing Pipeline ---
    fprintf('  [SSLIP Pipeline] 1. Preprocessing data...\n');
    [prepData, ebsdID] = preprocessSSLIP(ebsd, DeformationData, cfg.preprocess);

    % --- 2. Mathematical Solver Pipeline ---
    fprintf('  [SSLIP Pipeline] 2. Solving using IDMethod %d...\n', cfg.solver.IDMethod);
    
    switch cfg.solver.IDMethod
        case 1
            [slipIDcor, residualEeff] = solveSSLIP_Combined(sSLocal, prepData.Hxx, prepData.Hxy, prepData.Hyx, prepData.Hyy, cfg.solver);
        case 2
            [slipIDcor, residualEeff] = solveSSLIP_Constrained(sSLocal, prepData.Hxx, prepData.Hxy, prepData.Hyx, prepData.Hyy, cfg.solver);
        case 3
            [slipIDcor, residualEeff] = solveSSLIP_SingleSlip(sSLocal, prepData.Hxx, prepData.Hxy, prepData.Hyx, prepData.Hyy, cfg.solver);
        otherwise
            error('IDMethod %d is unknown. Expected 1, 2, or 3.', cfg.solver.IDMethod);
    end

    % --- 3. EBSD Object Mapping Pipeline ---
    fprintf('  [SSLIP Pipeline] 3. Mapping results back to MTEX EBSD object...\n');
    
    ebsdID.prop.residualEeff = residualEeff';
    
    % Safely extract the rotation pseudo-system if it was calculated
    if isfield(cfg.solver, 'enableRotation') && cfg.solver.enableRotation
        ebsdID.prop.rotationIDcor = slipIDcor(end, :)';
        slipIDcor = slipIDcor(1:end-1, :);
    end
    
    ebsdID.prop.slipIDcor = slipIDcor';
    optOut = cfg.solver;

end
