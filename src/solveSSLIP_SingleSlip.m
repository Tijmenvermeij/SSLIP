function [slipIDcor, residualEeff] = solveSSLIP_SingleSlip(sS, Hxx, Hxy, Hyx, Hyy, cfg_solver)
% solveSSLIP_SingleSlip Checks for each pixel if a SINGLE slip system fits the activity.
% It loops over each single slip system, calls solveSSLIP_Constrained, and aggregates.

    NoSs = 1:length(sS);

    % initialize matrices
    nPixels = length(Hxx(:));
    slipIDcor = zeros(length(NoSs), nPixels);
    residualEeff = zeros(length(NoSs), nPixels);

    % don't use constraints at single slip ID (since it has no benefit
    % usually, and is much faster without it)
    if isfield(cfg_solver, 'posConstr') && cfg_solver.posConstr == 1
        warning('option "posConstr" changed to 0, since it has no added value for single slip ID, and is much slower')
        cfg_solver.posConstr = 0;
    end

    % loop over single slip systems, and perform SSLIP (without minimization)
    for j = 1:length(NoSs)
        [slipIDcor(j,:), residualEeff(j,:)] = solveSSLIP_Constrained(sS(j), Hxx, Hxy, Hyx, Hyy, cfg_solver);
    end
   
    % Calculate effective strain for thresholding if needed
    Eeff = calcEffectiveE(Hxx, Hxy, Hyx, Hyy);
    

    
    % use residual to "filter" the fields, only leaving activities with low residual
    if ~isfield(cfg_solver,'threshResidualFraction')
        % "clean" slipID by using threshold on residual, based on single systems
        goodData = residualEeff < cfg_solver.threshResidual;
    else
        % "clean" slipID by using threshold on residual, based on factor of total effective shear
        threshResidual = Eeff(:) * cfg_solver.threshResidualFraction;
        threshResidual(threshResidual < cfg_solver.threshResidual) = cfg_solver.threshResidual;
        goodData = residualEeff < threshResidual';
    end

    % this option leaves only 1 slip system active at each datapoint, which has the lowest residual.
    if isfield(cfg_solver,'singleSlipPerPixel') && cfg_solver.singleSlipPerPixel
        [~, minThreshInd] = min(residualEeff, [], 1);
        slipIDcor_new = zeros(size(slipIDcor));
        for j=1:nPixels
            slipIDcor_new(minThreshInd(j), j) = slipIDcor(minThreshInd(j), j);
        end
        slipIDcor = slipIDcor_new;
    end
    
    % make all slip activities with residual above threshold zero.
    slipIDcor(~goodData) = 0;
end
