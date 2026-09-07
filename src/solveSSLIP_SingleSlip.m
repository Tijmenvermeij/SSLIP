function [slipIDcor, residualEeff, cfg_solver] = solveSSLIP_SingleSlip(sS, Hxx, Hxy, Hyx, Hyy, cfg_solver)
% Check each supplied slip system independently, then filter by residual.
% sS is already selected and ordered by SSLIP; cfg_solver.NoSs is not reapplied.
% Activities and residuals are systems-by-pixels. The optional third output
% reports the executed options, including the existing posConstr = 0 override.
% SSLIP supplies the defaults. This method retains its residual-based filter;
% it does not apply the combined solver's minEeff cutoff or rotation basis.
%
% Adapted from Philipp (PhilKro), PR #2:
% https://github.com/Tijmenvermeij/SSLIP/pull/2
% Original commit: d3ee1a54eb9f8248565e0063060e749f9f56b9ce.
% Retains the existing SSLIPConstr function and single-slip fit conventions.

if isfield(cfg_solver,'enableRotation') && cfg_solver.enableRotation
    error('SSLIP:RotationRequiresMethod1', ...
        'enableRotation is supported only with IDMethod = 1.');
end

nSystems = length(sS);
nPixels = numel(Hxx);
slipIDcor = zeros(nSystems,nPixels);
residualEeff = zeros(nSystems,nPixels);

% Preserve the fast unconstrained single-system fit and its reported options.
if isfield(cfg_solver,'posConstr') && cfg_solver.posConstr == 1
    warning('option "posConstr" changed to 0, since it has no added value for single slip ID, and is much slower')
    cfg_solver.posConstr = 0;
end

for j = 1:nSystems
    [slipIDcor(j,:),residualEeff(j,:)] = SSLIPConstr(sS(j),Hxx,Hxy,Hyx,Hyy,cfg_solver);
end

if ~isfield(cfg_solver,'threshResidualFraction')
    goodData = residualEeff < cfg_solver.threshResidual;
else
    Eeff = calcEffectiveE(Hxx,Hxy,Hyx,Hyy);
    threshResidual = Eeff(:) * cfg_solver.threshResidualFraction;
    threshResidual(threshResidual < cfg_solver.threshResidual) = cfg_solver.threshResidual;
    goodData = residualEeff < threshResidual';
end

% Keep the original amplitude of the system with the lowest residual.
if isfield(cfg_solver,'singleSlipPerPixel') && cfg_solver.singleSlipPerPixel
    [~,minThreshInd] = min(residualEeff,[],1);
    singleSlipID = zeros(size(slipIDcor));
    for j = 1:nPixels
        singleSlipID(minThreshInd(j),j) = slipIDcor(minThreshInd(j),j);
    end
    slipIDcor = singleSlipID;
end

% Rejected fits retain zero activity and their original per-system residual.
slipIDcor(~goodData) = 0;
end
