function plotSSLIP(slipIDcor,residualEeff,ebsdID,sSLocal,opt)
% Plot activity and optional residual fields using the supplied arrays.
% Activity-plot separation adapted from Philipp (PhilKro), PR #2:
% https://github.com/Tijmenvermeij/SSLIP/pull/2 (commit d3ee1a5).

% Keep the explicit arrays authoritative when the stored fields differ.
ebsdID.prop.slipIDcor = slipIDcor;
caxisMinMax = plotSSLIP_SlipActivities(ebsdID,sSLocal,opt);

if opt.plotResidual && (opt.IDMethod == 1 || opt.IDMethod == 2)
    ebsdID.prop.residualEeff = residualEeff;
    opt.caxisMinMax = caxisMinMax;
    plotSSLIP_Residual(ebsdID,opt);
end
end
