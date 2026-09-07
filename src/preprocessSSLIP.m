function [prepData, ebsdID] = preprocessSSLIP(ebsd, DeformationData, cfg_prep)
% Filter and coarse-grain displacements, then calculate their gradients.
% DeformationData contains U and V in the same shape and order as ebsd.
% cfg_prep contains filterSize and coarsegrain (defaults supplied by SSLIP).
% prepData contains U, V, Hxx, Hxy, Hyx, Hyy, and Eeff on the returned grid.
% This first integration supports displacement input only.
%
% Preprocessing separation adapted from Philipp (PhilKro), PR #2:
% https://github.com/Tijmenvermeij/SSLIP/pull/2
% Original commit: d3ee1a54eb9f8248565e0063060e749f9f56b9ce.
% Retains Tijmen Vermeij's existing grid mapping and displacement processing.

if ~isfield(DeformationData,'U') || ~isfield(DeformationData,'V')
    error('SSLIP:MissingDisplacements','DeformationData must contain U and V.');
end
U = DeformationData.U;
V = DeformationData.V;
if ~isequal(size(ebsd),size(U),size(V))
    error('SSLIP:DisplacementSizeMismatch','ebsd data must have the same size as U and V.');
end

% Carry displacements through gridify so coordinates and values stay aligned.
ebsd.prop.U = vector3d(U,V,zeros(size(U)));
ebsd = ebsd.gridify;
X = ebsd.x;
Y = ebsd.y;
U = ebsd.prop.U.x;
V = ebsd.prop.U.y;

if cfg_prep.filterSize ~= 0
    filt_opts.filt_std = cfg_prep.filterSize;
    % Preserve the existing convention for non-correlated DIC points:
    % zero displacement components become NaN only when filtering is enabled.
    U(U == 0) = NaN;
    V(V == 0) = NaN;
    data = filterDisplacements(U,V,filt_opts);
else
    data.U = U;
    data.V = V;
end

crs = coarsegrainDisp(data.U,data.V,X(1,:),Y(:,1)',cfg_prep.coarsegrain);
ebsdID = dummyEBSDSimple(ebsd.orientations(1),crs.X,crs.Y);
prepData.U = crs.f;
prepData.V = crs.g;
[Hxx,Hxy] = gradient(prepData.U,crs.pixelsize(1),crs.pixelsize(2));
[Hyx,Hyy] = gradient(prepData.V,crs.pixelsize(1),crs.pixelsize(2));
prepData.Hxx = Hxx;
prepData.Hxy = Hxy;
prepData.Hyx = Hyx;
prepData.Hyy = Hyy;
prepData.Eeff = calcEffectiveE(Hxx,Hxy,Hyx,Hyy);

ebsdID.prop.Eeff = prepData.Eeff;
ebsdID.prop.U = prepData.U;
ebsdID.prop.V = prepData.V;
ebsdID.prop.Hxx = Hxx;
ebsdID.prop.Hxy = Hxy;
ebsdID.prop.Hyx = Hyx;
ebsdID.prop.Hyy = Hyy;
end
