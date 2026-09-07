function [prepData, ebsdID] = preprocessSSLIP(ebsd, DeformationData, cfg_prep)
% Filter and coarse-grain displacements, then calculate their gradients.
% DeformationData contains U,V or Hxx,Hxy,Hyx,Hyy, shaped and ordered as ebsd.
% cfg_prep contains filterSize and coarsegrain (defaults supplied by SSLIP).
% prepData contains Hxx,Hxy,Hyx,Hyy,Eeff and, when supplied, processed U,V.
% Direct gradients are ready for fitting: filterSize=0 and coarsegrain=1.
% Their zeros are physical values; use NaN for missing data. Displacements
% are not reconstructed or filled in when only gradients are supplied.
%
% Preprocessing separation adapted from Philipp (PhilKro), PR #2:
% https://github.com/Tijmenvermeij/SSLIP/pull/2
% Original commit: d3ee1a54eb9f8248565e0063060e749f9f56b9ce.
% Retains Tijmen Vermeij's existing grid mapping and displacement processing.

if ~isstruct(DeformationData) || ~isscalar(DeformationData)
    error('SSLIP:InvalidDeformationData','DeformationData must be a scalar structure.');
end
gradientFields = {'Hxx','Hxy','Hyx','Hyy'};
hasGradientFields = isfield(DeformationData,gradientFields);
if any(hasGradientFields)
    if any(isfield(DeformationData,{'U','V'}))
        error('SSLIP:AmbiguousDeformationData','Supply displacements or gradients, not both.');
    end
    if ~all(hasGradientFields)
        error('SSLIP:MissingGradients','Supply all four fields: Hxx, Hxy, Hyx, Hyy.');
    end
    if cfg_prep.filterSize ~= 0 || cfg_prep.coarsegrain ~= 1
        error('SSLIP:GradientPreprocessing', ...
            'Supplied gradients must be ready for fitting: set filterSize=0 and coarsegrain=1.');
    end

    % A clean copy retains the grid/orientation metadata without stale fits.
    % Attach each component before gridify so values follow their coordinates.
    ebsdID = ebsd;
    ebsdID.prop = struct();
    for k = 1:numel(gradientFields)
        name = gradientFields{k};
        value = DeformationData.(name);
        if ~isfloat(value) || ~isreal(value) || any(isinf(value),'all')
            error('SSLIP:InvalidGradients', ...
                'Gradient components must be real floating-point arrays; use NaN for missing values.');
        end
        if ~isequal(size(value),size(ebsd))
            error('SSLIP:GradientSizeMismatch','Each gradient component must have the same size as ebsd.');
        end
        ebsdID.prop.(name) = double(value);
    end
    ebsdID = ebsdID.gridify;
    for k = 1:numel(gradientFields)
        name = gradientFields{k};
        prepData.(name) = ebsdID.prop.(name);
    end
    prepData.Eeff = calcEffectiveE(prepData.Hxx,prepData.Hxy,prepData.Hyx,prepData.Hyy);
    ebsdID.prop.Eeff = prepData.Eeff;
    return;
end

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
