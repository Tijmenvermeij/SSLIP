function [prepData, ebsdID] = preprocessSSLIP(ebsd, DeformationData, cfg_prep)
% preprocessSSLIP Handles data filtering, coarse-graining, and gradient extraction 
% before running the SSLIP optimizer.

if nargin < 3, cfg_prep = struct; end

X = ebsd.x;
Y = ebsd.y;

hasGradients = isfield(DeformationData, 'Hxx') && isfield(DeformationData, 'Hxy') && isfield(DeformationData, 'Hyx') && isfield(DeformationData, 'Hyy');
hasDisps = isfield(DeformationData, 'U') && isfield(DeformationData, 'V');

if ~hasGradients && ~hasDisps
    error('DeformationData must contain either U and V, or Hxx, Hxy, Hyx, Hyy');
end

% Initialize dummy structure for spatial coordinates
crs = struct();

% --- PROCESS DISPLACEMENTS ---
if hasDisps
    U = DeformationData.U;
    V = DeformationData.V;
    
    if cfg_prep.filterSize ~= 0
        filt_opts.filt_std = cfg_prep.filterSize;
        U(U == 0) = NaN;
        V(V == 0) = NaN;
        data = filterDisplacements(U, V, filt_opts);
    else
        data.U = U;
        data.V = V;
    end
    
    crs = coarsegrainDisp(data.U, data.V, X(1,:), Y(:,1)', cfg_prep.coarsegrain);
    
    if ~hasGradients
        % Calculate gradients numerically
        [Hxx, Hxy] = gradient(crs.f, crs.pixelsize(1), crs.pixelsize(2));
        [Hyx, Hyy] = gradient(crs.g, crs.pixelsize(1), crs.pixelsize(2));
    end
    
    prepData.U = crs.f;
    prepData.V = crs.g;
end

% --- PROCESS GRADIENTS ---
if hasGradients
    Hxx = DeformationData.Hxx; Hxy = DeformationData.Hxy;
    Hyx = DeformationData.Hyx; Hyy = DeformationData.Hyy;
    
    % If raw grid gradients were passed (same size as original EBSD)
    if isequal(size(Hxx), size(ebsd.x))
        if cfg_prep.filterSize ~= 0
            filt_opts.filt_std = cfg_prep.filterSize;
            fH1 = filterDisplacements(Hxx, Hxy, filt_opts);
            fH2 = filterDisplacements(Hyx, Hyy, filt_opts);
            Hxx = fH1.U; Hxy = fH1.V;
            Hyx = fH2.U; Hyy = fH2.V;
        end
        if cfg_prep.coarsegrain > 1
            crsH1 = coarsegrainDisp(Hxx, Hxy, X(1,:), Y(:,1)', cfg_prep.coarsegrain);
            crsH2 = coarsegrainDisp(Hyx, Hyy, X(1,:), Y(:,1)', cfg_prep.coarsegrain);
            Hxx = crsH1.f; Hxy = crsH1.g;
            Hyx = crsH2.f; Hyy = crsH2.g;
        end
    end
    
    if ~hasDisps
        % Need spatial grid info even if Displacements aren't provided
        if cfg_prep.coarsegrain > 1
            crs = coarsegrainDisp(zeros(size(ebsd.x)), zeros(size(ebsd.y)), X(1,:), Y(:,1)', cfg_prep.coarsegrain);
        else
            crs.X = X; crs.Y = Y;
        end
        prepData.U = zeros(size(Hxx));
        prepData.V = zeros(size(Hxx));
    end
end

% Create EBSD mapping object
ebsdID = dummyEBSDSimple(ebsd.orientations(1), crs.X, crs.Y);

prepData.Hxx = Hxx; prepData.Hxy = Hxy;
prepData.Hyx = Hyx; prepData.Hyy = Hyy;
prepData.Eeff = calcEffectiveE(Hxx, Hxy, Hyx, Hyy);

% Store preprocessed data in EBSD object
ebsdID.prop.U = prepData.U;
ebsdID.prop.V = prepData.V;
ebsdID.prop.Hxx = Hxx; ebsdID.prop.Hxy = Hxy;
ebsdID.prop.Hyx = Hyx; ebsdID.prop.Hyy = Hyy;
ebsdID.prop.Eeff = prepData.Eeff;

end
