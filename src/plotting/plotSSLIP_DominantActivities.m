function plotSSLIP_DominantActivities(ebsdID, sSLocal, opt)
% plotSSLIP_DominantActivities Plots the top active slip systems and bundles
% the remainder into a 'leftover activity' plot.
% NOTE: Supports context boundaries via plotSSLIP_Context (requires opt.plotGrainContext=1 and opt.grains)

if nargin < 3
    opt = struct;
end

% =========================================================================
% Parse Options & Set Defaults
% =========================================================================
% General Plotting Options
if ~isfield(opt, 'cmap'),               opt.cmap = parula(256); try opt.cmap = viridis(256); catch, end; end % Colormap for activity fields
if ~isfield(opt, 'clim_percentile'),    opt.clim_percentile = 99; end     % Percentile used to clip color axes (avoids extreme outliers dominating colorbar)
if ~isfield(opt, 'stress'),             opt.stress = stressTensor.uniaxial(xvector); end % Macroscopic stress tensor used for calculating Schmid Factors
if ~isfield(opt, 'posConstr'),          opt.posConstr = 0; end            % Flag (0/1): If 1, activities are assumed strictly positive
if ~isfield(opt, 'num_top_systems'),    opt.num_top_systems = 3; end      % Number of dominant slip systems (or bundles) to display before grouping remainder into 'Leftover'
if ~isfield(opt, 'sizeAdjust'),         opt.sizeAdjust = 1; end           % Scale factor for MTEX figure layout spacing
if ~isfield(opt, 'plotGrainContext'),   opt.plotGrainContext = 0; end     % Flag (0/1): Overlay grain boundaries (requires opt.grains to be set)
if ~isfield(opt, 'ori'),                opt.ori = []; end                 % Orientation object used to recover crystal coordinates from specimen coordinates (inv(ori)*sSLocal)

% Output / Export Options
if ~isfield(opt, 'saveFig'),            opt.saveFig = 0; end              % Flag (0/1): Save the resulting figure to disk
if ~isfield(opt, 'plotname'),           opt.plotname = 'SSLIP'; end       % Base string used for exported filenames
if ~isfield(opt, 'outputReport'),       opt.outputReport = 0; end         % Flag (0/1) or String: Export a CSV report of slip activities. If 1, auto-names the file. If string, uses as filepath.

% Similarity Bundling & Activity Thresholding
if ~isfield(opt, 'bundleSimilar'),      opt.bundleSimilar = 0; end        % Flag (0/1): Accumulate crystallographically and kinematically similar systems into common bundles
if ~isfield(opt, 'bundleStrainThresh'), opt.bundleStrainThresh = 0.05; end% Maximum allowable difference (norm) between theoretical 2D deformation gradients to bundle systems together
if ~isfield(opt, 'activityThreshold'),  opt.activityThreshold = 0.01; end % Minimum threshold for a slip system's activity metric to be considered 'active' (otherwise it goes to leftovers)
if ~isfield(opt, 'activityPercentile'), opt.activityPercentile = 85; end  % The statistical percentile of the pixel-wise field used to rank system activity (99 is highly sensitive to localized bands)
if ~isfield(opt, 'enableRotation'),     opt.enableRotation = 0; end       % Flag (0/1): Should be 1 if solving included pseudo-slip rigid body rotations
if ~isfield(opt, 'normalizeInplane'),   opt.normalizeInplane = 0; end     % Flag (0/1): Normalize the 2D deformation gradients prior to similarity checking

bundleAngleThresh = 1e-5; % Hardcoded strict threshold for crystallographic equivalence (angle between normal and burgers vectors)
% =========================================================================

% =========================================================================
% Precompute String and Vector Representations for all Systems
% =========================================================================
num_total_sys = length(sSLocal);
sys_nx = zeros(1, num_total_sys);
sys_ny = zeros(1, num_total_sys);
sys_nz = zeros(1, num_total_sys);
sys_bx = zeros(1, num_total_sys);
sys_by = zeros(1, num_total_sys);
sys_bz = zeros(1, num_total_sys);
sys_strings = cell(1, num_total_sys);

for i = 1:num_total_sys
    obj = sSLocal(i);
    if ~isempty(opt.ori)
        sys_rep = (inv(opt.ori) * obj);
        sys_rep_n = round(sys_rep.n);
        sys_rep_b = round(sys_rep.b, 'round3IndexDirection');
        sys_nx(i) = sys_rep_n.h; sys_ny(i) = sys_rep_n.k; sys_nz(i) = sys_rep_n.l;
        sys_bx(i) = sys_rep_b.u; sys_by(i) = sys_rep_b.v; sys_bz(i) = sys_rep_b.w;
        sys_strings{i} = sprintf('n=(%d,%d,%d) b=[%d,%d,%d]', sys_nx(i), sys_ny(i), sys_nz(i), sys_bx(i), sys_by(i), sys_bz(i));
    else
        sys_rep = obj;
        sys_nx(i) = sys_rep.n.x; sys_ny(i) = sys_rep.n.y; sys_nz(i) = sys_rep.n.z;
        sys_bx(i) = sys_rep.b.x; sys_by(i) = sys_rep.b.y; sys_bz(i) = sys_rep.b.z;
        sys_strings{i} = sprintf('n=(%.2f,%.2f,%.2f) b=[%.2f,%.2f,%.2f]', sys_nx(i), sys_ny(i), sys_nz(i), sys_bx(i), sys_by(i), sys_bz(i));
    end
end
% =========================================================================

% Extract 2D deformation gradient for all systems (used for similarity checking)
Hslip = sSLocal.deformationTensor.matrix;
if isfield(opt, 'enableRotation') && opt.enableRotation
    Hslip(:,:,1,end+1) = [0 -1 0; 1 0 0; 0 0 0];
    if opt.posConstr
        Hslip(:,:,1,end+1) = [0 1 0; -1 0 0; 0 0 0];
    end
end
Hslip11 = reshape(Hslip(1,1,:), 1, []);
Hslip12 = reshape(Hslip(1,2,:), 1, []);
Hslip21 = reshape(Hslip(2,1,:), 1, []);
Hslip22 = reshape(Hslip(2,2,:), 1, []);
A_mat = [Hslip11; Hslip12; Hslip21; Hslip22]; 
if isfield(opt, 'normalizeInplane') && opt.normalizeInplane
    A_mat = A_mat ./ sqrt(sum(A_mat.^2, 1));
end

% Rank systems by percentile activity (robust to outliers, sensitive to bands)
gamma_abs = abs(ebsdID.prop.slipIDcor);
num_sys = size(gamma_abs, 2);
metric_vals = zeros(1, num_sys);
for i = 1:num_sys
    valid_g = gamma_abs(:, i);
    valid_g = valid_g(~isnan(valid_g) & ~isinf(valid_g));
    if isempty(valid_g)
        metric_vals(i) = 0;
    else
        metric_vals(i) = prctile(valid_g, opt.activityPercentile);
    end
end

[sort_metric_acts, sort_sys_idx] = sort(metric_vals, 'descend');

% Form bundles
bundles = {}; % cell array of system indices
rep_sys_list = []; % Array of representative system indices for fast vectorization
leftover_idx = [];

if opt.bundleSimilar
    for i = 1:length(sort_sys_idx)
        sys_i = sort_sys_idx(i);
        
        % If activity is below threshold, it's leftover
        if metric_vals(sys_i) < opt.activityThreshold
            leftover_idx = [leftover_idx, sys_i];
            continue;
        end
        
        % Check similarity against existing bundles (Vectorized)
        matched_bundle = 0;
        if ~isempty(rep_sys_list) && sys_i <= length(sSLocal)
            % Filter out pseudo-slip systems from the comparison list
            valid_reps_mask = rep_sys_list <= length(sSLocal);
            valid_reps = rep_sys_list(valid_reps_mask);
            
            if ~isempty(valid_reps)
                % Vectorized Condition A: Crystallographic Equivalence
                obj_i = sSLocal(sys_i);
                obj_reps = sSLocal(valid_reps);
                ang_n = angle(obj_i.n, obj_reps.n);
                ang_b = angle(obj_i.b, obj_reps.b);
                cond_A = (ang_n(:)' < bundleAngleThresh) & (ang_b(:)' < bundleAngleThresh);
                
                % Vectorized Condition B: 2D Projection Similarity
                A_i = A_mat(:, sys_i);
                A_reps = A_mat(:, valid_reps);
                norm1 = sqrt(sum((A_reps - A_i).^2, 1));
                norm2 = sqrt(sum((A_reps + A_i).^2, 1));
                cond_B = (norm1 < opt.bundleStrainThresh) | (norm2 < opt.bundleStrainThresh);
                
                % Find first match
                match_local_idx = find(cond_A & cond_B, 1);
                if ~isempty(match_local_idx)
                    % Map back to original bundle index
                    valid_bundle_indices = find(valid_reps_mask);
                    matched_bundle = valid_bundle_indices(match_local_idx);
                end
            end
        end
        
        if matched_bundle > 0
            bundles{matched_bundle}(end+1) = sys_i;
        else
            % Limit number of top bundles, otherwise it falls to leftover
            if length(bundles) < opt.num_top_systems
                bundles{length(bundles)+1} = sys_i;
                rep_sys_list(end+1) = sys_i;
            else
                leftover_idx = [leftover_idx, sys_i];
            end
        end
    end
else
    % No bundling: just take top N active systems above threshold
    for i = 1:length(sort_sys_idx)
        sys_i = sort_sys_idx(i);
        if length(bundles) < opt.num_top_systems && metric_vals(sys_i) >= opt.activityThreshold
            bundles{end+1} = sys_i;
        else
            leftover_idx = [leftover_idx, sys_i];
        end
    end
end

n_top = length(bundles);

% Calculate Leftover activity map
leftover_activity = zeros(size(ebsdID.prop.slipIDcor, 1), 1);
for i = 1:length(leftover_idx)
    leftover_activity = leftover_activity + ebsdID.prop.slipIDcor(:, leftover_idx(i));
end

% Build combined top field maps
top_fields = zeros(size(ebsdID.prop.slipIDcor, 1), n_top);
for t_i = 1:n_top
    for mem = 1:length(bundles{t_i})
        sys = bundles{t_i}(mem);
        top_fields(:, t_i) = top_fields(:, t_i) + ebsdID.prop.slipIDcor(:, sys);
    end
end

% Determine common color limits
all_fields_to_plot = [top_fields, leftover_activity];
validData = all_fields_to_plot(~isinf(all_fields_to_plot) & ~isnan(all_fields_to_plot));
if isempty(validData)
    caxisMinMax = [0 0.1];
else
    caxisMinMax = [prctile(validData, 100 - opt.clim_percentile), prctile(validData, opt.clim_percentile)];
    if caxisMinMax(1) < 0 && opt.posConstr, caxisMinMax(1) = 0; end
    if caxisMinMax(2) <= caxisMinMax(1), caxisMinMax(2) = caxisMinMax(1) + 0.05; end
end

% Generate Output Report if requested
if isfield(opt, 'outputReport') && ~isempty(opt.outputReport) && (ischar(opt.outputReport) || isstring(opt.outputReport) || opt.outputReport == 1)
    reportData = struct('BundleID', {}, 'SystemID', {}, 'Activity', {}, 'SchmidFactor', {}, 'H11', {}, 'H12', {}, 'H21', {}, 'H22', {}, 'IsRepresentative', {}, 'n_x', {}, 'n_y', {}, 'n_z', {}, 'b_x', {}, 'b_y', {}, 'b_z', {});
    % Gather all systems to write to avoid duplicate assignment code
    rep_sys_to_write = [];
    rep_b_ids = string.empty;
    rep_is_rep = logical.empty;
    
    for t_i = 1:n_top
        for mem = 1:length(bundles{t_i})
            rep_sys_to_write(end+1) = bundles{t_i}(mem);
            rep_b_ids(end+1) = string(t_i);
            rep_is_rep(end+1) = (mem == 1);
        end
    end
    for i = 1:length(leftover_idx)
        rep_sys_to_write(end+1) = leftover_idx(i);
        rep_b_ids(end+1) = "Leftover";
        rep_is_rep(end+1) = false;
    end
    
    for idx = 1:length(rep_sys_to_write)
        sys = rep_sys_to_write(idx);
        if sys <= length(sSLocal)
            sf = abs(sSLocal(sys).SchmidFactor(opt.stress));
            nx = sys_nx(sys); ny = sys_ny(sys); nz = sys_nz(sys);
            bx = sys_bx(sys); by = sys_by(sys); bz = sys_bz(sys);
        else
            sf = NaN;
            nx = NaN; ny = NaN; nz = NaN; bx = NaN; by = NaN; bz = NaN;
        end
        
        reportData(idx).BundleID = rep_b_ids(idx);
        reportData(idx).SystemID = sys;
        reportData(idx).Activity = metric_vals(sys);
        reportData(idx).SchmidFactor = sf;
        reportData(idx).H11 = A_mat(1, sys);
        reportData(idx).H12 = A_mat(2, sys);
        reportData(idx).H21 = A_mat(3, sys);
        reportData(idx).H22 = A_mat(4, sys);
        reportData(idx).IsRepresentative = rep_is_rep(idx);
        reportData(idx).n_x = nx;
        reportData(idx).n_y = ny;
        reportData(idx).n_z = nz;
        reportData(idx).b_x = bx;
        reportData(idx).b_y = by;
        reportData(idx).b_z = bz;
    end
    
    rep_path = opt.outputReport;
    if ~ischar(rep_path) && ~isstring(rep_path)
        if isfield(opt, 'plotname')
            rep_path = [opt.plotname, '_DominantSlipReport.csv'];
        else
            rep_path = 'DominantSlipReport.csv';
        end
    end
    writeSSLIP_Report(rep_path, reportData);
end

figure;
num_plots = n_top + 1;
if num_plots <= 2
    layout_grid = [1, 2];
elseif num_plots <= 4
    layout_grid = [2, 2];
else
    layout_grid = [2, ceil(num_plots/2)];
end
mtexFig = newMtexFigure('layout', layout_grid);

for t_i = 1:n_top
    nextAxis;
    rep_sys = bundles{t_i}(1);
    num_sys_in_bundle = length(bundles{t_i});
    
    % Plot context if provided using centralized helper
    plotSSLIP_Context(opt);
    
    plot(ebsdID, top_fields(:, t_i), 'micronbar', 'off');
    
    % Use representative system for quiver and title
    if rep_sys <= length(sSLocal)
        sys_obj = sSLocal(rep_sys);
        sf_val = abs(sys_obj.SchmidFactor(opt.stress));
        
        % Quiver trace and slip direction
        sys_obj.b.z = sys_obj.b.z .* sign(sys_obj.b.z) * sign(plottingConvention.default.outOfScreen.z);
        ebsdTrace = ebsdID(round(size(ebsdID,1)/2), round(size(ebsdID,2)/2));
        quiver_scale = 0.2 * (max(ebsdID.x,[],'all') - min(ebsdID.x,[],'all'));
        hold on; quiver(ebsdTrace, quiver_scale * sys_obj.trace, 'color', 'r', 'linewidth', 1.5);
        hold on; quiver(ebsdTrace, quiver_scale * sys_obj.b.normalize, 'color', 'r', 'linewidth', 1.5);
        hold off;
        
        sys_str = sys_strings{rep_sys};
        
        if num_sys_in_bundle > 1
            title(sprintf('[Bundle %d] %s\n%d systems | SF=%.2f | P%d=%.4f', t_i, sys_str, num_sys_in_bundle, sf_val, opt.activityPercentile, metric_vals(rep_sys)));
        else
            title(sprintf('[Top %d] %s\nSF=%.2f | P%d=%.4f', t_i, sys_str, sf_val, opt.activityPercentile, metric_vals(rep_sys)));
        end
    else
        % It's a pseudo-slip rotation system
        if num_sys_in_bundle > 1
            title(sprintf('[Bundle %d] Rotation\n%d systems | P%d=%.4f', t_i, num_sys_in_bundle, opt.activityPercentile, metric_vals(rep_sys)));
        else
            title(sprintf('[Top %d] Rotation\nP%d=%.4f', t_i, opt.activityPercentile, metric_vals(rep_sys)));
        end
    end
    
    caxis(caxisMinMax);
    mtexColorbar('title', '\gamma');
    mtexColorMap(opt.cmap);
end

% Leftover activity
nextAxis;
plotSSLIP_Context(opt);
plot(ebsdID, leftover_activity, 'micronbar', 'off');
title(sprintf('Leftover Activity\nSum of %d systems', length(leftover_idx)));
caxis(caxisMinMax);
mtexColorbar('title', '\gamma');
mtexColorMap(opt.cmap);
hold off;

if isfield(opt, 'sizeAdjust')
    mtexFig.figSizeFactor = opt.sizeAdjust;
    mtexFig.drawNow;
end

if opt.saveFig && isfield(opt, 'plotname') && ~isempty(opt.plotname)
    exportScaledFigure(gcf, [opt.plotname, '_SlipActivities_Dominant.jpg'], 'target_px', 2000, 'dpi', 300, 'nodatestring');
end

end
