function plotSSLIP_SlipActivities(ebsdID, sSLocal, opt)
% plotSSLIP_SlipActivities Plots the physical slip system activities in a grid or singly.
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
if ~isfield(opt, 'plotSingle'),         opt.plotSingle = 0; end           % Flag (0/1): If 1, plot and save each slip system activity as its own individual figure
if ~isfield(opt, 'plotTraces'),         opt.plotTraces = 1; end           % Flag (0/1): Overlay quiver vectors showing the theoretical slip trace and burgers direction
if ~isfield(opt, 'logscale'),           opt.logscale = 0; end             % Flag (0/1): Plot color axes on a logarithmic scale
if ~isfield(opt, 'layout'),             opt.layout = []; end              % Explicit layout grid for the MTEX figure (e.g., [2 3])
if ~isfield(opt, 'NoSs'),               opt.NoSs = 1:length(sSLocal); end % Array of slip system indices to explicitly plot
if ~isfield(opt, 'sizeAdjust'),         opt.sizeAdjust = 1; end           % Scale factor for MTEX figure layout spacing
if ~isfield(opt, 'plotGrainContext'),   opt.plotGrainContext = 0; end     % Flag (0/1): Overlay grain boundaries (requires opt.grains to be set)
if ~isfield(opt, 'ori'),                opt.ori = []; end                 % Orientation object used to recover crystal coordinates from specimen coordinates (inv(ori)*sSLocal)

% Output / Export Options
if ~isfield(opt, 'saveFig'),            opt.saveFig = 0; end              % Flag (0/1): Save the resulting figure(s) to disk
if ~isfield(opt, 'plotname'),           opt.plotname = 'SSLIP'; end       % Base string used for exported filenames
if ~isfield(opt, 'saveExt'),            opt.saveExt = ''; end             % Secondary file extension to save (e.g., '.pdf')
if ~isfield(opt, 'outputReport'),       opt.outputReport = 0; end         % Flag (0/1) or String: Export a CSV report of slip activities. If 1, auto-names the file. If string, uses as filepath.

% Similarity Bundling
if ~isfield(opt, 'bundleSimilar'),      opt.bundleSimilar = 0; end        % Flag (0/1): Accumulate crystallographically and kinematically similar systems into common bundles
if ~isfield(opt, 'bundleStrainThresh'), opt.bundleStrainThresh = 0.05; end% Maximum allowable difference (norm) between theoretical 2D deformation gradients to bundle systems together
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
        sys_rep_b = round(sys_rep.b, 'uvw');
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
NoSs = opt.NoSs;
slipIDcor = ebsdID.prop.slipIDcor;

% calc SF, used in the title of each slip activity field
SF = abs(SchmidFactor(sSLocal(NoSs), opt.stress));
SF = round(SF*100)/100;

% make sure slip directions are plotted above the activity fields
sSLocal.b.z = sSLocal.b.z .* sign(sSLocal.b.z) * sign(plottingConvention.default.outOfScreen.z);

if opt.bundleSimilar
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

    bundles = {}; 
    rep_sys_list = []; % Array of representative system indices for fast vectorization
    for i = 1:length(NoSs)
        sys_i = NoSs(i);
        
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
            bundles{matched_bundle}(end+1) = i;
        else
            bundles{end+1} = i;
            rep_sys_list(end+1) = sys_i;
        end
    end

    % Build bundled arrays
    n_bundles = length(bundles);
    bundled_slipIDcor = zeros(size(slipIDcor, 1), n_bundles);
    bundled_NoSs = zeros(1, n_bundles);
    bundled_SF = zeros(1, n_bundles);
    bundled_count = zeros(1, n_bundles);
    
    for b = 1:n_bundles
        for mem = 1:length(bundles{b})
            idx = bundles{b}(mem);
            bundled_slipIDcor(:, b) = bundled_slipIDcor(:, b) + slipIDcor(:, idx);
        end
        rep_idx = bundles{b}(1);
        bundled_NoSs(b) = NoSs(rep_idx);
        bundled_SF(b) = SF(rep_idx);
        bundled_count(b) = length(bundles{b});
    end
    
    % Override variables for the rest of the script
    slipIDcor = bundled_slipIDcor;
    NoSs = bundled_NoSs;
    SF = bundled_SF;
else
    bundled_count = ones(1, length(NoSs));
    bundles = num2cell(1:length(NoSs));
end

% Generate Output Report if requested
if isfield(opt, 'outputReport') && ~isempty(opt.outputReport) && (ischar(opt.outputReport) || isstring(opt.outputReport) || opt.outputReport == 1)
    reportData = struct('BundleID', {}, 'SystemID', {}, 'Activity', {}, 'SchmidFactor', {}, 'H11', {}, 'H12', {}, 'H21', {}, 'H22', {}, 'IsRepresentative', {}, 'n_x', {}, 'n_y', {}, 'n_z', {}, 'b_x', {}, 'b_y', {}, 'b_z', {});
    idx_count = 1;
    % Original NoSs before overriding
    orig_NoSs = opt.NoSs;
    for b = 1:length(bundles)
        for mem = 1:length(bundles{b})
            idx = bundles{b}(mem);
            sys = orig_NoSs(idx);
            isRep = (mem == 1);
            
            % Compute activity metric if needed (nth percentile of raw field)
            valid_g = abs(ebsdID.prop.slipIDcor(:, idx));
            valid_g = valid_g(~isnan(valid_g) & ~isinf(valid_g));
            if isempty(valid_g), act_val = 0; else, act_val = prctile(valid_g, opt.clim_percentile); end
            
            if sys <= length(sSLocal)
                obj = sSLocal(sys);
                sf = abs(obj.SchmidFactor(opt.stress));
                
                nx = sys_nx(sys); ny = sys_ny(sys); nz = sys_nz(sys);
                bx = sys_bx(sys); by = sys_by(sys); bz = sys_bz(sys);
                
                % Compute H components for original system if A_mat not available
                Hslip = obj.deformationTensor.matrix;
                h11 = Hslip(1,1); h12 = Hslip(1,2); h21 = Hslip(2,1); h22 = Hslip(2,2);
            else
                sf = NaN;
                h11 = NaN; h12 = NaN; h21 = NaN; h22 = NaN;
                nx = NaN; ny = NaN; nz = NaN; bx = NaN; by = NaN; bz = NaN;
            end
            
            reportData(idx_count).BundleID = string(b);
            reportData(idx_count).SystemID = sys;
            reportData(idx_count).Activity = act_val;
            reportData(idx_count).SchmidFactor = sf;
            reportData(idx_count).H11 = h11;
            reportData(idx_count).H12 = h12;
            reportData(idx_count).H21 = h21;
            reportData(idx_count).H22 = h22;
            reportData(idx_count).IsRepresentative = isRep;
            reportData(idx_count).n_x = nx;
            reportData(idx_count).n_y = ny;
            reportData(idx_count).n_z = nz;
            reportData(idx_count).b_x = bx;
            reportData(idx_count).b_y = by;
            reportData(idx_count).b_z = bz;
            idx_count = idx_count + 1;
        end
    end
    rep_path = opt.outputReport;
    if ~ischar(rep_path) && ~isstring(rep_path)
        if isfield(opt, 'plotname')
            rep_path = [opt.plotname, '_SlipActivitiesReport.csv'];
        else
            rep_path = 'SlipActivitiesReport.csv';
        end
    end
    writeSSLIP_Report(rep_path, reportData);
end

if opt.logscale && ~opt.posConstr
    slipIDcor = abs(slipIDcor);
    warning('Since logaritmic plotting is required (opt.logscale==1), while pos. constr. was not used (opt.posConstr==0), absolute values of slip activities are plotted')
end

maxFields = 24;
caxisMinMax = [0 0];

if ~opt.plotSingle
    if isempty(opt.layout)
        figure; f1=newMtexFigure;
    else
        figure; f1=newMtexFigure('layout',opt.layout);
    end
    
    for i = 1:min([length(NoSs) , maxFields])
        nextAxis
        plotSSLIP_Context(opt);
        ha(i) = plot(ebsdID, slipIDcor(:,i),'micronbar','off' ); 
        sys_obj = sSLocal(NoSs(i));
        sys_str = sys_strings{NoSs(i)};
        
        if bundled_count(i) > 1
            title(sprintf('[Bundle %d] #%d %s\n%d systems | SF=%.2f', i, NoSs(i), sys_str, bundled_count(i), SF(i)));
        else
            title(sprintf('[Top %d] #%d %s | SF=%.2f', i, NoSs(i), sys_str, SF(i)));
        end
        
        if opt.plotTraces
            ebsdTrace = ebsdID(round(size(ebsdID,1)/2),round(size(ebsdID,2)/2));
            hold on
            quiver(ebsdTrace,0.2*(max(ebsdID.x,[],'all')-min(ebsdID.x,[],'all')) * sSLocal(NoSs(i)).trace,'color','r');
            quiver(ebsdTrace,0.2*(max(ebsdID.x,[],'all')-min(ebsdID.x,[],'all')) * sSLocal(NoSs(i)).b.normalize,'color','r');
            hold off
        end
    end
    
    validData = slipIDcor(~isinf(slipIDcor) & ~isnan(slipIDcor));
    
    if opt.logscale || isfield(opt,'logmin')
        caxisMinMax(1) = opt.logmin;
    elseif opt.posConstr == 1
        caxisMinMax(1) = 0;
    elseif isfield(opt,'maxE')
        caxisMinMax(1) = -opt.maxE;
    else
        caxisMinMax(1) = prctile(validData, 100 - opt.clim_percentile);
    end
    
    if isfield(opt,'maxE')
        caxisMinMax(2) = opt.maxE;
    else
        caxisMinMax(2) = prctile(validData, opt.clim_percentile);
    end
    
    if isempty(caxisMinMax(2)) || isnan(caxisMinMax(2)), caxisMinMax(2) = 0.05; end
    if caxisMinMax(2) <= caxisMinMax(1)
        if opt.logscale && caxisMinMax(1) > 0
            caxisMinMax(2) = caxisMinMax(1) * 10;
        else
            caxisMinMax(2) = caxisMinMax(1) + 0.05;
        end
    end

    for i=1:min([length(NoSs) , maxFields])
        ha(i).Parent.CLim = caxisMinMax;
        if opt.logscale
            ha(i).Parent.ColorScale = 'log';
        end
    end

    mtexColorMap(opt.cmap)
    if opt.logscale && ~opt.posConstr
        mtexColorbar('title','|\gamma|')
    else
        mtexColorbar('title','\gamma')
    end
    
    if isfield(opt,'sizeAdjust')
        f1.figSizeFactor = opt.sizeAdjust;
        f1.innerPlotSpacing = f1.innerPlotSpacing * opt.sizeAdjust;
        f1.drawNow;
    end
    
    if opt.saveFig && isfield(opt, 'plotname') && ~isempty(opt.plotname)
        saveFigure([opt.plotname, '_SlipActivities.png'])
        if isfield(opt,'saveExt')
            saveFigure([opt.plotname, '_SlipActivities', opt.saveExt])
        end
    end

    % If more than 24
    if length(NoSs) > maxFields
        if isempty(opt.layout)
            figure; f2=newMtexFigure;
        else
            figure; f2=newMtexFigure('layout',opt.layout);
        end
        
        for i = maxFields+1 : length(NoSs)
            nextAxis
            plotSSLIP_Context(opt);
            ha2(i-maxFields) = plot(ebsdID, slipIDcor(:,i),'micronbar','off' ); 
            sys_obj = sSLocal(NoSs(i));
            sys_str = sys_strings{NoSs(i)};
            
            if bundled_count(i) > 1
                title(sprintf('[Bundle %d] #%d %s\n%d systems | SF=%.2f', i, NoSs(i), sys_str, bundled_count(i), SF(i)));
            else
                title(sprintf('[Top %d] #%d %s | SF=%.2f', i, NoSs(i), sys_str, SF(i)));
            end
            if opt.plotTraces
                ebsdTrace = ebsdID(round(size(ebsdID,1)/2),round(size(ebsdID,2)/2));
                hold on
                quiver(ebsdTrace,0.2*(max(ebsdID.x,[],'all')-min(ebsdID.x,[],'all')) * sSLocal(NoSs(i)).trace,'color','r');
                quiver(ebsdTrace,0.2*(max(ebsdID.x,[],'all')-min(ebsdID.x,[],'all')) * sSLocal(NoSs(i)).b.normalize,'color','r');
                hold off
            end
        end
        
        for i=1:length(NoSs)-maxFields
            ha2(i).Parent.CLim = caxisMinMax;
            if opt.logscale
                ha2(i).Parent.ColorScale = 'log';
            end
        end

        mtexColorMap(opt.cmap)
        if opt.logscale && ~opt.posConstr
            mtexColorbar('title','|\gamma|')
        else
            mtexColorbar('title','\gamma')
        end
        
        if isfield(opt,'sizeAdjust')
            f2.figSizeFactor = opt.sizeAdjust;
            f2.innerPlotSpacing = f2.innerPlotSpacing * opt.sizeAdjust;
            f2.drawNow;
        end
        
        if opt.saveFig && isfield(opt, 'plotname') && ~isempty(opt.plotname)
            saveFigure([opt.plotname, '_SlipActivities_24plus.png'])
            if isfield(opt,'saveExt')
                saveFigure([opt.plotname, '_SlipActivities_24plus',opt.saveExt])
            end
        end
    end
else
    % Single plots
    if opt.saveFig && isfield(opt, 'plotname') && ~isempty(opt.plotname)
        if ~exist(opt.plotname,'dir')
            mkdir(opt.plotname)
        end
    end
    
    validSingle = slipIDcor(~isinf(slipIDcor) & ~isnan(slipIDcor));
    
    if opt.logscale || isfield(opt,'logmin')
        cmin = opt.logmin;
    elseif opt.posConstr == 1
        cmin = 0;
    elseif isfield(opt,'maxE')
        cmin = -opt.maxE;
    else
        cmin = prctile(validSingle(:), 100 - opt.clim_percentile);
    end
    
    if isfield(opt,'maxE')
        cmax = opt.maxE;
    else
        cmax = prctile(validSingle(:), opt.clim_percentile);
    end
    
    if isempty(cmax) || isnan(cmax), cmax = cmin + 0.05; end
    if cmax <= cmin
        if opt.logscale && cmin > 0
            cmax = cmin * 10;
        else
            cmax = cmin + 0.05;
        end
    end
    
    for i = 1:length(NoSs) 
        figure;
        f1=newMtexFigure;
        plotSSLIP_Context(opt);
        ha = plot(ebsdID, slipIDcor(:,i),'micronbar','off' ); 
        sys_obj = sSLocal(NoSs(i));
        sys_str = sys_strings{NoSs(i)};
        
        if bundled_count(i) > 1
            title(sprintf('[Bundle %d] #%d %s\n%d systems | SF=%.2f', i, NoSs(i), sys_str, bundled_count(i), SF(i)));
        else
            title(sprintf('[Top %d] #%d %s | SF=%.2f', i, NoSs(i), sys_str, SF(i)));
        end
        if opt.plotTraces
            ebsdTrace = ebsdID(round(size(ebsdID,1)/2),round(size(ebsdID,2)/2));
            hold on
            quiver(ebsdTrace,0.2*(ebsdID.xmax-ebsdID.xmin) * sSLocal(NoSs(i)).trace,'color','r','linewidth',5);
            quiver(ebsdTrace,0.2*(ebsdID.xmax-ebsdID.xmin) * sSLocal(NoSs(i)).b.normalize,'color','r','linewidth',5);
            hold off
        end
        
        caxis([cmin cmax])
        mtexColorMap(opt.cmap)
        if opt.logscale && ~opt.posConstr
            mtexColorbar('title','|\gamma|')
        else
            mtexColorbar('title','\gamma')
        end

        set(gca,'fontSize',50)
        f1.drawNow;
        
        if isfield(opt,'sizeAdjust')
            f1.figSizeFactor = opt.sizeAdjust;
            f1.drawNow;
        end
        
        if opt.saveFig && isfield(opt, 'plotname') && ~isempty(opt.plotname)
            saveFigure([opt.plotname,'/sS',num2str(NoSs(i)) '.png'])
            if isfield(opt,'saveExt')
                saveFigure([opt.plotname,'/sS',num2str(NoSs(i)),opt.saveExt])
            end
        end
    end
end
end
