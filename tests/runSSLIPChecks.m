function runSSLIPChecks
% Run focused regression checks with the currently initialized MTEX version.
% Requires MATLAB Optimization Toolbox and Parallel Computing Toolbox.
% Run in a separate MATLAB session: the plotting checks close figures.

root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
pathCleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(fullfile(root, 'src'));
oldVisibility = get(groot, 'defaultFigureVisible');
figureCleanup = onCleanup(@() set(groot, 'defaultFigureVisible', oldVisibility)); %#ok<NASGU>
set(groot, 'defaultFigureVisible', 'off');
ps = parallel.Settings; %#ok<NASGU> % Initialize the toolbox settings tree.
settingsRoot = settings;
poolSetting = settingsRoot.parallel.client.pool.AutoCreate;
hadTemporaryValue = hasTemporaryValue(poolSetting);
oldAutoCreate = poolSetting.ActiveValue;
poolCleanup = onCleanup(@() restorePoolSetting(poolSetting, hadTemporaryValue, oldAutoCreate)); %#ok<NASGU>
poolSetting.TemporaryValue = false;

fprintf('VALIDATION MTEX: %s\n', getMTEXpref('version'));
sS = slipSystem(xvector, yvector);
opt = struct('minEeff', .005, 'threshResidual', .01, ...
    'normalizeInplane', 0, 'posConstr', 0);
hxy = [.02 .2 -.2 .001 NaN];
z = zeros(size(hxy));
[gFixed,rFixed,flagsFixed] = SSLIPConeprogConstrMinAbs(sS,z,hxy,z,z,opt);
assert(max(abs(gFixed(1:3) - [.01 .19 -.19])) < 2e-6);
assert(gFixed(4) == 0 && rFixed(4) == 0);
assert(isnan(gFixed(5)) && isnan(rFixed(5)));
assert(all(flagsFixed(1:3) == 1));
assert(all(isnan(flagsFixed(4:5))));
opt.threshResidualFraction = .1;
[g,r] = SSLIPConeprogConstrMinAbs(sS,z,hxy,z,z,opt);
expected = [.01, .2-.02/sqrt(2), -(.2-.02/sqrt(2))];
assert(max(abs(g(1:3)-expected)) < 2e-6);
assert(max(abs(r(1:3)-[.01; .02/sqrt(2); .02/sqrt(2)])) < 2e-6);
assert(g(4) == 0 && isnan(g(5)));
opt.threshResidualFraction = 0;
[gZero,rZero] = SSLIPConeprogConstrMinAbs(sS,z,hxy,z,z,opt);
assert(max(abs(gZero(1:4)-gFixed(1:4))) < 2e-6);
assert(max(abs(rZero(1:4)-rFixed(1:4))) < 2e-6);
opt.threshResidualFraction = .1;
opt.posConstr = 1;
[gPos,rPos] = SSLIPConeprogConstrMinAbs(sS,z(1:2),hxy(1:2),z(1:2),z(1:2),opt);
assert(max(abs(gPos-expected(1:2))) < 2e-6);
assert(all(gPos >= 0));
sZero = sS;
sZero.CRSS = 0;
[gCRSS,~] = SSLIPConeprogConstrMinAbs(sZero,z(1:2),hxy(1:2),z(1:2),z(1:2),opt);
assert(max(abs(gCRSS-gPos)) < 2e-6);
[gInfeasible,rInfeasible,flagInfeasible] = SSLIPConeprogConstrMinAbs(sS,0,-.2,0,0,opt);
assert(isnan(gInfeasible) && isnan(rInfeasible) && flagInfeasible < 0);
fprintf('PASS: coneprog tolerances, signs, invalid pixels, CRSS independence, and per-pixel exit flags.\n');

% Exercise the full single-slip call path with affine displacement fields.
[X,Y] = meshgrid(0:3,0:3);
CS = crystalSymmetry('m-3m');
ori = orientation.byEuler(0,0,0,CS);
ebsd = dummyEBSDSimple(ori,X,Y);
opt.IDMethod = 3;
opt.posConstr = 0;
opt.filterSize = 0;
opt.coarsegrain = 1;
opt.plotSSLIP = 0;
opt.plotDefGrad = 0;
opt.saveFig = 0;
opt.threshResidualFraction = .001;
[ebsdFloor,~] = SSLIP(ebsd,.02*Y,.005*Y,sS,opt);
assert(max(abs(ebsdFloor.prop.slipIDcor(:)-.02)) < 1e-10);
opt.threshResidualFraction = .1;
[ebsdRelative,~] = SSLIP(ebsd,.2*Y,.015*Y,sS,opt);
assert(max(abs(ebsdRelative.prop.slipIDcor(:)-.2)) < 1e-10);
opt = rmfield(opt,'threshResidualFraction');
[ebsdReject,~] = SSLIP(ebsd,.2*Y,.015*Y,sS,opt);
assert(all(ebsdReject.prop.slipIDcor(:) == 0));
fprintf('PASS: single-slip floor, relative threshold acceptance, fixed-threshold rejection.\n');

% Selecting one system per pixel must retain its fitted amplitude.
opt.singleSlipPerPixel = 1;
opt.threshResidual = .2;
sTwo = [slipSystem(xvector,yvector); slipSystem(yvector,xvector)];
[ebsdFirst,~] = SSLIP(ebsd,.2*Y,.1*X,sTwo,opt);
assert(max(abs(ebsdFirst.prop.slipIDcor(1,:)-.2)) < 1e-10);
assert(all(ebsdFirst.prop.slipIDcor(2,:) == 0));
[ebsdSecond,~] = SSLIP(ebsd,.1*Y,.2*X,sTwo,opt);
assert(all(ebsdSecond.prop.slipIDcor(1,:) == 0));
assert(max(abs(ebsdSecond.prop.slipIDcor(2,:)-.2)) < 1e-10);
fprintf('PASS: singleSlipPerPixel preserves the best-fit activity.\n');

% Check real residual plots, including degenerate data ranges.
popt = struct('stress',stressTensor.uniaxial(xvector),'NoSs',1, ...
    'layout',[], 'plotTraces',0, 'logscale',0, 'posConstr',0, ...
    'maxE',1, 'saveFig',0, 'cmap',parula(256),'plotResidual',1,'IDMethod',1);
residualCases = {ones(16,1)*.02, zeros(16,1), NaN(16,1), ...
    [.01;.03;NaN;Inf;zeros(12,1)]};
expectedMax = [.02 1 1 .03];
for k=1:numel(residualCases)
    plotSSLIP(linspace(.01,.2,16),residualCases{k},ebsd,sS,popt);
    ax = findall(gcf,'Type','axes');
    assert(any(arrayfun(@(a) isequal(a.CLim,[0 expectedMax(k)]),ax)));
    close all;
end
fprintf('PASS: residual plots with constant, zero, missing, and nonfinite values.\n');

% Smoke-test every trace plotting branch, including >24 slip systems.
popt.plotResidual = 0;
popt.plotTraces = 1;
popt.traceScale = .5;
plotSSLIP(linspace(.01,.2,16),zeros(16,1),ebsd,sS,popt);
close all;
popt.plotSingle = 1;
plotSSLIP(linspace(.01,.2,16),zeros(16,1),ebsd,sS,popt);
close all;
popt.plotSingle = 0;
popt.NoSs = 1:25;
sMany = repmat(sS,25,1);
plotSSLIP(repmat(linspace(.01,.2,16),25,1),zeros(16,1),ebsd,sMany,popt);
close all;
fprintf('PASS: traces in combined, individual, and >24-system plots.\n');
fprintf('ALL SSLIP REVIEW CHECKS PASSED.\n');
end

function restorePoolSetting(poolSetting, hadTemporaryValue, oldValue)
if hadTemporaryValue
    poolSetting.TemporaryValue = oldValue;
else
    clearTemporaryValue(poolSetting);
end
end
