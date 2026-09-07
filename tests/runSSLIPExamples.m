function results = runSSLIPExamples(outputDirectory)
% Run both supplied examples in a temporary copy with the active MTEX version.
% Saves numeric results and the final activity plot for each example.
% Use a separate MATLAB session: the supplied scripts close figures.
if nargin < 1
    outputDirectory = tempname;
end
if ~exist(outputDirectory,'dir'), mkdir(outputDirectory); end
[ok,info] = fileattrib(outputDirectory);
assert(ok, 'Cannot access output directory.');
outputDirectory = info.Name;
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
pathCleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
oldVisibility = get(groot,'defaultFigureVisible');
figureCleanup = onCleanup(@() set(groot,'defaultFigureVisible',oldVisibility)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
oldRng = rng;
rngCleanup = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(0);
ps = parallel.Settings; %#ok<NASGU>
settingsRoot = settings;
poolSetting = settingsRoot.parallel.client.pool.AutoCreate;
hadTemporaryValue = hasTemporaryValue(poolSetting);
oldAutoCreate = poolSetting.ActiveValue;
poolCleanup = onCleanup(@() restorePoolSetting(poolSetting,hadTemporaryValue,oldAutoCreate)); %#ok<NASGU>
poolSetting.TemporaryValue = false;
% Keep example-generated files out of the source checkout.
work = tempname(outputDirectory);
mkdir(work);
copyfile(fullfile(root,'src'),fullfile(work,'src'));
copyfile(fullfile(root,'examples'),fullfile(work,'examples'));
mkdir(fullfile(work,'data'));
copyfile(fullfile(root,'data','NiSuperAloy_Aligned.mat'),fullfile(work,'data'));
results = struct;
for script = {'NiSuperAlloyExperiment','virtualExperimentHCP'}
    name = script{1};
    result = runExample(fullfile(work,'examples',[name '.m']));
    result.mtexVersion = getMTEXpref('version');
    results.(name) = result;
    save(fullfile(outputDirectory,[name '_numeric.mat']),'-struct','result');
    exportgraphics(gcf,fullfile(outputDirectory,[name '.png']),'Resolution',120);
    fprintf('PASS: %s (%d pixels, %d slip systems)\n', ...
        name,size(result.gamma,2),size(result.gamma,1));
    close all;
end
fprintf('Example results saved in %s\n',outputDirectory);
end

function result = runExample(scriptPath)
% Isolate the scripts' clear statements from the runner's workspace.
run(scriptPath);
result.gamma = ebsdID.prop.slipIDcor;
result.residual = ebsdID.prop.residualEeff;
result.solverExitFlag = ebsdID.prop.solverExitFlag;
result.coords = [ebsdID.x(:),ebsdID.y(:)];
result.H = [ebsdID.prop.Hxx(:),ebsdID.prop.Hxy(:),ebsdID.prop.Hyx(:),ebsdID.prop.Hyy(:)];
result.slipTensors = sSLocal.deformationTensor.matrix;
result.options = optOut;
assert(size(result.gamma,2) == numel(ebsdID.x));
assert(any(isfinite(result.gamma(:)) & abs(result.gamma(:)) > 1e-5));
assert(any(isfinite(result.residual(:))));
assert(~any(double(optOut.plotname)<32),'Example output name contains a control character.');
end

function restorePoolSetting(poolSetting,hadTemporaryValue,oldValue)
if hadTemporaryValue
    poolSetting.TemporaryValue = oldValue;
else
    clearTemporaryValue(poolSetting);
end
end
