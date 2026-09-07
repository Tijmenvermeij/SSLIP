function comparison = compareSSLIPExamples(referenceDirectory,candidateDirectory)
% Compare numeric outputs saved by runSSLIPExamples under two MTEX versions.
% Reports solver-status differences explicitly; does not assume unique L1 fits.
names = {'NiSuperAlloyExperiment','virtualExperimentHCP'};
rows = struct([]);
for k=1:numel(names)
    name = names{k};
    a = load(fullfile(referenceDirectory,[name '_numeric.mat']));
    b = load(fullfile(candidateDirectory,[name '_numeric.mat']));
    assert(isequal(size(a.gamma),size(b.gamma)), 'Activity dimensions differ.');
    assert(isequal(size(a.slipTensors),size(b.slipTensors)), 'Slip tensor dimensions differ.');
    coordinateDifference = maxCommonDifference(a.coords,b.coords);
    assert(coordinateDifference < 1e-8, 'Pixel grids differ; align them before comparing.');
    tensorDifference = maxCommonDifference(a.slipTensors,b.slipTensors);
    assert(tensorDifference < 1e-12, 'Slip tensors differ; check system ordering and orientation.');
    sameStatus = a.solverExitFlag == b.solverExitFlag | ...
        (isnan(a.solverExitFlag) & isnan(b.solverExitFlag));
    common = all(isfinite(a.gamma),1) & all(isfinite(b.gamma),1);
    A = reshape(a.slipTensors(1:2,1:2,:),4,[]);
    rows(k).example = string(name);
    rows(k).pixels = size(a.gamma,2);
    rows(k).changedSolverStatus = nnz(~sameStatus);
    rows(k).maxInputDifference = maxCommonDifference(a.H,b.H);
    rows(k).maxActivityDifference = maxCommonDifference(a.gamma(:,common),b.gamma(:,common));
    rows(k).maxResidualDifference = maxCommonDifference(a.residual,b.residual);
    rows(k).maxReconstructionDifference = maxCommonDifference(A*a.gamma(:,common),A*b.gamma(:,common));
    rows(k).maxObjectiveDifference = maxCommonDifference(sum(abs(a.gamma(:,common)),1),sum(abs(b.gamma(:,common)),1));
    fprintf('%s: %d changed solver outcomes; max activity difference %.3g; residual difference %.3g.\n', ...
        name,rows(k).changedSolverStatus,rows(k).maxActivityDifference,rows(k).maxResidualDifference);
end
comparison = struct2table(rows);
end

function delta = maxCommonDifference(x,y)
assert(isequal(size(x),size(y)), 'Compared arrays have different dimensions.');
common = isfinite(x) & isfinite(y);
if any(common(:))
    delta = max(abs(x(common)-y(common)));
else
    delta = NaN;
end
end
