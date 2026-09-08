function runSSLIPChecks
% Run focused regression checks with the currently initialized MTEX version.
% Requires Optimization, Parallel Computing, and Image Processing Toolboxes.
% Run in a separate MATLAB session: the plotting checks close figures.

root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
pathCleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(genpath(fullfile(root,'src')));
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

% Philipp's rotation basis, including both signs under positive slip bounds.
% These two slip tensors and the rotation tensor are linearly independent.
sRotation = [slipSystem(vector3d(1,0,1),yvector); ...
    slipSystem(vector3d(1,1,0),vector3d(1,-1,0))];
tensor = sRotation.deformationTensor.matrix;
physicalA = [reshape(tensor(1,1,:),1,[]); reshape(tensor(1,2,:),1,[]); ...
    reshape(tensor(2,1,:),1,[]); reshape(tensor(2,2,:),1,[])];
rotationBasis = [0;-1;1;0];
expectedSlip = [.2 0 0 .2 .1; .1 0 0 .1 .2];
expectedRotation = [0 .03 -.03 .02 -.02];
rotationOpt = struct('minEeff',0,'threshResidual',1e-6, ...
    'normalizeInplane',0,'posConstr',0,'enableRotation',1);
% coneprog permits constraint violation up to its convergence tolerance.
coneOptions = optimoptions('coneprog');
constraintTolerance = coneOptions.ConstraintTolerance;
for normalized = 0:1
    rotationOpt.normalizeInplane = normalized;
    modelA = physicalA;
    if normalized
        modelA = modelA ./ vecnorm(modelA);
    end
    inputH = modelA*expectedSlip + rotationBasis*expectedRotation;
    for positive = 0:1
        rotationOpt.posConstr = positive;
        [gRot,rRot,flagRot] = SSLIPConeprogConstrMinAbs(sRotation, ...
            inputH(1,:),inputH(2,:),inputH(3,:),inputH(4,:),rotationOpt);
        assert(isequal(size(gRot),[3 5]) && all(flagRot == 1));
        assert(max(abs(gRot(1:2,:)-expectedSlip),[],'all') < 2e-5);
        assert(max(abs(gRot(3,:)-expectedRotation)) < 2e-5);
        reconstructedH = modelA*gRot(1:2,:) + rotationBasis*gRot(3,:);
        assert(max(abs(vecnorm(reconstructedH-inputH)'-rRot)) < 1e-10);
        assert(all(rRot <= rotationOpt.threshResidual+constraintTolerance), ...
            'Maximum residual %.9g with threshold %.9g (normalized=%d, positive=%d).', ...
            max(rRot),rotationOpt.threshResidual,normalized,positive);

        % A diagonal perturbation cannot be represented by these bases.
        % Check the residual budget and fit accuracy for a noisy mixture.
        noisyH = inputH(:,4) + [1e-4;0;0;1e-4];
        rotationOpt.threshResidual = .002;
        [gNoise,rNoise,flagNoise] = SSLIPConeprogConstrMinAbs(sRotation, ...
            noisyH(1),noisyH(2),noisyH(3),noisyH(4),rotationOpt);
        assert(flagNoise == 1 && rNoise <= .002+constraintTolerance);
        assert(max(abs(gNoise-[expectedSlip(:,4);expectedRotation(4)])) < .01);
        assert(abs(norm([modelA rotationBasis]*gNoise-noisyH)-rNoise) < 1e-10);
        rotationOpt.threshResidual = 1e-6;
    end
end
rotationOpt.normalizeInplane = 0;
rotationOpt.posConstr = 0;
negativeH = physicalA*[-.2;.1] + rotationBasis*(-.03);
[gNegative,~,flagNegative] = SSLIPConeprogConstrMinAbs(sRotation, ...
    negativeH(1),negativeH(2),negativeH(3),negativeH(4),rotationOpt);
assert(flagNegative == 1 && max(abs(gNegative-[-.2;.1;-.03])) < 2e-5);
rotationOpt.minEeff = .005;
[gSkip,rSkip,flagSkip] = SSLIPConeprogConstrMinAbs(sRotation, ...
    [0 NaN],[-.03 0],[.03 0],[0 0],rotationOpt);
assert(all(gSkip(:,1) == 0) && rSkip(1) == 0 && isnan(flagSkip(1)));
assert(all(isnan(gSkip(:,2))) && isnan(rSkip(2)) && isnan(flagSkip(2)));
fprintf('PASS: rotation signs, radians after normalization, mixtures, noise, and skipped/invalid pixels.\n');

% Exercise the full single-slip call path with affine displacement fields.
[X,Y] = meshgrid(0:3,0:3);
CS = crystalSymmetry('m-3m');
ori = orientation.byEuler(0,0,0,CS);
ebsd = dummyEBSDSimple(ori,X,Y);

% Philipp's stress alignment must return the actual fitting basis, retaining
% the full list and original NoSs labels even for a reordered subset.
stressSystems = [sRotation(1); slipSystem(xvector,zvector); sRotation(2)];
stressSystems.CRSS = [0;2;3];
originalTensors = stressSystems.deformationTensor.matrix;
stressOpt = struct('IDMethod',1,'posConstr',1,'minEeff',0, ...
    'threshResidual',1e-6,'filterSize',0,'coarsegrain',1, ...
    'plotSSLIP',0,'plotDefGrad',0,'saveFig',0,'NoSs',[3 1], ...
    'stress',stressTensor([2 -1 .4;-1 -2 .2;.4 .2 .5]));
% Independently specified signs: resolved shear is negative for system 1,
% positive for systems 2 and 3. The fitted order is system 3 then system 1.
alignedA = [physicalA(:,2),-physicalA(:,1)];
stressGamma = [.12;.21];
for normalized = 0:1
    stressOpt.normalizeInplane = normalized;
    fitA = alignedA;
    if normalized, fitA = fitA ./ vecnorm(fitA); end
    for method = [1 2]
        stressOpt.IDMethod = method;
        stressOpt.enableRotation = method == 1;
        theta = -.025*stressOpt.enableRotation;
        stressH = fitA*stressGamma + rotationBasis*theta;
        stressData = struct('Hxx',ones(size(X))*stressH(1), ...
            'Hxy',ones(size(X))*stressH(2),'Hyx',ones(size(X))*stressH(3), ...
            'Hyy',ones(size(X))*stressH(4));
        [stressFit,stressOut,alignedSystems] = SSLIP(ebsd,stressData,stressSystems.',stressOpt);
        assert(isequal(size(alignedSystems),[3 1]) && isequal(stressOut.NoSs,[3 1]));
        assert(isequal(alignedSystems.CRSS,stressSystems.CRSS));
        expectedTensors = originalTensors;
        expectedTensors(:,:,1) = -expectedTensors(:,:,1);
        assert(isequal(alignedSystems.deformationTensor.matrix,expectedTensors));
        assert(isequal(stressSystems.deformationTensor.matrix,originalTensors));
        assert(max(abs(stressFit.prop.slipIDcor-stressGamma),[],'all') < 2e-5);
        fittedTensors = alignedSystems(stressOut.NoSs).deformationTensor.matrix;
        returnedA = [reshape(fittedTensors(1,1,:),1,[]); reshape(fittedTensors(1,2,:),1,[]); ...
            reshape(fittedTensors(2,1,:),1,[]); reshape(fittedTensors(2,2,:),1,[])];
        if normalized, returnedA = returnedA ./ vecnorm(returnedA); end
        reconstructed = returnedA*stressFit.prop.slipIDcor;
        if method == 1
            assert(all(stressFit.prop.solverExitFlag == 1));
            assert(max(abs(stressFit.prop.rotationIDcor-theta)) < 2e-5);
            reconstructed = reconstructed + rotationBasis*stressFit.prop.rotationIDcor';
        end
        assert(max(abs(reconstructed-stressH),[],'all') < 2e-5);
    end
end

% A vector3d is equivalent to uniaxial tension. Zero resolved shear must not
% erase a direction, and zero stress must leave every system unchanged.
stressOpt.IDMethod = 2; stressOpt.enableRotation = 0;
stressOpt.stress = vector3d(1,-2,0);
[vectorFit,~,vectorSystems] = SSLIP(ebsd,stressData,stressSystems,stressOpt);
stressOpt.stress = stressTensor.uniaxial(stressOpt.stress);
[tensorFit,~,tensorSystems] = SSLIP(ebsd,stressData,stressSystems,stressOpt);
assert(isequaln(vectorFit.prop,tensorFit.prop));
assert(isequal(vectorSystems.deformationTensor.matrix,tensorSystems.deformationTensor.matrix));
expectedTensors = originalTensors;
expectedTensors(:,:,[1 3]) = -expectedTensors(:,:,[1 3]);
assert(isequal(vectorSystems.deformationTensor.matrix,expectedTensors));
stressOpt.stress = stressTensor(zeros(3));
[~,~,zeroStressSystems] = SSLIP(ebsd,stressData,stressSystems,stressOpt);
assert(isequal(zeroStressSystems.deformationTensor.matrix,originalTensors));

% Signed methods must retain supplied signs even when a stress is present.
% Method 3 disables posConstr before solving, so alignment must not run there.
stressOpt.stress = stressTensor([2 -1 .4;-1 -2 .2;.4 .2 .5]);
stressOpt.IDMethod = 1; stressOpt.posConstr = 0;
[signedFit,~,signedSystems] = SSLIP(ebsd,stressData,stressSystems,stressOpt);
assert(isequal(signedSystems.deformationTensor.matrix,originalTensors));
assert(max(abs(signedFit.prop.slipIDcor-[.12;-.21]),[],'all') < 2e-5);
stressOpt.IDMethod = 3; stressOpt.posConstr = 1; stressOpt.threshResidual = 1;
[~,singleStressOut,singleStressSystems] = SSLIP(ebsd,stressData,stressSystems,stressOpt);
assert(singleStressOut.posConstr == 0);
assert(isequal(singleStressSystems.deformationTensor.matrix,originalTensors));

stressOpt.IDMethod = 2;
for invalidStress = {eye(3),vector3d(0,0,0),[xvector;yvector],stressTensor(cat(3,eye(3),eye(3)))}
    stressOpt.stress = invalidStress{1};
    caught = false;
    try
        SSLIP(ebsd,stressData,stressSystems,stressOpt);
    catch exception
        caught = any(strcmp(exception.identifier,{'SSLIP:InvalidStressType','SSLIP:InvalidStressValue'}));
    end
    assert(caught,'Invalid stress must be rejected for positive-constrained fits.');
end
stressOpt = rmfield(stressOpt,'stress');
lastwarn('');
[~,~,manualSystems] = SSLIP(ebsd,stressData,stressSystems,stressOpt);
[~,warningId] = lastwarn;
assert(strcmp(warningId,'SSLIP:MissingStress'));
assert(isequal(manualSystems.deformationTensor.matrix,originalTensors));
fprintf('PASS: stress alignment, returned basis, subsets, signed modes, rotation, zero shear, and input validation.\n');

% Preprocessing must preserve gradients and coordinate ordering on a shifted
% rectangular grid, including shuffled EBSD points and coarse-grained data.
[prepX,prepY] = meshgrid(2:2:20,-9:2:5);
prepU = 2 + .02*prepX - .03*prepY;
prepV = -1 + .04*prepX + .05*prepY;
prepGrid = dummyEBSDSimple(ori,prepX,prepY);
order = [2:2:numel(prepX),1:2:numel(prepX)]';
positions = vector3d(prepX(order),prepY(order),zeros(numel(order),1));
prepPoints = EBSD(positions,repmat(ori,numel(order),1), ...
    ones(numel(order),1),CS,struct());
for coarsegrain = [1 2]
    prepOpt = struct('filterSize',0,'coarsegrain',coarsegrain);
    [prepared,preparedGrid] = preprocessSSLIP(prepGrid, ...
        struct('U',prepU,'V',prepV),prepOpt);
    [shuffled,shuffledGrid] = preprocessSSLIP(prepPoints, ...
        struct('U',prepU(order),'V',prepV(order)),prepOpt);
    expectedSize = [8 10] / 2^(coarsegrain-1);
    assert(isequal(size(preparedGrid),expectedSize));
    assert(isequaln(prepared,shuffled));
    assert(isequal(preparedGrid.x,shuffledGrid.x) && isequal(preparedGrid.y,shuffledGrid.y));
    assert(max(abs(prepared.U-(2+.02*preparedGrid.x-.03*preparedGrid.y)),[],'all') < 1e-12);
    assert(max(abs(prepared.V-(-1+.04*preparedGrid.x+.05*preparedGrid.y)),[],'all') < 1e-12);
    assert(max(abs(prepared.Hxx-.02),[],'all') < 1e-12);
    assert(max(abs(prepared.Hxy+.03),[],'all') < 1e-12);
    assert(max(abs(prepared.Hyx-.04),[],'all') < 1e-12);
    assert(max(abs(prepared.Hyy-.05),[],'all') < 1e-12);
end

% Retain the existing zero/missing-displacement convention during filtering.
maskU = ones(size(prepX)); maskU(3,3) = 0;
maskV = ones(size(prepY)); maskV(4,4) = NaN;
prepOpt.coarsegrain = 1;
[unfiltered,~] = preprocessSSLIP(prepGrid,struct('U',maskU,'V',maskV),prepOpt);
assert(unfiltered.U(3,3) == 0 && isnan(unfiltered.V(4,4)));
prepOpt.filterSize = 1;
[filtered,~] = preprocessSSLIP(prepGrid,struct('U',maskU,'V',maskV),prepOpt);
assert(isequal(isnan(filtered.U),maskU == 0));
assert(isequal(isnan(filtered.V),isnan(maskV)));
assert(max(abs(filtered.U(isfinite(filtered.U))-1)) < 1e-12);
assert(max(abs(filtered.V(isfinite(filtered.V))-1)) < 1e-12);
% A mismatch in only one dimension must also be rejected.
caught = false;
try
    preprocessSSLIP(prepGrid,struct('U',prepU(:,1:7),'V',prepV),prepOpt);
catch exception
    caught = strcmp(exception.identifier,'SSLIP:DisplacementSizeMismatch');
end
assert(caught,'Displacement sizes must match the EBSD grid.');
fprintf('PASS: preprocessing grid order, affine gradients, coarse-graining, masks, and input sizes.\n');

% Ready-to-fit gradients follow EBSD point ordering and preserve physical
% zeros/NaNs without inventing displacement fields or reprocessing H.
gradientData = struct('Hxx',.001*prepX+.002*prepY,'Hxy',zeros(size(prepX)), ...
    'Hyx',ones(size(prepX))*.03,'Hyy',.002*prepX);
gradientData.Hxx(3,3) = NaN;
gradientFields = {'Hxx','Hxy','Hyx','Hyy'};
shuffledData = struct;
for k = 1:4
    name = gradientFields{k};
    shuffledData.(name) = gradientData.(name)(order);
end
readyOpt = struct('filterSize',0,'coarsegrain',1);
[ready,readyGrid] = preprocessSSLIP(prepGrid,gradientData,readyOpt);
[shuffledReady,shuffledReadyGrid] = preprocessSSLIP(prepPoints,shuffledData,readyOpt);
assert(isequaln(ready,shuffledReady));
assert(isequal(readyGrid.x,shuffledReadyGrid.x) && isequal(readyGrid.y,shuffledReadyGrid.y));
assert(~isfield(ready,'U') && ~isfield(readyGrid.prop,'V'));
assert(isequal(ready.Hxy,zeros(size(prepX))) && isnan(ready.Hxx(3,3)));
for k = 1:4
    name = gradientFields{k};
    assert(isequaln(ready.(name),gradientData.(name)));
end
badGradientData = {rmfield(gradientData,'Hyy'),gradientData,gradientData,gradientData,gradientData};
badGradientData{2}.Hyy = zeros(2,3);
badGradientData{3}.U = prepU;
badGradientData{4}.Hxx(1) = Inf;
badGradientData{5}.Hxy(1) = 1i;
expectedErrors = {'SSLIP:MissingGradients','SSLIP:GradientSizeMismatch', ...
    'SSLIP:AmbiguousDeformationData','SSLIP:InvalidGradients','SSLIP:InvalidGradients'};
for k = 1:numel(badGradientData)
    caught = false;
    try
        preprocessSSLIP(prepGrid,badGradientData{k},readyOpt);
    catch exception
        caught = strcmp(exception.identifier,expectedErrors{k});
    end
    assert(caught,'Invalid gradient input must be rejected.');
end
for settingsPair = [1 0;1 2]
    rejectedOpt = struct('filterSize',settingsPair(1),'coarsegrain',settingsPair(2));
    caught = false;
    try
        preprocessSSLIP(prepGrid,gradientData,rejectedOpt);
    catch exception
        caught = strcmp(exception.identifier,'SSLIP:GradientPreprocessing');
    end
    assert(caught,'Ready-to-fit gradients must not be filtered or coarse-grained again.');
end

% Four-input SSLIP uses the same physical/rotation solve and keeps the second
% output as options. Gradient-only plots contain five available fields.
gradientH = physicalA(:,2)*.2 + rotationBasis*(-.03);
fitData = struct;
for k = 1:4, fitData.(gradientFields{k}) = ones(size(X))*gradientH(k); end
gradientOpt = struct('NoSs',2,'minEeff',0,'threshResidual',1e-6, ...
    'enableRotation',1,'plotSSLIP',0,'plotDefGrad',1,'saveFig',0,'cmap',parula(256));
[gradientFit,gradientOut] = SSLIP(ebsd,fitData,sRotation,gradientOpt);
assert(gradientOut.filterSize==0 && gradientOut.coarsegrain==1 && gradientOut.NoSs==2);
assert(max(abs(gradientFit.prop.slipIDcor-.2),[],'all') < 2e-5);
assert(max(abs(gradientFit.prop.rotationIDcor+.03)) < 2e-5);
assert(all(gradientFit.prop.solverExitFlag==1) && ~isfield(gradientFit.prop,'U'));
assert(numel(findall(gcf,'Type','axes'))==5);
close all;
for k = 1:4, fitData.(gradientFields{k}) = zeros(size(X)); end
gradientOpt.enableRotation = 0;
gradientOpt.IDMethod = 3;
[zeroGradientFit,~] = SSLIP(gradientFit,fitData,sRotation,gradientOpt);
assert(all(zeroGradientFit.prop.slipIDcor==0,'all'));
assert(~isfield(zeroGradientFit.prop,'rotationIDcor') && ~isfield(zeroGradientFit.prop,'solverExitFlag'));
assert(numel(findall(gcf,'Type','axes'))==5);
close all;
fprintf('PASS: direct gradients, grid alignment, validation, rotation, missing displacements, and zero-field plots.\n');

opt.IDMethod = 3;
opt.posConstr = 0;
opt.filterSize = 0;
opt.coarsegrain = 1;
opt.plotSSLIP = 0;
opt.plotDefGrad = 0;
opt.saveFig = 0;
opt.threshResidualFraction = .001;
[ebsdFloor,floorOpt] = SSLIP(ebsd,.02*Y,.005*Y,sS,opt);
assert(max(abs(ebsdFloor.prop.slipIDcor(:)-.02)) < 1e-10);
[structuredFloor,structuredOpt] = SSLIP(ebsd,struct('U',.02*Y,'V',.005*Y),sS,opt);
assert(isequaln(structuredFloor.prop,ebsdFloor.prop) && isequaln(structuredOpt,floorOpt));
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

% The extracted single-slip solver uses the supplied system order and reports
% the positive-constraint override through the unchanged public SSLIP call.
sOrdered = [sTwo(1); sRotation(2); sTwo(2)];
singleOpt = opt;
singleOpt.NoSs = [3 1];
singleOpt.posConstr = 1;
singleOpt.singleSlipPerPixel = 0;
singleOpt.threshResidual = .3;
[ordered,orderedOpt] = SSLIP(ebsd,.2*Y,-.3*X,sOrdered,singleOpt);
assert(isequal(size(ordered.prop.slipIDcor),[2 16]));
assert(max(abs(ordered.prop.slipIDcor-[-.3;.2]),[],'all') < 1e-12);
assert(isequal(size(ordered.prop.residualEeff),[2 16]));
assert(isequal(orderedOpt.NoSs,[3 1]) && orderedOpt.posConstr == 0);
assert(singleOpt.posConstr == 1 && endsWith(orderedOpt.plotname,'_singleSlip'));
singleOpt.posConstr = 0;
singleOpt.singleSlipPerPixel = 1;
[tiedGamma,~] = solveSSLIP_SingleSlip(sOrdered(singleOpt.NoSs), ...
    zeros(2,3),ones(2,3)*.2,ones(2,3)*.2,zeros(2,3),singleOpt);
assert(max(abs(tiedGamma(1,:)-.2)) < 1e-12);
assert(all(tiedGamma(2,:) == 0));

% Normalized coefficients retain their original meaning. The helper counts
% all pixels of matrix input and does not apply method 1's minEeff cutoff.
singleCfg = struct('posConstr',0,'threshResidual',1e-9,'minEeff',1);
singleH = physicalA(:,1)*.2;
for normalized = 0:1
    singleCfg.normalizeInplane = normalized;
    [singleGamma,singleResidual,executedCfg] = solveSSLIP_SingleSlip(sRotation(1), ...
        repmat(singleH(1),2,3),repmat(singleH(2),2,3), ...
        repmat(singleH(3),2,3),repmat(singleH(4),2,3),singleCfg);
    expectedGamma = .2 * norm(physicalA(:,1))^normalized;
    assert(isequal(size(singleGamma),[1 6]) && isequal(size(singleResidual),[1 6]));
    assert(max(abs(singleGamma-expectedGamma)) < 1e-12);
    assert(max(abs(singleResidual)) < 1e-12 && isequal(executedCfg,singleCfg));
end
singleCfg.normalizeInplane = 0;
[lowOrMissing,resLowOrMissing] = solveSSLIP_SingleSlip(sS,[0 0],[.001 NaN],[0 0],[0 0],singleCfg);
assert(isequal(lowOrMissing,[.001 0]));
assert(resLowOrMissing(1) == 0 && isnan(resLowOrMissing(2)));
singleCfg.threshResidual = calcEffectiveE(1,0,0,0);
[boundaryGamma,boundaryResidual] = solveSSLIP_SingleSlip(sS,1,.2,0,0,singleCfg);
assert(boundaryGamma == 0 && boundaryResidual == singleCfg.threshResidual);
singleCfg.enableRotation = 1;
caught = false;
try
    solveSSLIP_SingleSlip(sS,0,.2,0,0,singleCfg);
catch exception
    caught = strcmp(exception.identifier,'SSLIP:RotationRequiresMethod1');
end
assert(caught,'The single-slip helper must reject unsupported rotation.');
fprintf('PASS: single-slip system order, options, ties, normalization, matrix inputs, and rejection semantics.\n');

% Rotation extraction must preserve physical channels and original NoSs.
rotationOpt.minEeff = 0;
rotationOpt.IDMethod = 1;
rotationOpt.filterSize = 0;
rotationOpt.coarsegrain = 1;
rotationOpt.plotSSLIP = 0;
rotationOpt.plotDefGrad = 0;
rotationOpt.saveFig = 0;
rotationOpt.NoSs = 2;
pipelineH = physicalA(:,2)*.2 + rotationBasis*(-.03);
[ebsdRotation,optRotation] = SSLIP(ebsd, ...
    pipelineH(1)*X+pipelineH(2)*Y, pipelineH(3)*X+pipelineH(4)*Y, ...
    sRotation,rotationOpt);
assert(isequal(size(ebsdRotation.prop.slipIDcor),[1 16]));
assert(max(abs(ebsdRotation.prop.slipIDcor(:)-.2)) < 2e-5);
assert(isequal(size(ebsdRotation.prop.rotationIDcor),[16 1]));
assert(max(abs(ebsdRotation.prop.rotationIDcor(:)+.03)) < 2e-5);
assert(all(ebsdRotation.prop.solverExitFlag == 1) && optRotation.NoSs == 2);
for unsupportedMethod = [2 3]
    rotationOpt.IDMethod = unsupportedMethod;
    caught = false;
    try
        SSLIP(ebsd,.2*Y,.1*X,sRotation,rotationOpt);
    catch exception
        caught = strcmp(exception.identifier,'SSLIP:RotationRequiresMethod1');
    end
    assert(caught,'Unsupported methods must reject rotation before fitting.');
end
assert(~isfield(ebsdFirst.prop,'rotationIDcor'));
fprintf('PASS: rotation output mapping, selected-system labels, and unsupported-method guards.\n');

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

% The residual helper must receive the activity limits and the explicit
% residual argument, even if the EBSD object contains an older residual.
staleResidual = ebsd;
staleResidual.prop.residualEeff = ones(16,1)*9;
popt.residualScaleSame = 1;
popt.logscale = 1;
popt.logmin = .01;
plotSSLIP(linspace(.01,.2,16),ones(16,1)*.03,staleResidual,sS,popt);
assert(isequal(gca().CLim,[.01 1]) && strcmp(gca().ColorScale,'log'));
assert(strcmp(gca().Title.String,'residual Eeff,mean=0.03'));
close all;
popt = rmfield(popt,{'residualScaleSame','logmin'});
popt.logscale = 0;

% The separated plotters can use saved fields without running a solver.
% Check custom deformation limits and the existing export filenames.
plotDirectory = tempname;
mkdir(plotDirectory);
oldDirectory = pwd;
directoryCleanup = onCleanup(@() cd(oldDirectory));
cd(plotDirectory);
dopt = struct('cmap',parula(256),'logscale',1,'logmin',.001, ...
    'maxE',.1,'DefGradLim',[-.07 .07],'saveFig',1,'saveExt','.jpg', ...
    'coarsegrain',1,'filterSize',0,'casename','check','comment','');
plotSSLIP_DeformationFields(preparedGrid,dopt);
axesList = findall(gcf,'Type','axes');
effectiveAxis = axesList(arrayfun(@(a) strcmp(a.Title.String,'E_{eff}'),axesList));
gradientAxes = axesList(arrayfun(@(a) startsWith(a.Title.String,'H_{'),axesList));
assert(isscalar(effectiveAxis) && isequal(effectiveAxis.CLim,[.001 .1]));
assert(strcmp(effectiveAxis.ColorScale,'log') && numel(gradientAxes)==4);
assert(all(arrayfun(@(a) isequal(a.CLim,[-.07 .07]),gradientAxes)));
assert(isfile('SSLIP_CGR_1_Filt_0_check___gradients.png'));
assert(isfile('ssAnalysis_CoarseGr_1_Filt_0_check___disp_grad_tensor.jpg'));
close all;

ropt = struct('cmap',parula(256),'saveFig',1,'saveExt','.jpg','plotname','check');
storedResult = ebsdRotation;
storedResult.prop.residualEeff = ones(16,1)*.02;
plotSSLIP_Residual(storedResult,ropt);
assert(isequal(gca().CLim,[0 .02]));
assert(isfile('check_Residual.png') && isfile('check_Residual.jpg'));
close all;
storedResult.prop.rotationIDcor = ones(16,1)*(-.03);
storedResult.prop.rotationIDcor(end) = NaN;
originalRotation = storedResult.prop.rotationIDcor;
rotationFigure = plotSSLIP_Rotation(storedResult,ropt);
assert(isgraphics(rotationFigure,'figure') && isequal(rotationFigure,gcf));
assert(max(abs(gca().CLim-[-.03 .03]/degree)) < 1e-12);
assert(isequal(gca().Colormap,blue2redColorMap));
assert(strcmp(gca().Title.String,'Inferred rotation correction'));
assert(isequaln(storedResult.prop.rotationIDcor,originalRotation));
assert(isfile('check_rotation.png'));
close all;
assert(isempty(plotSSLIP_Rotation(ebsd,ropt)));
assert(isempty(findall(groot,'Type','figure')));
cd(oldDirectory);
clear directoryCleanup;
fprintf('PASS: separate deformation, residual, and rotation plots, scales, units, and exports.\n');

% Plot the selected physical system from a rotation-corrected fit. The
% original system index is valid because SSLIP retains the full input list.
popt.NoSs = 2;
plotSSLIP(ebsdRotation.prop.slipIDcor,ebsdRotation.prop.residualEeff, ...
    ebsdRotation,sRotation,popt);
close all;
popt.NoSs = 1;

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
