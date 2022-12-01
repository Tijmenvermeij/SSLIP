function [slipID,residualEeff] = SSLIPConstr(sS,Hxx,Hxy,Hyx,Hyy,options)
    
%%% function to perform slip id analysis using def grad tensor component
%%% solving. Only uses constraints

%%% Tijmen Vermeij / TUe / t.vermeij@tue.nl

if nargin < 6
    options = struct;
end

%% set options if not given
% normalize the theoretucal disp grad tensor towards the "in-plane
% configuration"
if ~isfield(options,'normalizeInplane')
    options.normalizeInplane = 0;
end

% get theoretical disp grad tensor components from the slip systems:
Hslip = sS.deformationTensor.matrix;

% extract the theoretical 2D components
Hslip11 = reshape( Hslip(1,1,:) , 1, []);
Hslip12 = reshape( Hslip(1,2,:) , 1, []);
Hslip21 = reshape( Hslip(2,1,:) , 1, []);
Hslip22 = reshape( Hslip(2,2,:) , 1, []);

% put the theoretical components in matrix form for solving with backslash
% operator (A*x = B) (rows: different components, columns: different slip systems)
A = [Hslip11
    Hslip12
    Hslip21
    Hslip22];

% if necessary, normalize inplane
if options.normalizeInplane
    A = A./sqrt(sum(A.^2,1));
end

% put the experimental disp grad components in a matrix (rows: different
% components, cols: different positions in map (all points))
B = [Hxx(:)'
    Hxy(:)'
    Hyx(:)'
    Hyy(:)'];

% use a least square solver or backslash operator for solving
if options.posConstr
    %%% positive constraint of slip amplitudes, use lsqlin.
    
    % Convert A and B to a form that is suitable for "lsqlin" function (which can handle the positive constraint).
    % Lsqlin cannot solve multiple systems like the backslash operator (i.e. calculate slip at all points at once), so we need to make it into one system

    % number of points in map
    n = length(Hxx(:));
    % no of slipsystems
    N = length(sS);

    % repeat the A matrix (which is the same for all points) over the diagonal,
    % n times, such that the system can be solved at once (but still for every point separately)
    % Matrix will be very large, so use sparse matrix
    sparseA = sparse(repmat(A,1,n));
    sparseA = mat2cell(sparseA,size(A,1),repmat(size(A,2),1,n));
    sparseA = sparse(blkdiag(sparseA{:}));

    % change the B matrix such that it corresponds to the A matrix; i.e. stack
    % everything in one column
    
    % convert NaNs to 0 to avoid problems (will be replaced later)
    nanB = isnan(B(1,:));
    B(isnan(B)) = 0;
    
    colB=B(:);

    % solve with lsqlin

    % define positive constraints by setting the lower bound of the solution to zero
    constr = zeros(size(sparseA,2),1);
    lsqlinoptions = optimoptions('lsqlin','Algorithm','interior-point','Display','iter-detailed','OptimalityTolerance',1e-15);

    %% solve with lsqlin
    x = lsqlin(sparseA,colB,[],[],[],[],constr,[],[],lsqlinoptions);

    % reshape result to matrix for easier plotting etc
    x = reshape(x,N,n);
    
    % input NaNs again
    x(:,nanB) = NaN;
else
    % without constraints, use backslash operator. significantly faster
    x = A\B;
end

% calculate residual displacement gradient tensor components and from that calc an "effective
% shear strain residual"
HCalc = A*x;
residualH = B - HCalc;
residualEeff = calcEffectiveE(residualH(1,:),residualH(2,:),residualH(3,:),residualH(4,:));

% store solutions correctly
residualEeff = residualEeff';

slipID = x;

end

