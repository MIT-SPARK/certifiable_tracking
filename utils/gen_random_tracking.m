function problem = gen_random_tracking(problem)
%% Generate tracking problem from random data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Problem variables
N = problem.N;
K = problem.K;
L = problem.L;

intraRadius = problem.intraRadius;
velocityBound = problem.velocityBound;
translationBound = problem.translationBound;
noiseSigma = problem.noiseSigma;
dt = problem.dt;

accNoiseBound = problem.accelerationNoiseBound;
rotNoiseBound = problem.rotationNoiseBound;

% TODO: INCORPORATE MAX POSITION BOUNDS
% if (problem.velocityBound*L*dt >= problem.translationBound)
%     warning("Object may leave translationBound during motion. " + ...
%         "Please set velocityBound*L*dt <= translationBound.")
% end

% Weights!
% TODO: make this more random?
problem.weights = ones(N,L);
% no support for different x,y,z weights (repeat each weight 3 times)
problem.weights_position = ones(3*(L-1),1);
problem.weights_velocity = ones(3*(L-1),1);
problem.weights_rotation = ones(3*(L-1),1);
problem.weights_rotrate  = ones(3*(L-1),1);

%% generate a mean shape centered at zero
mean_shape = randn(3,N);
mean_shape = mean_shape - mean(mean_shape,2);

shapes = zeros(3,N,K);
for k = 1:K
    shapes(:,:,k) = mean_shape + intraRadius * randn(3,N);
end
B = reshape(shapes, 3*N, K);

%% ground truth c, v, p, R, dR
% allow override by specifying in problem struct
if ~isfield(problem,'c_gt')
    c_gt = rand(K,1);
    c_gt = c_gt/sum(c_gt); % sum to one
else
    c_gt = problem.c_gt;
end

if ~isfield(problem,'v_gt')
    v_0 = randn(3,1);
    v_0 = v_0/norm(v_0);
    v_0 = velocityBound*rand*v_0;
    v_gt = repmat(v_0,1,1,L);
else
    v_gt = problem.v_gt;
end

if ~isfield(problem,'p_gt')
    p_gt = zeros(3,1,L);
    p_0 = randn(3,1);
    p_0 = p_0/norm(p_0);
    p_0 = translationBound*rand*p_0;
    p_gt(:,:,1) = p_0;
else
    p_gt = problem.p_gt;
end

if ~isfield(problem,'R_gt')
    R_gt = zeros(3,3,L);
    R_gt(:,:,1) = rand_rotation;
else
    R_gt = problem.R_gt;
end
if ~isfield(problem,'dR_gt')
    dR_gt = rand_rotation;
    dR_gt = repmat(dR_gt,1,1,L);
else
    dR_gt = problem.dR_gt;
end

%% generate noise-y measurements
shape = reshape(B*c_gt,3,N);
y = zeros(3*N,L);

for l = 1:L
    R = R_gt(:,:,l);
    p = p_gt(:,:,l);

    scene = R * shape + p + noiseSigma * randn(3,N);
    y(:,l) = reshape(scene, 3*N, 1);

    if (l < L)
        % Apply random (bounded) acceleration (TODO: make actual bound)
        dR_gt(:,:,l) = dR_gt(:,:,l) * rand_rotation('RotationBound',rotNoiseBound);
        v_gt(:,:,l) = v_gt(:,:,l) + accNoiseBound*randn(3,1)*dt;

        R_gt(:,:,l+1) = R * dR_gt(:,:,l);
        p_gt(:,:,l+1) = p + v_gt(:,:,l) * dt;
    end
end

%% generate outliers
% nrOutliers = round(N*outlierRatio);
% if nrOutliers > 0
%     fprintf('Category registration: generate %d outliers.\n',nrOutliers);
%     outlierIDs = N-nrOutliers+1:N;
%     outliers = randn(3,nrOutliers);
%     scene(:,outlierIDs) = outliers;
% else
%     outlierIDs = [];
% end
% 
% theta_gt = ones(N,1);
% theta_gt(outlierIDs) = -1;

%% Save
problem.type = 'tracking';
problem.B = B;
problem.shapes = shapes;
problem.y = y;
% problem.nrOutliers = nrOutliers;
% problem.outlierIDs = outlierIDs;
% problem.theta_gt = theta_gt;
problem.c_gt = c_gt;

problem.p_gt = p_gt;
problem.v_gt = v_gt;
problem.R_gt = R_gt;
problem.dR_gt = dR_gt;

problem.cBound = 1.0;

noiseBoundSq = max(4e-2, noiseSigma^2 * chi2inv(0.99,3));
problem.noiseBoundSq = noiseBoundSq;
problem.noiseBound = sqrt(problem.noiseBoundSq);

% also save a "ground truth vector" that we can quickly compare
rh_gt = zeros(9*(L-1), 1);
s_gt = zeros(3*L,1);
for l = 1:L
    R_cur = problem.R_gt(:,:,l);
    dR_cur = problem.dR_gt(:,:,l);
    if (l < L)
        rh_gt((9*(l-1)+1):(9*l)) = reshape(R_cur * dR_cur,9,1);
    end
    s_gt((3*(l-1)+1):(3*l)) = R_cur' * problem.p_gt(:,:,l);
end
x_gt = [reshape(problem.R_gt, problem.L*9,1,1);
        reshape(problem.dR_gt,problem.L*9,1,1);
        rh_gt; 
        reshape(problem.p_gt,problem.L*3,1,1); 
        s_gt];
problem.x_gt = x_gt;

end

%% Helper function for rotations
function R = rand_rotation(varargin)
% generate random rotation matrix

params = inputParser;
params.CaseSensitive = false;

params.addParameter('RotationBound',2*pi,...
    @(x) 0.0<=x && x<=2*pi);

params.parse(varargin{:});

RotationBound = params.Results.RotationBound;

angle = RotationBound*rand - RotationBound/2;
axis  = randn(3,1);
axis  = axis / norm(axis);
R     = axang2rotm([axis' angle]);
end