function problem = gen_random_tracking_old(problem)
%% Generate tracking problem from random data
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

% TODO: INCORPORATE MAX POSITION BOUNDS
% if (problem.velocityBound*L*dt >= problem.translationBound)
%     warning("Object may leave translationBound during motion. " + ...
%         "Please set velocityBound*L*dt <= translationBound.")
% end

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
    v_gt = randn(3,1);
    v_gt = v_gt/norm(v_gt);
    v_gt = velocityBound*rand*v_gt;
else
    v_gt = problem.v_gt;
end

if ~isfield(problem,'p_gt')
    p_gt = randn(3,1);
    p_gt = p_gt/norm(p_gt);
    p_gt = translationBound*rand*p_gt;
else
    p_gt = problem.p_gt;
end

if ~isfield(problem,'R_gt')
    R_gt = rand_rotation;
else
    R_gt = problem.R_gt;
end
if ~isfield(problem,'dR_gt')
    dR_gt = rand_rotation;
else
    dR_gt = problem.dR_gt;
end

%% generate noise-y measurements
shape = reshape(B*c_gt,3,N);
y = zeros(3*N,L);

R = R_gt;
p = p_gt;
for l = 1:L
    scene = R * shape + p + noiseSigma * randn(3,N);
    y(:,l) = reshape(scene, 3*N, 1);

    R = R * dR_gt;
    p = p + v_gt * dt;
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