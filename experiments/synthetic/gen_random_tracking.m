function problem = gen_random_tracking(problem)
%% Generate tracking problem from random data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Problem variables
N = problem.N_VAR;
K = problem.K;
L = problem.L;

intraRadius = problem.intraRadius;
velocityBound = problem.velocityBound;
translationBound = problem.translationBound;
noiseSigmaSqrt = problem.noiseSigmaSqrt;
dt = problem.dt;

accNoiseSigmaSqrt = problem.accelerationNoiseBoundSqrt;
rotNoiseBound = problem.rotationNoiseBound;

% TODO: INCORPORATE MAX POSITION BOUNDS
% if (problem.velocityBound*L*dt >= problem.translationBound)
%     warning("Object may leave translationBound during motion. " + ...
%         "Please set velocityBound*L*dt <= translationBound.")
% end

% Weights!
noiseBoundSq = problem.noiseBound^2;
% problem.covar_measure = ones(N,L)*(noiseBoundSq/9);
% problem.covar_velocity = ones(L-2,1)*problem.covar_measure(1)*1;
% problem.kappa_rotrate  = ones(L-2,1)*(2/problem.covar_velocity(1));
% notes on covariances:
% - covar_measure is scaled by noise bound is the bound is 3 standard
%   deviations
% - covar_velocity can be adjusted based on how linear you believe the
%   system is. In simulated data we are fairly confident.
% - kappa_rotrate: see covar_velocity. Note that covar ~ 1/(2*kappa)
% In this case we scale both by covar_measure.

%% generate a mean shape centered at zero
% allow override by specifying in problem struct
if ~isfield(problem,'B')

    mean_shape = randn(3,N);
    mean_shape = mean_shape - mean(mean_shape,2);
    
    shapes = zeros(3,N,K);
    for k = 1:K
        shapes(:,:,k) = mean_shape + intraRadius * randn(3,N);
    end
    B = reshape(shapes, 3*N, K);

    if isfield(problem,'true_rand')
    shapes = randn(3,N,K);
    B = reshape(shapes, 3*N,K);
    end

else
    B = problem.B;
    shapes = problem.shapes;
end

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
    v_gt = repmat(v_0,1,1,L-1);
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
    dR_gt = repmat(dR_gt,1,1,L-1);
else
    dR_gt = problem.dR_gt;
end
%% generate noise-y measurements
shape = reshape(B*c_gt,3,N);
y = zeros(3*N,L);

for l = 1:L
    R = R_gt(:,:,l);
    p = p_gt(:,:,l);

    scene = R * shape + p + noiseSigmaSqrt * randn(3,N);
    y(:,l) = reshape(scene, 3*N, 1);

    if (l < L)
        % Apply random (bounded) acceleration (TODO: make actual bound)
        dR_gt(:,:,l) = dR_gt(:,:,l) * rand_rotation('RotationBound',rotNoiseBound);
        v_gt(:,:,l) = v_gt(:,:,l) + accNoiseSigmaSqrt*randn(3,1)*dt;

        if strcmp(problem.velprior, "body") || strcmp(problem.velprior, "body-sym")
            dR = dR_gt(:,:,l);
            v = v_gt(:,:,l);
            p = p_gt(:,:,l);
        
            % spiral to next
            pts = get_spiral_pts(R, dR, v, p, dt, 2);
            % pts = sim_dynamics(R, dR, v, p, dt, 20, false);
            
            p_gt(:,:,l+1) = pts(:,end);
            R_gt(:,:,l+1) = R * dR_gt(:,:,l);
        elseif strcmp(problem.velprior, "world")
            R_gt(:,:,l+1) = R * dR_gt(:,:,l);
            p_gt(:,:,l+1) = R_gt(:,:,l+1)*(R'*p + v_gt(:,:,l) * dt);
        elseif strcmp(problem.velprior, "grav-world")
            error("Selected prior is not implemented")
        else
            error("Selected prior is not implemented")
        end
    end
end

%% generate outliers
nrOutliers = round(N*L*problem.outlierRatio);
if nrOutliers > 0
    fprintf('Tracking: generate %d outliers.\n',nrOutliers);

    outlierIDs = randperm(N*L,nrOutliers);
    outliers = randn(3,nrOutliers);
    
    curOutIdx = 1;
    for id = outlierIDs
        l = floor((id-1)/N) + 1;
        i = id - (l-1)*N;
        y(ib3(i),l) = outliers(curOutIdx);
        curOutIdx = curOutIdx + 1;
    end
else
    outlierIDs = [];
end

theta_gt = ones(N*L,1);
theta_gt(outlierIDs) = -1;
outliers_gt = outlierIDs;
inliers_gt = [];
for i = 1:(N*L)
    if ~ismember(i,outliers_gt)
        inliers_gt(end+1) = i;
    end
end

%% Save
if strcmp(problem.velprior, "body") || strcmp(problem.velprior, "body-sym")
    % convert p to s, sh for saving
    s_gt  = zeros(3,1,L);
    sh_gt = zeros(3,1,L-1);
    for l = 1:L
        s_gt(:,:,l) = R_gt(:,:,l)'*p_gt(:,:,l);
        if (l > 1)
            sh_gt(:,:,l-1) = dR_gt(:,:,l-1)*s_gt(:,:,l);
        end
    end
    problem.s_gt = s_gt;
    problem.sh_gt = sh_gt;
end
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

problem.noiseBoundSq = noiseBoundSq;
problem.noiseBound = sqrt(problem.noiseBoundSq);

problem.theta_gt = theta_gt;
problem.inliers_gt = inliers_gt;
problem.outliers_gt = outliers_gt;

% also save a "ground truth vector" that we can quickly compare
if strcmp(problem.velprior, "body") || strcmp(problem.velprior, "body-sym")
    for l = 1:L
        R_cur = problem.R_gt(:,:,l);
        if (l < L)
            dR_cur = problem.dR_gt(:,:,l);
        end
    end
    x_gt = [reshape(problem.R_gt, problem.L*9,1,1);
            reshape(problem.dR_gt,(problem.L-1)*9,1,1); 
            reshape(problem.s_gt,problem.L*3,1,1); 
            reshape(problem.v_gt,(problem.L-1)*3,1,1)];
    problem.x_gt = x_gt;
elseif strcmp(problem.velprior, "world")
    rh_gt = zeros(9*(L-1), 1);
    s_gt = zeros(3*L,1);
    for l = 1:L
        R_cur = problem.R_gt(:,:,l);
        if (l < L)
            dR_cur = problem.dR_gt(:,:,l);
            rh_gt((9*(l-1)+1):(9*l)) = reshape(R_cur * dR_cur,9,1);
        end
        s_gt((3*(l-1)+1):(3*l)) = R_cur' * problem.p_gt(:,:,l);
    end
    x_gt = [reshape(problem.R_gt, problem.L*9,1,1);
            reshape(problem.dR_gt,(problem.L-1)*9,1,1);
            rh_gt; 
            reshape(problem.p_gt,problem.L*3,1,1); 
            s_gt];
    problem.x_gt = x_gt;
elseif strcmp(problem.velprior, "grav-world")
    error("Selected prior is not implemented")
else
    error("Selected prior is not implemented")
end

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