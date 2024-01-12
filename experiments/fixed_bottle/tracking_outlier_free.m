%% Dense SDP relaxation for certifiable tracking
%  Generic, tunable script to run one iteration of dense tracking.
%    Operates on random data with no outlier support.
%    Run setup.m once to set up paths.
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Generate random tracking problem
problem.bag = "2024-01-11-17-30-52.bag";
problem.batchsize = 10; % number of visible keyframes to solve over

% Set bounds based on problem setting
problem.dt = 1/30; % set based on camera rate
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.noiseBound = 0.2;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% regen if batch size changes.
problem.regen_sdp = true; % when in doubt, set to true

% add shape, measurements, outliers
load("cad_frame.mat");
problem.shapes = annotatedPointsWrtTarget'; % 3 x N x K
problem_list = bag2problem(problem);

%% Solve for each batch
for j = 1:length(problem_list)
curproblem = problem_list(j);

soln = solve_weighted_tracking(curproblem);

% soln_pace = pace_with_EKF(problem);

% soln = solve_full_tracking(problem,lambda);
% Ap = solve_nopos_tracking(problem);
pause
end

%% Check solutions
% eigenvalue plot
L = problem.L;
figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
hold on

% if strcmp(problem.velprior, "body")
%     slices = 1:(1+9*(3*L-2)+3*L);
%     Xopt_pRemoved = soln.raw.Xopt{1}(slices, slices);
%     bar([zeros(3*(L-1),1);eig(Xopt_pRemoved)]);
%     title("Eigenvalues of Relaxed Solution")
% elseif strcmp(problem.velprior, "world")
%     slices = [1:(1+9*(3*L-2)),(1+9*(3*L-2)+3*L+1):(9*(3*L-2)+6*L)];
%     Xopt_pRemoved = soln.raw.Xopt{1}(slices, slices);
%     bar([zeros(3*L+1,1);eig(Xopt_pRemoved)]);
%     title("Eigenvalues of Relaxed Solution")
% elseif strcmp(problem.velprior, "grav-world")
%     error("Selected prior is not implemented")
% else
%     error("Selected prior is not implemented")
% end
hold off

% raw error
% x_err = norm(problem.x_gt - soln.x_est);

% projected errors
R_err = zeros(L,1);
dR_err = zeros(L-1,1);
p_err = zeros(L,1);
v_err = zeros(L-1,1);
for l = 1:L
    % R
    R_err(l) = getAngularError(problem.R_gt(:,:,l), soln.R_est(:,:,l));
    % dR
    if (l < L)
        dR_err(l) = getAngularError(problem.dR_gt(:,:,l), soln.dR_est(:,:,l));
    end
    % p
    p_err(l) = norm(problem.p_gt(:,:,l) - soln.p_est(:,:,l));
    % v
    if (l < L)
        v_err(l) = norm(problem.v_gt(:,:,l) - soln.v_est(:,:,l));
    end
end

% shape error
c_err = norm(problem.c_gt - soln.c_est);

% PACE errors
norm(problem.p_gt - soln_pace.p_raw,'fro')
norm(problem.p_gt - soln_pace.p_smoothed,'fro')
norm(problem.p_gt - soln.p_est,'fro')

% Plot trajectory!
plot_trajectory(problem,soln)
