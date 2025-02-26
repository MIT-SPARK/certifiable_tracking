%% Certifiable Tracking, Outlier-Free Case
%  Generic, tunable script to run one iteration of dense tracking.
%    Operates on random data with no outlier support.
%    Run setup.m once to set up paths.
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
rng("default")

%% Generate random tracking problem
problem.N_VAR = 10; % nr of keypoints
problem.K = 3; % nr of shapes
problem.L = 8; % nr of keyframes in horizon

problem.outlierRatio = 0.0; % TODO: no support for outliers
problem.noiseSigmaSqrt = 0.1; % [m]
problem.covar_measure_base = 0.05^2;
problem.covar_velocity_base = 0.05^2;
problem.covar_rotrate_base = 0.05^2;

problem.noiseBound = 0.3; %chi2inv(0.95,3*problem.N_VAR*problem.L)*problem.noiseSigmaSqrt^2;
problem.processNoise = 0.01;

problem.intraRadius = 0.2; 
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity

problem.accelerationNoiseBoundSqrt = 0;%0.6;
problem.rotationNoiseBound = 0;%pi/32; % rad

% regen if pbound, vbound, N, L, K change.
problem.regen_sdp = true; % when in doubt, set to true

% Optional: use a specified velocity trajectory
% problem = make_trajectory(problem);
% problem.dR_gt = repmat(axang2rotm([0,0,1,pi/2]),[1,1,problem.L-1]);
% problem.R_gt = repmat(eye(3,3),[1,1,problem.L]);
% problem.dR_gt = repmat(axang2rotm([0,0,1,1]),[1,1,problem.L-1]);
% problem.v_gt = repmat([1;0;1],[1,1,problem.L-1]);

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 0;
problem.lambda = lambda;

%% Solve!
soln = solve_weighted_tracking(problem);

pace = pace_raw(problem);
% paceukf = pace_py_UKF(problem,pace);
% paceekf = pace_ekf2(problem,pace);


%% Check solutions
% eigenvalue plot
L = problem.L;
figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
title("Eigenvalues of Relaxed Solution")

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

% Plot trajectory!
plot_trajectory2(problem,soln)

compare(problem, soln, pace, paceukf, paceekf);

function compare(gt, ours, pace, paceukf, paceekf)
L = gt.L;
% compare position
epace.p = norm(gt.p_gt(:,:,end) - pace.p(:,:,end));
eukf.p = norm(gt.p_gt(:,:,end) - paceukf.p(:,:,end));
eekf.p = norm(gt.p_gt(:,:,end) - paceekf.p(:,:,end));
eours.p = norm(gt.p_gt(:,:,end) - ours.p_est(:,:,end));

% compare rotation
epace.R = zeros(L,1);
eukf.R = zeros(L,1);
eours.R = zeros(L,1);
for l = 1:L
    epace.R(l) = getAngularError(gt.R_gt(:,:,l), pace.R(:,:,l));
    eukf.R(l) = getAngularError(gt.R_gt(:,:,l), paceukf.R(:,:,l));
    eours.R(l) = getAngularError(gt.R_gt(:,:,l), ours.R_est(:,:,l));
end

fprintf("            PACE      +UKF      OURS      LEKF \n")
fprintf("Position: %.2e, %.2e, %.2e, %.2e\n",epace.p,eukf.p,eours.p, eekf.p);
fprintf("Rotation: %.2e, %.2e, %.2e\n",epace.R(end),eukf.R(end),eours.R(end));
end