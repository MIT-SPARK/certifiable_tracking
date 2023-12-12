%% Dense SDP relaxation for certifiable tracking
%  Generic, tunable script to run one iteration of dense tracking.
%    Operates on random data with no outlier support.
%
% Lorenzo Shaikewitz for SPARK Lab

% TODOS:
% 1) selector for PACE
% 2) selector for data checking (below)
% 3) selector for viz

clc; clear; close all; restoredefaultpath
% rng("default")

%% dependencies
% Change paths here
certifiablyrobustperceptionpath = "../CertifiablyRobustPerception";
mosekpath   = 'C:/Program Files/Mosek/10.1/toolbox/r2017a';
sdpnalpath  = '../SDPNALv1.0';

% add external paths
spotpath    = certifiablyrobustperceptionpath + '/spotless';
stridepath  = certifiablyrobustperceptionpath + '/STRIDE';
manoptpath  = certifiablyrobustperceptionpath + '/manopt';
addpath(certifiablyrobustperceptionpath + '/utils')
addpath(genpath(spotpath)) % Use spotless for defining polynomials
addpath(certifiablyrobustperceptionpath + '/SDPRelaxations') % implementations for SDP relaxation

% add internal paths
addpath('./solvers')
addpath('./visualization')
addpath('./utils')

path.stridepath = stridepath;
path.mosekpath  = mosekpath;
path.manoptpath = manoptpath;

%% Generate random tracking problem
problem.N_VAR = 11; % nr of keypoints
problem.K = 3; % nr of shapes
problem.L = 10; % nr of keyframes in horizon

problem.outlierRatio = 0.0; % TODO: no support for outliers
problem.noiseSigmaSqrt = 0.01; % [m]
problem.intraRadius = 0.2; 
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

% problem.velprior = "body";       % constant body frame velocity
problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

problem.accelerationNoiseBoundSqrt = 0.01;
problem.rotationNoiseBound = pi/32; % rad

% Optional: use a specified velocity trajectory
% problem = make_trajectory(problem);

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 0.0;
problem.lambda = lambda;

problem.mosekpath = mosekpath;

%% Solve!
soln = solve_weighted_tracking(problem);

soln_pace = pace_with_EKF(problem, path);

% soln = solve_full_tracking(problem,lambda);
% Ap = solve_nopos_tracking(problem);

%% Check solutions
% eigenvalue plot
L = problem.L;
figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
hold on

if strcmp(problem.velprior, "body")
    slices = 1:(1+9*(3*L-2)+3*L);
    Xopt_pRemoved = soln.raw.Xopt{1}(slices, slices);
    bar([zeros(3*(L-1),1);eig(Xopt_pRemoved)]);
    title("Eigenvalues of Relaxed Solution")
elseif strcmp(problem.velprior, "world")
    slices = [1:(1+9*(3*L-2)),(1+9*(3*L-2)+3*L+1):(9*(3*L-2)+6*L)];
    Xopt_pRemoved = soln.raw.Xopt{1}(slices, slices);
    bar([zeros(3*L+1,1);eig(Xopt_pRemoved)]);
    title("Eigenvalues of Relaxed Solution")
elseif strcmp(problem.velprior, "grav-world")
    error("Selected prior is not implemented")
else
    error("Selected prior is not implemented")
end
hold off

% raw error
x_err = norm(problem.x_gt - soln.x_est);

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
