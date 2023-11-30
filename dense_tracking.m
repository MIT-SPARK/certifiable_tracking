%% Dense SDP relaxation for certifiable tracking
% 
% PUMO:
% Pose estimation Using Multiple Observations
%
% Lorenzo Shaikewitz for SPARK Lab

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
problem.N_VAR = 11;
problem.K = 3;
problem.L = 10;

problem.outlierRatio = 0.0; % TODO: no support for outliers
problem.noiseSigma = 0.01;
problem.intraRadius = 0.2;
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 2.0;

problem.accelerationNoiseBound = 0.05;
problem.rotationNoiseBound = pi/32; % rad

% Override ground truths (for testing)
% problem.dR_gt = repmat(eye(3),1,1,problem.L);
% problem.R_gt = repmat(eye(3),1,1,problem.L);
% problem.p_gt = zeros(3,1,problem.L);
% problem.v_gt = zeros(3,1,problem.L);

% Optional: use a specified velocity trajectory
problem = make_trajectory(problem);

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 0.0;
problem.lambda = lambda;

problem.mosekpath = mosekpath;

%% Solve!
soln = solve_weighted_tracking(problem);

soln_pace = [];
for l = 1:problem.L
    pace_problem = problem;
    pace_problem.weights = ones(problem.N_VAR,1);
    pace_problem.scene = reshape(problem.y(:,l),[3,problem.N_VAR]);
    [R_est,t_est,c_est,out] = outlier_free_category_registration(pace_problem, path, 'lambda',lambda);
    s.R_est = R_est; s.p_est = t_est;
    s.c_est = c_est; s.out = out;
    soln_pace = [soln_pace; s];
end

%% Check solutions
% eigenvalue plot
L = problem.L;
figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
hold on
slices = [1:(1+9*(3*L-1)),(1+9*(3*L-1)+3*L+1):(9*(3*L-1)+6*L)];
Xopt_pRemoved = soln.raw.Xopt{1}(slices, slices);
bar([zeros(3*L+1,1);eig(Xopt_pRemoved)]);
title("Eigenvalues of Relaxed Solution")
hold off

% raw error
x_err = norm(problem.x_gt - soln.x_est);

% raw error excluding p
x_gt = problem.x_gt;
x_est = soln.x_est;
x_gt_no_p = [x_gt(1:(9*(3*L-1) + 1)); x_gt((9*(3*L-1) + 3*L + 1):end)];
x_est_no_p = [x_est(1:(9*(3*L-1) + 1)); x_est((9*(3*L-1) + 3*L + 1):end)];
x_err_no_p = norm(x_gt_no_p - x_est_no_p);

% projected errors
R_err = zeros(L,1);
dR_err = zeros(L,1);
p_err = zeros(L,1);
v_err = zeros(L,1);
p_err_bad = zeros(L,1);
for l = 1:L
    % R
    R_err(l) = getAngularError(problem.R_gt(:,:,l), soln.R_est(:,:,l));
    % dR
    dR_err(l) = getAngularError(problem.dR_gt(:,:,l), soln.dR_est(:,:,l));
    % p
    p_err(l) = norm(problem.p_gt(:,:,l) - soln.p_est(:,:,l));
    % v
    v_err(l) = norm(problem.v_gt(:,:,l) - soln.v_est(:,:,l));
    % bad p
    p_err_bad(l) = norm(problem.p_gt(:,:,l) - soln.p_est_raw(:,:,l));
end

% shape error
c_err = norm(problem.c_gt - soln.c_est);

% Plot trajectory!
plot_trajectory(problem,soln, soln_pace)
