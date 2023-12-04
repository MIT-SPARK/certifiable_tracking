%% Minimal Demo of Certifiable Tracking
% Requires no installation/setup
%
% Key variables:
% problem:
%   struct, contains all problem setup info and ground truths.
% 
% soln:
%   struct, contains raw MOSEK output (soln.raw) and derived outputs,
%   including estimated positions, rotations, etc.
%

clc; clear; close all; restoredefaultpath

%% Load variables
load("data/demo_data.mat")

% add internal paths
addpath('./visualization')
addpath('./utils')

%% Plot
L = problem.L;
% eigenvalue plot
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