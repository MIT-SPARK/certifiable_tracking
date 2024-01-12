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
problem.L = 10; % batch size

% Set bounds based on problem setting
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.noiseBound = 0.2;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% regen if batch size changes.
problem.regen_sdp = false; % when in doubt, set to true

% add shape, measurements, outliers
load("cad_frame.mat");
problem.shapes = annotatedPointsWrtTarget'; % 3 x N x K
problems = bag2problem(problem);

%% Solve for each batch
solns = [];
for j = 1:length(problems)
curproblem = problems(j);

soln = solve_weighted_tracking(curproblem);

% soln_pace = pace_with_EKF(problem);

% soln = solve_full_tracking(problem,lambda);
% Ap = solve_nopos_tracking(problem);
solns = [solns; soln];
break
if problem.regen_sdp
    break;
    disp("SDP data generated. Rerun with regen_sdp true for faster results.")
end
end

%% Check solutions

figure(1);
for j = 1:length(solns)

problem = problems(j);
soln = solns(j);

% eigenvalue plot
L = problem.L;
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

% Plot trajectory!
figure(1);
hold on
p_est = reshape(soln.p_est,[3,L,1]);
plot3(p_est(1,:),p_est(2,:),p_est(3,:),'x', 'MarkerSize',30,'LineWidth',2);
end
