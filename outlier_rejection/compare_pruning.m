%% This test script compares the time-dependent pruing method
% with ROBIN, its time-independent cousin.
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

n = 50;
counts = zeros(n,3); % robin, milp, clique
robin_time = 0;
milp_time = 0;
clique_time = 0;

for iii = 1:100

%% Generate random tracking problem
problem.N_VAR = 11; % nr of keypoints
problem.K = 3; % nr of shapes
problem.L = 10; % nr of keyframes in horizon

problem.outlierRatio = 0.5;
problem.noiseSigmaSqrt = 0.01; % [m]
problem.noiseBound = 0.05;
problem.intraRadius = 0.2;
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 0.5;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

problem.accelerationNoiseBoundSqrt = 0.0;%0.01;
problem.rotationNoiseBound = 0;%pi/32; % rad

% regen if pbound, vbound, N, L, K change.
problem.regen_sdp = false; % when in doubt, set to true

problem.N = problem.N_VAR*problem.L; % How many measurements this problem has (updated by ROBIN)

% Optional: use a specified velocity trajectory
% problem = make_trajectory(problem);

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 0.0;
problem.lambda = lambda;

%% Solve!
% prune outliers with ROBIN
tic
problem_robin = robin_prune(problem);
robin_time = robin_time + toc;

% prune outliers with MILP
tic
problem_lorenzo = lorenzo_prune(problem);
milp_time = milp_time + toc;

% prune outliers with max weighted clique
tic
problem_lorenzo_graph = lorenzo_graph_prune(problem);
clique_time = clique_time + toc;

%% Check solutions
disp(iii)
if isequal(problem.inliers_gt,problem_robin.priorinliers)
    counts(iii,1) = 1;
end
if isequal(problem.inliers_gt,problem_lorenzo.priorinliers)
    counts(iii,2) = 1;
end
if isequal(problem.inliers_gt,problem_lorenzo_graph.priorinliers)
    counts(iii,3) = 1;
end

end