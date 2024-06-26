%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through ROBIN+GNC
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Generate random tracking problem based on PASCAL
problem.category = "car";
problem.L = 10; % nr of keyframes in horizon

problem.outlierRatio = 0.3;
problem.noiseSigmaSqrt = 0.01; % [m]
problem.noiseBound = 3*problem.noiseSigmaSqrt;
problem.intraRadius = 0.2;
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

problem.accelerationNoiseBoundSqrt = 0.0;%0.01;
problem.rotationNoiseBound = 0;%pi/32; % rad

% regen if pbound, vbound, N, L, K change.
problem.regen_sdp = false; % when in doubt, set to true

% Optional: use a specified velocity trajectory
% problem = make_trajectory(problem);

% add shape, measurements, outliers
problem = gen_pascal_tracking(problem);
lambda = 0.0;
problem.lambda = lambda;

% data for GNC
problem.N = problem.N_VAR*problem.L; % How many measurements this problem has (updated by ROBIN)
problem.outliers = []; % outlier indicies
problem.priors = [];
problem.dof = 3;

%% Solve!
% prune outliers with ROBIN
problem = robin_prune(problem);

% run GNC
[inliers, info] = gnc(problem, @solver_for_gnc, 'NoiseBound', problem.noiseBound,'MaxIterations',100,'FixPriorOutliers',true);
% convert to true inliers
inliers = problem.priorinliers(inliers);

%% Check solutions
if isequal(problem.inliers_gt,inliers)
    disp("Correct inliers found after " + string(info.Iterations) + " iterations.");
else
    disp("Inliers not found after " + string(info.Iterations) + " iterations.");
end
