%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through lorenzo+GNC
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Generate random tracking problem
problem.category = "car";
problem.L = 11; % nr of keyframes in horizon

problem.outlierRatio = 0.5;
problem.noiseSigmaSqrt = 0.05; % [m]
problem.noiseBound = 3*problem.noiseSigmaSqrt;
problem.noiseBound_GNC = 0.15;
problem.noiseBound_GRAPH = 0.20;
problem.processNoise = 0.5;
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

% for GNC
problem.N = problem.N_VAR*problem.L; % How many measurements this problem has (updated by ROBIN)
problem.outliers = []; % outlier indicies
problem.priors = [];
problem.dof = 3;


%% Solve!
% pace first
% soln_pace = pace_py_UKF(problem,true,true);

% prune outliers with max weighted clique
problem = lorenzo_prune(problem);

% run GNC
[inliers, info] = gnc_custom(problem, @solver_for_gnc, 'NoiseBound', problem.noiseBound,'MaxIterations',100,'FixPriorOutliers',true);
% convert to true inliers
inliers = problem.priorinliers(inliers);

%% Check solutions
if isequal(problem.inliers_gt,inliers)
    disp("Correct inliers found after " + string(info.Iterations) + " iterations.");
else
    if (isempty(setdiff(inliers,problem.inliers_gt)))
        disp("Subset of correct inliers found after " + string(info.Iterations) + " iterations.");
    else
        disp("Inliers not found after " + string(info.Iterations) + " iterations.");
    end
end