%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through GNC
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Generate random tracking problem
problem.N_VAR = 11; % nr of keypoints
problem.K = 3; % nr of shapes
problem.L = 10; % nr of keyframes in horizon

problem.outlierRatio = 0.1;
problem.noiseSigmaSqrt = 0.1; % [m]
problem.noiseBound = 10*problem.noiseSigmaSqrt;
problem.processNoise = 0.15;
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
problem.regen_sdp = true; % when in doubt, set to true

problem.N = problem.N_VAR*problem.L; % How many measurements this problem has
problem.outliers = []; % outlier indicies
problem.priors = [];
problem.dof = 3;

% Optional: use a specified velocity trajectory
% problem = make_trajectory(problem);

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 0.0;
problem.lambda = lambda;

%% Solve!
[inliers, info] = gnc_custom(problem, @solver_for_gnc, 'NoiseBound', problem.noiseBound,'MaxIterations',100,'Debug',true);

soln_pace = pace_py_UKF(problem,true);

%% Check solutions
if isequal(problem.inliers_gt,inliers)
    disp("Correct inliers found after " + string(info.Iterations) + " iterations.");
else
    disp("Inliers not found after " + string(info.Iterations) + " iterations.");
end

% play done sound
load handel.mat
sound(y,2*Fs);
