%% Can some kind of calibration of noise do anything?
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Generate random tracking problem
problem.bag = "../datasets/racecar_fixed/2024-01-30-18-14-08.bag";
% problem.bag = "../datasets/racecar_fixed/2_2024-01-24-10-02-26.bag";
problem.L = 10; % batch size

% Set bounds based on problem setting
problem.translationBound = 5.0;
problem.velocityBound = 2.0;
problem.noiseBound = 0.02;
problem.covar_velocity_base = 0.5^2;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% regen if batch size changes.
problem.regen_sdp = true; % when in doubt, set to true

% add shape, measurements, outliers
load("racecar_cad.mat");
problem.shapes = racecar_cad' / 1000; % 3 x N x K

[problems, gt, sd] = bag2problem(problem, 15, 50.0);
problem = problems{1};
soln = solve_weighted_tracking(problem);
problem.regen_sdp = false;


noiseBoundRange = [0.01, 0.05];
covarVelRange = [0.03, 0.05];

% Initial guess
initial_guess = [mean(noiseBoundRange), mean(covarVelRange)];

fun = @(x) f(x(1),x(2),problem);

% Options for the optimization algorithm
options = optimset('Display', 'iter', 'MaxIter', 10, 'TolX', 1e-3, 'TolFun', 1e-4);

% Minimize the function using Nelder-Mead method (fminsearch)
[optimized_xy, optimized_fval] = fminsearch(fun, initial_guess, options);

function g = f(noiseBound, covarVelociy, problem)
    problem.noiseBound = noiseBound;
    problem.covar_velocity_base = covarVelociy^2;
    [problems, gt, sd] = bag2problem(problem, 15, 50.0);
    problem = problems{1};

    problem.type = "tracking";
    problem.N = problem.N_VAR*problem.L; % How many measurements this problem has (updated by ROBIN)
    problem.outliers = []; % outlier indicies
    problem.priors = [];
    problem.dof = 3;
    problem = lorenzo_prune(problem);
    [inliers, info] = gnc_custom(problem, @solver_for_gnc, 'NoiseBound', problem.noiseBound,'MaxIterations',100,'FixPriorOutliers',false);
    soln = info.f_info.soln;

    % soln = solve_weighted_tracking(problem);

    g = abs(soln.gap);
end