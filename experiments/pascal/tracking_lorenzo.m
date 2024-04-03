%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through lorenzo+GNC
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Generate random tracking problem
problem.category = "aeroplane";
problem.L = 8; % nr of keyframes in horizon

problem.outlierRatio = 0.05;
problem.outlierVariance = 1.6;

problem.noiseSigmaSqrt = 0.05*0.2; % [m]
problem.covar_measure_base = 1;
problem.covar_velocity_base = 1;
problem.covar_rotrate_base = 1;

problem.noiseBound = 0.15*0.2;
problem.noiseBound_GNC = 0.01;%0.15*0.3;
problem.noiseBound_GNC_residuals = 1;
problem.noiseBound_GRAPH = 0.15*0.2;
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
problem.regen_sdp = true; % when in doubt, set to true
problem.usecBound = false;

% Optional: use a specified velocity trajectory
% problem = make_trajectory(problem);

% add shape, measurements, outliers
problem = gen_pascal_tracking(problem);
lambda = 0.;
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
% problem = lorenzo_prune(problem);

% run GNC
[inliers, info] = gnc2(problem, @solver_for_gnc,'barc2',problem.noiseBound_GNC, 'ContinuationFactor', 1.8);%, 'NoiseBound', problem.noiseBound_GNC,'MaxIterations',100, 'Debug', true);
% convert to true inliers
% inliers = problem.priorinliers(inliers);

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

view_gnc(problem,info);

% % ground truth version
% theta_gt = problem.theta_gt;
% theta_gt(theta_gt<0) = 0;
% [fgt, infogt] = solver_for_gnc(problem, 'Weights',theta_gt);
% if fgt < info.soln.obj_est
%     fprintf("Converged to local minimum %3.2e > %3.2e.\n",info.soln.obj_est, fgt)
% else
%     fprintf("Global minimum: %3.2e <= %3.2e (diff: %1.1e).\n",info.soln.obj_est, fgt, fgt - info.soln.obj_est)
% end
