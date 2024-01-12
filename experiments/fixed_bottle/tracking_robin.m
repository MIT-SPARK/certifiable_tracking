%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through ROBIN+GNC
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Define settings for batch processing
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

% data for GNC
curproblem.type = "tracking";
curproblem.N = curproblem.N_VAR*curproblem.L; % How many measurements this problem has (updated by ROBIN)
curproblem.outliers = []; % outlier indicies
curproblem.priors = [];
curproblem.dof = 3;

% prune outliers with ROBIN
% curproblem = robin_prune(curproblem);

% run GNC
[inliers, info] = gnc(curproblem, @solver_for_gnc, 'NoiseBound', curproblem.noiseBound,'MaxIterations',100,'FixPriorOutliers',false);
% convert to true inliers
inliers = curproblem.priorinliers(inliers);

% Check solutions
if isequal(curproblem.inliers_gt,inliers)
    disp("Correct inliers found after " + string(info.Iterations) + " iterations.");
else
    disp("Inliers not found after " + string(info.Iterations) + " iterations.");
end

end
