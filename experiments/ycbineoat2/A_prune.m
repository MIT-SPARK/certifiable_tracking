%% Step 1: Prune Outliers
% Generate frame problems & prune outliers with MILP.
% This must be done sequentially due to python implementation.
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment Settings
videos = ["cracker_box_reorient", "cracker_box_yalehand0", ...
          "sugar_box1", "sugar_box_yalehand0", ...
          "mustard0", "mustard_easy_00_02", ...
          "bleach0", "bleach_hard_00_03_chaitanya", ...
          "tomato_soup_can_yalehand0"];

% parameters to change
video = videos(5);
maxL = 6;
skip = 1; % sequential frames or skip frames

% for pruning
noiseBound_GRAPH = 0.01;
% for GNC
noiseBound_GNC = 0.01;
% for solver
velocityBound = 1.5;
covar_measure_base = 0.01^2;
covar_velocity_base = 0.01^2;
covar_rotrate_base = 0.01^2;

savename = "ycbineoat_" + video;
jsondir = "../datasets/YCBInEOAT/";

%% Generate frame problems
% Define json, generate problem
problem = struct();
problem.json = jsondir + video + ".json";
problem.L = maxL; % batch size
problem.savefile = jsondir + video + "_ours.json";
problem.velprior = "body";

% Add shape, split into batches
[problems, gt, teaser, shapes] = json2frameproblem(problem, skip);
numProblemsToSolve = length(problems);

% precompute max/min distances for ROBIN
disp("Computing max/min distances...")
min_max_dists = robin_min_max_dists(problems{1}.shapes, true);
disp("Finished computing max/min distances")

%% Run pruning
% must be run sequentially due to COPT solver
disp("Pruning " + numProblemsToSolve + " problems...");
    
for batch = 1:numProblemsToSolve % no parfor!
    curproblem = problems{batch};
    % Pruning parameters
    curproblem.noiseBound_GRAPH = noiseBound_GRAPH;

    % stock data for GNC/pruning
    curproblem.type = "tracking";
    curproblem.N = curproblem.N_VAR*curproblem.L; % updated by pruning
    curproblem.outliers = [];
    curproblem.priors = [];
    curproblem.dof = 0;

    % run pruning!
    % no warmstart--does not help in this case
    curproblem = lorenzo_prune(curproblem, min_max_dists, false);

    % save pruned problem
    problems{batch} = curproblem;
    fprintf("%d/%d\n",batch,numProblemsToSolve);
end

%% Save pruned problems
save(savename+"_a","problems");
save(savename+"_gt","gt","teaser");