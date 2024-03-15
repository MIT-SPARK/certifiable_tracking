function ycbineoat(videoNumber, skip_pruning)
%% Function version of A_prune and B_solve_many
% Run D_visualize once data is collected
%
% Lorenzo Shaikewitz for SPARK Lab
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Prune Outliers
% Generate frame problems & prune outliers with MILP.
% This must be done sequentially due to python implementation.

%% Experiment Settings
videos = ["cracker_box_reorient", "cracker_box_yalehand0", ...
          "sugar_box1", "sugar_box_yalehand0", ...
          "mustard0", "mustard_easy_00_02", ...
          "bleach0", "bleach_hard_00_03_chaitanya", ...
          "tomato_soup_can_yalehand0"];

% parameters to change
video = videos(videoNumber);
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
jsondir = "../datasets/ycbineoat/";

if ~skip_pruning
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
else
    load(savename +"_a")
    numProblemsToSolve = length(problems);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: Solve With Many Keypoints (spherical)
% Load pruned problems and solve, parallelized.

%% Run
disp("Solving " + string(numProblemsToSolve) + " problems...")
solns = cell(numProblemsToSolve,1);

% L should change for the first problem.L - 2 problems
parfor batch = 1:min(maxL-2, numProblemsToSolve) % PARFOR
    curproblem = problems{batch};
    curproblem.sdp_filename = "sdpdata" + curproblem.L;
    curproblem.regen_sdp = true;

    curproblem.noiseBound_GNC = noiseBound_GNC;
    curproblem.velocityBound = velocityBound;
    curproblem.covar_measure_base = covar_measure_base;
    curproblem.covar_velocity_base = covar_velocity_base;
    curproblem.covar_rotrate_base = covar_rotrate_base;

    soln = solveBatch(curproblem);
    solns{batch} = soln;

    % report bad results
    if (soln.gap > 1e-2)
        fprintf("Batch %d failed: %.4e\n",batch,soln.gap)
    elseif (isnan(soln.gap))
        fprintf("Batch %d failed: NaN\n",batch)
    end
end

% Now that L is fixed, run through the remainder of the problems
parfor batch = min(maxL-2, numProblemsToSolve)+1:numProblemsToSolve
    curproblem = problems{batch};
    curproblem.sdp_filename = "sdpdata" + curproblem.L;
    curproblem.regen_sdp = false;

    curproblem.noiseBound_GNC = noiseBound_GNC;
    curproblem.velocityBound = velocityBound;
    curproblem.covar_measure_base = covar_measure_base;
    curproblem.covar_velocity_base = covar_velocity_base;
    curproblem.covar_rotrate_base = covar_rotrate_base;

    if (maxL ~= curproblem.L)
        curproblem.L = maxL;
        curproblem.y = curproblem.y(:,(end-maxL+1):end);
        problems{batch} = curproblem;
    end

    soln = solveBatch(curproblem);
    solns{batch} = soln;

    % report bad results
    if (soln.gap > 1e-2)
        fprintf("Batch %d failed: %.4e\n",batch,soln.gap)
    elseif (isnan(soln.gap))
        fprintf("Batch %d failed: NaN\n",batch)
    end
end

%% Save
save(savename+"_b","problems","solns")