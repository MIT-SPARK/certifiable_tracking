%% Step 3: Solve Reduced Version
% Load solved spherical problems and rerun without spherical
% Run Experiment Settings section before this
%
% Lorenzo Shaikewitz for SPARK Lab

load(savename + "_b") % problems, solns
numProblemsToSolve = length(problems);

%% Select the closest inlier from the spherical pool
% aggregate over all solved batches

% get unique shapes
[C, ia, ic] = unique(problems{1}.shapes(:,:,1)','stable','rows');
shapes = problems{1}.shapes(:,ia,:);

% placeholder for measurements
measurements_reduced = zeros(3*length(ia),numProblemsToSolve + 2); % 3N x L
gaps = ones(numProblemsToSolve+2,1);

% loop through solutions
for index = 1:length(solns)
    soln = solns{index};
    problem = problems{index};

    for l = 1:problem.L
        globalIdx = index+2 - (problem.L-l);
        % check the quality of this sample
        if soln.gap > gaps(globalIdx)
            % don't update
            continue
        end
        gaps(globalIdx) = soln.gap;

        yl = reshape(problem.y(:,l),[3,problem.N_VAR])';

        idxBase = problem.N_VAR*(l-1);
        for kpt = 1:length(ia)
            % if original keypoint is inlier, just keep it
            if sum(soln.inliers==(idxBase+kpt)) > 0
                measurements_reduced(ib3(kpt),globalIdx) = yl(kpt,:);
                continue
            end
            
            % Otherwise, pick the inlier keypoint closest to the original
            % use x-y since data is stored as pixels
            t = 1:length(ic);
            toCheck = t(ic == kpt);
            associatedInliers = intersect(idxBase + toCheck, soln.inliers);
            if isempty(associatedInliers)
                continue
            end
            yassoc = yl(associatedInliers - idxBase,:);

            dists = yassoc - [yl(kpt,1), yl(kpt,2), 0];
            dists = vecnorm(dists(:,1:2)');
            [~,minIdx] = min(dists);
            
            measurements_reduced(ib3(kpt),globalIdx) = yassoc(minIdx,:);
        end
    end
end

%% Construct reduced problems
problems_reduced = cell(numProblemsToSolve,1);
for index = 1:numProblemsToSolve
    problem = problems{index};
    N = length(ia);
    K = size(shapes,3);

    problem.shapes = shapes;
    problem.B = reshape(problem.shapes, 3*N, K);
    globalRan = (index+2 - problem.L + 1):(index+2);
    problem.y = measurements_reduced(:,globalRan);

    problem.N_VAR = N;
    problem.N = problem.N_VAR*problem.L;
    problem.prioroutliers = [];
    problem.priorinliers = [];

    problems_reduced{index} = problem;
end
problems = problems_reduced;

%% Prune outliers
% precompute max/min distances for ROBIN
disp("Computing max/min distances...")
min_max_dists = robin_min_max_dists(problems{1}.shapes);
disp("Finished computing max/min distances")

%% Run pruning
% must be run sequentially due to COPT solver
disp("Pruning " + numProblemsToSolve + " problems...");
    
for batch = 1:numProblemsToSolve
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
    % warmstart should help because problems are small
    curproblem = lorenzo_prune(curproblem, min_max_dists, true);

    % save pruned problem
    problems{batch} = curproblem;
    fprintf("%d/%d\n",batch,numProblemsToSolve);
end

%% Solve each batch!
disp("Solving " + string(numProblemsToSolve) + " problems...")
solns = cell(numProblemsToSolve,1);

% L should change for the first problem.L - 2 problems
parfor batch = 1:min(maxL-2, numProblemsToSolve)
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

%% Save results
save(savename+"_c","problems","solns")