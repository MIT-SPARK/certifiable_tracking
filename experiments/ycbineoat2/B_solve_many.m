%% Step 2: Solve With Many Keypoints (spherical)
% Load pruned problems and solve, parallelized.
% Run Experiment Settings section before this
%
% Lorenzo Shaikewitz for SPARK Lab

load(savename +"_a")
numProblemsToSolve = length(problems);

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