%% Helper function: solve each batch
function soln = solveBatch(problem)    
    % run GNC
    soln = struct();
    [inliers, info] = gnc2(problem, @solver_for_gnc,'barc2',problem.noiseBound_GNC, 'ContinuationFactor', 1.6);

    % report compact form
    soln.p = info.f_info.soln.p_est;
    soln.R = info.f_info.soln.R_est;
    soln.c = info.f_info.soln.c_est;
    soln.gap = info.f_info.soln.gap;
    soln.solvetime = info.f_info.soln.solvetime;
    soln.iterations = info.Iterations;
    soln.inliers = problem.priorinliers(inliers);
end