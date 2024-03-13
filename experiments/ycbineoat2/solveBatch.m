%% Helper function: solve each batch
function soln = solveBatch(problem)    
    % run GNC
    soln = struct();
    try
        [inliers, info] = gnc_custom(problem, @solver_for_gnc, 'NoiseBound', problem.noiseBound_GNC,'MaxIterations',100,'FixPriorOutliers',false);
    
        % report compact form
        soln.p = info.f_info.soln.p_est;
        soln.R = info.f_info.soln.R_est;
        soln.c = info.f_info.soln.c_est;
        soln.gap = info.f_info.soln.gap;
        soln.solvetime = info.f_info.soln.solvetime;
        soln.iterations = info.Iterations;
        soln.inliers = problem.priorinliers(inliers);
    catch
        % report NaNs
        soln.p = ones(3,1,problem.L)*NaN;
        soln.R = ones(3,3,problem.L)*NaN;
        soln.c = ones(problem.K,1)*NaN;
        soln.gap = NaN;
        soln.solvetime = NaN;
        soln.iterations = NaN;
        soln.inliers = NaN;
    end
end