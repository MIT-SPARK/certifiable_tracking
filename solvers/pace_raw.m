function soln = pace_raw(problem, gnc, robin)
% A naive tracking approach: just perform pose estimation at each time.
% 
% Lorenzo Shaikewitz for SPARK Lab

if nargin < 2
    gnc = false;
end
if nargin < 3
    robin = false;
end

%% Run PACE at each time step
soln_pace = [];
for l = 1:problem.L
    pace_problem = problem;
    
    % scale weights
    pace_problem.weights = ones(problem.N_VAR,1);
    noiseBoundSq = problem.noiseBound^2;
    pace_problem.weights = pace_problem.weights.*(9/noiseBoundSq);

    % convert y into native PACE form
    pace_problem.scene = reshape(problem.y(:,l),[3,problem.N_VAR]);
    if (gnc)
        pace_problem.N = pace_problem.N_VAR;
        pace_problem.L = 1;

        % prune with ROBIN if requested
        if (robin)
            pace_problem.y = problem.y(:,l);
            pace_problem = robin_prune(pace_problem);
            pace_problem.N_VAR = pace_problem.N;
        end

        % ignore prior outliers (bad keypoints or ROBIN)
        if isfield(pace_problem,'prioroutliers')
            pace_problem.scene(:,pace_problem.prioroutliers) = [];
            pace_problem.weights(pace_problem.prioroutliers) = [];
            pace_problem.shapes(:,pace_problem.prioroutliers,:) = [];
        end

        % Relax and solve
        SDP = relax_category_registration_v2(pace_problem,'lambda',problem.lambda,'checkMonomials',false);
        try
            out = gnc_category_registration(pace_problem,SDP,path,'lambda',problem.lambda);
            [R_est,t_est] = invert_transformation(out.R_est,out.t_est);
            c_est = out.c_est;
            gap = out.gap;
        catch
            % GNC may fail. Catch these and report NaN for pose estimate
            out = NaN;
            t_est = ones(3,1)*NaN;
            R_est = ones(3,3)*NaN;
            c_est = ones(problem.K,1)*NaN;
            gap = NaN;
            disp("PACE GNC Failed.")
        end
    else
        % Run without GNC
        [R_est,t_est,c_est,out] = outlier_free_category_registration(pace_problem, 'lambda',problem.lambda);
        gap = out.gap;
    end
    s.R_est = R_est; s.p_est = t_est;
    s.c_est = c_est; s.out = out;
    s.s_est = R_est'*t_est;
    s.gap = gap;
    soln_pace = [soln_pace; s];
end

%% Save
K = problem.K;
L = problem.L;
soln.p = reshape([soln_pace.p_est],[3,1,L]);
soln.R = reshape([soln_pace.R_est],[3,3,L]);
soln.c = reshape([soln_pace.c_est],[K,1,L]);

end