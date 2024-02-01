%% Runs batch-level tracking on racecar data with PACE
% mostly used for sanity checks.
% 
% Lorenzo Shaikewitz for SPARK Lab

problem.bag = "../datasets/racecar_fixed/2024-01-30-18-14-08.bag";
problem.L = 10; % batch size

% Set bounds based on problem setting
problem.translationBound = 5.0; % [m]
problem.velocityBound = 1.0; % [m/s]
problem.noiseBound_GNC = 0.05;
problem.noiseBound_GRAPH = 0.05;
problem.noiseBound = 0.05;
problem.covar_velocity_base = 0.001^2;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% add shape, measurements, outliers
load("racecar_cad.mat");
problem.shapes = racecar_cad' / 1000; % 3 x N x K [m]
[problems, gt, sd] = bag2problem(problem, 32, 40.0); % 15 -> 50?

gnc=true;
robin=true;
soln_pace = [];
disp("Solving " + string(length(problems)) + " problems...")
for j = 1:length(problems)
    problem = problems{j};
    for l = 1:problem.L
        pace_problem = problem;
        pace_problem.weights = ones(problem.N_VAR,1);
        pace_problem.scene = reshape(problem.y(:,l),[3,problem.N_VAR]);
        if (gnc)
            pace_problem.N = pace_problem.N_VAR;
            pace_problem.L = 1;
            if (robin)
                pace_problem.y = problem.y(:,l);
                pace_problem = robin_prune(pace_problem);
                pace_problem.N_VAR = pace_problem.N;
            end
            if isfield(pace_problem,'prioroutliers')
                pace_problem.scene(:,pace_problem.prioroutliers) = [];
                pace_problem.weights(pace_problem.prioroutliers) = [];
                pace_problem.shapes(:,pace_problem.prioroutliers,:) = [];
            end
            SDP = relax_category_registration_v2(pace_problem,'lambda',problem.lambda,'checkMonomials',false);
            try
                out = gnc_category_registration(pace_problem,SDP,path,'lambda',problem.lambda);
                [R_est,t_est] = invert_transformation(out.R_est,out.t_est);
                c_est = out.c_est;
            catch
                out = NaN;
                t_est = ones(3,1)*NaN;
                R_est = ones(3,3)*NaN;
                c_est = ones(problem.K,1)*NaN;
                disp("PACE GNC Failed.")
            end
        else
            [R_est,t_est,c_est,out] = outlier_free_category_registration(pace_problem, 'lambda',problem.lambda);
        end
        s.R_est = R_est; s.p_est = t_est;
        s.c_est = c_est; s.out = out;
        s.s_est = R_est'*t_est;
        soln_pace = [soln_pace; s];
    end
end

%% PLOT
figure
p = [soln_pace.p_est];
plot3(p(1,:),p(2,:),p(3,:),'.r', 'MarkerSize',10);