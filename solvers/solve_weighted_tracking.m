function soln = solve_weighted_tracking(problem)
%% Redirects to nonminimal solver based on priors
% 
% Lorenzo Shaikewitz for SPARK Lab

%% deal with prior outliers
L = problem.L;
if isfield(problem,'prioroutliers')
    % add ROBIN priors to weights
    w = ones(problem.N_VAR*L-length(problem.prioroutliers),1);

    for i = 1:length(problem.prioroutliers)
        o = problem.prioroutliers(i);
        w = [w(1:o-1),0.0,w(o:end)];
    end
    problem.covar_measure = reshape(w.^(-1),[problem.N_VAR, problem.L]);
end

%% Set weights
noiseBoundSq = problem.noiseBound^2;
problem.covar_measure = problem.covar_measure.*((noiseBoundSq/9));
if (isfield(problem,"covar_velocity_base"))
    problem.covar_velocity = ones(L-2,1)*problem.covar_velocity_base;
else
    problem.covar_velocity = ones(L-2,1)*problem.covar_measure(1);
end
if (isfield(problem,"kappa_rotrate_base"))
    problem.kappa_rotrate = ones(L-2,1)*problem.kappa_rotrate_base;
else
    problem.kappa_rotrate  = ones(L-2,1)*(2/problem.covar_velocity(1));
end

%% Run solver!
% default to body frame
if ~isfield(problem,'velprior')
    problem.velprior = 'body';
end

% redirect to appropriate solver
if strcmp(problem.velprior, "body")
    soln = solve_tracking_body(problem);
elseif strcmp(problem.velprior,"body-sym")
    soln = solve_tracking_body_sym(problem);
elseif strcmp(problem.velprior, "world")
    soln = solve_tracking_world(problem);
elseif strcmp(problem.velprior, "grav-world")
    soln = solve_tracking_grav_world(problem);
else
    error("Selected prior is not implemented")
end

end