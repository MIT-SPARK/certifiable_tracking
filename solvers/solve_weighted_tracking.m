function soln = solve_weighted_tracking(problem)
%% Redirects to nonminimal solver based on priors
% 
% Lorenzo Shaikewitz for SPARK Lab

%% deal with prior outliers
L = problem.L;
if isfield(problem,'prioroutliers')
    % add ROBIN priors to weights
    w = ones(1,problem.N_VAR*L-length(problem.prioroutliers));

    for i = 1:length(problem.prioroutliers)
        o = problem.prioroutliers(i);
        w = [w(1:o-1),0.0,w(o:end)];
    end
    problem.covar_measure = reshape(w.^(-1),[problem.N_VAR, problem.L]);
elseif isfield(problem,'covar_measure')
    disp("OVERRIDING COVAR MEASURE. YOU PROBABLY DON'T WANT TO DO THIS.")
    % use only for testing: skipping points or other tests
else
    problem.covar_measure = ones(problem.N_VAR,L);
end
if (isfield(problem,"covar_measure_base"))
    problem.covar_measure = problem.covar_measure*problem.covar_measure_base;
end

%% Set weights
% multiplying by a constant amount (noiseBound) does nothing to
% optimization problem.
% noiseBoundSq = problem.noiseBound^2;
% problem.covar_measure = problem.covar_measure.*((noiseBoundSq/9));
if (isfield(problem,"covar_velocity_base"))
    problem.covar_velocity = ones(L-2,1)*problem.covar_velocity_base;
else
    base = mean(problem.covar_measure(~isinf(problem.covar_measure)));
    problem.covar_velocity = ones(L-2,1)*base;
end
if (isfield(problem,"kappa_rotrate_base"))
    problem.kappa_rotrate = ones(L-2,1)*problem.kappa_rotrate_base;
elseif (isfield(problem,"covar_rotrate_base"))
    problem.kappa_rotrate = ones(L-2,1)*(1/problem.covar_rotrate_base*(1/2));
else
    base = mean(problem.covar_velocity(~isinf(problem.covar_velocity)));
    problem.kappa_rotrate  = ones(L-2,1)*(1/base*(1/2));
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