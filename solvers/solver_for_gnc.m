function [f_val, info] = solver_for_gnc(problem, varargin)
%% Redirects to nonminimal solver based on priors
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Parse GNC params
% varargin must have:
% - weights
default_weights = ones(1, problem.N);

params = inputParser;
params.CaseSensitive = false;
params.addParameter('Weights', default_weights, @(x) isnumeric(x));
params.parse(varargin{:});
% assert(numel(params.Results.Weights) == problem.N, 'The weights should be a 1xN vector');

%% Convert to weights--ignore ROBIN outliers
w = params.Results.Weights;
if isfield(problem,'prioroutliers')
    % add ROBIN priors to weights
    for i = 1:length(problem.prioroutliers)
        o = problem.prioroutliers(i);
        w = [w(1:o-1),0.0,w(o:end)];
    end
end
problem.covar_measure = reshape(w.^(-1),[problem.N_VAR, problem.L]);

%% Set weights
L = problem.L;
noiseBoundSq = problem.noiseBound^2;
problem.covar_measure = problem.covar_measure.*((noiseBoundSq/9));
if (isfield(problem,"covar_velocity_base"))
    problem.covar_velocity = ones(L-2,1)*problem.covar_velocity_base;
else
    base = mean(problem.covar_measure(~isinf(problem.covar_measure)));
    problem.covar_velocity = ones(L-2,1)*base;
end
if (isfield(problem,"kappa_rotrate_base"))
    problem.kappa_rotrate = ones(L-2,1)*problem.kappa_rotrate_base;
else
    base = mean(problem.covar_velocity(~isinf(problem.covar_velocity)));
    problem.kappa_rotrate  = ones(L-2,1)*(2/base);
end

%% Redirect to appropriate solver
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

%% Pull out relevant results
f_val = soln.raw.obj(1);
info.soln = soln;
info.residuals = reshape(soln.residuals, problem.N_VAR*problem.L,1);

% for ROBIN: remove ROBIN outliers from residuals
if isfield(problem,'prioroutliers')
    info.residuals(problem.prioroutliers) = [];
end

end