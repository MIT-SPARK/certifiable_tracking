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
assert(numel(params.Results.Weights) == problem.N, 'The weights should be a 1xN vector');

problem.covar_measure = reshape(params.Results.Weights.^(-1),[problem.N_VAR, problem.L]);

%% Redirect
% default to body frame
if ~isfield(problem,'velprior')
    problem.velprior = 'body';
end

% redirect to appropriate solver
if strcmp(problem.velprior, "body")
    soln = solve_tracking_body(problem);
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

ef = eig(soln.raw.Xopt{1});
if (ef(end-4) > 1e-6)
    disp("Not convergent")
end
% bar(eig(soln.raw.Xopt{1}));

% info must have:
% - residuals
info.residuals = reshape(soln.residuals, problem.N_VAR*problem.L,1);

end