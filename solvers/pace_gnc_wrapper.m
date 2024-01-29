function [f_val, info] = pace_gnc_wrapper(problem, varargin)
%% Wrapper for running PACE with GNC
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

%% Convert to weights--ignore ROBIN outliers
w = params.Results.Weights;
if isfield(problem,'prioroutliers')
    % Simply remove prior outliers from the list of measurements
    TODO
end
problem.covar_measure = reshape(w.^(-1),[problem.N_VAR, problem.L]);

% redirect to appropriate solver
soln = pace_py_UKF(problem);

%% Pull out relevant results
f_val = soln.raw.obj(1);
info.soln = soln;

% ef = eig(soln.raw.Xopt{1});
% if (ef(end-4) > 1e-6)
%     disp("Not convergent")
% end
% bar(eig(soln.raw.Xopt{1}));

% info must have:
% - residuals
info.residuals = reshape(soln.residuals, problem.N_VAR*problem.L,1);

% for ROBIN: remove ROBIN outliers from residuals
if isfield(problem,'prioroutliers')
    info.residuals(problem.prioroutliers) = [];
end

end