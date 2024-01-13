function soln = solve_weighted_tracking(problem)
%% Redirects to nonminimal solver based on priors
% 
% Lorenzo Shaikewitz for SPARK Lab

% deal with prior outliers
if isfield(problem,'prioroutliers')
    % add ROBIN priors to weights
    w = reshape(problem.covar_measure.^(-1),[problem.N_VAR*problem.L,1]);

    for i = 1:length(problem.prioroutliers)
        o = problem.prioroutliers(i);
        w = [w(1:o-1),0.0,w(o:end)];
    end
end
problem.covar_measure = reshape(w.^(-1),[problem.N_VAR, problem.L]);

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

end