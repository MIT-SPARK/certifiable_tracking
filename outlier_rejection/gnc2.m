function [inliers, info] = gnc2(problem, f, varargin)
%% GNC loop, based on Pasquale Antonante's implementation
% Key variables to watch:
% - residuals (comes from solver)
% - weights (passed to solver)
%
%
% Lorenzo Shaikewitz for SPARK Lab

%% Preliminaries
% stopping parameters
params = inputParser;
params.addParameter('MaxIterations', 1e2, @(x) isscalar(x));
params.addParameter('StopThreshold', 1e-6, @(x) isscalar(x));
params.addParameter('TightnessThreshold', 1e-6, @(x) isscalar(x));
params.addParameter('ContinuationFactor', 1.4, @(x) isscalar(x));
params.addParameter('barc2', 1, @(x) isscalar(x));
params.addParameter('FailGracefully', true, @(x) islogical(x)) % return first run if failure
params.addParameter('Debug', false, @(x) islogical(x));
params.parse(varargin{:});

maxSteps = params.Results.MaxIterations;
stopTh = params.Results.StopThreshold;
tightTh = params.Results.TightnessThreshold; % not used
debug = params.Results.Debug;
N = problem.N;

% GNC Parameters
barc2 = params.Results.barc2;
divFactor = params.Results.ContinuationFactor;

%% GNC Main Loop
weights = ones(N,1);
prevCost = 1e6;
costDiff = 1e6;
gap = 1e6;

% for logging
if debug
    history.weights = weights;
    history.fcost = prevCost;
    history.costDiff = 0;
    history.gap = gap;
    history.mu = [];
    history.residuals = zeros(N,1);
    history.solns = {};
end
history.infos = {};
failed = false;

for itr = 0:maxSteps
    % Termination conditions
    if costDiff < stopTh
        fprintf("GNC converged %3.2e < %3.2e.\n", costDiff, stopTh);
        break;
    end
    % if gap < tightTh
    %     % WARNING: this may cause issues!
    %     fprintf("GNC converged to tight solution %3.2e < %3.2e.\n",gap,tightTh);
    %     break;
    % end
    if (max(abs(weights)) < 1e-6)
        fprintf("GNC encounters numerical issues. The solution is likely wrong.\n");
        failed = true;
        break;
    end

    % Run nonminimal solver
    [~, info] = f(problem,'Weights',weights);
    residuals = info.residuals*problem.covar_measure_base;
    gap = info.soln.gap;
    f_cost = residuals(:)'*weights(:) + barc2*sum(weights==0) + problem.lambda*norm(info.soln.c_est);

    % Initialize mu
    if itr < 1
        maxResidual = max(residuals);
        mu = barc2 / (2*maxResidual - barc2);
        history.mu = mu;
        if ~problem.usecBound
            problem.regen_sdp = false;
        end
        if mu < 0
            break;
        end
    end

    % Weights update
    th1 = (mu+1)/mu * barc2;
    th2 = mu/(mu+1) * barc2;
    for i = 1:N
        if residuals(i) - th1 >= 0
            weights(i) = 0;
        elseif residuals(i) - th2 <= 0
            weights(i) = 1;
        else
            weights(i) = sqrt(barc2*mu*(mu+1)/residuals(i)) - mu;
        end
    end

    % Update cost diff
    costDiff = abs(f_cost - prevCost);
    prevCost = f_cost;

    % Increase mu
    mu = mu*divFactor;

    % save history
    if debug
        history.weights = [history.weights, weights];
        history.fcost = [history.fcost, f_cost];
        history.costDiff = [history.costDiff, costDiff];
        history.gap = [history.gap, gap];
        history.mu = [history.mu, mu];
        history.residuals = [history.residuals, residuals];
        history.solns{end+1} = info.soln;
    end
    history.infos{end+1} = info;
end

if (failed) && (params.Results.FailGracefully)
    info = history.info{1};
    info.failed = true;
end

%% Convert to output
theta_est = zeros(N,1);
theta_est(weights>0.5) = 1;
theta_est(weights<0.5) = -1;

allPoints = 1:N;
inliers = allPoints(theta_est > 0);
info.f_info = info;
info.Iterations = itr;

%% Debug: print history
if debug
    history.residuals = history.residuals(:,[2:end,1]);
    info.history = history;

    % replay history
    fprintf(' itr |    obj        delobj  |    mu    |   sumw   |  sumout  |   maxres   |   gap   |\n')
    for itr = 1:size(info.history.weights,2)
        f_cost = history.fcost(itr);
        costDiff = history.costDiff(itr);
        mu = history.mu(itr);
        weights = history.weights(:,itr);
        residuals = history.residuals(:,itr);
        gap = history.gap(itr);
        fprintf('%4d | %3.4e %3.4e | %3.2e | %3.2e | %3.2e | %3.4e | %1.1e |\n',...
            itr-1,f_cost,costDiff,mu,sum(weights),N-sum(weights),max(residuals),gap);
    end
end
