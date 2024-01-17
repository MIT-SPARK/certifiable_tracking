function problem = lorenzo_prune(problem, min_max_dists)
% Uses maximum weighted clique calculation to reject inliers efficiently.
%
% TODO:
% - update to accomodate measuring only a few keypoints!!!
% 
% Lorenzo Shaikewitz for SPARK Lab

N = problem.N_VAR;

% If not passed in: compute CAD min and max distances
if nargin == 1
    % this part can be precomputed based on CAD database
    min_max_dists = robin_min_max_dists(problem.shapes);
end
cdmin = min_max_dists{1};
cdmax = min_max_dists{2};

out = py.outlier_rejection.prune_outliers_weighted.prune_outliers(py.numpy.array(problem.y), cdmin, cdmax, problem.noiseBound, 2*problem.noiseBound);
priorinliers = sort(double(out))+1;
prioroutliers = setdiff(1:problem.N_VAR*problem.L,priorinliers);

% save
if isfield(problem,'prioroutliers')
    prioroutliers = [problem.prioroutliers, prioroutliers];
end
problem.prioroutliers = sort(prioroutliers);
problem.priorinliers = priorinliers;
problem.N = problem.N - length(prioroutliers);

end