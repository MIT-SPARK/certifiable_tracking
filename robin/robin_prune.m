function problem = robin_prune(problem, min_max_dists)
% Uses ROBIN to compute largest set of compatible inliers at each time step
% adds field 'problem.prioroutliers' to problem for use with GNC
%
% TODO:
% - Add more advanced pruning
% - update to accomodate measuring only a few keypoints
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

outlierpriors = [];
for l = 1:problem.L
    % convert measurements to python
    yl = reshape(problem.y(:,l),[3,N]);
    yl_np = py.numpy.array(yl);
    
    % prune outliers with ROBIN
    out = py.prune_outliers.robin_prune_outliers(yl_np, cdmin, cdmax, problem.noiseBound, 'maxclique');
    out = cell(out);
    inlier_indicies = double(out{1}) + 1;

    % convert to outlier index list
    outliers = setdiff(1:N, inlier_indicies);
    outlierpriors = [outlierpriors, (l-1)*N + outliers];
end

% save
problem.outlierpriors = outlierpriors;

end