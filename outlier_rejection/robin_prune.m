function problem = robin_prune(problem, min_max_dists)
% Uses ROBIN to compute largest set of compatible inliers at each time step
% adds field 'problem.prioroutliers' to problem for use with GNC
% also updates problem.N
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

prioroutliers = cell(1,problem.L);
if isfield(problem,'prioroutliers')
    prioroutliers_temp = sort(problem.prioroutliers)-1;
    for ii = 1:length(prioroutliers_temp)
        l = floor((prioroutliers_temp(ii))/problem.N_VAR)+1;
        out = prioroutliers_temp(ii) - (l-1)*problem.N_VAR;
        if (l > length(prioroutliers))
            prioroutliers{l} = {out};
        else
            prioroutliers{l}{end+1} = out;
        end
    end
end

prioroutliers_new = [];
for l = 1:problem.L
    % convert measurements to python
    yl = reshape(problem.y(:,l),[3,N]);
    yl_np = py.numpy.array(yl);
    
    % prune outliers with ROBIN
    out = py.outlier_rejection.prune_outliers.robin_prune_outliers(yl_np, cdmin, cdmax, problem.noiseBound, prioroutliers{l}, 'maxclique');
    out = cell(out);
    inlier_indicies = double(out{1}) + 1;

    % convert to outlier index list
    outliers = setdiff(1:N, inlier_indicies);
    prioroutliers_new = [prioroutliers_new, (l-1)*N + outliers];
end

if isfield(problem,'prioroutliers')
    % setdiff(problem.prioroutliers,prioroutliers)
    prioroutliers = union(problem.prioroutliers, prioroutliers_new);
else
    prioroutliers = prioroutliers_new;
end

% save
problem.prioroutliers = prioroutliers;
problem.priorinliers = setdiff(1:problem.N_VAR*problem.L,problem.prioroutliers);
problem.N = problem.N - length(prioroutliers);

end