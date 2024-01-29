function problem = lorenzo_prune(problem, min_max_dists)
% Uses maximum weighted clique calculation to reject inliers efficiently.
% 
% Lorenzo Shaikewitz for SPARK Lab

% If not passed in: compute CAD min and max distances
if nargin == 1
    % this part can be precomputed based on CAD database
    min_max_dists = robin_min_max_dists(problem.shapes);
end
cdmin = min_max_dists{1};
cdmax = min_max_dists{2};

prioroutliers = {};
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

out = py.outlier_rejection.prune_outliers_weighted.prune_outliers(py.numpy.array(problem.y), cdmin, cdmax, problem.noiseBound, 2*problem.noiseBound, py.list(prioroutliers));
priorinliers = sort(double(out))+1;
prioroutliers = setdiff(1:problem.N_VAR*problem.L,priorinliers);

if isfield(problem,'prioroutliers')
    % setdiff(problem.prioroutliers,prioroutliers)
    prioroutliers = union(problem.prioroutliers, prioroutliers);
end

% save
problem.prioroutliers = sort(prioroutliers);
problem.priorinliers = setdiff(1:problem.N_VAR*problem.L,prioroutliers);
problem.N = problem.N - length(prioroutliers);

end