N = problem.N_VAR;

% get cad_dist_min and max
shapes = permute(problem.shapes,[3,1,2]); % 3 x K x N -> K x 3 x N
shapes_np = py.numpy.array(shapes);
out = cell(py.prune_outliers.compute_min_max_distances(shapes_np));
cdmin = out{1};
cdmax = out{2};

% convert yi to python
i = 1;
yi = reshape(problem.y(:,i),[3,N]);
yi_np = py.numpy.array(yi);

% prune outliers with ROBIN
out = py.prune_outliers.robin_prune_outliers(yi_np, cdmin, cdmax, problem.noiseBound, 'maxclique');
out = cell(out);
inlier_indicies = double(out{1}) + 1;

disp("Found " + length(inlier_indicies) + " inliers.")