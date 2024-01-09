function problem = clipper_reject(problem)
%% Rejects outliers using weighted maximum clique.
% - at each time step, use shape-informed invarients
%   to build consistency graph.
% - weight edges of consistency graph by time consistency
%   of edge (using rigid body/constant shape assumption)
% - return pruned problem.

NOISE_BOUND = problem.noiseBound;
RIGID_NOISE_BOUND = problem.noiseBound;

%% Build consistency graph
minmax_dist = robin_min_max_dists(problem.shapes);
C_list = shape_consistency(problem, minmax_dist, NOISE_BOUND);

%% Weight edges
A_list = time_consistency(problem, C_list, RIGID_NOISE_BOUND);

%% Solve with clipper at each time step
for l = 1:problem.L
    M = A_list(:,:,l) + eye(problem.N_VAR);
    C = C_list(:,:,l);

    u = py.outlier_rejection.run_clipper.run_clipper(py.numpy.array(M), py.numpy.array(C));
    idx = double(u)+1

    % TODO: SAVE INLIER + OUTLIER LIST
end

end