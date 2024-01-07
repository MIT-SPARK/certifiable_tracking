function C_list = shape_consistency(problem, minmax_dist, NOISE_BOUND)
% For each edge, use shape invarients to compute
% whether the edge is feasible given the shape library.
% 
% Lorenzo Shaikewitz for SPARK Lab

% setup
N = problem.N_VAR;
L = problem.L;

cdmin = minmax_dist{1};
cdmax = minmax_dist{2};

C_list = zeros(N,N,L);

for l = 1:L
    % convert measurements to python
    yl = reshape(problem.y(:,l),[3,N]);
    yl_np = py.numpy.array(yl);
    
    % prune outliers with ROBIN
    out = py.robin.prune_outliers.make_graph(yl_np, cdmin, cdmax, NOISE_BOUND);
    out = double(out);

    % save
    C_list(:,:,l) = out;
end

end