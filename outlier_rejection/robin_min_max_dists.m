function out = robin_min_max_dists(shapes)
% Precomputes min and max distances for ROBIN from CAD library
%   This can be called once and reused if the shape library is unchanged.
%
%   shapes: 3 x K x N array
%   Example usage:
%   out = robin_min_max_dists(problem.shapes)
% 
% Lorenzo Shaikewitz for SPARK Lab

shapes = permute(shapes,[3,1,2]); % 3 x K x N -> K x 3 x N
shapes_np = py.numpy.array(shapes);
out = cell(py.robin.prune_outliers.compute_min_max_distances(shapes_np));

end

