function out = robin_min_max_dists(shapes, interp)
% Precomputes min and max distances for ROBIN from CAD library
%   This can be called once and reused if the shape library is unchanged.
%
%   shapes: 3 x K x N array
%   Example usage:
%   out = robin_min_max_dists(problem.shapes)
% 
% Lorenzo Shaikewitz for SPARK Lab

if nargin < 2
    interp = false;
end
% remove repeated keypoints
if interp
    % assume keypoints are repeated uniformly across shape lib
    [C, ia, ic] = unique(shapes(:,:,1)','stable','rows');
    shapes = shapes(:,ia,:);
end

shapes = permute(shapes,[3,1,2]); % 3 x N x K -> K x 3 x N
shapes_np = py.numpy.array(shapes);
out = cell(py.outlier_rejection.prune_outliers.compute_min_max_distances(shapes_np));

if interp
    % convert to full matrix for easier editing
    key = triu(ones(length(ia)),1);
    L = zeros(length(ia));
    U = zeros(length(ia));
    L(key==1) = double(out{1});
    L = L'+L;
    L = L(ic,ic);
    U(key==1) = double(out{2});
    U = U'+U;
    U = U(ic,ic);

    key = triu(ones(length(ic)),1);
    out{1} = py.numpy.array(L(key==1));
    out{2} = py.numpy.array(U(key==1));
end

end

