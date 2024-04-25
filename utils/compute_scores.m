function [add,adds] = compute_scores(gt, est, pcfile_gt, pcfile_est, max_dist)
%% Compute the ADD, ADD-S score.
% 
% Lorenzo Shaikewitz for SPARK Lab

% preliminaries
if (isfield(est,"p"))
L = length(est.p);
T_gt = repmat(eye(4),[1,1,L]);
T_est = repmat(eye(4),[1,1,L]);
for l = 1:L
    T_gt(1:3,1:3,l) = gt.R(:,:,l);
    T_gt(1:3,4,l) = gt.p(:,:,l);

    T_est(1:3,1:3,l) = est.R(:,:,l);
    T_est(1:3,4,l) = est.p(:,:,l);
end
else
    T_gt = gt;
    T_est = est;
    L = length(T_est);
end

% Read point clouds
if isstring(pcfile_gt)
pc_gt = pcread(pcfile_gt);
pc_est = pcread(pcfile_est);
N_pts = length(pc_gt.Location);

% homogeneous coords
pc_gt = [double(pc_gt.Location)';ones(1,N_pts)];
pc_est = [double(pc_est.Location)';ones(1,N_pts)];
else

pc_gt = pcfile_gt;
pc_est = pcfile_est;
N_pts = size(pc_est,2);
% homogeneous coords
pc_gt = [pc_gt;ones(1,N_pts)];
pc_est = [pc_est;ones(1,N_pts)];
end

%% ADD score: compute mean distance between points
all_err = zeros(N_pts,L);
add = zeros(L,1);
for l = 1:L
    all_err(:,l) = vecnorm( ...
        T_est(:, :, l)*pc_est - ...
        T_gt( :, :, l)*pc_gt);
    add(l) = mean(all_err(:,l),"omitmissing");
end
add = add(~isnan(add));
add = VOCap(add, max_dist)*100;

%% ADD-S score
adds = zeros(L,1);
for l = 1:L
    e = T_est(:,:,l)*pc_est;
    e = e(1:3,:,:)'; % ours
    g = T_gt( :, :, l)*pc_gt;
    g = g(1:3,:,:)'; % ground truth

    % create kd tree using ground truth data
    model = createns(g);
    % find the ground truth point nearest each measurement
    idx_g = knnsearch(model,e);

    % Compute the score (TODO:CHECK)
    adds(l) = mean(vecnorm((e - g(idx_g,:))'));

    % e = T_est(:, :, l)*pc_est;
    % e = e(1:3,:,:);
    % g = T_gt( :, :, l)*pc_gt;
    % g = g(1:3,:,:);
    % nn_index = KDTreeSearcher(e');
    % [nn_dists, ~] = knnsearch(nn_index, g', 'K', 1);
    % 
    % adds(l) = mean(vecnorm(e(:,nn_dists) - g(:,:,:)));
end
adds = VOCap(adds, max_dist)*100;
end

function ap = VOCap(rec, max_dist)
    rec = sort(rec);
    n = length(rec);
    if n == 0
        ap = 0;
        return;
    end
    
    prec = (1:n)' / n;
    if nargin >= 2
        index = find(rec < max_dist);
    else
        index = 1:length(rec);
    end
    
    if isempty(index)
        ap = 0;
        return;
    end
    
    rec = rec(index);
    prec = prec(index);

    rec = [0; rec; 0.1];
    prec = [0; prec; prec(end)];

    for i = 2:length(prec)
        prec(i) = max(prec(i), prec(i-1));
    end

    diff_rec = [0; diff(rec)];
    ap = sum(diff_rec .* prec) * 10;
end