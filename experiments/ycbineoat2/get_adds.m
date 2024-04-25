function [add, adds] = get_adds(gt, est, pcfile_gt)

L = length(est.p);
fprintf("\n")

%% Step 0: load point clouds
shape_gt = pcread(pcfile_gt);
shape_gt = double(shape_gt.Location)';
shape_est = shape_gt; % for now, just use ground truth (TODO: use est pc from shape coefficient)

%% step 1.A: compute ADD score for each pose estimate
% this score is mean distance between predicted and gt point clouds
% (including predicted/gt transforms)
add = zeros(L,1);
for l = 1:L
    pc_pred = est.R(:,:,l)*shape_est + est.p(:,:,l);
    pc_gt = gt.R(:,:,l)*shape_gt + gt.p(:,:,l);

    add(l) = mean(vecnorm(pc_pred - pc_gt));
end

%% step 1.B: compute ADD-S score for each pose estimate
% this score is mean distance computed using closest point distance
adds = zeros(L,1);
for l = 1:L
    pc_pred = est.R(:,:,l)*shape_est + est.p(:,:,l);
    pc_gt = gt.R(:,:,l)*shape_gt + gt.p(:,:,l);
    
    T = delaunayn(pc_gt');
    % min distance between each predicted point and the gt point cloud
    [~, dist] = dsearchn(pc_gt', T, pc_pred');
    adds(l) = mean(dist);

    if ~mod(l,50)
        fprintf("%d%%\n",round(l/L*100))
    end
end

end