function [add,adds] = compute_scores(gt, est, pcfile_gt, pcfile_est, threshold)
%% Compute the ADD, ADD-S score.
% 
% Lorenzo Shaikewitz for SPARK Lab

% preliminaries
T_gt = gt;
T_est = est;
L = length(T_est);

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

%% step 1.A: compute ADD score for each pose estimate
% this score is mean distance between predicted and gt point clouds
% (including predicted/gt transforms)
add = zeros(L,1);
for l = 1:L
    pts_pred = T_est(:,:,l)*pc_est;
    pts_gt = T_gt(:,:,l)*pc_gt;

    add(l) = mean(vecnorm(pts_pred - pts_gt));
end

%% step 1.B: compute ADD-S score for each pose estimate
% this score is mean distance computed using closest point distance
adds = zeros(L,1);
for l = 1:L
    pts_pred = T_est(:,:,l)*pc_est;
    pts_pred = pts_pred(1:3,:);
    pts_gt = T_gt(:,:,l)*pc_gt;
    pts_gt = pts_gt(1:3,:);
    
    T = delaunayn(pts_gt');
    % min distance between each predicted point and the gt point cloud
    [~, dist] = dsearchn(pts_gt', T, pts_pred');
    adds(l) = mean(dist);

    if ~mod(l,50)
        fprintf("%d%%\n",round(l/L*100))
    end
end

%% Get AUC
add = get_auc(add, threshold);
adds = get_auc(adds, threshold);

end

function score = get_auc(metric, threshold)
%% step 2: compute area under curve (AUC)
% curve in question is accuracy-threshold curve
% see pose-cnn figure 8
% generate curve
thresh = linspace(0,threshold,100); % x-axis
accuracy = zeros(length(thresh),1); % y-axis
for t = 1:length(thresh)
    accuracy(t) = sum(metric < thresh(t))/length(metric);
end
% figure
plot(thresh,accuracy);

% area under curve!
max_score = threshold*1;
score = trapz(thresh, accuracy) / max_score;

end