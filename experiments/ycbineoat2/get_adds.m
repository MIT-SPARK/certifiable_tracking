function [add, adds] = get_adds(gt, est, pcfiles, calc_adds)

L = length(est.p);
fprintf("\n")

%% Step 0: load point clouds
shapes = {};
for s = 1:length(pcfiles)
    shape = pcread(pcfiles(s));
    shapes{s} = double(shape.Location)';
end

%% step 1.A: compute ADD score for each pose estimate
% this score is mean distance between predicted and gt point clouds
% (including predicted/gt transforms)
add = zeros(L,1);
for l = 1:L
    if isnan(est.c(l))
        add(l) = 100;
        continue
    else
        % pc_pred = est.R(:,:,l)*shapes{est.c(l)} + est.p(:,:,l);
        pc_pred = est.R(:,:,l)*shapes{gt.c} + est.p(:,:,l);
    end
    pc_gt = gt.R(:,:,l)*shapes{gt.c} + gt.p(:,:,l);

    add(l) = mean(vecnorm(pc_pred - pc_gt));
end

%% step 1.B: compute ADD-S score for each pose estimate
% this score is mean distance computed using closest point distance
adds = zeros(L,1);
if (calc_adds)
pc_gt = shapes{gt.c};
% T = delaunayn(pc_gt');

for l = 1:L
    if isnan(est.c(l))
        adds(l) = 100;
        continue
    else
        % pc_pred = est.R(:,:,l)*shapes{est.c(l)} + est.p(:,:,l);
        pc_pred = est.R(:,:,l)*shapes{gt.c} + est.p(:,:,l);
        % this makes it much faster
        pc_pred = gt.R(:,:,l)'*(pc_pred - gt.p(:,:,l));
    end
    % pc_gt = gt.R(:,:,l)*shapes{gt.c} + gt.p(:,:,l);
    % 
    % T = delaunayn(pc_gt');
    % min distance between each predicted point and the gt point cloud
    % [~, dist] = dsearchn(pc_gt', T, pc_pred');
    [~, dist] = knnsearch(pc_gt', pc_pred');
    adds(l) = mean(dist);

    if ~mod(l,50)
        fprintf("%d%%\n",round(l/L*100))
    end
end
end

end