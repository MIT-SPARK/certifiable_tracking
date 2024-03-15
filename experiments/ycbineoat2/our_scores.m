%% Compute ADD, ADD-S scores for US
% assumes bundlesdf data is in folder labeled "out_*"
% also assumes gt pose in same folder
%
% Lorenzo Shaikewitz for SPARK Lab

clear; close all; clc;

% directories
out_dir = "~/tracking/datasets/ycbineoat/results/"; % where we put everything out
gt_dir = "~/tracking/datasets/ycbineoat/";
models_dir = "~/tracking/datasets/ycbineoat/models/";

% Load data
videos = ["cracker_box_reorient", "cracker_box_yalehand0", ...
          "sugar_box1", "sugar_box_yalehand0", ...
          "mustard0", "mustard_easy_00_02", ...
          "bleach0", "bleach_hard_00_03_chaitanya", ...
          "tomato_soup_can_yalehand0"];

% parameters to change
objects = ["bleach", "mustard", "cracker", "sugar", "soup"];

for v = videos
    % load poses
    savename = "ycbineoat_" + v;
    load(savename + "_b")
    % load(savename + "_c") % problems, solns
    load(savename + "_gt");
    
    est_legit.p = zeros(3,1,length(solns))*NaN;
    est_legit.R = zeros(3,3,length(solns))*NaN;
    est_ignoringbad.p = zeros(3,1,length(solns))*NaN;
    est_ignoringbad.R = zeros(3,3,length(solns))*NaN;
    est_bestrun.p = zeros(3,1,length(solns))*NaN;
    est_bestrun.R = zeros(3,3,length(solns))*NaN;
    gaps = [];
    for j = 1:length(solns)
        problem = problems{j};
        soln = solns{j};
        L_cur = problem.L;
        idx = problem.startIdx:(problem.startIdx + L_cur);
        % compute errors
        % 1) legit: use latest estimate no matter what
        est_legit.p(:,:,idx(end)) = soln.p(:,:,end);
        est_legit.R(:,:,idx(end)) = soln.R(:,:,end);
        % 2) ignoringbad: throw away bad
        if (soln.gap < 0.3)
            est_ignoringbad.p(:,:,idx(end)) = soln.p(:,:,end);
            est_ignoringbad.R(:,:,idx(end)) = soln.R(:,:,end);
        end
        % 3) bestrun: use estimate that has best run (by gap)
        gaps = [gaps; abs(soln.gap)];
        for l = 1:L_cur-1
            % for each pose estimate, excluding the current one
            testrange = length(gaps) - (1:(L_cur - l));
            newgood = true;
            for l_test = testrange
                if (l_test < 1)
                    break
                end
                % for each other gap we can compare with
                if (abs(soln.gap) > gaps(l_test))
                    % fails if any of them have better gaps
                    newgood = false;
                    break
                end
            end
            if (newgood)
                % update!
                est_bestrun.p(:,:,idx(l)) = soln.p(:,:,l);
                est_bestrun.R(:,:,idx(l)) = soln.R(:,:,l);
            end
        end
        est_bestrun.p(:,:,idx(end)) = soln.p(:,:,end);
        est_bestrun.R(:,:,idx(end)) = soln.R(:,:,end);
    end
    est_ignoringbad.p = est_ignoringbad.p(:,:,1:end-1);
    est_legit.p = est_legit.p(:,:,1:end-1);
    est_bestrun.p = est_bestrun.p(:,:,1:end-1);
    est_ignoringbad.R = est_ignoringbad.R(:,:,1:end-1);
    est_legit.R = est_legit.R(:,:,1:end-1);
    est_bestrun.R = est_bestrun.R(:,:,1:end-1);

    % align first frame to gt
    % invthing = (bundlePoses(:,:,1)\gtPoses(:,:,1));
    % for i = 1:length(poseTxt)
    %     bundlePoses(:,:,i) = bundlePoses(:,:,i)*invthing;
    % end
    
    % get model directory
    for obj = objects
        if (contains(savename,obj))
            o = obj;
            break
        end
    end
    pcfile = models_dir + obj + ".ply";

    % compute scores!
    [add_ours_legit, adds_ours_legit] = compute_scores(gt, est_legit, pcfile, pcfile, 0.1);
    [add_ours_ignoringbad, adds_ours_ignoringbad] = compute_scores(gt, est_ignoringbad, pcfile, pcfile, 0.1);
    [add_ours_bestrun, adds_ours_bestrun] = compute_scores(gt, est_bestrun, pcfile, pcfile, 0.1);
    
    [add_teaser, adds_teaser] = compute_scores(gt, teaser, pcfile, pcfile, 0.1);
    
    scores.(v).add_ours_legit = add_ours_legit;
    scores.(v).adds_ours_legit = adds_ours_legit;
    scores.(v).add_ours_ignoringbad = add_ours_ignoringbad;
    scores.(v).adds_ours_ignoringbad = adds_ours_ignoringbad;
    scores.(v).add_ours_bestrun = add_ours_bestrun;
    scores.(v).adds_ours_bestrun = adds_ours_bestrun;
    scores.(v).add_teaser = add_teaser;
    scores.(v).adds_teaser = adds_teaser;
end

save("our_scores.mat","scores")

%% Mean
fn = fieldnames(scores);
add = 0;
adds = 0;
for k = 1:numel(fn)
    add = add + scores.(videos{k}).add_ours_legit;
    adds = adds + scores.(videos{k}).adds_ours_legit;
end
add / numel(fn)
adds / numel(fn)

%% Reformat into table
fields = fieldnames(scores);
T = table('Size',[2,length(fields)],'VariableTypes',repelem(["double"],length(fields)),'VariableNames',fields);
for vi = 1:length(fields)
    v = fields(vi);
    v = v{1};
    if ~isfield(scores,v)
        continue
    end
    if (isempty(T.Properties.RowNames))
        % add row names
        T.Properties.RowNames = fieldnames(scores.(v));
    end
    c = struct2cell(scores.(v));
    T.(v) = [c{:}]';
end