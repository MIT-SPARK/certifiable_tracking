%% Compute ADD, ADD-S scores for BundleSDF
% assumes bundlesdf data is in folder labeled "out_*"
% also assumes gt pose in same folder
%
% Lorenzo Shaikewitz for SPARK Lab

% directories
out_dir = "~/tracking/debug/"; % where BundleSDF puts everything out
gt_dir = "~/tracking/datasets/ycbineoat/";
models_dir = "~/tracking/datasets/ycbineoat/models/";

% Load BundleSDF data
folders = dir(out_dir + "out_*");
folders = {folders([folders.isdir]).name};
objects = ["bleach", "mustard", "cracker", "sugar", "soup"];

for f = folders
    % load BundleTrack poses
    d = out_dir + f + "/ob_in_cam/";
    dg = gt_dir + f{1}(5:end) + "/annotated_poses/"; % remove out_
    poseTxt = {dir(d + "*.txt").name};
    gtposeTxt = {dir(dg + "*.txt").name};
    
    bundlePoses = zeros(4,4,length(poseTxt));
    gtPoses = zeros(4,4,length(poseTxt));
    for i = 1:length(poseTxt)
        p = poseTxt(i);
        bundlePoses(:,:,i) = readmatrix(d + p);
        pg = gtposeTxt(i);
        gtPoses(:,:,i) = readmatrix(dg + pg);
    end
    

    % align first frame of bundle to gt
    invthing = (bundlePoses(:,:,1)\gtPoses(:,:,1));
    for i = 1:length(poseTxt)
        bundlePoses(:,:,i) = bundlePoses(:,:,i)*invthing;
    end
    
    % get model directory
    for obj = objects
        if (contains(f,obj))
            o = obj;
            break
        end
    end
    pcfile = models_dir + obj + ".ply";

    % compute scores!
    [add,adds] = compute_scores(gtPoses, bundlePoses, pcfile, pcfile);
    scores.(f{1}).add = add;
    scores.(f{1}).adds = adds;
end

save('bundlesdf_scores.mat',scores)

%% Mean
fn = fieldnames(scores);
add = 0;
adds = 0;
for k = 1:numel(fn)
    add = add + scores.(fn{k}).add;
    adds = adds + scores.(fn{k}).adds;
end
add / numel(fn)
adds / numel(fn)