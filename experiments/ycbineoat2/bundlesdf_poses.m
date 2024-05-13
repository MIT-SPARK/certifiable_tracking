%% Compute ADD, ADD-S scores for BundleSDF
% assumes bundlesdf data is in folder labeled "out_*"
% also assumes gt pose in same folder
%
% Lorenzo Shaikewitz for SPARK Lab

% directories
out_dir = "~/tracking/datasets/BundleSDF_data/"; % where BundleSDF puts everything out
gt_dir = "~/tracking/datasets/YCBInEOAT/";

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

    % save as metadata
    json = gt_dir + "TODO";
    jsonout = json + "_bundle";
    save2json(json, jsonout, bundlePoses)
end

function save2json(jsonin, jsonout, T_est)
    L_big = length(T_est);
    for l = 1:L_big
        % convert to mm
        T_est(1:3,4,l) = T_est(1:3,4,l)*1000.0;
    end
    
    fid = fopen(jsonin);
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid);
    data = jsondecode(str);
    
    for l = 1:length(T_est)
        data(l).bundlesdf_pose = T_est(:,:,l);
    end
    
    cocoString = jsonencode(data, "PrettyPrint",true);
    fid = fopen(jsonout, 'w');
    fprintf(fid, '%s', cocoString);
    fclose(fid);
    end