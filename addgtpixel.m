
% files to open
jsonfile = "/home/lorenzo/research/tracking/datasets/YCBInEOAT/mustard_easy_00_02.json";
gtfile = "/home/lorenzo/research/tracking/keypoint_detection/SoftDrone-ROS/softdrone_target_pose_estimator/notebooks/metrics/mustard_easy_gt.json";

% open
fid = fopen(jsonfile);
raw = fread(fid,Inf);
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);

fid = fopen(gtfile);
raw = fread(fid,Inf);
str = char(raw'); 
fclose(fid); 
data_gt = jsondecode(str);

% modify
for l = 1:length(data)
    data(l).gt_pixel_est_depth_world_keypoints = data_gt(l).gt_pixel_est_depth_world_keypoints;
end
return

%% save
cocoString = jsonencode(data, "PrettyPrint",true);
fid = fopen(jsonfile, 'w');
fprintf(fid, '%s', cocoString);
fclose(fid);