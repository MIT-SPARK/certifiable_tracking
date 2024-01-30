%% Quick script to correct a faulty CAD frame
% 
% Lorenzo Shaikewitz for SPARK Lab

jsonfile = "../datasets/ycbineoat/cracker_box_yalehand0_metrics.json";
% Load shapes
load("../datasets/ycbineoat/cheese.mat");
shapes = annotatedPoints'; % 3 x N x K [mm]
shapes_homo = [shapes; ones(1,size(shapes,2))];

% Load data
fid = fopen(jsonfile);
raw = fread(fid,Inf);
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);

% Use gt data to fix CAD frame
gtpose = data(1).est_teaser_pose;
