%% Topic Definitions.

% Change this to point to rosbag
bagName = "../datasets/racecar_fixed/1_2024-01-30-14-00-09.bag";

% Pose estimate from Teaser + Fixed Lag Smoother (ie. softdrone algo)
sdPoseTopic = "/gtsam_tracker_node/target_global_odom_estimate";

% Ground truth target pose in the global space
gtPoseTopic = "/gtsam_tracker_node_secondary/target_global_odom_estimate";

% 3D keypoints (xyz) in the global space
keypointsTopic = "/reproject_keypoints_node/keypoints_3d_global_out";

%% Load Data.
bag = rosbag(bagName);

%% Load Messages
sdPoseSelect = select(bag, 'Topic', sdPoseTopic);
sdPoseMessages = readMessages(sdPoseSelect, 'DataFormat', 'struct');

gtPoseSelect = select(bag, 'Topic', gtPoseTopic);
gtPoseMessages = readMessages(gtPoseSelect, 'DataFormat', 'struct');

keypointsSelect = select(bag, 'Topic', keypointsTopic);
keypointsMessages = readMessages(keypointsSelect, 'DataFormat', 'struct');
