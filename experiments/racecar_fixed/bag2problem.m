function problem_list = bag2problem(problem)
%% generates a problem from rosbag data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Define shape data
% shapes is 3 x N x K
problem.N_VAR = size(problem.shapes,2);
problem.K = size(problem.shapes,3);
N = problem.N_VAR;
K = problem.K;

problem.B = reshape(problem.shapes, 3*N, K);
B = problem.B;

% set lambda if K > N
if K > N
    problem.lambda = 1.0;
    disp("Set lambda to " + string(problem.lambda));
else
    problem.lambda = 0.0;
end

%% Load and parse bag data into problem format
[stamps, measurements, gt, sd] = parseBag(problem.bag);
% TODO: PROCESS GROUND TRUTH CORRECTLY
tot_L = length(stamps);

L = problem.L;

problem_list = {};

% define batch problems
for batch = 1:floor(tot_L/L)
    idxrange = ((batch-1)*L+1):(batch*L);
    curproblem = problem;
     
    % interpolate between measurements
    t = stamps(idxrange);
    t = t - t(1);
    m = measurements(:,:,idxrange); % 3 x N x L
    m(m==0) = NaN;
    t_even = linspace(0,t(end),L);
    dt = t_even(2) - t_even(1);
    
    y = zeros(3*N,L);
    for i = 1:N
        xtemp = interp1(t,reshape(m(1,i,:),[L,1,1]),t_even);
        ytemp = interp1(t,reshape(m(2,i,:),[L,1,1]),t_even);
        ztemp = interp1(t,reshape(m(3,i,:),[L,1,1]),t_even);
        y(ib3(i),:) = [xtemp;ytemp;ztemp];
    end
    % change weights to ignore nans
    prioroutliers = [];
    for l = 1:L
        yl = y(:,l);
        i3_nan = strfind(isnan(yl)',true(1,3));
        for j = 1:length(i3_nan)
            i3 = i3_nan(j);
            i = (i3-1)/3 + 1;
            prioroutliers(end+1) = i + N*(l-1);
        end
        curproblem.prioroutliers = prioroutliers;
    end

    % y cannot have nans in it
    y(isnan(y)) = 0.0;

    % set covariances
    noiseBoundSq = problem.noiseBound^2;
    weights = ones(N*L-length(prioroutliers),1)*((noiseBoundSq/9).^(-1));
    covar_velocity = ones(L-2,1)*weights(1)*1;
    kappa_rotrate  = ones(L-2,1)*(2/covar_velocity(1));

    % Save gt and sd poses
    t_gt = gt.stamps - t(1);
    t_sd = sd.stamps - t(1);
    p_gt = zeros(3,1,L);
    p_sd = zeros(3,1,L);
    for i = 1:3
        p_gt(i,:,:) = interp1(t_gt,squeeze(gt.p_gt(i,:,:)),t_even);
        p_sd(i,:,:) = interp1(t_sd,squeeze(sd.p_sd(i,:,:)),t_even);
    end

    % save
    curproblem.y = y;
    curproblem.weights = weights;
    curproblem.covar_velocity = covar_velocity;
    curproblem.kappa_rotrate = kappa_rotrate;
    curproblem.dt = dt;
    curproblem.p_gt = p_gt;
    curproblem.p_sd = p_sd;
    problem_list{end+1} = curproblem;
end

end

function [stamps, measurements, gt, sd] = parseBag(bagfile)

%% Keypoint Topics
% 3D keypoints (xyz) in world space
worldKeypointsTopic = "/reproject_keypoints_node/keypoints_3d_global_out";

% Pose estimate from Teaser + Fixed Lag Smoother (ie. softdrone algo)
sdPoseTopic = "/gtsam_tracker_node/target_global_odom_estimate";

% Ground truth target pose in the global space
gtPoseTopic = "/gtsam_tracker_node_secondary/target_global_odom_estimate";

%% Load Data.
bag = rosbag(bagfile);

% Load Keypoints
worldKeypointsSelect = select(bag, 'Topic', worldKeypointsTopic);
worldKeypointsMessages = readMessages(worldKeypointsSelect, 'DataFormat', 'struct');

sdPoseSelect = select(bag, 'Topic', sdPoseTopic);
sdPoseMessages = readMessages(sdPoseSelect, 'DataFormat', 'struct');

gtPoseSelect = select(bag, 'Topic', gtPoseTopic);
gtPoseMessages = readMessages(gtPoseSelect, 'DataFormat', 'struct');

%% Format Keypoint structs into arrays.

% numPixelKeypointMessages and numWorldKeypointMessages should be the same,
% but numPixelKeypointMessages may be 1 greater if the rosbag recording
% is cancelled at the right time. So we will use the world's length.
numKeypointsMessages = length(worldKeypointsMessages);

% This is how many keypoints are annotated in the CAD frame.
numKeypoints = length(worldKeypointsMessages{1}.Keypoints3D_);

worldKeypoints = zeros(3,numKeypoints, numKeypointsMessages);
worldKeypointsStamps = zeros(numKeypointsMessages, 1);

startSec = worldKeypointsMessages{1}.Header.Stamp.Sec;
startNsec = double(worldKeypointsMessages{1}.Header.Stamp.Nsec)*1e-9;
for i = 1:numKeypointsMessages
    currentWorldKeypointsMessage = worldKeypointsMessages{i};
    stamp = currentWorldKeypointsMessage.Header.Stamp;
    worldKeypointsStamps(i) = double(stamp.Sec-startSec) + double(stamp.Nsec) * 1e-9 - startNsec;
    X = [currentWorldKeypointsMessage.Keypoints3D_.X];
    Y = [currentWorldKeypointsMessage.Keypoints3D_.Y];
    Z = [currentWorldKeypointsMessage.Keypoints3D_.Z];
    for j = 1:numKeypoints
        worldKeypoints(1,j,i) = X(j);
        worldKeypoints(2,j,i) = Y(j);
        worldKeypoints(3,j,i) = Z(j);
    end
end

stamps = worldKeypointsStamps; % L x 1
measurements = worldKeypoints; % 3 x N x L

% Repeat for gt poses
N = length(gtPoseMessages);
p_gt = zeros(3,1,N);
R_gt = zeros(3,3,N);
stamps_gt = zeros(N,1);

for i = 1:N
    cur = gtPoseMessages{i};
    stamp = cur.Header.Stamp;
    stamps_gt(i) = double(stamp.Sec-startSec) + double(stamp.Nsec) * 1e-9 - startNsec;
    p_gt(1,:,i) = cur.Pose.Pose.Position.X;
    p_gt(2,:,i) = cur.Pose.Pose.Position.Y;
    p_gt(3,:,i) = cur.Pose.Pose.Position.Z;

    q = zeros(4,1);
    q(1) = cur.Pose.Pose.Orientation.W;
    q(2) = cur.Pose.Pose.Orientation.X;
    q(3) = cur.Pose.Pose.Orientation.Y;
    q(4) = cur.Pose.Pose.Orientation.Z;
    R_gt(:,:,i) = quat2rotm(q');
end
gt.p_gt = p_gt;
gt.R_gt = R_gt;
gt.stamps = stamps_gt;

% Repeat for sd poses
N = length(sdPoseMessages);
p_sd = zeros(3,1,N);
R_sd = zeros(3,3,N);
stamps_sd = zeros(N,1);

for i = 1:N
    cur = sdPoseMessages{i};
    stamp = cur.Header.Stamp;
    stamps_sd(i) = double(stamp.Sec-startSec) + double(stamp.Nsec) * 1e-9 - startNsec;
    p_sd(1,:,i) = cur.Pose.Pose.Position.X;
    p_sd(2,:,i) = cur.Pose.Pose.Position.Y;
    p_sd(3,:,i) = cur.Pose.Pose.Position.Z;

    q = zeros(4,1);
    q(1) = cur.Pose.Pose.Orientation.W;
    q(2) = cur.Pose.Pose.Orientation.X;
    q(3) = cur.Pose.Pose.Orientation.Y;
    q(4) = cur.Pose.Pose.Orientation.Z;
    R_sd(:,:,i) = quat2rotm(q');
end
sd.p_sd = p_sd;
sd.R_sd = R_sd;
sd.stamps = stamps_sd;

end