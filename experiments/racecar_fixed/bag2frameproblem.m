function [problem_list, gt, sd] = bag2frameproblem(problem, startTime, endTime)
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
[stamps, measurements, gt, sd] = parseBag(problem.bag, startTime, endTime);
% TODO: PROCESS GROUND TRUTH CORRECTLY
tot_L = length(stamps);

L = problem.L;

problem_list = {};

% define batch problems
for batch = 1:tot_L
    idxrange = max(1,batch-L):batch;
    L_cur = length(idxrange);
    if (L_cur < 3)
        continue
    end
    curproblem = problem;
    curproblem.startIdx = idxrange(1);
     
    % interpolate between measurements
    t_init = stamps(idxrange);
    t = t_init - t_init(1);
    m = measurements(:,:,idxrange); % 3 x N x L
    m(m==0) = NaN;
    t_even = linspace(0,t(end),L_cur);
    dt = t_even(2) - t_even(1);
    
    y = zeros(3*N,L_cur);
    for i = 1:N
        xtemp = interp1(t,reshape(m(1,i,:),[L_cur,1,1]),t_even);%,'spline');
        ytemp = interp1(t,reshape(m(2,i,:),[L_cur,1,1]),t_even);%,'spline');
        ztemp = interp1(t,reshape(m(3,i,:),[L_cur,1,1]),t_even);%,'spline');

        % FOR NO INTERPOLATION:
        % xtemp = reshape(m(1,i,:),[L,1,1])';
        % ytemp = reshape(m(2,i,:),[L,1,1])';
        % ztemp = reshape(m(3,i,:),[L,1,1])';

        y(ib3(i),:) = [xtemp;ytemp;ztemp];
    end
    % change weights to ignore nans
    prioroutliers = [];
    for l = 1:L_cur
        yl = y(:,l);
        i3_nan = strfind(isnan(yl)',true(1,3));
        for j = 1:length(i3_nan)
            i3 = i3_nan(j);
            i = (i3-1)/3 + 1;
            if (i == round(i))
                prioroutliers(end+1) = i + N*(l-1);
            end
        end
        curproblem.prioroutliers = prioroutliers;
    end
    if (length(prioroutliers) >= N*L_cur)
        disp("Batch " + batch + " failed.")
        continue
    end

    % y cannot have nans in it
    y(isnan(y)) = 0.0;

    % set covariances
    noiseBoundSq = problem.noiseBound^2;
    weights = ones(1,N*L_cur-length(prioroutliers))*((noiseBoundSq/9).^(-1));
    if (~isfield(curproblem,"covar_velocity_base"))
        covar_velocity = ones(L_cur-2,1)*weights(1).^(-1);
    else
        covar_velocity = ones(L_cur-2,1)*curproblem.covar_velocity_base;
    end
    kappa_rotrate  = ones(L_cur-2,1)*(2/covar_velocity(1));

    % Save gt and sd poses
    t_gt = gt.stamps - t_init(1);
    t_sd = sd.stamps - t_init(1);
    p_gt = zeros(3,1,L_cur);
    p_sd = zeros(3,1,L_cur);
    for i = 1:3
        p_gt(i,:,:) = interp1(t_gt,squeeze(gt.p(i,:,:)),t_even);
        p_sd(i,:,:) = interp1(t_sd,squeeze(sd.p(i,:,:)),t_even);
    end

    % save
    curproblem.L = L_cur;
    curproblem.noiseBoundSq = curproblem.noiseBound^2;
    curproblem.cBound = 1.0;
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

function [stamps, measurements, gt, sd] = parseBag(bagfile, startTime, endTime)

%% Keypoint Topics
% 3D keypoints (xyz) in world space
worldKeypointsTopic = "/reproject_keypoints_node/keypoints_3d_global_out";

% Pose estimate from Teaser + Fixed Lag Smoother (ie. softdrone algo)
gtPoseTopic = "/gtsam_tracker_node/target_global_odom_estimate";

% Ground truth target pose in the global space
sdPoseTopic = "/gtsam_tracker_node_secondary/target_global_odom_estimate"; % target_global_odom_estimate or teaser_global

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

worldKeypoints = zeros(3,numKeypoints, numKeypointsMessages)*NaN;
worldKeypointsStamps = zeros(numKeypointsMessages, 1);
stamps = [];

startSec = worldKeypointsMessages{1}.Header.Stamp.Sec;
startNsec = double(worldKeypointsMessages{1}.Header.Stamp.Nsec)*1e-9;
for i = 1:numKeypointsMessages
    currentWorldKeypointsMessage = worldKeypointsMessages{i};
    stamp = currentWorldKeypointsMessage.Header.Stamp;
    worldKeypointsStamps(i) = double(stamp.Sec-startSec) + double(stamp.Nsec) * 1e-9 - startNsec;

    % crop by start or end time
    if (worldKeypointsStamps(i) < startTime) || (worldKeypointsStamps(i) > endTime)
        continue
    end
    stamps(end+1) = worldKeypointsStamps(i);

    X = [currentWorldKeypointsMessage.Keypoints3D_.X];
    Y = [currentWorldKeypointsMessage.Keypoints3D_.Y];
    Z = [currentWorldKeypointsMessage.Keypoints3D_.Z];
    for j = 1:numKeypoints
        worldKeypoints(1,j,i) = X(j);
        worldKeypoints(2,j,i) = Y(j);
        worldKeypoints(3,j,i) = Z(j);
    end
end

% stamps = worldKeypointsStamps; % L x 1
measurements = worldKeypoints; % 3 x N x L
measurements = reshape(measurements,[3*numKeypoints*numKeypointsMessages,1]);
measurements = reshape(measurements(~isnan(measurements)),3,7,[]);

% Repeat for gt poses
N = length(gtPoseMessages);
p_gt = zeros(3,1,N)*NaN;
R_gt = zeros(3,3,N)*NaN;
stamps_gt = [];

for i = 1:N
    cur = gtPoseMessages{i};
    stamp = cur.Header.Stamp;
    s = double(stamp.Sec-startSec) + double(stamp.Nsec) * 1e-9 - startNsec;
    
    % crop by start or end time
    if (s < startTime) || (s > endTime)
        continue
    end
    stamps_gt(end+1) = s;

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
gt.p = reshape(p_gt,[3*N,1]);
gt.p = reshape(gt.p(~isnan(gt.p)),3,1,[]);
gt.R = reshape(R_gt,[3*3*N,1]);
gt.R = reshape(gt.R(~isnan(gt.R)),3,3,[]);
gt.stamps = stamps_gt;

% Repeat for sd poses
N = length(sdPoseMessages);
p_sd = zeros(3,1,N)*NaN;
R_sd = zeros(3,3,N)*NaN;
stamps_sd = [];

for i = 1:N
    cur = sdPoseMessages{i};
    stamp = cur.Header.Stamp;
    s = double(stamp.Sec-startSec) + double(stamp.Nsec) * 1e-9 - startNsec;
    
    % crop by start or end time
    if (s < startTime) || (s > endTime)
        continue
    end
    stamps_sd(end+1) = s;

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
sd.p = reshape(p_sd,[3*N,1]);
sd.p = reshape(sd.p(~isnan(sd.p)),3,1,[]);
sd.R = reshape(R_sd,[3*3*N,1]);
sd.R = reshape(sd.R(~isnan(sd.R)),3,3,[]);
sd.stamps = stamps_sd;

end