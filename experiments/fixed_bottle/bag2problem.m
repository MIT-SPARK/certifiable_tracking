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
[stamps, measurements] = parseBag(problem.bag);
tot_L = length(stamps);

batchsize = problem.batchsize;
dt = problem.dt;

problem_list = [];

for batch = 1:floor(tot_L/batchsize)
    idxrange = ((batch-1)*batchsize+1):(batch*batchsize + 1);
    curproblem = problem;
    % figure out L
    t = stamps(idxrange);
    L = round((t(end) - t(1))/dt) + 1;
    curproblem.L = L;

    tint = round((t - t(1))/dt) + 1;

    % set weights and create y
    y = zeros(3*N,L);
    covar_measure = ones(N,L)*Inf;
    covar_velocity = ones(L-2,1)*Inf; % l = 2, ..., L-1
    kappa_rotrate  = ones(L-2,1)*0.0; % l = 2, ..., L-1
    
    for j = 1:length(tint)
        l = tint(j);
        lmeasure = idxrange(j);
        measure = reshape(measurements(:,:,lmeasure)',[3*N,1]);
        y(:,l) = measure;
        covar_measure(:,l) = 1.0;
        if (l > 1 && l < L)
            covar_velocity(l-1) = 1.0;
            kappa_rotrate(l-1) = 1.0;
        end
    end
    % save
    curproblem.y = y;
    curproblem.covar_measure = covar_measure;
    curproblem.covar_velocity = covar_velocity;
    curproblem.kappa_rotrate = kappa_rotrate;
    problem_list = [problem_list; curproblem];
end

end

function [stamps, measurements] = parseBag(bagfile)

% Keypoint Topics
% 2D keypoints in the pixel space
pixelKeypointsTopic = "/keypoint_detector_node/keypoints_out";
% 3D keypoints (xyz) in world space
worldKeypointsTopic = "/reproject_keypoints_node/keypoints_3d_out";

%% Load Data.
bag = rosbag(bagfile);

% Load Keypoints
pixelKeypointsSelect = select(bag, 'Topic', pixelKeypointsTopic);
pixelKeypointsMessages = readMessages(pixelKeypointsSelect, 'DataFormat', 'struct');

worldKeypointsSelect = select(bag, 'Topic', worldKeypointsTopic);
worldKeypointsMessages = readMessages(worldKeypointsSelect, 'DataFormat', 'struct');

%% Format Keypoint structs into arrays.

% numPixelKeypointMessages and numWorldKeypointMessages should be the same,
% but numPixelKeypointMessages may be 1 greater if the rosbag recording
% is cancelled at the right time. So we will use the world's length.
numKeypointsMessages = length(worldKeypointsMessages);

% This is how many keypoints are annotated in the CAD frame.
numKeypoints = length(pixelKeypointsMessages{1}.Keypoints2D_);

pixelKeypoints = zeros(numKeypoints, 2, numKeypointsMessages);
pixelKeypointsStamps = zeros(numKeypointsMessages, 1);

worldKeypoints = zeros(numKeypoints, 3, numKeypointsMessages);
worldKeypointsStamps = zeros(numKeypointsMessages, 1);

startSec = pixelKeypointsMessages{1}.Header.Stamp.Sec;
startNsec = double(pixelKeypointsMessages{1}.Header.Stamp.Nsec)*1e-9;
for i = 1:numKeypointsMessages

    currentPixelKeypointsMessage = pixelKeypointsMessages{i};
    stamp = currentPixelKeypointsMessage.Header.Stamp;
    pixelKeypointsStamps(i) = stamp.Sec + stamp.Nsec * 1e-9;
    X = [currentPixelKeypointsMessage.Keypoints2D_.X];
    Y = [currentPixelKeypointsMessage.Keypoints2D_.Y];
    for j = 1:numKeypoints
        pixelKeypoints(j,1,i) = X(j);
        pixelKeypoints(j,2,i) = Y(j);
    end

    currentWorldKeypointsMessage = worldKeypointsMessages{i};
    stamp = currentPixelKeypointsMessage.Header.Stamp;
    worldKeypointsStamps(i) = double(stamp.Sec-startSec) + double(stamp.Nsec) * 1e-9 - startNsec;
    X = [currentWorldKeypointsMessage.Keypoints3D_.X];
    Y = [currentWorldKeypointsMessage.Keypoints3D_.Y];
    Z = [currentWorldKeypointsMessage.Keypoints3D_.Z];
    for j = 1:numKeypoints
        worldKeypoints(j,1,i) = X(j);
        worldKeypoints(j,2,i) = Y(j);
        worldKeypoints(j,3,i) = Z(j);

    end
end

stamps = worldKeypointsStamps;
measurements = worldKeypoints / 1000;

end