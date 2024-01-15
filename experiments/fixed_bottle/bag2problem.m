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

    % save
    curproblem.y = y;
    curproblem.weights = weights;
    curproblem.covar_velocity = covar_velocity;
    curproblem.kappa_rotrate = kappa_rotrate;
    curproblem.dt = dt;
    problem_list{end+1} = curproblem;
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

worldKeypoints = zeros(3,numKeypoints, numKeypointsMessages);
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
        worldKeypoints(1,j,i) = X(j);
        worldKeypoints(2,j,i) = Y(j);
        worldKeypoints(3,j,i) = Z(j);

    end
end

stamps = worldKeypointsStamps;
measurements = worldKeypoints / 1000.0;

end