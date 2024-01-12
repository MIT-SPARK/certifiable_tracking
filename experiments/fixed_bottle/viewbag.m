%% Topic Definitions.

% Change this to point to rosbag
bagName = "2024-01-11-17-30-52.bag";

% Image Topics
rgbImageTopic = "/target_cam/color/image_raw";
% This topic has keypoints annotated on it (green) directly from ros.
annotatedRgbImageTopic = "/keypoint_detector_node/annotated_img_out";
depthImageTopic = "/target_cam/depth/image_rect_raw/compressed";

% Keypoint Topics
% 2D keypoints in the pixel space
pixelKeypointsTopic = "/keypoint_detector_node/keypoints_out";
% 3D keypoints (xyz) in world space
worldKeypointsTopic = "/reproject_keypoints_node/keypoints_3d_out";

%% Load Data.

bag = rosbag(bagName);

% Load Images
rgbImageSelect = select(bag, 'Topic', rgbImageTopic);
rgbImageMessages = readMessages(rgbImageSelect, 'DataFormat', 'struct');

annotatedRgbImageSelect = select(bag, 'Topic', annotatedRgbImageTopic);
annotatedRgbImageMessages = readMessages(annotatedRgbImageSelect, 'DataFormat', 'struct');

depthImageSelect = select(bag, 'Topic', depthImageTopic);
depthImageMessages = readMessages(depthImageSelect, 'DataFormat', 'struct');

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


%% Loop through RGB images and look at keypoints

fig = figure;
axes = gca;

for i = 1:numKeypointsMessages

    % I found that the rgb images are recorded one frame behind (ahead?) of
    % keypoint messages so adding 1 here fixes sync issues.
    currentRgbImage = rosReadImage(annotatedRgbImageMessages{i + 1});

    % See if keypoints messages are close to the pre-annotated ones.
    currentPixelKeypointsMessage = pixelKeypointsMessages{i};
    for j = 1:numKeypoints
        currentRgbImage = insertShape(currentRgbImage, "FilledCircle", ...
                                      [pixelKeypoints(j,1,i), pixelKeypoints(j,2,i), 3], ...
                                      "Color", "red", "Opacity", 1);
    end

    % disp(currentPixelKeypoints);
    imshow(currentRgbImage, 'Parent', axes);
    drawnow;
end
%% Generate a ground truth trajectory.

% We know in reality the object is traveling in a circle and its x-axis is tangent to
% the circle.

% First, lets plot the trajectory of the centroid from the 3D world
% keypoints.

worldKeypoints(worldKeypoints==0) = NaN;
centroids = mean(worldKeypoints, 1, 'omitnan');

fig = figure;
plot3(squeeze(centroids(:,1,:)), ...
      squeeze(centroids(:,2,:)), ...
      squeeze(centroids(:,3,:)));
axis equal;

% The pose should be tangent to this ellipse




