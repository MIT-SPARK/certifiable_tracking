function [problem_list, gt] = gtkeypoints2problem(problem, path)
%% generates a problem from ground truth generated keypoints
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Define shape data
load(path.cad);
problem.shapes = cad_keypoints / 1000; % 3 x N x K [m]

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

%% Define ground truth
problem.c_gt = 1; %(TODO)

%% Load and parse into problem format
[gt, stamps, depths, intrinsics] = importFromFiles(path); % slow
measurements = genMeasurements(problem, gt, depths, intrinsics);

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
            if (i == round(i))
                prioroutliers(end+1) = i + N*(l-1);
            end
        end
        curproblem.prioroutliers = prioroutliers;
    end
    if (length(prioroutliers) >= N*L)
        disp("Batch " + batch + " failed.")
        continue
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

function [gt, stamps, depths, intrinsics] = importFromFiles(path)
% Intrinsics
intrinsics = readmatrix(path.K);

% poses & depth & rgb
poseFiles = dir(path.poses);
depthFiles = dir(path.depth);
% rgbFiles = dir(path.rgb);
L = length(poseFiles);

poses = {};
depths = {};
% rgb = {};
gt.p = zeros(3,1,L);
gt.R = zeros(3,3,L);
for l = 1:L
    p = strcat(poseFiles(l).folder, '/', poseFiles(l).name);
    d = strcat(depthFiles(l).folder, '/', depthFiles(l).name);
    % r = strcat(rgbFiles(l).folder, '/', rgbFiles(l).name);
    
    % this is slow
    poses{l} = readmatrix(p);
    depths{l} = cast(imread(d), 'double');
    % rgb{l} = imread(r);

    % convert to m
    poses{l}(1:3,4) = poses{l}(1:3,4);
    depths{l} = depths{l} / 1000.0;

    % save the poses in gt format
    gt.p(:,:,l) = poses{l}(1:3,4);
    gt.R(:,:,l) = poses{l}(1:3,1:3);
end

% just make up the stamps
HZ = 1/30;
stamps = 0:HZ:(HZ*(L-1));

end

%% Generate measurements
% TODO: incorporate problem noise info.
function measurements = genMeasurements(problem, gt, depths, intrinsics)
    occludedThresh = 0.025; % [m]
    N = problem.N_VAR;
    L = length(depths);

    cad_keypoints_tilde = [problem.shapes*problem.c_gt; ones(1,N)];
    
    % convert CAD keypoints into camera frame
    cam_keypoints = zeros(4,N,L);
    pixel_keypoints = zeros(2,N,L);
    occludedPoints = ones(N,L);
    keypoints_simulated = zeros(3,N,L);
    for l = 1:L
        pose = zeros(4);
        pose(4,4) = 1;
        pose(1:3,1:3) = gt.R(:,:,l);
        pose(1:3,4) = gt.p(:,:,l);
        cam_keypoints(:,:,l) = pose*cad_keypoints_tilde;
    
        % convert camera frame keypoints to pixels
        pixel_keypoints_temp = [intrinsics, zeros(3,1)] * cam_keypoints(:,:,l);
        pixel_keypoints(:,:,l) = rdivide(pixel_keypoints_temp(1:2,:),pixel_keypoints_temp(3,:));
    
        % check for occulded points and export
        for i = 1:N
            px = round(pixel_keypoints(1,i,l));
            py = round(pixel_keypoints(2,i,l));
    
            occluded = (px <= 0) || (px > size(depths{l},2));
            occluded = occluded || (py <= 0) || (py > size(depths{l},1));
    
            % use depth for in-image outliers
            if ~(occluded)
                z_image = depths{l}(py, px);
                occluded = (abs(z_image - cam_keypoints(3,i,l)) > occludedThresh);
            end
            occludedPoints(i,l) = occluded;
    
            % export
            if ~(occluded)
                keypoints_simulated(:,i,l) = transformPixelsToCameraFrame(px,py,z_image,intrinsics);
            end
        end
    end
    measurements = keypoints_simulated;
end

function pointWrtCam = transformPixelsToCameraFrame(px, py, z, intrinsics)
    K = intrinsics;
    fx = K(1,1);
    fy = K(2,2);
    cx = K(3,1);
    cy = K(3,2);
    x = (px - cx) * z / fx;
    y = (py - cy) * z / fy; 
    pointWrtCam = [x, y, z];
end