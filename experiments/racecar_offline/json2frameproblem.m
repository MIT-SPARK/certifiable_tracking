function [problem_list, gt, teaser] = json2frameproblem(problem)
%% generates a problem from json metadata file
% frame-based: optimize each frame by the last L frames.
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Load and parse bag data
[stamps, measurements, gt, teaser, shapes] = parseJson(problem.json);
if (~isnan(shapes))
    problem.shapes = shapes;
end

%% Define shape data
% shapes is 3 x N x K
problem.N_VAR = size(problem.shapes,2);
problem.K = size(problem.shapes,3);
N = problem.N_VAR;
K = problem.K;

% set lambda if K > N
if K > N
    problem.lambda = 1.0;
    disp("Set lambda to " + string(problem.lambda));
else
    problem.lambda = 0.0;
end
problem.B = reshape(problem.shapes, 3*N, K);

%% Parse into problem format
tot_L = length(stamps);

L = problem.L;

problem_list = {};

% define batch problems
for batch = 1:tot_L
    
    idxrange = max(1,batch-L+1):batch;
    L_cur = length(idxrange);
    if (L_cur < 3)
        continue
    end
    curproblem = problem;
    curproblem.startIdx = idxrange(1);
     
    % interpolate between measurements
    t = stamps(idxrange);
    t = t - t(1);
    m = measurements(:,:,idxrange); % 3 x N x L
    m(m==0) = NaN;
    t_even = linspace(0,t(end),L_cur);
    dt = t_even(2) - t_even(1);
    
    y = zeros(3*N,L_cur);
    for i = 1:N
        xtemp = interp1(t,reshape(m(1,i,:),[L_cur,1,1]),t_even);
        ytemp = interp1(t,reshape(m(2,i,:),[L_cur,1,1]),t_even);
        ztemp = interp1(t,reshape(m(3,i,:),[L_cur,1,1]),t_even);
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

    % save
    curproblem.L = L_cur;
    curproblem.y = y;
    % curproblem.weights = weights;
    % curproblem.covar_velocity = covar_velocity;
    % curproblem.kappa_rotrate = kappa_rotrate;
    % curproblem.dt = dt;
    curproblem.dt = 1.0; % USE DT=1 for GREATER ROBUSTNESS

    problem_list{end+1} = curproblem;
end

end

function [stamps, keypoints, gt, teaser, shapes] = parseJson(jsonfile)

% Open the json file
fid = fopen(jsonfile);
raw = fread(fid,Inf);
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);

% get CAD keypoints if there
if isfield(data, "interp_cad_keypoints")
    shapes = data(1).interp_cad_keypoints' / 1000.0;
    shapes = shapes(:,1:7); % not interp
    N = size(shapes,2);
    % shapes = NaN;
else
    shapes = NaN;
end

keypoints = [data.est_world_keypoints];
% keypoints = [data.est_interp_world_keypoints];
keypoints = reshape(keypoints,[N,3,size(data,1)]);
keypoints = permute(keypoints,[2,1,3]) / 1000.0; % [m]

bigL = size(data,1);

% transform keypoints
cam_wrt_world = [data.cam_wrt_world];
cam_wrt_world = reshape(cam_wrt_world, [4,4,size(data,1)]);
cam_wrt_world(1:3,4,:) = cam_wrt_world(1:3,4,:) / 1000.0; % [m]
for l = 1:bigL
    kpts = keypoints(:,:,l);
    T = cam_wrt_world(:,:,l);
    kpts = [kpts; ones(1,N)];
    kpts = T * kpts;
    keypoints(:,:,l) = kpts(1:3,:);
end

% gt poses
poses = reshape([data.gt_teaser_pose],[4,4,bigL]);
poses(1:3,4,:) = poses(1:3,4,:) / 1000.0; % [m]
for l = 1:bigL
    T = cam_wrt_world(:,:,l);
    poses(:,:,l) = T*poses(:,:,l);
end
p_gt = poses(1:3,4,:); % [m]
R_gt = poses(1:3,1:3,:);
gt.p = p_gt;
gt.R = R_gt;

% Teaser poses
poses = reshape([data.est_teaser_pose],[4,4,bigL]);
poses(1:3,4,:) = poses(1:3,4,:) / 1000.0; % [m]
for l = 1:bigL
    T = cam_wrt_world(:,:,l);
    poses(:,:,l) = T*poses(:,:,l);
end
p_teaser = poses(1:3,4,:); % [m]
R_teaser = poses(1:3,1:3,:);
teaser.p = p_teaser;
teaser.R = R_teaser;

% make up stamps
% stamps = 0:(1/30):(1/30*(bigL-1));
stamps = [data.timestamp];
stamps = stamps - stamps(1);

% gt poses to velocity estimate
% v_gt = zeros(3,1,bigL);
% for l = 1:bigL-1
%     v_gt(:,:,l) = (p_gt(:,:,l+1) - p_gt(:,:,l)) / (stamps(l+1) - stamps(l));
% end

end