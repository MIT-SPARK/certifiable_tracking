function [problem_list, gt, teaser, shapes] = json2frameproblem(problem, skip)
%% generates a problem from json metadata file at the frame level
% horizon-based: for each frame, optimize with last L frames
% 
% Lorenzo Shaikewitz for SPARK Lab

if nargin < 2
    skip = 1;
end

%% Load and parse bag data
[stamps, measurements, gt, teaser, shapes] = parseJson(problem.json, problem, skip);
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

    % set covariances
    % noiseBoundSq = problem.noiseBound^2;
    % weights = ones(1,N*L_cur)*((noiseBoundSq/9).^(-1));
    % if (isfield(problem,"covar_velocity_base"))
    %     covar_velocity = ones(L_cur-2,1)*problem.covar_velocity_base;
    % else
    %     covar_velocity = ones(L_cur-2,1)*weights(1).^(-1);
    % end
    % if (isfield(problem,"kappa_rotrate_base"))
    %     kappa_rotrate = ones(L-2,1)*problem.kappa_rotrate_base;
    % else
    %     kappa_rotrate  = ones(L-2,1)*(2/covar_velocity(1));
    % end

    % save
    curproblem.L = L_cur;
    curproblem.y = y;
    % curproblem.weights = weights;
    % curproblem.covar_velocity = covar_velocity;
    % curproblem.kappa_rotrate = kappa_rotrate;
    curproblem.dt = 1;%dt;

    problem_list{end+1} = curproblem;
end

end

function [stamps, keypoints, gt, teaser, shapes] = parseJson(jsonfile, problem, skip)

% Open the json file
fid = fopen(jsonfile);
raw = fread(fid,Inf);
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);

% get CAD keypoints if there
if (problem.object == "cracker") || (problem.object == "sugar")
    load("../datasets/YCBInEOAT/shapes/shapes_box.mat","shapes");
    shapes = shapes(:,:,1:2); % TEMP
elseif (problem.object == "mustard") || (problem.object == "bleach")
    load("../datasets/YCBInEOAT/shapes/shapes_bottle.mat","shapes");
elseif (problem.object == "tomato")
    load("../datasets/YCBInEOAT/shapes/shapes_can.mat","shapes");
end
N = size(shapes,2);
field = "est_world_keypoints";


% if isfield(data, "interp_cad_keypoints")
%     shapes = data(1).interp_cad_keypoints' / 1000.0;
%     % REMOVE INTERP:
%     shapes = shapes(:,1:52);
% 
%     % [C, ia, ic] = unique(shapes','stable','rows');
%     % % remove nonspherical keypoints/shapes
%     % % keep = [6,7,8,10,12,14,15,16];
%     % keep = [6,8,14,16];
%     % keepVec = ia;
%     % for n = 1:length(ia)
%     %     keepVec = [keepVec; length(ia) + 16*(n-1) + keep'];
%     % end
%     % shapes = shapes(:,keepVec,:);
% 
%     N = size(shapes,2);
%     % field = "est_interp_world_keypoints";
%     field = "est_world_keypoints";
% else
%     shapes = NaN;
%     N = size(problem.shapes,2);
%     field = "est_world_keypoints";
% end

keypoints = [data.(field)];
% if isfield(data,"interp_cad_keypoints")
%     keypoints = reshape(keypoints,[length(ic),3,size(data,1)]);
%     keypoints = keypoints(keepVec,:,:);
% end

keypoints = reshape(keypoints,[N,3,size(data,1)]);
keypoints = permute(keypoints,[2,1,3]) / 1000.0; % [m]

keypoints = keypoints(:,:,1:skip:end);

bigL = size(data,1);

% gt poses
poses = reshape([data.gt_teaser_pose],[4,4,bigL]);
p_gt = poses(1:3,4,:) / 1000; % [m]
R_gt = poses(1:3,1:3,:);
gt.p = p_gt(:,:,1:skip:end);
gt.R = R_gt(:,:,1:skip:end);

% Teaser poses
poses = reshape([data.est_teaser_pose],[4,4,bigL]);
p_teaser = poses(1:3,4,:) / 1000; % [m]
R_teaser = poses(1:3,1:3,:);
teaser.p = p_teaser(:,:,1:skip:end);
teaser.R = R_teaser(:,:,1:skip:end);

% make up stamps
stamps = 0:(1/30):(1/30*(bigL-1));
stamps = stamps(1:skip:end);

end