function [problem_list, gt, teaser] = json2frameproblem(problem)
%% generates a problem from json metadata file
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
[stamps, measurements, gt, teaser] = parseJson(problem.json, N);
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
    noiseBoundSq = problem.noiseBound^2;
    weights = ones(N*L_cur-length(prioroutliers),1)*((noiseBoundSq/9).^(-1));
    covar_velocity = ones(L_cur-2,1)*weights(1)*1;
    kappa_rotrate  = ones(L_cur-2,1)*(2/covar_velocity(1));

    % save
    curproblem.L = L_cur;
    curproblem.y = y;
    curproblem.weights = weights;
    curproblem.covar_velocity = covar_velocity;
    curproblem.kappa_rotrate = kappa_rotrate;
    curproblem.dt = dt;

    problem_list{end+1} = curproblem;
end

end

function [stamps, keypoints, gt, teaser] = parseJson(jsonfile, N)

% Open the json file
fid = fopen(jsonfile);
raw = fread(fid,Inf);
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);
keypoints = [data.est_world_keypoints];
keypoints = reshape(keypoints,[N,3,size(data,1)]);
keypoints = permute(keypoints,[2,1,3]) / 1000.0; % [m]

bigL = size(data,1);

% gt poses
poses = reshape([data.gt_teaser_pose],[4,4,bigL]);
p_gt = poses(1:3,4,:) / 1000; % [m]
R_gt = poses(1:3,1:3,:);
gt.p = p_gt;
gt.R = R_gt;

% Teaser poses
poses = reshape([data.est_teaser_pose],[4,4,bigL]);
p_teaser = poses(1:3,4,:) / 1000; % [m]
R_teaser = poses(1:3,1:3,:);
teaser.p = p_teaser;
teaser.R = R_teaser;

% make up stamps
stamps = 0:(1/30):(1/30*(bigL-1));

end