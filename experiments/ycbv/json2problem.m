function [problem_list, gt] = json2problem(problem)
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
[stamps, measurements, gt] = parseJson(problem.json, N);
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

function [stamps, measurements, gt] = parseJson(jsonfile, N)

% Open the json file
fid = fopen(jsonfile);
raw = fread(fid,Inf);
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);
annotations = data.annotations;

stamps = [annotations.timestamp]'; % L x 1
stamps = stamps - stamps(1);
bigL = length(stamps);

measurements = reshape([annotations.ground_truth_keypoints],[N,3,bigL]);
measurements = permute(measurements,[2,1,3]); % 3 x N x L
measurements = measurements / 1000; % [m]

% gt poses
p_gt = reshape([annotations.ground_truth_position],[3,1,bigL]) / 1000;
quat_gt = reshape([annotations.ground_truth_rotation],[4,1,bigL]);
R_gt = zeros(3,3,bigL);
for l = 1:bigL
    R_gt(:,:,l) = quat2rotm(quat_gt(:,:,l)');
end
gt.p = p_gt;
gt.R = R_gt;

end