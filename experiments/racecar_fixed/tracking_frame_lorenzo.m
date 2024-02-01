%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through Lorenzo+GNC, frame level
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Define settings for batch processing
problem.bag = "../datasets/racecar_fixed/2024-01-30-18-14-08.bag";
problem.L = 10; % batch size
times = [32, 40];

% Set bounds based on problem setting
problem.translationBound = 5.0; % [m]
problem.velocityBound = 1.0; % [m/s]
problem.noiseBound_GNC = 0.05;
problem.noiseBound_GRAPH = 0.05;
problem.noiseBound = 0.05;
problem.covar_velocity_base = 0.001^2;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% add shape, measurements, outliers
load("racecar_cad.mat");
problem.shapes = racecar_cad' / 1000; % 3 x N x K [m]
min_max_dists = robin_min_max_dists(problem.shapes);
[problems, gt, sd] = bag2frameproblem(problem, times(1), times(2)); % 15 -> 50?

%% Solve for each batch
solns = [];
solns_pace = [];
disp("Solving " + string(length(problems)) + " problems...")
last_L = 0;
for j = 1:length(problems)
% regen if batch size changes.

curproblem = problems{j};
curproblem.regen_sdp = (curproblem.L ~= last_L); % regen only first time
last_L = curproblem.L;

% data for GNC
curproblem.type = "tracking";
curproblem.N = curproblem.N_VAR*curproblem.L; % How many measurements this problem has (updated by ROBIN)
curproblem.outliers = []; % outlier indicies
curproblem.priors = [];
curproblem.dof = 3;

% soln_pace = pace_py_UKF(curproblem,true,true);

curproblem = lorenzo_prune(curproblem, min_max_dists);

% preprocess inliers
% if isfield(curproblem,'prioroutliers')
%     curproblem.prioroutliers = sort(curproblem.prioroutliers);
%     curproblem.N = curproblem.N - length(curproblem.prioroutliers);
% end

% run GNC
try
    [inliers, info] = gnc_custom(curproblem, @solver_for_gnc, 'NoiseBound', curproblem.noiseBound_GNC,'MaxIterations',100,'FixPriorOutliers',false);
    disp("GNC finished " + string(j))

    soln = info.f_info.soln;
    ef = eig(soln.raw.Xopt{1});
    if (ef(end-4) > 1e-4)
        disp("**Not convergent: " + string(soln.gap_nov))
    end
catch
    f = fieldnames(solns(1))';
    f{2,1} = {NaN};
    soln = struct(f{:});
    soln.p_est = ones(3,1,curproblem.L)*NaN;
    soln.R_est = ones(3,3,curproblem.L)*NaN;
    info.f_info.soln = soln;
    disp("GNC failed " + string(j))
end

solns = [solns; soln];
% solns_pace = [solns_pace; solns_pace];

if (mod(j,5) == 0)
    disp(j);
end

end

%% Check solutions
L = problem.L;
N = problem.N_VAR;

p_err = zeros(L*N,L)*NaN;
R_err = zeros(L*N,L)*NaN;

figure(1);
for j = 1:length(solns)

problem = problems{j};
soln = solns(j);
L_cur = problem.L;

if (soln.gap_nov > 0.1)
    continue
end

idx = problem.startIdx:(problem.startIdx + L_cur);
for l = 1:L_cur
    p_err(idx(l),l) = norm(soln.p_est(:,:,l) - problem.p_gt(:,:,l));
    % R_err(idx(l)) = getAngularError(problem.R_gt(:,:,l),soln.R_est(:,:,l));
end

% eigenvalue plot
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

% Plot trajectory!
figure(1);
p_est = reshape(soln.p_est,[3,L_cur,1]);
plot3(p_est(1,:),p_est(2,:),p_est(3,:),'.k', 'MarkerSize',10);
hold on
axis equal

p_quiv = repelem(p_est,3,1);
R_est = soln.R_est;
% quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,1,:)),squeeze(R_est(2,1,:)),squeeze(R_est(3,1,:)),'r');
% quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,2,:)),squeeze(R_est(2,2,:)),squeeze(R_est(3,2,:)),'g');
% quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,3,:)),squeeze(R_est(2,3,:)),squeeze(R_est(3,3,:)),'b');

end
title("OURS")
view(0,90)


%% Plot Ground Truth

figure
p_gt = reshape(gt.p,[3,size(gt.p,3),1]);
plot3(p_gt(1,:),p_gt(2,:),p_gt(3,:),'.k', 'MarkerSize',10);
hold on
axis equal

% R_est = soln.R_est;
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(gt.R(1,1,:)),squeeze(gt.R(2,1,:)),squeeze(gt.R(3,1,:)),'r');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(gt.R(1,2,:)),squeeze(gt.R(2,2,:)),squeeze(gt.R(3,2,:)),'g');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(gt.R(1,3,:)),squeeze(gt.R(2,3,:)),squeeze(gt.R(3,3,:)),'b');
title("Ground Truth")
view(0,90)

% and TEASER
figure
p_sd = reshape(sd.p,[3,size(sd.p,3),1]);
plot3(p_sd(1,:),p_sd(2,:),p_sd(3,:),'.k', 'MarkerSize',10);
hold on
axis equal

% R_est = soln.R_est;
quiver3(p_sd(1,:)',p_sd(2,:)',p_sd(3,:)',squeeze(sd.R(1,1,:)),squeeze(sd.R(2,1,:)),squeeze(sd.R(3,1,:)),'r');
quiver3(p_sd(1,:)',p_sd(2,:)',p_sd(3,:)',squeeze(sd.R(1,2,:)),squeeze(sd.R(2,2,:)),squeeze(sd.R(3,2,:)),'g');
quiver3(p_sd(1,:)',p_sd(2,:)',p_sd(3,:)',squeeze(sd.R(1,3,:)),squeeze(sd.R(2,3,:)),squeeze(sd.R(3,3,:)),'b');
title("Soft Drone")
view(0,90)

%% Quick test to see if keypoints fit
% figure
% for l = 1:L
%     for i = 1:N
%         pt = problem.y(ib3(i),l);
%         plot3(pt(1),pt(2),pt(3),'xk');
%         hold on
%     end
% end

% figure
% for j = 1:length(problems)
%     problem = problems{j};
%     keypoints = reshape(problem.y,[3,N,problem.L]);
%     for l = 1:problem.L
%         k = keypoints(:,:,l);
%         m = mean(k,2);
%         plot3(m(1),m(2),m(3),'rx','LineWidth',2)
%         hold on
%     end
% end

%% SAVE
save("drone.mat","solns")