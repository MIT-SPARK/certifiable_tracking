%% Batch-Level Script for Outlier-Free Tracking on YCBInEOAT Keypoint Data
%    Operates on YCBInEOAT keypoints with no outlier support.
%    Run setup.m once to set up paths.
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Define settings for batch processing
problem.json = "../datasets/ycbineoat/cracker_box_yalehand0_metrics.json";
problem.L = 10; % batch size

% Set bounds based on problem setting
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.noiseBound = 0.01;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% add shape, measurements, outliers
load("../datasets/ycbineoat/cheese.mat");
problem.shapes = annotatedPoints' / 1000; % 3 x N x K [m]
[problems, gt] = json2batchproblem(problem);

%% Solve for each batch
solns = [];
disp("Solving " + string(length(problems)) + " problems...")
for j = 1:length(problems)
curproblem = problems{j};
curproblem.regen_sdp = (j == 1);

soln = solve_weighted_tracking(curproblem);

% soln_pace = pace_with_EKF(problem);

ef = eig(soln.raw.Xopt{1});
if (ef(end-4) > 1e-6)
    disp("soln " + string(j) + " not convergent.")
end
solns = [solns; soln];

if (mod(j,5) == 0)
    disp(j);
end
end

%% Check solutions
L = problem.L;
N = curproblem.N_VAR;

p_err = zeros(L*length(solns),1);
R_err = zeros(L*length(solns),1);

figure(1);
for j = 1:length(solns)

problem = problems{j};
soln = solns(j);

% eigenvalue plot
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

% Plot trajectory!
figure(1);
axis equal
p_est = reshape(soln.p_est,[3,L,1]);
plot3(p_est(1,:),p_est(2,:),p_est(3,:),'.k', 'MarkerSize',10);
hold on

R_est = soln.R_est;
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,1,:)),squeeze(R_est(2,1,:)),squeeze(R_est(3,1,:)),'r');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,2,:)),squeeze(R_est(2,2,:)),squeeze(R_est(3,2,:)),'g');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,3,:)),squeeze(R_est(2,3,:)),squeeze(R_est(3,3,:)),'b');

idx = ((j-1)*L + 1):j*L;
for l = 1:L
    p_err(idx(l)) = norm(soln.p_est(:,:,l) - gt.p(:,:,idx(l)));
    R_err(idx(l)) = getAngularError(gt.R(:,:,idx(l)),soln.R_est(:,:,l));
end

end

%% Plot Ground Truth

figure
p_gt = reshape(gt.p,[3,size(gt.p,3),1]);
plot3(p_gt(1,:),p_gt(2,:),p_gt(3,:),'.k', 'MarkerSize',10);
hold on
axis equal

R_est = soln.R_est;
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(gt.R(1,1,:)),squeeze(gt.R(2,1,:)),squeeze(gt.R(3,1,:)),'r');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(gt.R(1,2,:)),squeeze(gt.R(2,2,:)),squeeze(gt.R(3,2,:)),'g');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(gt.R(1,3,:)),squeeze(gt.R(2,3,:)),squeeze(gt.R(3,3,:)),'b');

