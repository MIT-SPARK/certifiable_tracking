%% Dense SDP relaxation for certifiable tracking
%  Generic, tunable script to run one iteration of dense tracking.
%    Operates on racecar data with no outlier support.
%    Run setup.m once to set up paths.
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Generate random tracking problem
problem.bag = "../datasets/racecar_fixed/2024-01-30-18-14-08.bag";
% problem.bag = "../datasets/racecar_fixed/2_2024-01-24-10-02-26.bag";
problem.L = 10; % batch size

% Set bounds based on problem setting
problem.translationBound = 5.0;
problem.velocityBound = 2.0;
problem.noiseBound = 0.03;
problem.covar_velocity_base = 0.27^2;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% add shape, measurements, outliers
load("racecar_cad.mat");
problem.shapes = racecar_cad' / 1000; % 3 x N x K
[problems, gt, sd] = bag2gtproblem(problem, 15, 50.0);

%% Solve for each batch
solns = [];
disp("Solving " + string(length(problems)) + " problems...")
for j = 1:length(problems)
curproblem = problems{j};
curproblem.regen_sdp = (j == 1);

soln = solve_weighted_tracking(curproblem);

% soln_pace = pace_py_UKF(curproblem);

ef = eig(soln.raw.Xopt{1});
if (ef(end-4) > 1e-6)
    disp("soln " + string(j) + " not convergent: " + string(soln.gap))
end
solns = [solns; soln];

if (mod(j,5) == 0)
    disp(j);
end
end

%% Check solutions

L = problem.L;
p_err = zeros(L*length(solns),1);
R_err = zeros(L*length(solns),1);

figure(1);
for j = 1:length(solns)

problem = problems{j};
soln = solns(j);

% eigenvalue plot
L = problem.L;
N = problem.N_VAR;
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

idx = ((j-1)*L + 1):j*L;
for l = 1:L
    p_err(idx(l)) = norm(soln.p_est(:,:,l) - gt.p(:,:,idx(l)));
    R_err(idx(l)) = getAngularError(gt.R(:,:,idx(l)),soln.R_est(:,:,l));
end

% Plot trajectory!
figure(1);
p_est = reshape(soln.p_est,[3,L,1]);
plot3(p_est(1,:),p_est(2,:),p_est(3,:),'.k', 'MarkerSize',10);
hold on
axis equal

p_quiv = repelem(p_est,3,1);
R_est = soln.R_est;
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,1,:)),squeeze(R_est(2,1,:)),squeeze(R_est(3,1,:)),'r');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,2,:)),squeeze(R_est(2,2,:)),squeeze(R_est(3,2,:)),'g');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,3,:)),squeeze(R_est(2,3,:)),squeeze(R_est(3,3,:)),'b');

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
title("TEASER")
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
