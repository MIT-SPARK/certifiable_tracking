%% Dense SDP relaxation for certifiable tracking
%  Generic, tunable script to run one iteration of dense tracking.
%    Operates on random data with no outlier support.
%    Run setup.m once to set up paths.
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Generate random tracking problem
problem.bag = "../datasets/racecar_fixed/2_2024-01-24-10-02-26.bag";
problem.L = 10; % batch size

% Set bounds based on problem setting
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.noiseBound = 0.2;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% regen if batch size changes.
problem.regen_sdp = false; % when in doubt, set to true

% add shape, measurements, outliers
load("racecar_cad.mat");
problem.shapes = racecar_cad' / 1000; % 3 x N x K
problems = bag2problem(problem);

%% Solve for each batch
solns = [];
disp("Solving " + string(length(problems)) + " problems...")
for j = 1:length(problems)
curproblem = problems{j};

soln = solve_weighted_tracking(curproblem);

% soln_pace = pace_with_EKF(problem);

ef = eig(soln.raw.Xopt{1});
if (ef(end-4) > 1e-6)
    disp("soln " + string(j) + " not convergent.")
end
solns = [solns; soln];

if problem.regen_sdp
    break;
    disp("SDP data generated. Rerun with regen_sdp true for faster results.")
end
if (mod(j,5) == 0)
    disp(j);
end
end

%% Check solutions

figure(1);
for j = 1:length(solns)

problem = problems{j};
soln = solns(j);

% eigenvalue plot
L = problem.L;
N = problem.N_VAR;
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

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

%% Quick test to see if keypoints fit
% figure
% for l = 1:L
%     for i = 1:N
%         pt = problem.y(ib3(i),l);
%         plot3(pt(1),pt(2),pt(3),'xk');
%         hold on
%     end
% end
