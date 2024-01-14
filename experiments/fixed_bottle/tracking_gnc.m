%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through ROBIN+GNC
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Define settings for batch processing
problem.bag = "2024-01-11-17-30-52.bag";
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
load("cad_frame.mat");
problem.shapes = annotatedPointsWrtTarget' / 1000; % 3 x N x K
problems = bag2problem(problem);

%% Solve for each batch
solns = {};
for j = 1:length(problems)
curproblem = problems{j};

% data for GNC
curproblem.type = "tracking";
curproblem.N = curproblem.N_VAR*curproblem.L; % How many measurements this problem has (updated by ROBIN)
curproblem.outliers = []; % outlier indicies
curproblem.priors = [];
curproblem.dof = 3;

% preprocess inliers
if isfield(curproblem,'prioroutliers')
    curproblem.prioroutliers = sort(curproblem.prioroutliers);
    curproblem.N = curproblem.N - length(curproblem.prioroutliers);
end

% run GNC
[inliers, info] = gnc_custom(curproblem, @solver_for_gnc, 'NoiseBound', curproblem.noiseBound,'MaxIterations',100,'FixPriorOutliers',false);

solns{end+1} = info.f_info.soln;

disp("GNC finished " + string(j))

end

%% Check solutions
figure(1);
for j = 1:length(solns)

problem = problems{j};
soln = solns{j};

% eigenvalue plot
L = problem.L;
N = problem.N_VAR;
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

% Plot trajectory!
figure(1);
axis equal
p_est = reshape(soln.p_est,[3,L,1]);
plot3(p_est(1,:),p_est(2,:),p_est(3,:),'.k', 'MarkerSize',10);

p_quiv = repelem(p_est,3,1);
R_est = soln.R_est;
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,1,:)),squeeze(R_est(2,1,:)),squeeze(R_est(3,1,:)),'r');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,2,:)),squeeze(R_est(2,2,:)),squeeze(R_est(3,2,:)),'g');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,3,:)),squeeze(R_est(2,3,:)),squeeze(R_est(3,3,:)),'b');
hold on

end