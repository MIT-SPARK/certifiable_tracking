%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through Lorenzo+GNC, frame level
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Define settings for batch processing
problem.json = "../datasets/racecar_offline/racecar_fast2.json";
problem.L = 4; % batch size
problem.savefile = "../datasets/racecar_offline/racecar_fullsize_test_ours.json";

% Set bounds based on problem setting
problem.translationBound = 10.0; % [m]
problem.velocityBound = 5.0; % [m/s]
problem.noiseBound_GNC = 0.05;
problem.noiseBound_GNC_residuals = 1;
problem.noiseBound_GRAPH = 0.01;
problem.noiseBound = 0.05;

problem.covar_measure_base = 1;
problem.covar_velocity_base = 10;
problem.covar_rotrate_base = 10;

problem.velprior = "body";       % constant body frame velocity
problem.usecBound = false;

% add shape, measurements, outliers
% load("racecar_cad.mat");
% problem.shapes = racecar_cad' / 1000; % 3 x N x K [m]
[problems, gt, teaser] = json2frameproblem(problem);
min_max_dists = robin_min_max_dists(problems{1}.shapes);

%% Solve for each batch
solns = [];
last_L = 0;
disp("Solving " + string(length(problems)) + " problems...")
for j = 1:length(problems)
% regen if batch size changes.

curproblem = problems{j};
curproblem.regen_sdp = (curproblem.L~=last_L);
last_L = curproblem.L;

% data for GNC
curproblem.type = "tracking";
curproblem.N = curproblem.N_VAR*curproblem.L; % How many measurements this problem has (updated by ROBIN)
curproblem.outliers = []; % outlier indicies
curproblem.priors = [];
curproblem.dof = 3;

curproblem = lorenzo_prune(curproblem, min_max_dists);

% preprocess inliers
% if isfield(curproblem,'prioroutliers')
%     curproblem.prioroutliers = sort(curproblem.prioroutliers);
%     curproblem.N = curproblem.N - length(curproblem.prioroutliers);
% end

% run GNC
[inliers, info] = gnc2(curproblem, @solver_for_gnc, 'barc2', curproblem.noiseBound_GNC);
disp("GNC finished " + string(j) + " (" + info.Iterations + " iterations)")

soln = info.f_info.soln;
ef = eig(soln.raw.Xopt{1});
if (ef(end-1) > 1e-4)
    disp("**Not convergent: " + string(soln.gap_stable))
end

% view_gnc(curproblem,info);

solns = [solns; soln];

if (mod(j,5) == 0)
    disp(j);
end

end
save("racecar_L4.mat","solns");

%% Check solutions
est = struct();
est.p = zeros(3,1,length(solns)+2);
est.R = zeros(3,3,length(solns)+2);
est.gap = zeros(length(solns)+2);

for j = 1:length(solns)
    problem = problems{j};
    soln = solns(j);

    % true horizon-level estimate
    if j == 1
        est.p(:,:,1:3) = soln.p_est;
        est.R(:,:,1:3) = soln.R_est;
        est.gap(1:3) = soln.gap_stable;
    else
        est.p(:,:,j+2) = soln.p_est(:,:,end);
        est.R(:,:,j+2) = soln.R_est(:,:,end);
        est.gap(j+2) = soln.gap_stable;
    end
end

est.p(:,:,est.gap > 0.01) = NaN;
est.R(:,:,est.gap > 0.01) = NaN;

%% Plot solutions
options = {est, gt, teaser};
titles = ["OURS", "Ground Truth", "Teaser"];

for i = 1:length(options)
    figure
    p = reshape(options{i}.p,[3,size(options{i}.p,3),1]);
    % p_gt = p_gt(:,12*8:15*8)
    R = options{i}.R;
    plot3(p(1,:),p(2,:),p(3,:),'.k', 'MarkerSize',10);
    hold on
    axis equal
    
    % R_est = soln.R_est;
    quiver3(p(1,:)',p(2,:)',p(3,:)',squeeze(R(1,1,:)),squeeze(R(2,1,:)),squeeze(R(3,1,:)),'r');
    quiver3(p(1,:)',p(2,:)',p(3,:)',squeeze(R(1,2,:)),squeeze(R(2,2,:)),squeeze(R(3,2,:)),'g');
    quiver3(p(1,:)',p(2,:)',p(3,:)',squeeze(R(1,3,:)),squeeze(R(2,3,:)),squeeze(R(3,3,:)),'b');
    title(titles(i))
end

%% Plot Together
p_gt = gt.p;
t = 1:length(est.p);
figure
subplot(3,1,1)
plot(t,p_gt(1,1:length(est.p)),'DisplayName','Ground Truth')
hold on
plot(t,est.p(1,:),'DisplayName','Estimate')
ylabel("x")

legend('Location','ne')
title("Explict Comparison of Evaluated Trajectories")

subplot(3,1,2)
plot(t,p_gt(2,1:length(est.p)),'DisplayName','Ground Truth')
hold on
plot(t,est.p(2,:),'DisplayName','Estimate')
ylabel("y")
subplot(3,1,3)
plot(t,p_gt(3,1:length(est.p)),'DisplayName','Ground Truth')
hold on
plot(t,est.p(3,:),'DisplayName','Estimate')
xlabel("time")
ylabel("z")


%% Save Poses into JSON
L_big = length(est.p);
T_est = repmat(eye(4),[1,1,L_big]);
for l = 1:L_big
    T_est(1:3,1:3,l) = est.R(:,:,l);
    T_est(1:3,4,l) = est.p(:,:,l)*1000.0;
end

fid = fopen(problem.json); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid);
data = jsondecode(str);

for l = 1:length(T_est)
    data(l).cast_pose = T_est(:,:,l);
end

cocoString = jsonencode(data, "PrettyPrint",true);
fid = fopen(problem.savefile, 'w');
fprintf(fid, '%s', cocoString);
fclose(fid);

%% Try to estimate R_err
% figure
% axisdif = zeros(3,L_big);
% angdif = zeros(1,L_big);
% for l = 1:L_big
%     axang = rotm2axang(est.R(:,:,l)'*gt.R(:,:,l));
%     axisdif(:,l) = axang(1:3);
%     angdif(:,l) = axang(4)*180/pi;
%     quiver3(0,0,0,axang(1),axang(2),axang(3),'r');
%     hold on
% end
