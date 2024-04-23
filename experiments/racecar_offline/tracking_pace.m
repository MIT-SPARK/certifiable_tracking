%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through Lorenzo+GNC, batch level
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Define settings for batch processing
problem.json = "../datasets/racecar_offline/racecar_fast2.json";
problem.L = 3; % batch size
problem.savefile = "../datasets/racecar_offline/racecar_fullsize_test_ours.json";

% Set bounds based on problem setting
problem.translationBound = 10.0; % [m]
problem.velocityBound = 5.0; % [m/s]
problem.noiseBound_GNC = 0.05;
problem.noiseBound_GNC_residuals = 1;
problem.noiseBound_GRAPH = 0.01;
problem.noiseBound = 0.05;

problem.cBound = 1;
problem.noiseBoundSq = problem.noiseBound^2;

problem.covar_measure_base = 1;
problem.covar_velocity_base = 10;
problem.covar_rotrate_base = 10;

problem.velprior = "body";       % constant body frame velocity
problem.usecBound = false;

% add shape, measurements, outliers
% load("racecar_cad.mat");
% problem.shapes = racecar_cad' / 1000; % 3 x N x K [m]
[problems, gt, teaser] = json2batchproblem(problem);
min_max_dists = robin_min_max_dists(problems{1}.shapes);

%% Solve for each batch
solns = [];
disp("Solving " + string(length(problems)) + " problems...")
for j = 1:length(problems)
% regen if batch size changes.

curproblem = problems{j};
curproblem.regen_sdp = (j==1); % when in doubt, set to true

% run pace with GNC + ROBIN
t = tic;
pace = pace_raw(curproblem,true,true);
pace.fulltime = toc(t);

solns = [solns; pace];

if (mod(j,5) == 0)
    disp(j);
end

end

%% Check solutions
L = problem.L;
N = curproblem.N_VAR;

p_err = zeros(L*length(solns),1);
R_err = zeros(L*length(solns),1);
est.p = zeros(3,1,L*length(solns));
est.R = zeros(3,3,L*length(solns));

figure(1);
for j = 1:length(solns)

problem = problems{j};
soln = solns(j);

idx = ((j-1)*L + 1):j*L;
for l = 1:L
    p_err(idx(l)) = norm(soln.p(:,:,l) - gt.p(:,:,idx(l)));
    R_err(idx(l)) = getAngularError(gt.R(:,:,idx(l)),soln.R(:,:,l));

    % if (p_err(idx(l)) > 10)
    %     % don't plot
    %     soln.gap = 1;
    % end
end

% if (soln.gap > 0.5)
%     % don't plot
%     est.p(:,:,idx) = NaN;
%     est.R(:,:,idx) = NaN;
%     continue
% end

est.p(:,:,idx) = soln.p;
est.R(:,:,idx) = soln.R;

% eigenvalue plot
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

% Plot trajectory!
figure(1);
axis equal
p_est = reshape(soln.p,[3,L,1]);
plot3(p_est(1,:),p_est(2,:),p_est(3,:),'.k', 'MarkerSize',10);
hold on

R_est = soln.R;
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,1,:)),squeeze(R_est(2,1,:)),squeeze(R_est(3,1,:)),'r');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,2,:)),squeeze(R_est(2,2,:)),squeeze(R_est(3,2,:)),'g');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,3,:)),squeeze(R_est(2,3,:)),squeeze(R_est(3,3,:)),'b');

end

%% Plot Ground Truth
plotgt = gt;

figure
p_gt = reshape(plotgt.p,[3,size(plotgt.p,3),1]);
% p_gt = p_gt(:,12*8:15*8)
R_gt = plotgt.R;
plot3(p_gt(1,:),p_gt(2,:),p_gt(3,:),'.k', 'MarkerSize',10);
hold on
axis equal

% R_est = soln.R_est;
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(R_gt(1,1,:)),squeeze(R_gt(2,1,:)),squeeze(R_gt(3,1,:)),'r');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(R_gt(1,2,:)),squeeze(R_gt(2,2,:)),squeeze(R_gt(3,2,:)),'g');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(R_gt(1,3,:)),squeeze(R_gt(2,3,:)),squeeze(R_gt(3,3,:)),'b');

%% Plot Together
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

%% Print error metrics
gt.p = gt.p(:,:,1:length(est.p));
gt.R = gt.R(:,:,1:length(est.R));
% degcm
[est.degcm, est.p_err, est.R_err] = compute_degcm(gt(1:length(est)),est);
% [teaser.degcm, teaser.p_err, teaser.R_err] = compute_degcm(gt,teaser); % should remove 0s--those are where TEASER failed
degcm_10_5 = compute_degcm(gt,est,'degThreshold',10);
% c error
cerr = compute_cerr(solns, est, problem.shapes, problem.shapes(:,:,end));


function cerr = compute_cerr(solns, est, shapes, shape_gt)
    N = size(shapes,2);
    K = size(shapes,3);
    L = length(solns);
    
    B = reshape(shapes, [3*N,K]);

    cerr = zeros(3*L,N);
    for l = 1:L
        for i = 1:3
            b_est = reshape(B*solns(l).c(:,:,i),[3,N]);
            cerr(3*(l-1)+i,:) = vecnorm(b_est - shape_gt);
        end
    end

end