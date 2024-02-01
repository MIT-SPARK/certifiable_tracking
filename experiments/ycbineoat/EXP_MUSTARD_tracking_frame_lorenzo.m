%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through Lorenzo+GNC, frame level
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")

%% Define settings for batch processing
problem.json = "../datasets/ycbineoat/mustard0_metrics.json";
problem.L = 10; % batch size
problem.savefile = "../datasets/ycbineoat/mustard0_metrics.json";

% Set bounds based on problem setting
problem.translationBound = 5.0;
problem.velocityBound = 1.5;
problem.noiseBound_GNC = 0.01;
problem.noiseBound_GRAPH = 0.05;
problem.noiseBound = 0.01;
problem.covar_velocity_base = 0.05^2;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% add shape, measurements, outliers
load("../datasets/ycbineoat/mustard_syn.mat");
problem.shapes = annotatedPoints' / 1000; % 3 x N x K [m]
min_max_dists = robin_min_max_dists(problem.shapes);
[problems, gt, teaser] = json2frameproblem(problem);

%% Solve for each frame
solns = [];

disp("Solving " + string(length(problems)) + " problems...")
last_L = 0;
for j = 1:length(problems)
curproblem = problems{j};
curproblem.regen_sdp = (curproblem.L ~= last_L); % regen only first time
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
try
    [inliers, info] = gnc_custom(curproblem, @solver_for_gnc, 'NoiseBound', curproblem.noiseBound,'MaxIterations',100,'FixPriorOutliers',false);
    disp("GNC finished " + string(j))

    soln = info.f_info.soln;
    ef = eig(soln.raw.Xopt{1});
    if (ef(end-4) > 1e-4)
        disp("**Not convergent**")
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

if (mod(j,5) == 0)
    disp(j);
end

end

%% Check solutions
% TODO: UPDATE
L = problem.L;
N = curproblem.N_VAR;

p_err = zeros(L*N,L)*NaN;
R_err = zeros(L*N,L)*NaN;

est_legit.p = zeros(3,1,length(solns))*NaN;
est_legit.R = zeros(3,3,length(solns))*NaN;
est_ignoringbad.p = zeros(3,1,length(solns))*NaN;
est_ignoringbad.R = zeros(3,3,length(solns))*NaN;
est_bestrun.p = zeros(3,1,length(solns))*NaN;
est_bestrun.R = zeros(3,3,length(solns))*NaN;
gaps = [];

figure(1);
for j = 1:length(solns)

problem = problems{j};
soln = solns(j);
L_cur = problem.L;

% eigenvalue plot
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

% Plot trajectory!
% figure(1);
% axis equal
% p_est = reshape(soln.p_est,[3,L_cur,1]);
% plot3(p_est(1,:),p_est(2,:),p_est(3,:),'.k', 'MarkerSize',10);
% hold on
% 
% R_est = soln.R_est;
% quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,1,:)),squeeze(R_est(2,1,:)),squeeze(R_est(3,1,:)),'r');
% quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,2,:)),squeeze(R_est(2,2,:)),squeeze(R_est(3,2,:)),'g');
% quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,3,:)),squeeze(R_est(2,3,:)),squeeze(R_est(3,3,:)),'b');

idx = problem.startIdx:(problem.startIdx + L_cur);
for l = 1:L_cur
    p_err(idx(l),l) = norm(soln.p_est(:,:,l) - gt.p(:,:,idx(l)));
    R_err(idx(l),l) = getAngularError(gt.R(:,:,idx(l)),soln.R_est(:,:,l));
end


% compute errors
% 1) legit: use latest estimate no matter what
est_legit.p(:,:,idx(end)) = soln.p_est(:,:,end);
est_legit.R(:,:,idx(end)) = soln.R_est(:,:,end);
% 2) ignoringbad: throw away bad
if (soln.gap_nov < 0.3)
    est_ignoringbad.p(:,:,idx(end)) = soln.p_est(:,:,end);
    est_ignoringbad.R(:,:,idx(end)) = soln.R_est(:,:,end);
end
% 3) bestrun: use estimate that has best run (by gap)
gaps = [gaps; abs(soln.gap_nov)];
for l = 1:L_cur-1
    % for each pose estimate, excluding the current one
    testrange = length(gaps) - (1:(L_cur - l));
    newgood = true;
    for l_test = testrange
        if (l_test < 1)
            break
        end
        % for each other gap we can compare with
        if (abs(soln.gap_nov) > gaps(l_test))
            % fails if any of them have better gaps
            newgood = false;
            break
        end
    end
    if (newgood)
        % update!
        est_bestrun.p(:,:,idx(l)) = soln.p_est(:,:,l);
        est_bestrun.R(:,:,idx(l)) = soln.R_est(:,:,l);
    end
end
est_bestrun.p(:,:,idx(end)) = soln.p_est(:,:,end);
est_bestrun.R(:,:,idx(end)) = soln.R_est(:,:,end);


end

est_ignoringbad.p = est_ignoringbad.p(:,:,1:end-1);
est_legit.p = est_legit.p(:,:,1:end-1);
est_bestrun.p = est_bestrun.p(:,:,1:end-1);
est_ignoringbad.R = est_ignoringbad.R(:,:,1:end-1);
est_legit.R = est_legit.R(:,:,1:end-1);
est_bestrun.R = est_bestrun.R(:,:,1:end-1);

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

%% Save everything we need
% soln
save("ycb_mustard.mat","solns");

%% ADD and ADD-S scores
pcfile_gt = "~/Downloads/models/006_mustard_bottle/google_16k/nontextured.ply";
pcfile_est = "~/Downloads/models/006_mustard_bottle/google_16k/nontextured.ply";

% Compute scores!
[add_ours_legit, adds_ours_legit] = compute_scores(gt, est_legit, pcfile_gt, pcfile_est);
[add_ours_ignoringbad, adds_ours_ignoringbad] = compute_scores(gt, est_ignoringbad, pcfile_gt, pcfile_est);
[add_ours_bestrun, adds_ours_bestrun] = compute_scores(gt, est_bestrun, pcfile_gt, pcfile_est);

[add_teaser, adds_teaser] = compute_scores(gt, teaser, pcfile_gt, pcfile_est);

scores.add_ours_legit = add_ours_legit;
scores.adds_ours_legit = adds_ours_legit;
scores.add_ours_ignoringbad = add_ours_ignoringbad;
scores.adds_ours_ignoringbad = adds_ours_ignoringbad;
scores.add_ours_bestrun = add_ours_bestrun;
scores.adds_ours_bestrun = adds_ours_bestrun;
scores.add_teaser = add_teaser;
scores.adds_teaser = adds_teaser;

%% Save everything we need
% soln
save("ycb_mustard.mat","solns","scores");
