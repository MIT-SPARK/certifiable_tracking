%% IROS Experiment: Horizon-Based Performance on YCBInEOAT
% Dataset: YCBInEOAT
% Constants: K, N, L (max)
% Independent variable: choice of video
% Dependent variables: runtime, accuracy (p, R), ADD, ADD-S
%
% Lorenzo Shaikewitz for SPARK Lab

% TODO:
% 1) json2frameproblem needs to be updated to load the shape
% 2) results processing isn't smart, doesn't support multiple domains

clc; clear; close all

%% Experiment Settings
indepVar = "video"; % name of independent variable
savename = "ycbineoat_" + indepVar;
domain = ["cracker_box_reorient", "cracker_box_yalehand0", ...
          "sugar_box1", "sugar_box_yalehand0", ...
          "mustard0", "mustard_easy_00_02", ...
          "bleach0", "bleach_hard_00_03_chaitanya", ...
          "tomato_soup_can_yalehand0"];
domain = domain(5); % just do one for now
% directories
jsondir = "../datasets/ycbineoat/";

%% Loop
solns = cell(length(domain),1);
for iv = domain
    % Define json, generate problem
    problem = struct();
    problem.json = jsondir + iv + ".json";
    problem.L = 10; % batch size
    problem.savefile = jsondir + iv + "_ours.json";
    problem.velprior = "body";

    % Parameters
    problem.translationBound = 5.0;
    problem.velocityBound = 1.5;
    problem.noiseBound_GNC = 0.01; % TODO: tune?
    problem.noiseBound_GRAPH = 0.05; % TODO: tune?
    problem.noiseBound = 0.01; % TODO: tune?
    problem.covar_velocity_base = 0.05^2; % TODO: tune?
    skip = 1; % TODO: experiment with different skip?

    % Add shape, split into batches
    [problems, gt, teaser, shapes] = json2frameproblem(problem, skip);
    disp("Computing max/min distances...")
    min_max_dists = robin_min_max_dists(problems{1}.shapes);
    disp("Finished computing max/min distances")

    % For testing: limit to small set of problems
    numProblemsToSolve = length(problems);

    % define solution variable
    solnsIV = cell(numProblemsToSolve,1);

    % now we can start!
    disp("Starting " + string(iv));
    
    % we need to prune using serial processing
    % disp("Pruning " + string(numProblemsToSolve) + " problems...")
    for batch = 1:numProblemsToSolve
        curproblem = problems{batch};
        % data for GNC/pruning
        curproblem.type = "tracking";
        curproblem.N = curproblem.N_VAR*curproblem.L; % How many measurements this problem has (updated by ROBIN)
        curproblem.outliers = []; % outlier indicies
        curproblem.priors = [];
        curproblem.dof = 0;
        % run pruning!
        curproblem = lorenzo_prune(curproblem, min_max_dists);
        % save pruned problem
        problems{batch} = curproblem;
        disp("Finished " + string(batch));
    end

    disp("Solving " + string(numProblemsToSolve) + " problems...")

    % L should change for the first problem.L - 2 problems
    parfor batch = 1:min(problem.L-2, numProblemsToSolve)
        curproblem = problems{batch};
        curproblem.sdp_filename = "sdpdata" + curproblem.L;
        curproblem.regen_sdp = true;

        soln = solveBatch(curproblem);
        solnsIV{batch} = soln;

        % report bad results
        if (soln.gap_nov > 1e-2)
            fprintf("Batch %d failed: %.4e\n",batch,soln.gap_nov)
        elseif (isnan(soln.gap))
            fprintf("Batch %d failed: NaN\n",batch)
        end
    end
    
    % Now that L is fixed, run through the remainder of the problems
    parfor batch = min(problem.L-2, numProblemsToSolve):numProblemsToSolve
        curproblem = problems{batch};
        curproblem.sdp_filename = "sdpdata" + curproblem.L;
        curproblem.regen_sdp = false;

        soln = solveBatch(curproblem);
        solnsIV{batch} = soln;

        % report bad results
        if (soln.gap_nov > 1e-2)
            fprintf("Batch %d failed: %.4e\n",batch,soln.gap_nov)
        elseif (isnan(soln.gap))
            fprintf("Batch %d failed: NaN\n",batch)
        end
    end

    % save in solns
    solns{domain == iv} = [solnsIV{:}];
end

%% Check solutions
solns = solnsIV;
L = problem.L;
N = problems{1}.N_VAR;

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
soln = solns{j};
L_cur = problem.L;

% eigenvalue plot
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

% Plot trajectory!
% figure(1);
% axis equal
% p_est = reshape(soln.p,[3,L_cur,1]);
% plot3(p_est(1,:),p_est(2,:),p_est(3,:),'.k', 'MarkerSize',10);
% hold on
% 
% R_est = soln.R;
% quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,1,:)),squeeze(R_est(2,1,:)),squeeze(R_est(3,1,:)),'r');
% quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,2,:)),squeeze(R_est(2,2,:)),squeeze(R_est(3,2,:)),'g');
% quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,3,:)),squeeze(R_est(2,3,:)),squeeze(R_est(3,3,:)),'b');

idx = problem.startIdx:(problem.startIdx + L_cur);
for l = 1:L_cur
    p_err(idx(l),l) = norm(soln.p(:,:,l) - gt.p(:,:,idx(l)));
    R_err(idx(l),l) = getAngularError(gt.R(:,:,idx(l)),soln.R(:,:,l));
end


% compute errors
% 1) legit: use latest estimate no matter what
est_legit.p(:,:,idx(end)) = soln.p(:,:,end);
est_legit.R(:,:,idx(end)) = soln.R(:,:,end);
% 2) ignoringbad: throw away bad
if (soln.gap_nov < 0.3)
    est_ignoringbad.p(:,:,idx(end)) = soln.p(:,:,end);
    est_ignoringbad.R(:,:,idx(end)) = soln.R(:,:,end);
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
        est_bestrun.p(:,:,idx(l)) = soln.p(:,:,l);
        est_bestrun.R(:,:,idx(l)) = soln.R(:,:,l);
    end
end
est_bestrun.p(:,:,idx(end)) = soln.p(:,:,end);
est_bestrun.R(:,:,idx(end)) = soln.R(:,:,end);


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

R_est = soln.R;
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

%% Helper function: solve each batch
function soln = solveBatch(problem)
    % data for GNC
    % problem.type = "tracking";
    % problem.N = problem.N_VAR*problem.L; % How many measurements this problem has (updated by ROBIN)
    % problem.outliers = []; % outlier indicies
    % problem.priors = [];
    % problem.dof = 0;
    % 
    % problem = lorenzo_prune(problem, min_max_dists);
    
    % preprocess inliers (use only if pruning off)
    % if isfield(problem,'prioroutliers')
    %     problem.prioroutliers = sort(problem.prioroutliers);
    %     problem.N = problem.N - length(problem.prioroutliers);
    % end
    
    % run GNC
    soln = struct();
    try
        [~, info] = gnc_custom(problem, @solver_for_gnc, 'NoiseBound', problem.noiseBound_GNC,'MaxIterations',100,'FixPriorOutliers',false);
    
        % report compact form
        soln.p = info.f_info.soln.p_est;
        soln.R = info.f_info.soln.R_est;
        soln.c = info.f_info.soln.c_est;
        soln.gap = info.f_info.soln.gap;
        soln.gap_nov = info.f_info.soln.gap_nov;
        soln.solvetime = info.f_info.soln.solvetime;
        soln.iterations = info.Iterations;
    catch
        % report NaNs
        soln.p = ones(3,1,problem.L)*NaN;
        soln.R = ones(3,3,problem.L)*NaN;
        soln.c = ones(problem.K,1)*NaN;
        soln.gap = NaN;
        soln.gap_nov = NaN;
        soln.solvetime = NaN;
        soln.iterations = NaN;
    end
end
