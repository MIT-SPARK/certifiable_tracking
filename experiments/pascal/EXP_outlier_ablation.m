%% RAL Experiment: Outlier Robustness
% Dataset: pascal+car
% Constants: K, N, L, measurement noise, spiral motion
% Independent variable: #outliers, pruning type
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment settings
indepVar = "outlierratio"; % name of independent variable
savename = "pascalaeroplane_mle_" + indepVar;
lengthScale = 0.2; % smallest dimension
domain = [0.05:0.05:0.75];
num_repeats = 500; % 3: 500

%% Loop
results = cell(length(domain),1);
all_problems = {};
all_pruned_problems = {};
for index = 1:length(domain)
iv = domain(index);

problems = {};
problems_pruned = {};
for j = 1:num_repeats
problem = struct();
problem.category = "aeroplane";
problem.L = 8;

problem.outlierRatio = iv;
problem.outlierVariance = 1.0*lengthScale; % 3: 1.0 (no lengthScale)

problem.noiseSigmaSqrt = 0.01*lengthScale; % [m]
problem.accelerationNoiseBoundSqrt = 0.05*0.2;
problem.rotationKappa = 1/(0.05*0.2)^2*(1/2);

problem.covar_measure_base = problem.noiseSigmaSqrt^2;
problem.covar_velocity_base = problem.accelerationNoiseBoundSqrt^2;
problem.kappa_rotrate_base = problem.rotationKappa;

problem.noiseBound = 0.05*lengthScale;
problem.noiseBound_GNC = 0.01*lengthScale;
problem.noiseBound_GNC_residuals = 1;
problem.noiseBound_GRAPH = 0.01*lengthScale;

problem.translationBound = 10.0;
problem.velocityBound = 5.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity
problem.usecBound = false;

% problem.accelerationNoiseBoundSqrt = 0;
% problem.rotationNoiseBound = 0; % rad

% add shape, measurements, outliers
problem = gen_pascal_tracking(problem);
lambda = 0;
problem.lambda = lambda;
problem.useusecBound = false;

% for GNC
problem.N = problem.N_VAR*problem.L; % How many measurements this problem has (updated by ROBIN)
problem.outliers = []; % outlier indicies
problem.priors = [];
problem.dof = 3;

% pre-prune problems
problem_milp = lorenzo_prune(problem);

% save
problems{j} = problem;
problems_pruned{j} = problem_milp;
end
% save
all_problems{index} = problems;
all_pruned_problems{index} = problems_pruned;
end

%% can be parfor TODO
parfor index = 1:length(domain)
iv = domain(index);
problems = all_problems{index};
problems_pruned = all_pruned_problems{index};

resultsIV = struct();
resultsIV.(indepVar) = iv;
resultsIV.R_err_ours = zeros(num_repeats,1);
resultsIV.R_err_gnc = zeros(num_repeats,1);
resultsIV.R_err_milp = zeros(num_repeats,1);
resultsIV.p_err_ours = zeros(num_repeats,1);
resultsIV.p_err_gnc = zeros(num_repeats,1);
resultsIV.p_err_milp = zeros(num_repeats,1);
resultsIV.c_err_ours = zeros(num_repeats,1);
resultsIV.c_err_gnc = zeros(num_repeats,1);
resultsIV.c_err_milp = zeros(num_repeats,1);
resultsIV.iter_ours = zeros(num_repeats,1);
resultsIV.iter_gnc = zeros(num_repeats,1);
resultsIV.iter_milp = zeros(num_repeats,1);
resultsIV.time_ours = zeros(num_repeats,1);
resultsIV.time_gnc = zeros(num_repeats,1);
resultsIV.time_milp = zeros(num_repeats,1);
disp("Starting " + indepVar + "=" + string(iv));

for j = 1:num_repeats
problem = problems{j};
problem_milp = problems_pruned{j};

% regen only first time
problem.regen_sdp = false;
problem_milp.regen_sdp = (j==1);
% try
    t = tic;
    [inliers, info] = gnc2(problem_milp, @solver_for_gnc,'barc2',problem.noiseBound_GNC, 'ContinuationFactor', 1.6);
    soln_ours = info.f_info.soln;
    soln_ours.iters = info.Iterations;
    soln_ours.time = toc(t) + problem_milp.milptime;
% catch
%     soln_ours.p_est = ones(3,1,1)*NaN;
%     soln_ours.R_est = ones(3,3,1)*NaN;
%     soln_ours.c_est = ones(problem.K,1)*NaN;
%     soln_ours.iters = NaN;
%     soln_ours.time = NaN;
% end

% GNC only
% try
    t = tic;
    [inliers, info] = gnc2(problem, @solver_for_gnc,'barc2',problem.noiseBound_GNC, 'ContinuationFactor', 1.6);%,'Debug',true);
    soln_gnc = info.f_info.soln;
    soln_gnc.iters = info.Iterations;
    soln_gnc.time = toc(t);
% catch
%     soln_gnc.p_est = ones(3,1,1)*NaN;
%     soln_gnc.R_est = ones(3,3,1)*NaN;
%     soln_gnc.c_est = ones(problem.K,1)*NaN;
%     soln_gnc.iters = NaN;
%     soln_gnc.time = NaN;
% end

% MILP only
% try
    t = tic;
    soln = solve_weighted_tracking(problem_milp);
    soln_milp = soln;
    soln_milp.iters = 1;
    soln_milp.time = toc(t) + problem_milp.milptime;
% catch
%     soln_milp.p_est = ones(3,1,1)*NaN;
%     soln_milp.R_est = ones(3,3,1)*NaN;
%     soln_milp.c_est = ones(problem.K,1)*NaN;
%     soln_milp.iters = NaN;
%     soln_milp.time = NaN;
% end
disp("Finished iv="+iv+" ("+j+")");

R_err_ours = getAngularError(problem.R_gt(:,:,end), soln_ours.R_est(:,:,end));
p_err_ours = norm(problem.p_gt(:,:,end) - soln_ours.p_est(:,:,end));
c_err_ours = norm(problem.c_gt - soln_ours.c_est);
iters_ours = soln_ours.iters;
time_ours = soln_ours.time;

R_err_gnc = getAngularError(problem.R_gt(:,:,end), soln_gnc.R_est(:,:,end));
p_err_gnc = norm(problem.p_gt(:,:,end) - soln_gnc.p_est(:,:,end));
c_err_gnc = norm(problem.c_gt - soln_gnc.c_est);
iters_gnc = soln_gnc.iters;
time_gnc = soln_gnc.time;

R_err_milp = getAngularError(problem.R_gt(:,:,end), soln_milp.R_est(:,:,end));
p_err_milp = norm(problem.p_gt(:,:,end) - soln_milp.p_est(:,:,end));
c_err_milp = norm(problem.c_gt - soln_milp.c_est);
iters_milp = soln_milp.iters;
time_milp = soln_milp.time;

% save
resultsIV.R_err_ours(j) = R_err_ours;
resultsIV.R_err_gnc(j)  = R_err_gnc;
resultsIV.R_err_milp(j) = R_err_milp;
resultsIV.p_err_ours(j) = p_err_ours;
resultsIV.p_err_gnc(j)  = p_err_gnc;
resultsIV.p_err_milp(j) = p_err_milp;
resultsIV.c_err_ours(j) = c_err_ours;
resultsIV.c_err_gnc(j) = c_err_gnc;
resultsIV.c_err_milp(j) = c_err_milp;
resultsIV.c_err_ours(j) = c_err_ours;
resultsIV.c_err_gnc(j) = c_err_gnc;
resultsIV.c_err_milp(j) = c_err_milp;
resultsIV.iter_ours(j) = iters_ours;
resultsIV.iter_gnc(j) = iters_gnc;
resultsIV.iter_milp(j) = iters_milp;
resultsIV.time_ours(j) = time_ours;
resultsIV.time_gnc(j) = time_gnc;
resultsIV.time_milp(j) = time_milp;
end

results{index} = resultsIV;
end
results = [results{:}];
% save
save("../datasets/results/" + savename + ".mat","results")

%%
load("../datasets/results/" + savename + ".mat","results")

%% Display Results
% process into displayable form
settings.OURS = {'x-','DisplayName', 'OURS', 'Color', "#005b97",'LineWidth',2.5};
settings.PACEEKF = {'DisplayName', 'PACE-EKF', 'Color', "#D95319"};
settings.PACERAW = {'DisplayName', 'PACE-RAW', 'Color', "#EDB120"};
settings.GNC  = {'x:','DisplayName', 'OURS-GNC', 'Color', "#9A6324",'LineWidth',2}; % TODO: change colors
settings.MILP = {'x-.','DisplayName', 'OURS-MILP', 'Color', "#36B649",'LineWidth',2};
figure
tiledlayout(1,5);
domain = [results.(indepVar)];

% position figure
nexttile
errorshade([results.(indepVar)],[results.p_err_gnc]/lengthScale,hex2rgb(settings.GNC{5})); hold on;
errorshade([results.(indepVar)],[results.p_err_ours]/lengthScale,hex2rgb(settings.OURS{5}));
errorshade([results.(indepVar)],[results.p_err_milp]/lengthScale,hex2rgb(settings.MILP{5}));
b=loglog([results.(indepVar)],median([results.p_err_milp],"omitmissing")/lengthScale,settings.MILP{:}); hold on;
c=plot([results.(indepVar)],median([results.p_err_gnc],"omitmissing")/lengthScale,settings.GNC{:}); hold on;
a=plot([results.(indepVar)],median([results.p_err_ours],"omitmissing")/lengthScale,settings.OURS{:});
yscale log; %xscale log
xlabel(indepVar); ylabel("Position Error (normalized)");
xlim([domain(1), domain(end)])
ylim tight
title("Position Errors")

% Rotation figure
nexttile
errorshade([results.(indepVar)],[results.R_err_gnc],hex2rgb(settings.GNC{5})); hold on;
errorshade([results.(indepVar)],[results.R_err_ours],hex2rgb(settings.OURS{5}));
errorshade([results.(indepVar)],[results.R_err_milp],hex2rgb(settings.MILP{5}));
c=loglog([results.(indepVar)],median([results.R_err_gnc],"omitmissing"),settings.GNC{:}); hold on;
b=plot([results.(indepVar)],median([results.R_err_milp],"omitmissing"),settings.MILP{:});
a=plot([results.(indepVar)],median([results.R_err_ours],"omitmissing"),settings.OURS{:});

yscale log;% xscale log
xlabel(indepVar); ylabel("Rotation Error (deg)");
xlim([domain(1), domain(end)])
ylim tight
title("Rotation Errors")

lg = legend('Orientation','horizontal');
lg.Layout.Tile = 'south';

% shape figure
nexttile
errorshade([results.(indepVar)],[results.c_err_gnc],hex2rgb(settings.GNC{5})); hold on;
errorshade([results.(indepVar)],[results.c_err_milp],hex2rgb(settings.MILP{5}));
errorshade([results.(indepVar)],[results.c_err_ours],hex2rgb(settings.OURS{5}));
b=loglog([results.(indepVar)],median([results.c_err_gnc],"omitmissing"),settings.GNC{:});
b=loglog([results.(indepVar)],median([results.c_err_milp],"omitmissing"),settings.MILP{:});
a=plot([results.(indepVar)],median([results.c_err_ours],"omitmissing"),settings.OURS{:});
yscale log; %xscale log
xlabel(indepVar); ylabel("Shape Error");
xlim([domain(1), domain(end)])
ylim tight
title("Shape Errors")

% Iterations figure
nexttile
b=plot([results.(indepVar)],abs(median([results.iter_gnc])),settings.GNC{:}); hold on;
a=plot([results.(indepVar)],abs(median([results.iter_ours])),settings.OURS{:});
% b=plot([results.(indepVar)],abs(median([results.iter_milp])),settings.MILP{:});
errorshade([results.(indepVar)],abs([results.iter_ours]),hex2rgb(settings.OURS{5}));
errorshade([results.(indepVar)],abs([results.iter_gnc]),hex2rgb(settings.GNC{5}));
% errorshade([results.(indepVar)],abs([results.iter_milp]),hex2rgb(settings.MILP{5}));
% yscale log; xscale log
xlabel(indepVar); ylabel("Iterations");
xlim([domain(1), domain(end)])
ylim tight
title("GNC Iterations")

% time figure
nexttile
hold on
c=plot([results.(indepVar)],median([results.time_milp],"omitmissing"),settings.MILP{:});
b=plot([results.(indepVar)],median([results.time_gnc],"omitmissing"),settings.GNC{:});
a=plot([results.(indepVar)],median([results.time_ours],"omitmissing"),settings.OURS{:});
errorshade([results.(indepVar)],[results.time_ours],hex2rgb(settings.OURS{5}));
errorshade([results.(indepVar)],[results.time_gnc],hex2rgb(settings.GNC{5}));
errorshade([results.(indepVar)],[results.time_milp],hex2rgb(settings.MILP{5}));
xlabel(indepVar); ylabel("Time (s)");
xlim([domain(1), domain(end)])
ylim tight
title("Solve Time")
yscale log;

%% Change results format
% for i = 1:length(results)
%     results(i).R_err_ekf = results(i).R_err_ekf';
%     results(i).R_err_pace = results(i).R_err_pace';
%     results(i).p_err_ekf = results(i).p_err_ekf';
%     results(i).p_err_pace = results(i).p_err_pace';
%     results(i).c_err_pace = results(i).c_err_pace';
% end