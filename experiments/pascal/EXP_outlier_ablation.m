%% RAL Experiment: Outlier Robustness
% Dataset: pascal+car
% Constants: K, N, L, measurement noise, spiral motion
% Independent variable: #outliers, pruning type
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
%
% Lorenzo Shaikewitz for SPARK Lab

% problem: much too slow with PACE

clc; clear; close all

%% Experiment settings
indepVar = "outlierratio"; % name of independent variable
savename = "pascalcar_" + indepVar;
domain = [0.05:0.025:0.95];
num_repeats = 50;

%% Loop
results = cell(length(domain),1);
for index = 1:length(domain)
iv = domain(index);
resultsIV = struct();
resultsIV.(indepVar) = iv;
resultsIV.R_err_ours = zeros(num_repeats,1);
resultsIV.R_err_ekf = zeros(num_repeats,1);
resultsIV.R_err_pace = zeros(num_repeats,1);
resultsIV.p_err_ours = zeros(num_repeats,1);
resultsIV.p_err_ekf = zeros(num_repeats,1);
resultsIV.p_err_pace = zeros(num_repeats,1);
resultsIV.c_err_ours = zeros(num_repeats,1);
resultsIV.c_err_pace = zeros(num_repeats,1);
resultsIV.gap_ours = zeros(num_repeats,1);
resultsIV.gap_pace = zeros(num_repeats,1);
resultsIV.time_ours = zeros(num_repeats,1);
resultsIV.time_pace = zeros(num_repeats,1);
disp("Starting " + indepVar + "=" + string(iv));
for j = 1:num_repeats
problem = struct();
problem.category = "car";
problem.L = 8;

problem.outlierRatio = iv;
problem.noiseSigmaSqrt = 0.1*0.3; % [m] (0.3 is length scale for car)
problem.velocity_weight_multiplier = 1;
problem.rotrate_kappa_multiplier = 1;
problem.noiseBound = 0.1;
problem.noiseBound_GNC = 0.1; % TODO: set!
problem.noiseBound_GRAPH = 0.11;
problem.processNoise = 0.1;

problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity

problem.accelerationNoiseBoundSqrt = 0;
problem.rotationNoiseBound = 0; % rad

% add shape, measurements, outliers
problem = gen_pascal_tracking(problem);
lambda = 0;
problem.lambda = lambda;

% for GNC
problem.N = problem.N_VAR*problem.L; % How many measurements this problem has (updated by ROBIN)
problem.outliers = []; % outlier indicies
problem.priors = [];
problem.dof = 3;

% Solve!
pace = pace_raw(problem, true, true); % GNC+ROBIN
paceekf = pace_ekf(problem,pace);

% Save solutions: only use last error
% rotation error
R_err_ekf = getAngularError(problem.R_gt(:,:,end), paceekf.R(:,:,end));
R_err_pace = getAngularError(problem.R_gt(:,:,end), pace.R(:,:,end));
% position error
p_err_ekf = norm(problem.p_gt(:,:,end) - paceekf.p(:,:,end));
p_err_pace = norm(problem.p_gt(:,:,end) - pace.p(:,:,end));
% shape error
c_err_pace = norm(problem.c_gt - pace.c(:,:,end));
% time and gap
gap_pace = pace.gaps(end);
time_pace = pace.times(end);

% regen only first time
problem.regen_sdp = (j == 1);
problem = lorenzo_prune(problem);
try
    [inliers, info] = gnc_custom(problem, @solver_for_gnc, 'NoiseBound', problem.noiseBound_GNC,'MaxIterations',100,'FixPriorOutliers',false);
    disp("GNC finished " + string(j))

    soln = info.f_info.soln;
catch
    soln.p_est = ones(3,1,1)*NaN;
    soln.R_est = ones(3,3,1)*NaN;
    soln.c_est = ones(problem.K,1)*NaN;
    soln.gap = NaN;
    soln.solvetime = NaN;
    disp("\nGNC failed " + string(j))
end

R_err_ours = getAngularError(problem.R_gt(:,:,end), soln.R_est(:,:,end));
p_err_ours = norm(problem.p_gt(:,:,end) - soln.p_est(:,:,end));
c_err_ours = norm(problem.c_gt - soln.c_est);
gap_ours = soln.gap;
time_ours = soln.solvetime;

% save
resultsIV.R_err_ours(j) = R_err_ours;
resultsIV.R_err_ekf(j)  = R_err_ekf;
resultsIV.R_err_pace(j) = R_err_pace;
resultsIV.p_err_ours(j) = p_err_ours;
resultsIV.p_err_ekf(j)  = p_err_ekf;
resultsIV.p_err_pace(j) = p_err_pace;
resultsIV.c_err_ours(j) = c_err_ours;
resultsIV.c_err_pace(j) = c_err_pace;
resultsIV.gap_ours(j) = gap_ours;
resultsIV.gap_pace(j) = gap_pace;
resultsIV.time_ours(j) = time_ours;
resultsIV.time_pace(j) = time_pace;

end
results{index} = resultsIV;
end
results = [results{:}];
% save
save("../datasets/results/" + savename + ".mat","results")

%% Display Results
% process into displayable form
% settings.OURS = {'DisplayName', 'OURS','LineWidth',3};
ours_colors = ["#000000", "#002e4c", "#00395f", "#004471", "#005084",...
               "#005b97", "#0067aa","#0072bd", "#1a80c4", "#338eca",...
               "#4d9cd1","#66aad7","#80b9de"];
settings.PACEEKF = {'DisplayName', 'PACE-EKF', 'Color', "#D95319"};
settings.PACERAW = {'DisplayName', 'PACE-RAW', 'Color', "#EDB120"};
figure
tiledlayout(2,2);
set(0,'DefaultLineLineWidth',2)

display_range = 1:length(domain);
L = 10;

% Rotation figure
nexttile
hold on
c=plot([results.(indepVar)],median([results.R_err_pace]),'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.R_err_pace],get(c,'Color'));
% b=plot([results.(indepVar)],median([results.R_err_ekf]),'x-',settings.PACEEKF{:});
% errorshade([results.(indepVar)],[results.R_err_ekf],get(b,'Color'));
res = [results.R_err_ours];
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(1)};
a=plot([results.(indepVar)],median(res),'x-',plotsettings{:});
errorshade([results.(indepVar)],res,get(a,'Color'));

xlabel(indepVar); ylabel("Rotation Error (deg)");
title("Rotation Errors")

% position figure
nexttile
hold on
b=plot([results.(indepVar)],median([results.p_err_ekf])/0.3,'x-',settings.PACEEKF{:});
errorshade([results.(indepVar)],[results.p_err_ekf]/0.3,get(b,'Color'));
c=plot([results.(indepVar)],median([results.p_err_pace])/0.3,'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.p_err_pace]/0.3,get(c,'Color'));
res = [results.p_err_ours];
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(1)};
a=plot([results.(indepVar)],median(res),'x-',plotsettings{:});
errorshade([results.(indepVar)],res,get(a,'Color'));

xlabel(indepVar); ylabel("Position Error (normalized)");
title("Position Errors")

lg = legend('Orientation','horizontal');
lg.Layout.Tile = 'south';

% shape figure
nexttile
hold on
b=plot([results.(indepVar)],median([results.c_err_pace]),'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.c_err_pace],get(b,'Color'));
res = [results.c_err_ours];
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(1)};
a=plot([results.(indepVar)],median(res),'x-',plotsettings{:});
errorshade([results.(indepVar)],res,get(a,'Color'));

xlabel(indepVar); ylabel("Shape Error (normalized)");
title("Shape Errors")

% gap figure
nexttile
hold on
res = abs([results.gap_ours]);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(1)};
a=plot([results.(indepVar)],median(res),'x-',plotsettings{:});
errorshade([results.(indepVar)],res,get(a,'Color'));

xlabel(indepVar); ylabel("Gap");
title("Suboptimality Gaps")

% time figure
% nexttile
% hold on
% res = [results.time_ours];
% plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(1)};
% a=plot([results.(indepVar)],median(res),'x-',plotsettings{:});
% errorshade([results.(indepVar)],res,get(a,'Color'));
% 
% xlabel(indepVar); ylabel("Time (s)");
% title("Solve Time")

%% Change results format
% for i = 1:length(results)
%     results(i).R_err_ekf = results(i).R_err_ekf';
%     results(i).R_err_pace = results(i).R_err_pace';
%     results(i).p_err_ekf = results(i).p_err_ekf';
%     results(i).p_err_pace = results(i).p_err_pace';
%     results(i).c_err_pace = results(i).c_err_pace';
% end