%% RAL Experiment: How does the shape library affect performance?
% Dataset: synthetic
% Constants: N, L, measurement noise, spiral trajectory
% Independent variable: K
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment settings
indepVar = "K";
savename = "synthetic2_" + indepVar;
domain = [1, 5, 10, 25, 50, 100, 250, 500, 1000, 2000];
num_repeats = 100;

%% Loop
results = cell(length(domain),1);
parfor index = 1:length(domain)
iv = domain(index)
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
problem.N_VAR = 10; % nr of keypoints
problem.K = iv; % nr of shapes
problem.L = 10; % nr of keyframes in horizon
L = problem.L;

problem.outlierRatio = 0.0;
problem.noiseSigmaSqrt = 0.1; % [m]
problem.velocity_weight_multiplier = 1;
problem.rotrate_kappa_multiplier = 1;

problem.noiseBound = 0.3; %chi2inv(0.95,3*problem.N_VAR*problem.L)*problem.noiseSigmaSqrt^2;
problem.processNoise = 0.1;

problem.intraRadius = 0.2; 
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity

problem.accelerationNoiseBoundSqrt = 0;
problem.rotationNoiseBound = 0; % rad

% regen only first time
problem.regen_sdp = (j == 1);

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 250; %0.5*iv;
problem.lambda = lambda;

% Solve!
soln = solve_weighted_tracking(problem);
pace = pace_raw(problem);
paceekf = pace_ekf(problem,pace);

% Save solutions: only use last error
% rotation error
R_err_ours = getAngularError(problem.R_gt(:,:,L), soln.R_est(:,:,L));
R_err_ekf = getAngularError(problem.R_gt(:,:,L), paceekf.R(:,:,L));
R_err_pace = getAngularError(problem.R_gt(:,:,L), pace.R(:,:,L));
% shape error
c_err_ours = norm(problem.c_gt - soln.c_est);
c_err_pace = norm(problem.c_gt - pace.c(:,:,L));
% time and gap
gap_pace = pace.gaps(end);
time_pace = pace.times(end);

% save
resultsIV.R_err_ours(j) = R_err_ours;
resultsIV.R_err_ekf(j)  = R_err_ekf;
resultsIV.R_err_pace(j) = R_err_pace;
resultsIV.p_err_ours(j) = norm(problem.p_gt(:,:,L) - soln.p_est(:,:,L));
resultsIV.p_err_ekf(j)  = norm(problem.p_gt(:,:,L) - paceekf.p(:,:,L));
resultsIV.p_err_pace(j) = norm(problem.p_gt(:,:,L) - pace.p(:,:,L));
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
settings.OURS = {'DisplayName', 'OURS', 'Color', "#0072BD",'LineWidth',3};
settings.PACEEKF = {'DisplayName', 'PACE-EKF', 'Color', "#D95319"};
settings.PACERAW = {'DisplayName', 'PACE-RAW', 'Color', "#EDB120"};
figure
tiledlayout(2,2);
set(0,'DefaultLineLineWidth',2)

% Rotation figure
nexttile
hold on
c=plot([results.(indepVar)],median([results.R_err_pace]),'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.R_err_pace],get(c,'Color'));
% b=plot([results.(indepVar)],median([results.R_err_ekf]),'x-',settings.PACEEKF{:});
% errorshade([results.(indepVar)],[results.R_err_ekf],get(b,'Color'));
a=plot([results.(indepVar)],median([results.R_err_ours]),'x-',settings.OURS{:});
errorshade([results.(indepVar)],[results.R_err_ours],get(a,'Color'));

xlabel(indepVar); ylabel("Rotation Error (deg)");
title("Rotation Errors")

% position figure
nexttile
hold on
b=plot([results.(indepVar)],median([results.p_err_ekf]),'x-',settings.PACEEKF{:});
errorshade([results.(indepVar)],[results.p_err_ekf],get(b,'Color'));
c=plot([results.(indepVar)],median([results.p_err_pace]),'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.p_err_pace],get(c,'Color'));
a=plot([results.(indepVar)],median([results.p_err_ours]),'x-',settings.OURS{:});
errorshade([results.(indepVar)],[results.p_err_ours],get(a,'Color'));
xlabel(indepVar); ylabel("Position Error (normalized)");
title("Position Errors")

lg = legend('Orientation','horizontal');
lg.Layout.Tile = 'south';

% shape figure
nexttile
hold on
b=plot([results.(indepVar)],median([results.c_err_pace]),'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.c_err_pace],get(b,'Color'));
a=plot([results.(indepVar)],median([results.c_err_ours]),'x-',settings.OURS{:});
errorshade([results.(indepVar)],[results.c_err_ours],get(a,'Color'));
xlabel(indepVar); ylabel("Shape Error (normalized)");
title("Shape Errors")

% gap figure
nexttile
hold on
a=semilogy([results.(indepVar)],abs(median([results.gap_ours])),'x-',settings.OURS{:});
errorshade([results.(indepVar)],abs([results.gap_ours]),get(a,'Color'));
xlabel(indepVar); ylabel("Gap");
title("Suboptimality Gaps")

% time figure
% nexttile
% hold on
% a=plot([results.(indepVar)],median([results.time_ours]),'x-',settings.OURS{:});
% errorshade([results.(indepVar)],[results.time_ours],get(a,'Color'));
% xlabel(indepVar); ylabel("Time (s)");
% title("Solve Time")