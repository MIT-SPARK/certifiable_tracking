%% RSS Experiment: How does noise affect performance?
% Dataset: synthetic
% Constants: K, N, L, NO nonlinearities in gt
% Independent variable: noiseSigma
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment settings
indepVar = "noiseSigmaSqrt"; % name of independent variable
savename = "syn_" + indepVar;
domain = [0.01,0.02,0.05,0.1];
num_repeats = 50;
% SET INDEPENDENT VARIABLE, DEPENDENT VARS CORRECTLY IN LOOP

%% Loop
for iv = domain
resultsIV.(indepVar) = iv;
resultsIV.R_err_ours = zeros(num_repeats,1);
resultsIV.R_err_ukf = zeros(num_repeats,1);
resultsIV.R_err_pace = zeros(num_repeats,1);
resultsIV.p_err_ours = zeros(num_repeats,1);
resultsIV.p_err_ukf = zeros(num_repeats,1);
resultsIV.p_err_pace = zeros(num_repeats,1);
resultsIV.c_err_ours = zeros(num_repeats,1);
resultsIV.gap_ours = zeros(num_repeats,1);
resultsIV.time_ours = zeros(num_repeats,1);
disp("Starting " + indepVar + "=" + string(iv));
for j = 1:num_repeats

% Generate random tracking problem
problem.N_VAR = 11; % nr of keypoints
problem.K = 3; % nr of shapes

problem.L = 10; % nr of keyframes in horizon
L = problem.L;

problem.outlierRatio = 0.0;
problem.noiseSigmaSqrt = iv; % [m]
problem.noiseBoundSqrt = 3*iv;
problem.intraRadius = 0.2;
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

problem.accelerationNoiseBoundSqrt = 0;%0.01;
problem.rotationNoiseBound = 0;%pi/32; % rad

% regen if pbound, vbound, N, L, K change.
problem.regen_sdp = (j == 1); % regen only first time
problem.dR_gt = repmat(eye(3),[1,1,L-1]);

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 0.0;
problem.lambda = lambda;

% Solve!
soln = solve_weighted_tracking(problem);

soln_pace = pace_py_UKF(problem);

% Save solutions
% projected errors
R_err_ours = zeros(L,1);
R_err_ukf = zeros(L,1);
R_err_pace = zeros(L,1);
for l = 1:L
    R_err_ours(l) = getAngularError(problem.R_gt(:,:,l), soln.R_est(:,:,l));
    R_err_ukf(l) = getAngularError(problem.R_gt(:,:,l), soln_pace.R_smoothed(:,:,l));
    R_err_pace(l) = getAngularError(problem.R_gt(:,:,l), soln_pace.R_raw(:,:,l));
end

% shape error
c_err = norm(problem.c_gt - soln.c_est);

% save
resultsIV.R_err_ours(j) = norm(R_err_ours)/L;
resultsIV.R_err_ukf(j)  = norm(R_err_ukf)/L;
resultsIV.R_err_pace(j) = norm(R_err_pace)/L;
resultsIV.p_err_ours(j) = norm(problem.p_gt - soln.p_est,'fro')/L;
resultsIV.p_err_ukf(j)  = norm(problem.p_gt - soln_pace.p_smoothed,'fro')/L;
resultsIV.p_err_pace(j) = norm(problem.p_gt - soln_pace.p_raw,'fro')/L;
resultsIV.c_err_ours(j) = c_err;
resultsIV.gap_ours(j) = soln.gap;
resultsIV.time_ours(j) = soln.solvetime;
clear problem;
end
results(domain == iv) = resultsIV;
end
% save
save("../datasets/results/" + savename + ".mat","results")

%% Display Results
% process into displayable form
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};

% Rotation figure
figure
set(0,'DefaultLineLineWidth',2)
a=plot([results.(indepVar)],median([results.R_err_ours]),'x-','DisplayName','OURS');
hold on
b=plot([results.(indepVar)],median([results.R_err_ukf]),'x-','DisplayName','PACE-UKF');
c=plot([results.(indepVar)],median([results.R_err_pace]),'x-','DisplayName','PACE-RAW');

errorshade([results.(indepVar)],[results.R_err_ours],get(a,'Color'));
errorshade([results.(indepVar)],[results.R_err_ukf],get(b,'Color'));
errorshade([results.(indepVar)],[results.R_err_pace],get(c,'Color'));
legend
xlabel(indepVar); ylabel("Rotation Error (deg)");
title("Rotation Errors")

% position figure
figure
plot([results.(indepVar)],mean([results.p_err_ours]),'x-','DisplayName','OURS');
hold on
plot([results.(indepVar)],mean([results.p_err_ukf]),'x-','DisplayName','PACE-UKF');
plot([results.(indepVar)],mean([results.p_err_pace]),'x-','DisplayName','PACE-RAW');
legend
xlabel(indepVar); ylabel("Position Error (m)");
title("Position Errors")

% gap figure
figure
semilogy([results.(indepVar)],abs(mean([results.gap_ours])),'x-');
xlabel(indepVar); ylabel("Gap");
title("Suboptimality Gaps")

% time figure
figure
plot([results.(indepVar)],mean([results.time_ours]),'x-');
xlabel(indepVar); ylabel("Time (s)");
title("Solve Time")

%% Plotting Helper
function errorshade(x,y,color)

curve1 = prctile(y,75);
curve2 = prctile(y,25);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
obj = fill(x2, inBetween,color,'FaceAlpha',0.2,'EdgeColor','none');
obj.Annotation.LegendInformation.IconDisplayStyle = "off";

% TF = isoutlier(y);
% x_rep = repmat(x,[size(y,1),1]);
% plot(x_rep(TF),y(TF),"x",'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i},'MarkerSize',10);
end