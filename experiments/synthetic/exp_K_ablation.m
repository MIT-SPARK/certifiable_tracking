%% RSS Experiment: How does number of shapes affect performance?
% Dataset: synthetic
% Constants: K, N, noiseSigma, NO nonlinearities in gt
% Independent variable: K
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment settings
indepVar = "L"; % name of independent variable
savename = "syn_" + indepVar;
domain = [1:2:20, 30:10:100];
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
for j = 1:num_repeats

% Generate random tracking problem
problem.N_VAR = 11; % nr of keypoints
problem.K = iv; % nr of shapes

problem.L = 10; % nr of keyframes in horizon
L = problem.L;

problem.outlierRatio = 0.0;
problem.noiseSigmaSqrt = 0.1; % [m]
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

% Rotation figure
figure
set(0,'DefaultLineLineWidth',2)
plot([results.(indepVar)],mean([results.R_err_ours]),'x-','DisplayName','OURS');
hold on
plot([results.(indepVar)],mean([results.R_err_ukf]),'x-','DisplayName','PACE-UKF');
plot([results.(indepVar)],mean([results.R_err_pace]),'x-','DisplayName','PACE-RAW');
legend
xlabel("L"); ylabel("Rotation Error (deg)");
title("Rotation Errors")

% position figure
figure
plot([results.(indepVar)],mean([results.p_err_ours]),'x-','DisplayName','OURS');
hold on
plot([results.(indepVar)],mean([results.p_err_ukf]),'x-','DisplayName','PACE-UKF');
plot([results.(indepVar)],mean([results.p_err_pace]),'x-','DisplayName','PACE-RAW');
legend
xlabel("L"); ylabel("Position Error (m)");
title("Position Errors")

% gap figure
figure
semilogy([results.(indepVar)],mean([results.gap_ours]),'x-');
xlabel("L"); ylabel("Gap");
title("Suboptimality Gaps")

% time figure
figure
plot([results.(indepVar)],mean([results.time_ours]),'x-');
xlabel("L"); ylabel("Time (s)");
title("Solve Time")