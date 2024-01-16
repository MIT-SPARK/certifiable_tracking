%% RSS Experiment: What Time Step Should We Use?
% Dataset: synthetic
% Constants: K, N
% Independent variable: L
% Dependent variables: runtime, duality gap, accuracy
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment settings
Ls = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30];
num_repeats = 50;

%% Loop
for L = Ls
resultsL.L = L;
resultsL.R_err_ours = zeros(num_repeats,1);
resultsL.R_err_ukf = zeros(num_repeats,1);
resultsL.R_err_pace = zeros(num_repeats,1);
resultsL.p_err_ours = zeros(num_repeats,1);
resultsL.p_err_ukf = zeros(num_repeats,1);
resultsL.p_err_pace = zeros(num_repeats,1);
resultsL.c_err_ours = zeros(num_repeats,1);
resultsL.gap_ours = zeros(num_repeats,1);
resultsL.time_ours = zeros(num_repeats,1);
for j = 1:num_repeats

% Generate random tracking problem
problem.N_VAR = 11; % nr of keypoints
problem.K = 3; % nr of shapes

problem.L = L; % nr of keyframes in horizon

problem.outlierRatio = 0.0;
problem.noiseSigmaSqrt = 0.02; % [m]
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

soln_pace = pace_with_UKF(problem);

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
resultsL.R_err_ours(j) = norm(R_err_ours)/L;
resultsL.R_err_ukf(j)  = norm(R_err_ukf)/L;
resultsL.R_err_pace(j) = norm(R_err_pace)/L;
resultsL.p_err_ours(j) = norm(problem.p_gt - soln.p_est,'fro')/L;
resultsL.p_err_ukf(j)  = norm(problem.p_gt - soln_pace.p_smoothed,'fro')/L;
resultsL.p_err_pace(j) = norm(problem.p_gt - soln_pace.p_raw,'fro')/L;
resultsL.c_err_ours(j) = c_err;
resultsL.gap_ours(j) = soln.gap;
resultsL.time_ours(j) = soln.solvetime;
clear problem;
end
results(Ls == L) = resultsL;
end
% save
save("../datasets/results/results_L.mat","results")

%% Display Results
% process into displayable form

% Rotation figure
figure
set(0,'DefaultLineLineWidth',2)
plot([results.L],mean([results.R_err_ours]),'DisplayName','OURS');
hold on
plot([results.L],mean([results.R_err_ukf]),'DisplayName','PACE-UKF');
plot([results.L],mean([results.R_err_pace]),'DisplayName','PACE-RAW');
legend
xlabel("L"); ylabel("Rotation Error (deg)");
title("Rotation Errors")

% position figure
figure
plot([results.L],mean([results.p_err_ours]),'DisplayName','OURS');
hold on
plot([results.L],mean([results.p_err_ukf]),'DisplayName','PACE-UKF');
plot([results.L],mean([results.p_err_pace]),'DisplayName','PACE-RAW');
legend
xlabel("L"); ylabel("Position Error (m)");
title("Position Errors")

% gap figure
figure
semilogy([results.L],mean([results.gap_ours]));
xlabel("L"); ylabel("Gap");
title("Suboptimality Gaps")

% time figure
figure
plot([results.L],mean([results.time_ours]));
xlabel("L"); ylabel("Time (s)");
title("Solve Time")