%% Experiment: How does process noise affect UKF performance?
% Dataset: synthetic
% Constants: K, N, noiseSigma, NO nonlinearities in gt
% Independent variable: process noise
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment settings
indepVar = "processNoise"; % name of independent variable
savename = "pascalcar_" + indepVar;
domain = [1e-8, 1e-7, 1e-6, 0.00001, 0.0001, 0.0005, 0.00075, 0.001, 0.0025 0.005]; % for quick results
num_repeats = 150;
% SET INDEPENDENT VARIABLE, DEPENDENT VARS CORRECTLY IN LOOP

%% Loop
results = cell(length(domain),1);
parfor index = 1:length(domain)
iv = domain(index);
resultsIV = struct();
resultsIV.(indepVar) = iv;
R_err_ours = zeros(num_repeats,1);
R_err_ukf = zeros(num_repeats,1);
R_err_pace = zeros(num_repeats,1);
p_err_ours = zeros(num_repeats,1);
p_err_ukf = zeros(num_repeats,1);
p_err_pace = zeros(num_repeats,1);
c_err_ours = zeros(num_repeats,1);
gap_ours = zeros(num_repeats,1);
time_ours = zeros(num_repeats,1);
disp("Starting " + indepVar + "=" + string(iv));
for j = 1:num_repeats

% Generate random tracking problem
problem = struct();
problem.L = 6; % nr of keyframes in horizon
L = problem.L;
problem.category = "aeroplane";

problem.outlierRatio = 0.0;
problem.noiseSigmaSqrt = 0.05*0.02; % [m]
problem.noiseBound = 0.15*0.02;

problem.covar_measure_base = 0.01;
problem.covar_velocity_base = 0.001;
problem.covar_rotrate_base = 0.001;

problem.processNoise = iv;

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
% problem.sdp_filename = "sdpdata" + indepVar;

% add shape, measurements, outliers
problem = gen_pascal_tracking(problem);
lambda = 0;
problem.lambda = lambda;

% Solve!
soln = solve_weighted_tracking(problem);
pace = pace_raw(problem);
% paceukf = pace_py_UKF(problem,pace);
paceukf = pace_ekf(problem,pace);

% Save solutions
% projected errors
R_err_ours2 = zeros(L,1);
R_err_ukf2 = zeros(L,1);
R_err_pace2 = zeros(L,1);
for l = 1:L
    R_err_ours2(l) = getAngularError(problem.R_gt(:,:,l), soln.R_est(:,:,l));
    R_err_ukf2(l) = getAngularError(problem.R_gt(:,:,l), paceukf.R(:,:,l));
    R_err_pace2(l) = getAngularError(problem.R_gt(:,:,l), pace.R(:,:,l));
end

% shape error
c_err = norm(problem.c_gt - soln.c_est);

% save
R_err_ours(j) = R_err_ours2(end);
R_err_ukf(j)  = R_err_ukf2(end);
R_err_pace(j) = R_err_pace2(end);
p_err_ours(j) = norm(problem.p_gt(:,:,end) - soln.p_est(:,:,end));
p_err_ukf(j)  = norm(problem.p_gt(:,:,end) - paceukf.p(:,:,end));
p_err_pace(j) = norm(problem.p_gt(:,:,end) - pace.p(:,:,end));
c_err_ours(j) = c_err;
gap_ours(j) = soln.gap;
time_ours(j) = soln.solvetime;
end
resultsIV.R_err_ours = R_err_ours;
resultsIV.R_err_ukf  = R_err_ukf;
resultsIV.R_err_pace = R_err_pace;
resultsIV.p_err_ours = p_err_ours;
resultsIV.p_err_ukf  = p_err_ukf;
resultsIV.p_err_pace = p_err_pace;
resultsIV.c_err_ours = c_err_ours;
resultsIV.gap_ours   = gap_ours;
resultsIV.time_ours  = time_ours;
results{index} = resultsIV;
end
results = [results{:}];
% save
save("../datasets/results/" + savename + ".mat","results")

%% Display Results
% process into displayable form

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
a=plot([results.(indepVar)],median([results.p_err_ours]),'x-','DisplayName','OURS');
hold on
b=plot([results.(indepVar)],median([results.p_err_ukf]),'x-','DisplayName','PACE-UKF');
c=plot([results.(indepVar)],median([results.p_err_pace]),'x-','DisplayName','PACE-RAW');

errorshade([results.(indepVar)],[results.p_err_ours],get(a,'Color'));
errorshade([results.(indepVar)],[results.p_err_ukf],get(b,'Color'));
errorshade([results.(indepVar)],[results.p_err_pace],get(c,'Color'));
legend
xlabel(indepVar); ylabel("Position Error (m)");
title("Position Errors")
xscale log

% gap figure
figure
semilogy([results.(indepVar)],abs(median([results.gap_ours])),'x-');
xlabel(indepVar); ylabel("Gap");
title("Suboptimality Gaps")

% time figure
figure
plot([results.(indepVar)],median([results.time_ours]),'x-');
xlabel(indepVar); ylabel("Time (s)");
title("Solve Time")