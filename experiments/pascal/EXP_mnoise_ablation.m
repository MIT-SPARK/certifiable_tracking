%% IROS Experiment: How does noise affect performance?
% Dataset: pascal + car
% Constants: K, N, L, NO nonlinearities in gt
% Independent variable: noiseSigma (measurement)
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
% STATUS: NEEDS PACE UKF
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment settings
indepVar = "noiseSigmaSqrt"; % name of independent variable
savename = "pascalcar2_" + indepVar;
domain = [0.01:0.005:0.1];
num_repeats = 50;
% SET INDEPENDENT VARIABLE, DEPENDENT VARS CORRECTLY IN LOOP

%% Loop
results = cell(length(domain),1);
parfor index = 1:length(domain)
iv = domain(index)
resultsIV = struct();
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

problem = struct();
problem.L = 10; % nr of keyframes in horizon
L = problem.L;
problem.category = "car";

problem.outlierRatio = 0.0;
problem.noiseSigmaSqrt = iv; % [m]
problem.noiseBound = 0.1;
problem.processNoise = 0.12;
problem.covar_velocity_base = 0.001;
problem.kappa_rotrate_base = 500;

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
problem = gen_pascal_tracking(problem);
lambda = 0.0;
problem.lambda = lambda;

% Solve!
soln = solve_weighted_tracking(problem);
pace = pace_raw(problem);
paceukf = pace_py_UKF(problem,pace);

% Save solutions
% projected errors
R_err_ours = zeros(L,1);
R_err_ukf = zeros(L,1);
R_err_pace = zeros(L,1);
for l = 1:L
    R_err_ours(l) = getAngularError(problem.R_gt(:,:,l), soln.R_est(:,:,l));
    R_err_ukf(l) = getAngularError(problem.R_gt(:,:,l), paceukf.R(:,:,l));
    R_err_pace(l) = getAngularError(problem.R_gt(:,:,l), pace.R(:,:,l));
end

% shape error
c_err = norm(problem.c_gt - soln.c_est);

% save
resultsIV.R_err_ours(j) = norm(R_err_ours)/L;
resultsIV.R_err_ukf(j)  = norm(R_err_ukf)/L;
resultsIV.R_err_pace(j) = norm(R_err_pace)/L;
resultsIV.p_err_ours(j) = norm(problem.p_gt - soln.p_est,'fro')/L;
resultsIV.p_err_ukf(j)  = norm(problem.p_gt - paceukf.p,'fro')/L;
resultsIV.p_err_pace(j) = norm(problem.p_gt - pace.p,'fro')/L;
resultsIV.c_err_ours(j) = c_err;
resultsIV.gap_ours(j) = soln.gap;
resultsIV.time_ours(j) = soln.solvetime;
end
results{index} = resultsIV;
end
results = [results{:}];
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