%% RSS Experiment: What Time Step Should We Use?
% Dataset: pascal + car
% Constants: K, N, noiseSigma, NO nonlinearities in gt
% Independent variable: L
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
%
% Lorenzo Shaikewitz for SPARK Lab

% BROKEN!!!!!!

clc; clear; close all

%% Experiment settings
indepVar = "outlierratio"; % name of independent variable
savename = "pascalcar_" + indepVar;
domain = [0.05:0.05:0.95];
% domain = 0.5;
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

problem.L = 10; % nr of keyframes in horizon
L = problem.L;
problem.category = "car";

problem.outlierRatio = iv;
problem.noiseSigmaSqrt = 0.05; % [m]
problem.noiseBound = 3*problem.noiseSigmaSqrt;
problem.noiseBound_GNC = 0.15;
problem.noiseBound_GRAPH = 0.2;
problem.processNoise = 0.5;
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

% for GNC
problem.N = problem.N_VAR*problem.L; % How many measurements this problem has (updated by ROBIN)
problem.outliers = []; % outlier indicies
problem.priors = [];
problem.dof = 3;

% Solve!
soln_pace = pace_py_UKF(problem,true,false);
soln_pace_robin = pace_py_UKF(problem,true,true);

problem = lorenzo_prune(problem);

try
    [inliers, info] = gnc_custom(problem, @solver_for_gnc, 'NoiseBound', problem.noiseBound_GNC,'MaxIterations',100,'FixPriorOutliers',false);
    disp("GNC finished " + string(j))

    soln = info.f_info.soln;
    ef = eig(soln.raw.Xopt{1});
    if (ef(end-4) > 1e-4)
        disp("**Not convergent: " + string(soln.gap_nov))
    end
catch
    soln.p_est = ones(3,1,problem.L)*NaN;
    soln.R_est = ones(3,3,problem.L)*NaN;
    soln.c_est = ones(problem.K,1)*NaN;
    soln.gap = NaN;
    soln.solvetime = NaN;
    disp("GNC failed " + string(j))
end

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