%% Dense SDP relaxation for certifiable tracking
%  Generic, tunable script to run one iteration of dense tracking.
%    Operates on random data from PASCAL shapes with no outlier support.
%    Run setup.m once to set up paths.false
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all
% restoredefaultpath
% rng("default")
% rng(100)

%% Generate random tracking problem
problem.category = "aeroplane";
problem.L = 8; % nr of keyframes in horizon

problem.outlierRatio = 0.0;
problem.noiseSigmaSqrt = 0.1*0.2; % [m]
problem.noiseBound = 0.15*0.2;
problem.processNoise = 5e-2;

% MLE parameters
problem.accelerationNoiseBoundSqrt = 0.05*0.2;
problem.rotationKappa = 1/(0.05*0.2)^2*1/2;

problem.covar_measure_base = problem.noiseSigmaSqrt^2;
problem.covar_velocity_base = problem.accelerationNoiseBoundSqrt^2;
problem.kappa_rotrate_base = problem.rotationKappa;

pace_numbers = load("../datasets/results/pascalaeroplane_mle3_noiseSigmaSqrt.mat","results");
index = 4;
cur = pace_numbers.results(index);
problem.covar_measure_position = ones(1,3)*1e-10;%1*ones(1,3)*(mean(cur.p_err_pace.^2));
problem.covar_measure_rotation = ones(1,3)*1e-10;%1*ones(1,3)*(mean((cur.R_err_pace*pi/180.).^2));

problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

% regen if pbound, vbound, N, L, K change.
problem.regen_sdp = true; % when in doubt, set to true
problem.usecBound = false;

% Optional: use a specified velocity trajectory
% problem = make_trajectory(problem);
% problem.dR_gt = repmat(eye(3,3),[1,1,problem.L-1]);

% add shape, measurements, outliers
problem = gen_pascal_tracking(problem);
lambda = 0.0;
problem.lambda = lambda;

% problem.mosekpath = mosekpath;

%% Solve!
% soln = solve_weighted_tracking(problem);
pace = pace_raw(problem);

% temp
pace.p = problem.p_gt + 1*sqrt(problem.covar_measure_position(1))*randn(size(pace.p));
pace.R = problem.R_gt + 1*sqrt(problem.covar_measure_rotation(1))*randn(size(pace.R));

paceukf = pace_py_UKF(problem,pace);
% paceekf = pace_ekf(problem,pace);

mean(vecnorm(problem.p_gt(:,:,:) - pace.p(:,:,:)))
mean(vecnorm(problem.p_gt(:,:,:) - paceukf.p(:,:,:)))

figure
tiledlayout(3,1)
nexttile
plot(squeeze(pace.p(1,:,:)));
hold on
plot(squeeze(paceukf.p(1,:,:)));
plot(squeeze(problem.p_gt(1,:,:)),'k');

nexttile
plot(squeeze(pace.p(2,:,:)));
hold on
plot(squeeze(paceukf.p(2,:,:)));
plot(squeeze(problem.p_gt(2,:,:)),'k');

nexttile
plot(squeeze(pace.p(3,:,:)));
hold on
plot(squeeze(paceukf.p(3,:,:)));
plot(squeeze(problem.p_gt(3,:,:)),'k');

L = length(pace.R);
w_pace = zeros(3,L);
w_gt = zeros(3,L);
w_ukf = zeros(3,L);
for l = 1:length(pace.R)
    axang = rotm2axang(pace.R(:,:,l));
    w_pace(:,l) = axang(1:3)*axang(4);
    axang = rotm2axang(problem.R_gt(:,:,l));
    w_gt(:,l) = axang(1:3)*axang(4);
    axang = rotm2axang(paceukf.R(:,:,l));
    w_ukf(:,l) = axang(1:3)*axang(4);
end
figure
tiledlayout(3,1)
nexttile
plot(w_pace(1,:));
hold on
plot(w_ukf(1,:));
plot(w_gt(1,:),'k');

nexttile
plot(w_pace(2,:));
hold on
plot(w_ukf(2,:));
plot(w_gt(2,:),'k');

nexttile
plot(w_pace(3,:));
hold on
plot(w_ukf(3,:));
plot(w_gt(3,:),'k');
return


%% Check solutions
% eigenvalue plot
L = problem.L;
figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight

% if strcmp(problem.velprior, "body")
%     slices = 1:(1+9*(3*L-2)+3*L);
%     Xopt_pRemoved = soln.raw.Xopt{1}(slices, slices);
%     bar([zeros(3*(L-1),1);eig(Xopt_pRemoved)]);
%     title("Eigenvalues of Relaxed Solution")
% elseif strcmp(problem.velprior, "world")
%     slices = [1:(1+9*(3*L-2)),(1+9*(3*L-2)+3*L+1):(9*(3*L-2)+6*L)];
%     Xopt_pRemoved = soln.raw.Xopt{1}(slices, slices);
%     bar([zeros(3*L+1,1);eig(Xopt_pRemoved)]);
%     title("Eigenvalues of Relaxed Solution")
% elseif strcmp(problem.velprior, "grav-world")
%     error("Selected prior is not implemented")
% else
%     error("Selected prior is not implemented")
% end
% hold off

% raw error
% x_err = norm(problem.x_gt - soln.x_est);

% projected errors
R_err = zeros(L,1);
dR_err = zeros(L-1,1);
p_err = zeros(L,1);
v_err = zeros(L-1,1);
for l = 1:L
    % R
    R_err(l) = getAngularError(problem.R_gt(:,:,l), soln.R_est(:,:,l));
    % dR
    if (l < L)
        dR_err(l) = getAngularError(problem.dR_gt(:,:,l), soln.dR_est(:,:,l));
    end
    % p
    p_err(l) = norm(problem.p_gt(:,:,l) - soln.p_est(:,:,l));
    % v
    if (l < L)
        v_err(l) = norm(problem.v_gt(:,:,l) - soln.v_est(:,:,l));
    end
end

% shape error
c_err = norm(problem.c_gt - soln.c_est);

% Plot trajectory!
plot_trajectory2(problem,soln)

% temp for testing
% soln2.p = soln.p2_est;
% soln2.R = soln.R2_est;

compare(problem, soln, pace, paceukf, paceekf);

soln.gap_stable
soln.gap
% soln.gap2

function compare(gt, ours, pace, paceukf, paceekf)
L = gt.L;
% compare position
epace.p = vecnorm(gt.p_gt - pace.p);
eukf.p = vecnorm(gt.p_gt - paceukf.p);
eekf.p = vecnorm(gt.p_gt - paceekf.p);
eours.p = vecnorm(gt.p_gt - ours.p_est);

% compare rotation
epace.R = zeros(L,1);
eukf.R = zeros(L,1);
eours.R = zeros(L,1);
for l = 1:L
    epace.R(l) = getAngularError(gt.R_gt(:,:,l), pace.R(:,:,l));
    eukf.R(l) = getAngularError(gt.R_gt(:,:,l), paceukf.R(:,:,l));
    eours.R(l) = getAngularError(gt.R_gt(:,:,l), ours.R_est(:,:,l));
end

fprintf("           PACE    OUR2    OURS    LEKF \n")
fprintf("Position: %.4f, %.4f, %.4f, %.4f\n",epace.p(end),eukf.p(end),eours.p(end), eekf.p(end));
fprintf("Rotation: %.4f, %.4f, %.4f\n",mean(epace.R),mean(eukf.R),mean(eours.R));
end