function soln = pace_ekf(problem, pace)
% A "classic" approach to tracking using PACE and an EKF for smoothing.
%    Estimate the position of the object using PACE at each time step and
%    use CTRV (constant turn and velocity
% 
% INPUTS:
% - problem (struct): populated problem data
%
% RETURNS:
% - soln (struct): solution data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Define covariances and noise
% process noise (noise added to const. vel. model)
covar_velocity = problem.covar_velocity_base;
if (isfield(problem,"covar_rotrate_base"))
    covar_rotrate = problem.covar_rotrate_base;
else
    covar_rotrate = 1/(2*problem.kappa_rotrate_base);
end
processNoise = [0., covar_velocity, 0, covar_velocity, covar_rotrate, 0, covar_velocity];
processNoise = diag(processNoise);

% Measurement noise (TODO: wrong)
measureNoise = repmat(problem.covar_measure_base, [1,3]);
measureNoise = 0.1*diag(measureNoise);

% state covariance: initially the same as measurement noise (TODO: wrong)
covarState = repmat(problem.covar_measure_base,[1,7]);
covarState = 0.1*diag(covarState);

%% Convert PACE data into EKF form
L = problem.L;
dt = problem.dt;
% pace raw data
p = squeeze(pace.p);
R = pace.R;
% p = squeeze(problem.p_gt);
% R = problem.R_gt;

% rotate position to align with gt axis of rotation
gtaxang = rotm2axang(problem.dR_gt(:,:,1));
dR0 = axang2rotm(vrrotvec(gtaxang(1:3),[0,0,1]));
R0 = problem.R_gt(:,:,1)';
R_to_axis = dR0*R0;
p_aligned = R_to_axis*(p - p(:,1));

% estimate velocity and rotation rate to create initial state
% [x, vx, y, vy, ω, z, vz]
state0 = zeros(1,7);
state0(1) = p_aligned(1,1);
state0(3) = p_aligned(2,1);
state0(5) = gtaxang(4)*180/pi;
state0(6) = p_aligned(3,1);
state0(7) = 0;

state0(2) = 0;
state0(4) = 0;

%% Fuse with EKF
% Create EKF
EKF = trackingEKF(@constturn,@ctmeas,state0, ...
    'StateTransitionJacobianFcn',@constturnjac, ...
    'MeasurementJacobianFcn',@ctmeasjac,...
    'MeasurementNoise', measureNoise,...
    'StateCovariance', covarState,...
    'ProcessNoise', processNoise);

% Smooth measurements with EKF
p_smoothed = zeros(3,1,L);
p_smoothed(:,:,1) = p(:,1);
R_smoothed = zeros(3,3,L);
R_smoothed(:,:,1) = R(:,:,1);
for l = 2:L
    % prediction
    [xpred, Ppred] = predict(EKF,problem.dt);
    % correction
    [xcorr, Pcorr] = correct(EKF,p_aligned(:,l));

    % convert to p, R
    p_est = [xcorr(1); xcorr(3); xcorr(6)];
    p_smoothed(:,:,l) = R_to_axis'*p_est + p(:,1);
end

%% Save
soln.p = p_smoothed;
soln.R = R_smoothed;

end