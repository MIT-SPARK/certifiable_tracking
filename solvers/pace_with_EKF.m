function soln = pace_with_EKF(problem, path)
% A "classic" approach to tracking using PACE and an EKF for smoothing.
%    Estimate the position of the object using PACE at each time step and
%    incorporate linear priors/smoothing using an EKF.
% 
% INPUTS:
% - problem (struct): populated problem data
%
% RETURNS:
% - soln (struct): solution data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Define covariances and noise
% state covariance
covar_position = problem.covar_position';
covar_velocity = problem.covar_velocity';
covar_state = reshape([covar_position(1:3); covar_velocity(1:3)], 6,1);
covar_state = diag(covar_state);

% process noise (noise added to const. vel. model)
% TODO: this maybe makes less sense
processNoise = repmat([0,0.01^2],1,3);
processNoise = diag(processNoise);

% Measurement noise (todo: this isn't quite right)
measureNoise = repmat(problem.noiseSigmaSqrt^2, 1,3);
measureNoise = diag(measureNoise);


%% Run PACE at each time step
soln_pace = [];
for l = 1:problem.L
    pace_problem = problem;
    pace_problem.weights = ones(problem.N_VAR,1);
    pace_problem.scene = reshape(problem.y(:,l),[3,problem.N_VAR]);
    [R_est,t_est,c_est,out] = outlier_free_category_registration(pace_problem, path, 'lambda',problem.lambda);
    s.R_est = R_est; s.p_est = t_est;
    s.c_est = c_est; s.out = out;
    soln_pace = [soln_pace; s];
end

%% Pull out position, orientation
p = [soln_pace.p_est];
L = size(p,2);

%% Fuse with EKF
% Create EKF
% Prior: initial position/velocity
prior_pos = p(:,2);
prior_vel = (p(:,2)-p(:,1))/problem.dt;

prior_state = reshape([prior_pos'; prior_vel'],6,1);

% Define EKF
% EKF = trackingEKF(@constvel,@cvmeas,prior_state, ...
%     'StateTransitionJacobianFcn',@constveljac, ...
%     'MeasurementJacobianFcn',@cvmeasjac);

EKF = trackingEKF(State=prior_state,StateCovariance=covar_state, ...
        StateTransitionFcn=@constvw,ProcessNoise=processNoise, ...
        MeasurementFcn=@measureModel,MeasurementNoise=measureNoise);

% Smooth measurements with EKF
p_smoothed = zeros(3,1,L);
p_smoothed(:,:,1) = p(:,1);
p_smoothed(:,:,2) = p(:,2);
v_smoothed = zeros(3,1,L);
v_smoothed(:,:,2) = prior_vel;
for l = 3:L
    % prediction
    [xpred, Ppred] = predict(EKF,problem.dt);
    % correction
    [xcorr, Pcorr] = correct(EKF,p(:,l));

    % save
    xcorr = reshape(xcorr,2,3);
    p_smoothed(:,:,l) = xcorr(1,:)';
    v_smoothed(:,:,l) = xcorr(2,:)';
end

%% Save
soln.p_smoothed = p_smoothed;
soln.v_smoothed = v_smoothed;
soln.p_raw = reshape(p,[3,1,L]);

end

%% EKF State Functions
% compute next state from current state
function stateNext = constvw(state,dt)
    % No noise constant velocity and turn rate
    % State: [x, vx, y, vy, z]
    F = [1 dt 0  0 0  0; 
            0  1 0  0 0  0;
            0  0 1 dt 0  0;
            0  0 0  1 0  0;
            0  0 0  0 1 dt;
            0  0 0  0 0  1;];
    stateNext = F*state;
end

% Compute nominal measurement from state
function z = measureModel(state)
    % Measurement is just x, y, z position
    z = [state(1), state(3), state(5)];
end