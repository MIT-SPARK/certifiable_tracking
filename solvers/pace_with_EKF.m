function soln = pace_with_EKF(problem)
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

% TODO: SPIRAL VERSION!!

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

% FOR FULL DISTRIBUTION (TODO: FIX ALL THESE BEFORE TESTING!!)
% state covariance
covar_position = problem.covar_position(1:3)';
covar_velocity = problem.covar_velocity(1:3)';
kappa_rotation = problem.kappa_rotation(1:3)';
kappa_rotrate = problem.kappa_rotrate(1:3)';
covar_state_full = [covar_position, covar_velocity, kappa_rotation.^(-1), kappa_rotrate.^(-1)];
covar_state_full = diag(covar_state_full);

% process noise (noise added to const. vel. model)
% TODO: this maybe makes less sense
processNoise_full = repmat([0,0.01^2],1,6);
processNoise_full = diag(processNoise_full);

% Measurement noise (todo: this isn't quite right)
measureNoise_full = repmat(problem.noiseSigmaSqrt^2, 1,6);
measureNoise_full = diag(measureNoise_full);


%% Run PACE at each time step
soln_pace = [];
for l = 1:problem.L
    pace_problem = problem;
    pace_problem.weights = ones(problem.N_VAR,1);
    pace_problem.scene = reshape(problem.y(:,l),[3,problem.N_VAR]);
    [R_est,t_est,c_est,out] = outlier_free_category_registration(pace_problem, 'lambda',problem.lambda);
    s.R_est = R_est; s.p_est = t_est;
    s.c_est = c_est; s.out = out;
    s.s_est = R_est'*t_est;
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
        StateTransitionFcn=@constv,ProcessNoise=processNoise, ...
        MeasurementFcn=@measureModel_v,MeasurementNoise=measureNoise);

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

%% constant rotation rate version
s = [soln_pace.s_est]; % 3 x L
R = reshape([soln_pace.R_est],[3,3,L]); % 3 x 3 x L (CHECK)
L = size(p,2);

% prior: initial position/velocity, rotation/rate
prior_pos = s(:,2);
prior_r = R2r(R(:,:,2));

prior_dR = R(:,:,1)'*R(:,:,2);
prior_axang = rotm2axang(prior_dR);
prior_w = prior_axang(1:3)*prior_axang(4);

prior_v = (prior_dR*s(:,2) - s(:,1))/problem.dt;

prior_state = [prior_pos',prior_v',prior_r,prior_w]';

% define EKF
EKF2 = trackingEKF(State=prior_state,StateCovariance=covar_state_full, ...
        StateTransitionFcn=@constvw,ProcessNoise=processNoise_full, ...
        MeasurementFcn=@measureModel_vw,MeasurementNoise=measureNoise_full);

% Smooth measurements with EKF
s_smoothed = zeros(3,1,L);
s_smoothed(:,:,1) = s(:,1);
s_smoothed(:,:,2) = s(:,2);
v_smoothed2 = zeros(3,1,L);
v_smoothed2(:,:,2) = prior_v;
R_smoothed = zeros(3,3,L);
R_smoothed(:,:,1) = R(:,:,1);
R_smoothed(:,:,2) = R(:,:,2);
dR_smoothed = zeros(3,3,L);
dR_smoothed(:,:,2) = prior_dR;

p_smoothed2 = p_smoothed;
for l = 3:L
    % prediction
    [xpred, Ppred] = predict(EKF2,problem.dt);
    % correction
    [xcorr, Pcorr] = correct(EKF2,[s(:,l)',R2r(R(:,:,l))]');

    % save
    s_smoothed(:,:,l) = xcorr(1:3);
    v_smoothed2(:,:,l) = xcorr(4:6);
    R_smoothed(:,:,l) = r2R(xcorr(7:9));
    dR_smoothed(:,:,l) = axang2rotm([xcorr(10:12)'/norm(xcorr(10:12)), norm(xcorr(10:12))]);

    p_smoothed2(:,:,l) = R_smoothed(:,:,l)*s_smoothed(:,:,l);
end


%% Save
soln.p_smoothed = p_smoothed;
soln.v_smoothed = v_smoothed;
soln.p_raw = reshape(p,[3,1,L]);

soln.p_smoothed2 = p_smoothed2;

end

%% EKF State Functions--linear only
% compute next state from current state
function stateNext = constv(state,dt)
    % No noise constant velocity
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
function z = measureModel_v(state)
    % Measurement is just x, y, z position
    z = [state(1), state(3), state(5)];
end

%% EKF State Functions--Linear and Angular
% compute next state from current state
function stateNext = constvw(state,dt)
    % No noise constant velocity and turn rate
    % State: [x, y, z, vx, vy, vz,...
    %         r1, r2, r3, w1, w2, w3]
    % r: rodriguez parameters
    % see https://vanhunteradams.com/Estimation/MUKF.html
    s = state(1:3);
    v = state(4:6);
    r = state(7:9);
    w = state(10:12);

    % write angular velocity as dR
    dR = axang2rotm([w'/norm(w), norm(w*dt)]);

    % write r as rot matrix
    R = r2R(r);

    % compute next vals
    snext = dR'*(s + v*dt);
    vnext = v;
    Rnext = R*dR;
    wnext = w;

    % convert Rnext to rodriguez param
    rnext = R2r(Rnext)';
   
    stateNext = [snext; vnext; rnext; wnext];
end

% Compute nominal measurement from state
function z = measureModel_vw(state)
    % Measure x, y, z position & r1,r2,r3 angles
    z = [state(1:3); state(7:9)];
end

% convert R to rodriguez
function r = R2r(R)
    q = rotm2quat(R);
    r = 4*q(2:4) / (1 + q(1));
end
% convert rodriguez to R
function R = r2R(r)
    r = reshape(r,[length(r),1]);
    qw = (-r'*r + 16) / (r'*r + 16);
    qxyz = 1/4*(1+qw)*r;
    R = quat2rotm([qw, qxyz']);
end