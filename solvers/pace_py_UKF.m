function soln = pace_py_UKF(problem)
% A "classic" approach to tracking using PACE and an UKF for smoothing.
%    Estimate the position of the object using PACE at each time step and
%    incorporate linear priors/smoothing using a UKF.
% 
% INPUTS:
% - problem (struct): populated problem data
%
% RETURNS:
% - soln (struct): solution data
% 
% Lorenzo Shaikewitz for SPARK Lab

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

%% Define covariances and noise
% state covariance
covar_velocity = repmat(problem.covar_velocity(1),[1,3]);
covar_rotrate = repmat((2*problem.kappa_rotrate(1)).^(-1),[1,3]); % See SE-Sync

covar_position = covar_velocity.*covar_rotrate*problem.dt;
covar_rotation = covar_rotrate;

covar_state_full = [covar_position, covar_velocity, covar_rotation, covar_rotrate];
covar_state_full = diag(covar_state_full);
P = py.numpy.array(covar_state_full);

% process noise (noise added to const. vel. model)
% TODO: tune this?
processNoise_full = repmat([0.1^2,0.05^2],1,6);
processNoise_full = diag(processNoise_full);
Q = py.numpy.array(processNoise_full);

% Measurement noise (TODO: this is kinda cheating)
% measureNoise_full = repmat(problem.noiseSigmaSqrt^2, 1,6);
measureNoise_full = repmat(0.2^2,1,6);
measureNoise_full = diag(measureNoise_full);
R_covar = py.numpy.array(measureNoise_full);

%% Fuse with UKF
p = [soln_pace.p_est];
L = size(p,2);
R = reshape([soln_pace.R_est],[3,3,L]); % 3 x 3 x L (CHECK)

% remove the first element of p, R for UKF
p_meas = permute(p,[1,3,2]);
p_meas = py.numpy.array(p_meas(:,:,2:end));
R_meas = py.numpy.array(R(:,:,2:end));

% compute initial v, w from first two measurements.
v_init = R(:,:,1)'*(p(:,2) - p(:,1))/problem.dt;
dR_init = rotm2axang(R(:,:,1)'*R(:,:,2));
w_init = dR_init(1:3)'*dR_init(4);

% fuse!
out = py.solvers.ukf.run_ukf(problem.dt, L-1, p_meas, R_meas, v_init, w_init, Q, R_covar, P);
p_smoothed = double(out{1});
R_smoothed = double(out{2});

%% Save
soln.p_smoothed = p_smoothed;
soln.R_smoothed = R_smoothed;
soln.p_raw = reshape(p,[3,1,L]);
soln.R_raw = R;

end

%% UKF State Functions--Linear and Angular
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