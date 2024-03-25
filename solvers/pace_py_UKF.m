function soln = pace_py_UKF(problem, pace)
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

%% Define covariances and noise
rotBase = 1; % use 0.1% for rotation due to EKF operating in radians
% state covariance
if (isfield(problem,"covar_velocity_base"))
    covar_velocity = ones(1,3)*problem.covar_velocity_base;
else
    covar_velocity = ones(1,3)*(problem.noiseBound/3)^2;
end
if (isfield(problem,"kappa_rotrate_base"))
    covar_rotrate = ones(1,3)*(1/problem.kappa_rotrate_base*1/2)*rotBase;
elseif (isfield(problem,"covar_rotrate_base"))
    covar_rotrate = ones(1,3)*problem.covar_rotrate_base;
else
    covar_rotrate  = ones(1,3)*(problem.noiseBound/3)^2*rotBase;
end

covar_position = covar_velocity.*covar_rotrate*problem.dt;
covar_rotation = covar_rotrate;

covar_state_full = [covar_position, covar_velocity, covar_rotation, covar_rotrate];
covar_state_full = diag(covar_state_full);
P = py.numpy.array(covar_state_full);

% process noise (noise added to const. vel. model)
if (isfield(problem,"processNoise"))
    pn = problem.processNoise;
else
    pn = 0.05;
end
processNoise_full = [pn, pn, pn, rotBase*pn, rotBase*pn, rotBase*pn];
processNoise_full = diag(processNoise_full);
Q = py.numpy.array(processNoise_full);

% Measurement noise
measureNoise_full = repmat(problem.noiseBound^2/9,1,6);
measureNoise_full(4:end) = rotBase*measureNoise_full(4:end);
measureNoise_full = diag(measureNoise_full);
R_covar = py.numpy.array(measureNoise_full);

%% Fuse with UKF
p = pace.p;
L = size(p,3);
R = pace.R;

% remove the first element of p, R for UKF
p_meas = py.numpy.array(p(:,:,2:end));
R_meas = py.numpy.array(R(:,:,2:end));

% compute initial v, w from first two measurements.
v_init = R(:,:,1)'*(p(:,2) - p(:,1))/problem.dt;
dR_init = rotm2axang(R(:,:,1)'*R(:,:,2));
w_init = dR_init(1:3)'*dR_init(4) / problem.dt;

% fuse!
out = py.solvers.ukf.run_ukf(problem.dt, L-1, p_meas, R_meas, v_init, w_init, Q, R_covar, P);
p_smoothed = double(out{1});
R_smoothed = double(out{2});

p_smoothed_full = zeros(3,1,L);
p_smoothed_full(:,:,1:2) = reshape(p(:,1:2),[3,1,2]);
p_smoothed_full(:,:,3:end) = p_smoothed;
R_smoothed_full = zeros(3,3,L);
R_smoothed_full(:,:,1:2) = R(:,:,1:2);
R_smoothed_full(:,:,3:end) = R_smoothed;

%% Save
soln.p = p_smoothed_full;
soln.R = R_smoothed_full;

end