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

%% Scale covariances
L = problem.L;
if isfield(problem,'prioroutliers')
    % add ROBIN priors to weights
    w = ones(1,problem.N_VAR*L-length(problem.prioroutliers));

    for i = 1:length(problem.prioroutliers)
        o = problem.prioroutliers(i);
        w = [w(1:o-1),0.0,w(o:end)];
    end
    problem.covar_measure = reshape(w.^(-1),[problem.N_VAR, problem.L]);
elseif isfield(problem,'covar_measure')
    disp("OVERRIDING COVAR MEASURE. YOU PROBABLY DON'T WANT TO DO THIS.")
    % use only for testing: skipping points or other tests
else
    problem.covar_measure = ones(problem.N_VAR,L);
end
noiseBoundSq = problem.noiseBound^2;
problem.covar_measure = problem.covar_measure.*((noiseBoundSq/9));
if (isfield(problem,"covar_velocity_base"))
    problem.covar_velocity = ones(L-2,1)*problem.covar_velocity_base;
else
    base = mean(problem.covar_measure(~isinf(problem.covar_measure)));
    problem.covar_velocity = ones(L-2,1)*base;
end
if (isfield(problem,"kappa_rotrate_base"))
    problem.kappa_rotrate = ones(L-2,1)*problem.kappa_rotrate_base;
else
    base = mean(problem.covar_velocity(~isinf(problem.covar_velocity)));
    problem.kappa_rotrate  = ones(L-2,1)*(2/base);
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
if (isfield(problem,"processNoise"))
    pn = problem.processNoise;
else
    pn = 0.05;
end
processNoise_full = repmat(pn^2,1,12);
processNoise_full = diag(processNoise_full);
Q = py.numpy.array(processNoise_full);

% Measurement noise (TODO: this is kinda cheating)
% measureNoise_full = repmat(problem.noiseSigmaSqrt^2, 1,6);
measureNoise_full = repmat(problem.noiseBound^2,1,6);
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

p_smoothed_full = zeros(3,1,L);
p_smoothed_full(:,:,1:2) = reshape(p(:,1:2),[3,1,2]);
p_smoothed_full(:,:,3:end) = p_smoothed;
R_smoothed_full = zeros(3,3,L);
R_smoothed_full(:,:,1:2) = R(:,:,1:2);
R_smoothed_full(:,:,3:end) = R_smoothed;

%% Save
soln.p_smoothed = p_smoothed_full;
soln.R_smoothed = R_smoothed_full;

end