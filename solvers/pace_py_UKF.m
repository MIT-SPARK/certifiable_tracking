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
% process covariance (Q)
if (isfield(problem,"covar_velocity_base"))
    covar_velocity = ones(1,3)*problem.covar_velocity_base;
end
if (isfield(problem,"kappa_rotrate_base"))
    covar_rotrate = ones(1,3)*(1/problem.kappa_rotrate_base*1/2);
elseif (isfield(problem,"covar_rotrate_base"))
    covar_rotrate = ones(1,3)*problem.covar_rotrate_base;
end

covar_process = [covar_velocity, covar_rotrate];
Q = py.numpy.array(diag(covar_process));

% Measurement covariance (R)
% read this from the chart (3*mean err)
covar_measure_position = problem.covar_measure_position;
covar_measure_rotation = problem.covar_measure_rotation;
R_covar = 1*py.numpy.array(diag([covar_measure_position, covar_measure_rotation]));

% (initial) Error Covariance (P)
covar_error = [covar_measure_position, covar_measure_rotation, covar_velocity, covar_rotrate];
P = py.numpy.array(diag(covar_error));

%% Fuse with UKF
p = pace.p;
L = size(p,3);
R = pace.R;

% remove the first element of p, R for UKF
p_meas = py.numpy.array(p(:,:,2:end));
R_meas = py.numpy.array(R(:,:,2:end));

% compute initial v, w from first two measurements.
v_init = R(:,:,1)'*(p(:,2) - p(:,1))/problem.dt;
dR_init = py.numpy.array(R(:,:,1)'*R(:,:,2));

% fuse!
out = py.solvers.ukf2.run_ukf(problem.dt, L-1, p_meas, R_meas, v_init, dR_init, Q, R_covar, P);
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