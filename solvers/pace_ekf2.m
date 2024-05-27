function soln = pace_ekf2(problem, pace)
% A "classic" approach to tracking using PACE and an EKF for smoothing.
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
processNoise = [ones(1,3)*covar_velocity, 1e-2*ones(1,3)*covar_rotrate];
processNoise = 1*diag(processNoise);

% Measurement noise (TODO: wrong)
measureNoise = [problem.covar_measure_position, problem.covar_measure_rotation];
measureNoise = 1*diag(measureNoise);

% state covariance: initially the same as measurement noise (TODO: wrong)
covarState = 1*eye(12,12);
% covarState = diag(covarState);

%% Initialize
L = problem.L;
% pace raw data
p = pace.p;
R = pace.R;
s = p*0;
for l = 1:L
    s(:,:,l) = R(:,:,l)'*p(:,:,l);
end

state0 = [problem.p_gt(:,:,1); so3log(problem.R_gt(:,:,1)); problem.v_gt(:,:,1); so3log(problem.dR_gt(:,:,1))];

%% Fuse with EKF
% Create EKF
EKF = trackingEKF(@consttwist,@pacemeas,state0, ...
    'HasAdditiveProcessNoise', false, 'HasAdditiveMeasurementNoise', false,...
    'MeasurementNoise', measureNoise,...
    'StateCovariance', covarState,...
    'ProcessNoise', processNoise);

% Smooth measurements with EKF
p_smoothed = zeros(3,1,L);
s_smoothed = zeros(3,1,L);
R_smoothed = zeros(3,3,L);
for l = 1:L
    % prediction
    [xpred, Ppred] = predict(EKF,problem.dt);

    % correction
    meas = [p(:,:,l); so3log(R(:,:,l))];
    [xcorr, Pcorr] = correct(EKF,meas);

    % convert to p, R
    p_smoothed(:,:,l) = xcorr(1:3);
    R_smoothed(:,:,l) = so3exp(xcorr(4:6));
    
    % p_smoothed(:,:,l) = R_smoothed(:,:,l)*s_smoothed(:,:,l);
end

%% Save
soln.s = s_smoothed;
soln.p = p_smoothed;
soln.R = R_smoothed;

end

function statenew = consttwist(state, w, dt)
% state: [s, log(R), v, log(dR)]
p = state(1:3);
R = so3exp(state(4:6));
v = state(7:9);
dR = so3exp(state(10:12));

% update p, v
p2 = (p + R*v*dt);
R2 = R*dR;

% update v, dR (with noise)
v2 = v + w(1:3);
dR2 = dR*so3exp(w(4:6));

statenew = [p2; so3log(R2); v2; so3log(dR2)];
end

function measurement = pacemeas(state, v)
p = state(1:3);
R = so3exp(state(4:6));

% add noise
p_m = p + v(1:3);
R_m = R*so3exp(v(4:6));

measurement = [p_m; so3log(R_m)];
end

function R = so3exp(twist)
    twist = twist*1e0;
    if (norm(twist) < 1e-12)
        axang = [1,1,1,0];
    else
        axang = [twist/norm(twist); norm(twist)]';
    end
    R = axang2rotm(axang);
end
function twist = so3log(R)
    axang = rotm2axang(R);
    twist = axang(1:3)'*axang(4);
    twist = twist*1e-0;
end