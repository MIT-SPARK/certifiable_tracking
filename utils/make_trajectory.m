function problem = make_trajectory(problem)
%% Define a noise-free linear and angular velocity trajectory
%     many options: see below
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Setup
% problem variables
L = problem.L;
velocityBound = problem.velocityBound;

% generate random starting velocities
v_0 = randn(3,1);
v_0 = v_0/norm(v_0);
v_0 = velocityBound*rand*v_0;
v_gt = repmat(v_0,1,1,L);

dR_gt = rand_rotation;
dR_gt = repmat(dR_gt,1,1,L);

%% Jump in Linear Velocity
% for  l = 1:L
%     if (l > L/2)
%         v_gt(1,:,l) = -v_gt(1,:,l);
%         v_gt(2,:,l) = 2*v_gt(2,:,l);
%     end
% end

%% Quadratic motion (linear velocity)
% fix acceleration
% NOTE: this breaks it for large accelerations (e.g. 1.0)
acc = [0.;0.;-0.1];
for  l = 1:L
    v_gt(:,:,l) = v_gt(:,:,l) + acc*problem.dt*(l-1);
end

%% Save
problem.dR_gt = dR_gt;
problem.v_gt = v_gt;

end