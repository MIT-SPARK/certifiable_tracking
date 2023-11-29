function plot_trajectory(problem,soln,soln_pace)
% Plots the true and estimated trajectory, lin. interp. btwn time steps.
%   Assumes v, dR holds constant between time steps (for now).
%   Note that there will be discontinuities in the estimated data, since we
%   do not strictly enforce that the trajectory is continuous.
%
% INPUTS:
% - problem (struct): populated problem data
% - soln (struct): populated solution data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Process inputs
if nargin == 2
    % don't plot PACE
    plot_pace = false;
else
    plot_pace = true;
end

%% Variables and domain
L = problem.L;
dt = problem.dt;

interval_pts = 50.0;
t = linspace(0,dt*L,interval_pts*L+1);

%% Extract ground truth trajectory
traj_gt = zeros(3,L*interval_pts);

for l = 1:L
    ran = (interval_pts*(l-1) + 1):(interval_pts*l);
    if l == L
        % add last point to last interval
        ran = [ran,ran(end)+1];
    end

    t_cur = t(ran);
    v_cur = problem.v_gt(:,:,l) * (t_cur - t_cur(1));
    p_cur = problem.p_gt(:,:,l) + v_cur;

    traj_gt(:,ran) = p_cur;
end

%% Extract estimated trajectory
traj_est = zeros(3,L*interval_pts+1);

for l = 1:L
    ran = (interval_pts*(l-1) + 1):(interval_pts*l);
    if l == L
        % add final point to last interval
        ran = [ran,ran(end)+1];
    end

    t_cur = t(ran);
    v_cur = soln.v_est(:,:,l) * (t_cur - t_cur(1));
    p_cur = soln.p_est(:,:,l) + v_cur;

    traj_est(:,ran) = p_cur;
end

%% Plot!
% 3D Line Plot
figure
plot3(traj_gt(1,:),traj_gt(2,:),traj_gt(3,:),'DisplayName','Ground Truth')
hold on
plot3(traj_est(1,:),traj_est(2,:),traj_est(3,:), 'DisplayName','Estimate')

if (plot_pace)
    t_pace = linspace(0,dt*(L-1), L);
    p_est_pace = [soln_pace.p_est];
    plot3(p_est_pace(1,:),p_est_pace(2,:),p_est_pace(3,:),'x','DisplayName','PACE')
end

title("3D Trajectories")
legend('Location','ne')

% Explicit x, y, z comparison
figure
subplot(3,1,1)
plot(t,traj_gt(1,:),'DisplayName','Ground Truth')
hold on
plot(t,traj_est(1,:),'DisplayName','Estimate')
if (plot_pace); plot(t_pace,p_est_pace(1,:),'x','DisplayName','PACE'); end
ylabel("x")

legend('Location','ne')
title("Explict Comparison of Evaluated Trajectories")

subplot(3,1,2)
plot(t,traj_gt(2,:),'DisplayName','Ground Truth')
hold on
plot(t,traj_est(2,:),'DisplayName','Estimate')
if (plot_pace); plot(t_pace,p_est_pace(2,:),'x','DisplayName','PACE'); end
ylabel("y")
subplot(3,1,3)
plot(t,traj_gt(3,:),'DisplayName','Ground Truth')
hold on
plot(t,traj_est(3,:),'DisplayName','Estimate')
if (plot_pace); plot(t_pace,p_est_pace(3,:),'x','DisplayName','PACE'); end
xlabel("time")
ylabel("z")


end