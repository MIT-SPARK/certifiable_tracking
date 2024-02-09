function plot_trajectory2(problem,soln, varargin)
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

%% Preliminaries
L = problem.L;
intervalPoints = 50.0;
dt = problem.dt;

parser = inputParser;
addParameter(parser,'PlotAxes',true,@isboolean);
parse(parser,varargin{:});

%% Plot ground truth poses
figure
p_gt = reshape(problem.p_gt,[3,L,1]);
plot3(p_gt(1,:),p_gt(2,:),p_gt(3,:),'.k', 'MarkerSize',10);
hold on
axis equal

% Plot axes
if (parser.Results.PlotAxes)
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)', ...
    squeeze(problem.R_gt(1,1,:)),squeeze(problem.R_gt(2,1,:)),squeeze(problem.R_gt(3,1,:)),0.25,'r');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)', ...
    squeeze(problem.R_gt(1,2,:)),squeeze(problem.R_gt(2,2,:)),squeeze(problem.R_gt(3,2,:)),0.25,'g');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)', ...
    squeeze(problem.R_gt(1,3,:)),squeeze(problem.R_gt(2,3,:)),squeeze(problem.R_gt(3,3,:)),0.25,'b');
end

% Plot trajectory
traj_gt = zeros(3,(L-1)*intervalPoints);
for l = 1:L-1
    ran = (intervalPoints*(l-1) + 1):(intervalPoints*l);
    R = problem.R_gt(:,:,l);
    dR = problem.dR_gt(:,:,l);
    v = problem.v_gt(:,:,l);
    p = problem.p_gt(:,:,l);
    % generate spiral
    pts = get_spiral_pts(R, dR, v, p, dt/intervalPoints, 1+intervalPoints);
    traj_gt(:, ran) = pts(:,1:end-1);
    if (l == (L-1))
        traj_gt(:,end+1) = pts(:,end);
    end
end
plot3(traj_gt(1,:),traj_gt(2,:),traj_gt(3,:),'DisplayName','Ground Truth','Color',[0.4660 0.6740 0.1880])

%% Plot estimated trajectory
p_est = reshape(soln.p_est,[3,L,1]);
plot3(p_est(1,:),p_est(2,:),p_est(3,:),'xb', 'MarkerSize',10);

% Plot axes
if (parser.Results.PlotAxes)
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)', ...
    squeeze(soln.R_est(1,1,:)),squeeze(soln.R_est(2,1,:)),squeeze(soln.R_est(3,1,:)),0.25,'r');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)', ...
    squeeze(soln.R_est(1,2,:)),squeeze(soln.R_est(2,2,:)),squeeze(soln.R_est(3,2,:)),0.25,'g');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)', ...
    squeeze(soln.R_est(1,3,:)),squeeze(soln.R_est(2,3,:)),squeeze(soln.R_est(3,3,:)),0.25,'b');
end

% Plot trajectory
traj_est = zeros(3,(L-1)*intervalPoints);
for l = 1:L-1
    ran = (intervalPoints*(l-1) + 1):(intervalPoints*l);
    R = soln.R_est(:,:,l);
    dR = soln.dR_est(:,:,l);
    v = soln.v_est_corrected(:,:,l);
    p = soln.p_est(:,:,l);
    % generate spiral
    t = 0:dt/intervalPoints:dt;
    pts = p + R*v*t;
    traj_est(:, ran) = pts(:,1:end-1);
end
plot3(traj_est(1,:),traj_est(2,:),traj_est(3,:),'DisplayName','Estimate','Color',[0.8500 0.3250 0.0980])


end