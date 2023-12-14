function pts = sim_dynamics(R, dR, v, p0, dt, num_pts, use_v_true)
%% Simulate the dynamics of motion using a small time approximation.
% use_v_true: generates a trajectory based on the true velocity of the
% particle or instead based on the estimated velocity, which has straight
% lines between points.

%% Preliminaries
% Extract angular speed
w = rotm2axang(dR);
rot_axis = w(1:3)';
rot_speed = w(4);

% Set number of intermediate points
tot_num_pts = 1+100*num_pts + 1;
del_t = dt / (tot_num_pts-1); % time between intervals

pts = zeros(3,1,tot_num_pts);

%% Approx path
pts(:,:,1) = p0;
last_R = R;
for i = 2:tot_num_pts
    % estimate lin and ang distance travelled
    pv_cur = v*del_t;
    dR_cur = axang2rotm([rot_axis', rot_speed*del_t]);

    % get next point
    pts(:,:,i) = pts(:,:,i-1) + last_R*pv_cur;

    if (use_v_true)
        last_R = last_R*dR_cur;
    end
end

%% Finalize
pts = pts(:,:,2:100:end);

end