%% Test with Different R, dR
% L = 10;

%  Rs = repmat(eye(3), 1,1,L);
% dRs = repmat(axang2rotm([0,0,1,1]),1,1,L-1);
% for l = 2:L
%     Rs(:,:,l) = Rs(:,:,l-1)*dRs(:,:,l-1);
% end
vs = repmat([69.8303267605827;-11.8892471262851;-32.8564132917302],1,1,L);
% vs = repmat([-166.016018430323;-1047.64157595129;28.7619028834964],1,1,L);
% vs = repmat([221.800517521738;-21.8419553751813;481.804384694389],1,1,L);
% vs = repmat([-0.995837938256655;0.244557022344240;0.298766772260149],1,1,L);
% ps = zeros(3,1,L);

% vs = soln.v_est;
% ps = soln.p_est;
% Rs = soln.R_est;
% dRs = soln.dR_est;
% vs = v_to_v(vs,dRs,problem.dt,L);

vs = problem.v_gt;
ps = problem.p_gt;
Rs = problem.R_gt;
dRs = problem.dR_gt;

% vs = problem.v_gt;
% ps = soln.p_est;
% Rs = soln.R_est;
% dRs = soln.dR_est;

dt = 1;%problem.dt;

num_intermediate = 200;

plot_pts = zeros(3,1,(L-1)*num_intermediate);
for l = 1:L-1
    % current values
    R = Rs(:,:,l);
    dR = dRs(:,:,l);
    v = vs(:,:,l);
    p = ps(:,:,l);

    % spiral to next
    % pts = get_spiral_pts(R, dR, v, p, dt/(num_intermediate), num_intermediate + 1);
    pts = sim_dynamics(R, dR, v, p, dt, num_intermediate, false);

    ps(:,:,l+1) = pts(:,end);
    plot_pts(:,:,((l-1)*num_intermediate+1):(l*num_intermediate)) = pts(:,:,1:end-1);
    if (l == (l-1))
        plot_pts(:,:,end+1) = pts(:,end);
    end
end

% figure
plot3(plot_pts(1,:), plot_pts(2,:), plot_pts(3,:))
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
hold on