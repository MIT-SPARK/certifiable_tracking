function pts = get_spiral_pts(R, dR, v, p0, dt, num_pts)
    %% Produces points in a spiral shape.
    %    NOTE: VERIFIED USING SIM DYNAMICS!
    %
    %    R = w^R_l (rotation of frame l in world frame)
    %   dR = (l)^R_(l+1) (rotation of next frame in prev frame)
    %   v  = body velocity, in frame l
    %   p0 = initial position, in frame w
    %   dt = time between points
    %   num_pts = number of points to generate
    %
    % Lorenzo Shaikewitz for SPARK Lab
    
    %% Extract angular speed
    w = rotm2axang(dR);
    rot_axis = w(1:3)';
    rot_speed = w(4);
    
    t = 0:dt:((num_pts-1)*dt);
    
    % Case: no angular speed
    if (rot_speed == 0)
        % just use linear velocity
        pts = p0 + R*v*t;
        return
    end
    
    %% Compute linear part of motion
    % v_ll = proj_axis(v)
    v_ll = (v'*rot_axis)*rot_axis;
    
    % use axis of rotation to get R_from_z
    R_axis_to_z = axang2rotm(vrrotvec(rot_axis',[0,0,1]));
    if (~all(v_ll==0))
        % may be positive or negative
        lin_part = (R_axis_to_z*v_ll)*t;
    else
        lin_part = [0;0;0]*t;
    end
    z_move = lin_part(3,:);
    
    %% Compute circular part of motion
    v_perp = v - v_ll;
    % rotate into x-y plane
    x_y_vel = R_axis_to_z*v_perp;
    % rotate to x direction
    R_to_x = axang2rotm(vrrotvec(x_y_vel',[1,0,0]));
    x_vel = R_to_x*x_y_vel;
    
    % get radii (v = w*r)
    r = x_vel(1)/rot_speed;
    
    % compute circle (xdot(0) = w*r, ydot(0) = 0)
    x_circle = r*sin(rot_speed*t);
    y_circle = -r*cos(rot_speed*t) + r;
    
    % *start* at 0
    % x_circle = x_circle - x_circle(1);
    % y_circle = y_circle - y_circle(1);
    
    %% Rotate back into frame l
    spiral_zeroframe = [x_circle; y_circle; z_move];
    spiral_lframe = R_axis_to_z'*R_to_x'*spiral_zeroframe;
    spiral_wframe_offset = R*spiral_lframe;
    spiral_wframe = p0 + spiral_wframe_offset;% - spiral_wframe_offset(:,1);
    
    pts = spiral_wframe;
    
    end