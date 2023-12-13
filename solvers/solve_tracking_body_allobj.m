function soln = solve_tracking_body_allobj(problem)
% Solves const. vel. (body frame) optimization exactly via SDP
%   Assume the *body frame* velocity is constant and object is spinning.
%   Result: spiral trajectory.
%   Analytically remove velocity & shape. SDP variables are
%   * rotated position (s)
%   * last rotated position (sh)
%   * rotation (R)
%   * rotation change (dR)
%   * predicted rotation (Rh)
%
% INPUTS:
% - problem (struct): populated problem data
%
% RETURNS:
% - soln (struct): solution data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Process inputs
mosekpath = problem.mosekpath;

N = problem.N_VAR;
K = problem.K;
L = problem.L;

B = problem.B; % 3*N x K matrix of b_i(k)
y = problem.y; % 3*N x L matrix of y_i(t_l)
dt = problem.dt;

% Weights
% TODO: scale weights/covars by noisebound?
W = problem.covar_measure.^(-1); % N x L matrix of w_il
lambda = problem.lambda; % scalar
wp = problem.covar_position.^(-1); % L-1 vector
Wp = diag(reshape(repmat(wp',3,1),3*(L-1),1));
wv = problem.covar_velocity.^(-1); % L-2 vector
Wv = diag(reshape(repmat(wv',3,1),3*(L-2),1));
wr = problem.kappa_rotation; % L-1 vector
wd = problem.kappa_rotrate;  % L-1 vector

pBound = problem.translationBound;
vBound = problem.velocityBound;

% check lambda
if ((K >= 3*N) && (lambda == 0.0))
    error("lambda must be positive when there are more shapes (K) " + ...
        "than measurement points (3*N)");
end

%% Define objective
% optimization vector
d = 9*(3*L - 2) + 3*L + 3*(L-1); % 3L - 2 rotations, 3L rotated positions, 3L-1 time-varying rotated positions
% 3L - 2 rotations: L rotations, L-1 delta rotations, L-1 redundant rotations
x = msspoly('x',d);

% pull out individual variables
r  = x(1:(9*L));
dr = x((9*L + 1):(9*L + 9*(L-1)));
rh = x((18*L - 9 + 1):(18*L - 9 + 9*(L-1)));
s  = x((27*L - 18 + 1):(27*L - 18 + 3*L));
sh = x((30*L - 18 + 1):(30*L - 18 + 3*(L-1)));

% convert to useful form
R  = reshape(r ,3,3*L)';
dR = reshape(dr,3,3*(L-1))';
Rh = reshape(rh,3,3*(L-1))';
for l = 1:L
    R(ib3(l),:) =  R(ib3(l),:)';
    if (l < L)
       dR(ib3(l),:) = dR(ib3(l),:)';
        Rh(ib3(l),:) = Rh(ib3(l),:)';
    end
end

% define eliminated variables
% SHAPE
cbar = ones(K,1) / K;
sumh = msspoly(zeros(3*N,1));
sumW = zeros(N,1);
for l = 1:L
    for i = 1:N
        sumh(ib3(i), :) = sumh(ib3(i), :) + W(i,l)*( ...
                          R(ib3(l),:)'*y(ib3(i),l) ...
                          - s(ib3(l)));

        sumW(i) = sumW(i) + W(i,l);
    end
end
sumW = diag(reshape(repmat(sumW,1,3)',3*N,1));
H = [2*((B'*sumW*B) + lambda*eye(K)), ones(K,1);
     ones(1,K), 0];
b = [2*B'*sumh + 2*lambda*cbar; 1];
invH = H \ b;

c = invH(1:(end-1),:);

% VELOCITY
eye3LL = [zeros(3*(L-1),3), eye(3*(L-1))];
eye3LR = [eye(3*(L-1)), zeros(3*(L-1),3)];

eye3Lm3 = eye(3*(L-1));
eye3Lm3L = [zeros(3*(L-2),3), eye(3*(L-2))];
eye3Lm3R = [eye(3*(L-2)), zeros(3*(L-2),3)];

Av = dt*dt*(eye3Lm3'*Wp*eye3Lm3) + ((eye3Lm3L-eye3Lm3R)'*Wv*(eye3Lm3L-eye3Lm3R));

v = Av \ (dt*eye3Lm3'*Wp*(sh - eye3LR*s));

% MAIN OPTIMIZATION
prob_obj = 0;
for l = 1:L
    for i = 1:N
        obj2 = R(ib3(l),:)'*y(ib3(i),l) - B(ib3(i),:)*c - s(ib3(l));
        prob_obj = prob_obj + W(i,l) * (obj2' * obj2);
    end
end
% c regularization
prob_obj = prob_obj + lambda*((c - cbar)'*(c - cbar));
% p regularization
% prob_obj = prob_obj + lambda_p*(p(ib3(1))'*p(ib3(1)) + p(ib3(L))'*p(ib3(L)));
for l = 2:L
    % delta "p"
    delp = sh(ib3(l-1)) - (s(ib3(l-1)) + v(ib3(l-1))*dt);
    prob_obj = prob_obj + wp(l-1)*(delp'*delp);
    % delta v
    if (l < L)
        delv = v(ib3(l)) - v(ib3(l-1));
        prob_obj = prob_obj + wv(l-1)*(delv'*delv);
    end
    % delta R
    delR = reshape(R(ib3(l),:) - Rh(ib3(l-1),:),9,1);
    prob_obj = prob_obj + wr(l-1)*(delR'*delR);
    % dR
    if (l < L)
        deldR = reshape(dR(ib3(l),:) - dR(ib3(l-1),:),9,1);
        prob_obj = prob_obj + wd(l-1)*(deldR'*deldR);
    end
end

%% Define constraints
% EQUALITY
h = [];

% SO(3) constraints
for l = 1:L
    c1 = so3_constraints( R(ib3(l),:));
    if (l < L)
        c2 = so3_constraints(dR(ib3(l),:));
        c3 = so3_constraints(Rh(ib3(l),:));
        h = [h; c1; c2; c3];
    else
        h = [h; c1];
    end
end

% Rh = R dR constraint
for l = 2:L
    h = [h; reshape(Rh(ib3(l-1),:) - R(ib3(l-1),:)*dR(ib3(l-1),:),9,1)];

    % the following constraints do not help
%     h = [h; reshape(Rh(ib3(l-1),:)*dR(ib3(l-1),:)' - R(ib3(l-1),:),9,1)];
%     h = [h; reshape(R(ib3(l-1),:)'*Rh(ib3(l-1),:) - dR(ib3(l-1),:),9,1)];
end

% sh(l) = R'(l-1)*R(l)*s(l) constraint
% sh(l) = dR(l-1)*s(l)
for l = 2:L
    % dR version
    h = [h; sh(ib3(l-1)) - dR(ib3(l-1),:)*s(ib3(l))];
    h = [h; dR(ib3(l-1),:)'*sh(ib3(l-1)) - s(ib3(l))];
    % Long version (adding significantly slows solve time)
    % h = [h; R(ib3(l-1),:)*sh(ib3(l-1)) - R(ib3(l),:)*s(ib3(l))];
end

% INEQUALITY
% p,s in range for just first time (p'*p<=pBoundSq)
% TODO: enforce in range for ALL time steps - LC: agreed
% TODO: add in sh
pBoundSq = pBound^2;
g_s_first = pBoundSq*L - s(ib3(1))'*s(ib3(1));

% v bound (v'*v<=vBoundSq)
vBoundSq = vBound^2;
g_v = vBoundSq*L - v'*v;

% c bound (0<=c<=1)
cBoundSq = 1.0; % should just be 1
g_c = [cBoundSq - c'*c;c];

% g = [g_c];
% g = [g_c;g_v];
g = [g_s_first; g_v; g_c];

%% Complete problem definition
problem.vars = x;
problem.objective = prob_obj;
problem.equality = h; % equality
problem.inequality = g; % inequality

%% Relax!
kappa = 1; % relaxation order
[SDP,info] = dense_sdp_relax(problem,kappa);
% [SDP,info] = simple_sdp_relax(problem,kappa);

%% Solve using MOSEK
tic
prob = convert_sedumi2mosek(SDP.sedumi.At,...
                            SDP.sedumi.b,...
                            SDP.sedumi.c,...
                            SDP.sedumi.K);
addpath(genpath(mosekpath))
[~,res] = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
rmpath(genpath(mosekpath))
soln.solvetime = toc;

% figure; bar(eig(Xopt{1})); % if rank = 1, then relaxation is exact/tight

%% Compare to ground truth
% clip to first eigenvalue
[eigvecs, ~] = eig(Xopt{1});
vecmax = eigvecs(:,end);
% re-normalize so first element is 1
vecmax_normed = vecmax / vecmax(1);
x_est = vecmax_normed(2:end);

% Project to SO(3) and extract results
rs = x_est(1:(9*L));
Rs = projectRList(rs);
drs = x_est((9*L+1):(18*L-9));
dRs = projectRList(drs);
rhs = x_est((18*L-9+1):(27*L-18));
Rhs = projectRList(rhs);

s_est = reshape(full(dmsubs(s,x,x_est)),[3,1,L]);
v_est_raw = reshape(full(dmsubs(v,x,x_est)),[3,1,L-1]);
c_est = full(dmsubs(c,x,x_est));
sh_est_raw = reshape(full(dmsubs(sh,x,x_est)),[3,1,L-1]);

% estimate p from s
p_est = zeros(3,1,L);
for l = 1:L
    p_est(:,:,l) = Rs(:,:,l)*s_est(:,:,l);
end

% estimate sh from s
sh_est = zeros(3,1,L-1);
for l = 2:L
    % sh_est(:,:,l-1) = Rs(:,:,l-1)'*Rs(:,:,l)*s_est(:,:,l);
    sh_est(:,:,l-1) = dRs(:,:,l-1)*s_est(:,:,l);
end

% suboptimality gap
x_proj = [];
for l = 1:L
    r_temp = reshape(Rs(:,:,l),9,1);
    x_proj = [x_proj; r_temp];
end
for l = 1:L-1
    r_temp = reshape(dRs(:,:,l),9,1);
    x_proj = [x_proj; r_temp];
end
for l = 1:L-1
    % r_temp = reshape(Rhs(:,:,l),9,1);
    r_temp = reshape(Rs(:,:,l)*dRs(:,:,l),9,1);
    x_proj = [x_proj; r_temp];
end
x_proj = [x_proj; reshape(s_est,[3*L,1,1]); reshape(sh_est,[3*L-3,1,1])];

% compute gap
obj_est = dmsubs(prob_obj,x,x_proj);
gap = (obj_est - obj(1)) / obj_est;

% estimate v from projected version
v_est = reshape(full(dmsubs(v,x,x_proj)),[3,1,L-1]);

%% Pack into struct
% raw SDP/MOSEK data
soln.raw.Xopt = Xopt;
soln.raw.yopt = yopt;
soln.raw.Sopt = Sopt;
soln.raw.obj = obj;
soln.raw.relax_info = info;

% save estimates
soln.x_est = x_est;

soln.c_est = c_est;
soln.p_est = p_est;
soln.v_est = v_est;
soln.v_est_raw = v_est_raw;
soln.s_est = s_est;
soln.sh_est = sh_est;
soln.sh_est_raw = sh_est_raw;

soln.R_est = Rs;
soln.dR_est = dRs;
soln.Rh_est = Rhs;

soln.gap = gap;
soln.x_proj = x_proj;
soln.obj_est = obj_est;

end