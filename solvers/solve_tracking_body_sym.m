function soln = solve_tracking_body_sym(problem)
% Solves const. vel. (body frame) optimization exactly via SDP
%   Assume the *body frame* velocity is constant and object is spinning.
%   Result: spiral trajectory.
%   Analytically remove velocity & shape. SDP variables are
%   * rotated position (s)
%   * body velocity (v)
%   * rotation (R)
%   * rotation change (dR)
%   VERSION WITH RH AND P MOVED INTO CONSTRAINTS, FULLY SYMBOLIC
%
% INPUTS:
% - problem (struct): populated problem data
%
% RETURNS:
% - soln (struct): solution data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Process inputs
% mosekpath = problem.mosekpath;

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
wv = problem.covar_velocity.^(-1); % L-2 vector
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
d = 9*(2*L - 1) + 3*L + 3*(L-1); % 2L - 1 rotations, 3L rotated positions, 3L-1 body velocities
% 2L - 1 rotations: L rotations, L-1 delta rotations
x = msspoly('x',d);

% pull out individual variables
r  = x(1:(9*L));
dr = x((9*L + 1):(9*L + 9*(L-1)));
s  = x((18*L - 9 + 1):(18*L - 9 + 3*L));
v = x((21*L - 9 + 1):(21*L - 9 + 3*(L-1)));

% convert to useful form
R  = reshape(r ,3,3*L)';
dR = reshape(dr,3,3*(L-1))';
for l = 1:L
    R(ib3(l),:) =  R(ib3(l),:)';
    if (l < L)
       dR(ib3(l),:) = dR(ib3(l),:)';
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
% s regularization (DO NOT USE)
prob_obj = prob_obj + 0.001*(s(ib3(3))'*s(ib3(3)));
% for l = 2:L-1
%     % delta v
%     delv = v(ib3(l)) - v(ib3(l-1));
%     prob_obj = prob_obj + wv(l-1)*(delv'*delv);
%     % dR
%     deldR = reshape(dR(ib3(l),:) - dR(ib3(l-1),:),9,1);
%     prob_obj = prob_obj + wd(l-1)*(deldR'*deldR);
% end

%% Define constraints
% EQUALITY
h = [];

% SO(3) constraints
for l = 1:L
    c1 = so3_constraints( R(ib3(l),:));
    if (l < L)
        c2 = so3_constraints(dR(ib3(l),:));
        h = [h; c1; c2];
    else
        h = [h; c1];
    end
end

% R(t+1) = R(t) dR(t) constraint
for l = 2:L
    h = [h; reshape(R(ib3(l),:) - R(ib3(l-1),:)*dR(ib3(l-1),:),9,1)];
end

% sh(l) = s(l-1) + v(l-1)*dt constraint
for l = 2:L
    % dR version
    h = [h; dR(ib3(l-1),:)*s(ib3(l)) - s(ib3(l-1)) - v(ib3(l-1))*dt];
    h = [h; s(ib3(l)) - dR(ib3(l-1),:)'*s(ib3(l-1)) - dR(ib3(l-1),:)'*v(ib3(l-1))*dt];
end

% TEMP: R, v constraints
% for l = 2:L-1
%     % delta v
%     delv = v(ib3(l)) - v(ib3(l-1));
%     h = [h; wv(l-1)*(delv'*delv)];
%     % dR
%     deldR = reshape(dR(ib3(l),:) - dR(ib3(l-1),:),9,1);
%     h = [h; wd(l-1)*(deldR'*deldR)];
% end

% constraint on v(t1) as a function of v(t2)
% TODO: this may help solve time?
h = [h; dR(ib3(1),:)'*v(ib3(1))*dt - dR(ib3(2),:)*s(ib3(3)) + dR(ib3(1),:)'*s(ib3(1)) + v(ib3(2))*dt];

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
% addpath(genpath(mosekpath))
param = struct();
param.MSK_IPAR_LOG = 0;
[~,res] = mosekopt('minimize info',prob,param);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
% rmpath(genpath(mosekpath))
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

s_est = reshape(full(dmsubs(s,x,x_est)),[3,1,L]);
v_est = reshape(full(dmsubs(v,x,x_est)),[3,1,L-1]);
c_est = full(dmsubs(c,x,x_est));

% estimate p from s
p_est = zeros(3,1,L);
for l = 1:L
    p_est(:,:,l) = Rs(:,:,l)*s_est(:,:,l);
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
x_proj = [x_proj; reshape(s_est,[3*L,1,1]); reshape(v_est,[3*L-3,1,1])];

% compute gap (THIS TAKES FOREVER)
% obj_est = dmsubs(prob_obj,x,x_proj);
% gap = (obj_est - obj(1)) / obj_est;

% compute residuals (TODO: USE WHOLE OBJECTIVE?)
residuals = zeros(N, L);
for i = 1:N
    for l = 1:L
        residue = Rs(:,:,l)'*y(ib3(i),l) - B(ib3(i),:)*c_est - s_est(:,:,l);
        residuals(i,l) = residue'*residue;
    end
end
% residuals = residuals / (problem.noiseSigmaSqrt.^2);

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
soln.s_est = s_est;

soln.R_est = Rs;
soln.dR_est = dRs;

% soln.gap = gap;
% soln.x_proj = x_proj;
% soln.obj_est = obj_est;

soln.residuals = residuals;

end