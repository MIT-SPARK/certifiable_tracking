function soln = solve_nopos_tracking(problem)
% Solves const. vel. optimization exactly via SDP.
%   Analytically remove velocity & position & shape. SDP variables are
%   rotated position (s), rotation (R), rotation change (dR), and predicted
%   rotation (Rh) for each time step.
%
%   WARNING: YOU JUST GET QUARTIC
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
% TODO: SCALING MAY BE KEY TO TIGHTNESS (ALSO WHY DOES IT MESS v UP??)
W = problem.weights;% / problem.noiseBoundSq; % N x L matrix of w_il
lambda = problem.lambda; % scalar
Wp = problem.weights_position; % 3*(L-1) vector
Wv = problem.weights_velocity; % 3*(L-1) vector
Wr = problem.weights_rotation; % 3*(L-1) vector
Wd = problem.weights_rotrate;  % 3*(L-1) vector
lambda_p = 1.0; % TODO: FIX
% TODO: SCALE WEIGHTS by noisebound??

pBound = problem.translationBound;
vBound = problem.velocityBound;

% check lambda
if ((K >= 3*N) && (lambda == 0.0))
    error("lambda must be positive when there are more shapes (K) " + ...
        "than measurement points (3*N)");
end

%% Define objective
% optimization vector
d = 9*(3*L - 1) + 3*L; % 3L - 1 rotations, 3L rotated positions, 3L positions
x = msspoly('x',d);

% pull out individual variables
r  = x(1:(9*L));
dr = x((9*L + 1):(9*L + 9*L));
rh = x((18*L + 1):(18*L + 9*(L-1)));
s = x((18*L + 9*(L-1) + 1):(27*L - 9 + 3*L));

% convert to useful form
R  = reshape(r ,3,3*L)';
dR = reshape(dr,3,3*L)';
Rh = reshape(rh,3,3*(L-1))';
for l = 1:L
    R(ib3(l),:) =  R(ib3(l),:)';
   dR(ib3(l),:) = dR(ib3(l),:)';
    if (l < L)
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
Av = dt*dt*(eye3LR'*diag(Wp)*eye3LR) + ((eye3LL-eye3LR)'*diag(Wv)*(eye3LL-eye3LR));

% v = Av \ (dt*eye3LR'*diag(Wp)*(eye3LL-eye3LR) * p);
v_coeff = Av \ (dt*eye3LR'*diag(Wp)*(eye3LL-eye3LR));

% POSITION
Ap = 2*((eye3LL-eye3LR)'*diag(Wp)*(eye3LL-eye3LR)) ...
    -2*dt*(eye3LL-eye3LR)'*diag(Wp)*eye3LR*v_coeff ...
    +0;

soln = Ap;
return;


%     +2*(lambda_p*eye(3*L));
diagR = msspoly(zeros(3*L,3*L));
for l = 1:L
    diagR(ib3(l),ib3(l)) = R(ib3(l),:);
end

% M = [Ap, diagR';
%      diagR, zeros(3*L,3*L)];
% invert M using the Schur complement rule
% Sinv = -diagR*Ap*diagR'; % M/Ap
% ApInv = inv(Ap);
% invM = [eye(3*L), -ApInv*diagR';
%         zeros(3*L,3*L), eye(3*L)] * ...
%        [ApInv, zeros(3*L,3*L);
%         zeros(3*L,3*L), Sinv] * ...
%        [eye(3*L), zeros(3*L,3*L);
%         -diagR*ApInv, eye(3*L)];
% b = [zeros(3*L,1); s];
% p_full = invM * b;
% p = p_full(1:(end-3*L),:);
% p = -ApInv*diagR'*Sinv*s;
% p = ApInv*Ap*diagR'*s;
p = diagR'*s;
% Great.

% FINALIZE VELOCITY
v = v_coeff * p;


% MAIN OPTIMIZATION
prob_obj = 0;
for l = 1:L
    for i = 1:N
        obj2 = R(ib3(l),:)'*y(ib3(i),l) - B(ib3(i),:)*c - s(ib3(l));
        prob_obj = prob_obj + W(i,l) * (obj2' * obj2);
    end
end
prob_obj = prob_obj + lambda*((c - cbar)'*(c - cbar));
for l = 2:L
    % delta p
    delp = p(ib3(l)) - (p(ib3(l-1)) + v(ib3(l-1))*dt);
    prob_obj = prob_obj + Wp(3*(l-2)+1)*(delp'*delp);
    % delta v
    delv = v(ib3(l)) - v(ib3(l-1));
    prob_obj = prob_obj + Wv(3*(l-2)+1)*(delv'*delv);
    % delta R
    delR = reshape(R(ib3(l),:) - Rh(ib3(l-1),:),9,1);
    prob_obj = prob_obj + Wr(3*(l-2)+1)*(delR'*delR);
    % dR
    deldR = reshape(dR(ib3(l),:) - dR(ib3(l-1),:),9,1);
    prob_obj = prob_obj + Wd(3*(l-2)+1)*(deldR'*deldR);
end

%% Define constraints
% EQUALITY
h = [];

% SO(3) constraints
for l = 1:L
    c1 = so3_constraints( R(ib3(l),:));
    c2 = so3_constraints(dR(ib3(l),:));
    if (l < L)
        c3 = so3_constraints(Rh(ib3(l),:));
        h = [h; c1; c2; c3];
    else
        h = [h; c1; c2];
    end
end

% Rh = R dR constraint
for l = 2:L
    h = [h; reshape(Rh(ib3(l-1),:) - R(ib3(l-1),:)*dR(ib3(l-1),:),9,1)];

    % the following constraints do not help
%     h = [h; reshape(Rh(ib3(l-1),:)*dR(ib3(l-1),:)' - R(ib3(l-1),:),9,1)];
%     h = [h; reshape(R(ib3(l-1),:)'*Rh(ib3(l-1),:) - dR(ib3(l-1),:),9,1)];
end

% s = R' p constraint
for l = 1:L
    % explictly saying R' = R^{-1} helps the solution converge
    h = [h; s(ib3(l)) - R(ib3(l),:)'*p(ib3(l))];
    h = [h; R(ib3(l),:)*s(ib3(l)) - p(ib3(l))]; % this constraint more important
end

% INEQUALITY
% p,s in range for just first time (p'*p<=pBoundSq)
% TODO: enforce in range for ALL time steps
pBoundSq = pBound^2;
g_p_first = pBoundSq*L - p(ib3(1))'*p(ib3(1));
g_s_first = pBoundSq*L - s(ib3(1))'*s(ib3(1));

% v bound (v'*v<=vBoundSq)
vBoundSq = vBound^2;
g_v = vBoundSq*L - v'*v;

% c bound (0<=c<=1)
cBoundSq = 1.0; % should just be 1
g_c = [cBoundSq - c'*c;c];

% g = [g_c];
% g = [g_c;g_v];
g = [g_s_first; g_p_first; g_v; g_c];

%% Complete problem definition
problem.vars = x;
problem.objective = prob_obj;
problem.equality = h; % equality
problem.inequality = g; % inequsality

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
drs = x_est((9*L+1):(18*L));
dRs = projectRList(drs);
rhs = x_est((18*L+1):(27*L-9));
Rhs = projectRList(rhs);

s_est = reshape(full(dmsubs(s,x,x_est)),[3,1,L]);
v_est = reshape(full(dmsubs(v,x,x_est)),[3,1,L]);
p_est_raw = reshape(full(dmsubs(p,x,x_est)),[3,1,L]);
c_est = full(dmsubs(c,x,x_est));

% estimate p from s
p_est = zeros(3,1,L);
for l = 1:L
    p_est(:,:,l) = Rs(:,:,l)*s_est(:,:,l);
end

% duality gap
% obj_est = dmsubs(prob_obj,x,x_est);
% obj_gt = dmsubs(prob_obj,x,problem.x_gt);
% gap = abs(obj_est - obj_gt)/(1+abs(obj_est)+abs(obj_gt));

% FOR TESTING
% x_gt = x_gt(2:end);
% c_gt = dmsubs(c,x,x_gt); % will not match gt with noise
% p_gt = dmsubs(p,x,x_gt);
% v_gt = dmsubs(v,x,x_gt);

% can also test objective->0 in noise-free case
% obj_est = dmsubs(prob_obj,x,x_est)
% obj_gt = dmsubs(prob_obj,x,problem.x_gt)

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
soln.p_est_raw = p_est_raw;

soln.R_est = Rs;
soln.dR_est = dRs;
soln.Rh_est = Rhs;

% soln.gap = gap;

end

function Rs = projectRList(r)
% project list of Rs into SO(3) conveniently

N = length(r)/9;
temp_R  = reshape(r ,3, 3*N)';
Rs = zeros(3,3,N);
for i = 1:N
    Rs(:,:,i) = project2SO3(temp_R(ib3(i),:)');
end
end