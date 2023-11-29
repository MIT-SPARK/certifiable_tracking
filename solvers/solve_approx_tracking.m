function soln = solve_approx_tracking(problem,lambda)
% Solves const. vel. optimization with small ang. vel. assumption via SDP.
%   Analytically remove velocity, shape, & position.
%   SDP variables are rotation (R), rotation change (dR), and predicted
%   rotation (Rh) for each time step.
%   APPROXIMATION: use body frame position and assume R(t) ~ R(t+1) for
%   velocity term only.
%   TODO: THERE MAY STILL BE A BUG HERE!
%
% INPUTS:
% - problem (struct): populated problem data
% - lambda (float): value of lambda to use
%
% RETURNS:
% - soln (struct): solution data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Process inputs
mosekpath = problem.mosekpath;

N = problem.N;
K = problem.K;
L = problem.L;

B = problem.B; % 3*N x K matrix of b_i(k)
y = problem.y; % 3*N x L matrix of y_i(t_l)
dt = problem.dt;

pBound = problem.translationBound;
vBound = problem.velocityBound;

% check lambda
if ((K >= 3*N) && (lambda == 0.0))
    error("lambda must be positive when there are more shapes (K) " + ...
        "than measurement points (3*N)");
end

%% Define objective
% optimization vector
d = 9*(3*L - 1); % 3L - 1 rotations
x = msspoly('x',d);

% pull out individual variables
r  = x(1:(9*L));
dr = x((9*L + 1):(9*L + 9*L));
rh = x((18*L + 1):(18*L + 9*(L-1)));

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
eye3LL = [zeros(3*(L-1),3), eye(3*(L-1))];
eye3LR = [eye(3*(L-1)), zeros(3*(L-1),3)];

% VELOCITY COEFFICIENTS
Av = dt*dt*(eye3LR'*eye3LR) + ((eye3LL-eye3LR)'*(eye3LL-eye3LR));
Cv = Av \ (dt*eye3LR'*(eye3LL-eye3LR));
% v = Cv * p

% POSITION COEFFICIENTS
Ap = eye(3*L) + (Cv'*(eye3LL-eye3LR)'*(eye3LL-eye3LR)*Cv) + ...
     ((eye3LL - eye3LR - eye3LR*Cv*dt)'*(eye3LL - eye3LR - eye3LR*Cv*dt));
% p = Ap^{-1} * sum_i(h_i)

% SHAPE
Bsum = zeros(3,K);
for i = 1:N
    Bsum = Bsum + B(ib3(i),:);
end
Bsum = repmat(Bsum,L,1);

ySum = msspoly(zeros(3*L,1));
for l = 1:L
    for i = 1:N
        ySum(ib3(l)) = ySum(ib3(l)) + ...
            R(ib3(l),:)' * y(ib3(i),l);
    end
end

pbar = Ap \ ySum;
sumh = msspoly(zeros(3*N,1));
for l = 1:L
    for i = 1:N
        sumh(ib3(i), :) = sumh(ib3(i), :) + ...
                          R(ib3(l),:)'*y(ib3(i),l) + ...
                         -pbar(ib3(l));
    end
end

Ccp = Ap \ Bsum;
sumCcp = zeros(3,K);
for l = 1:L
    sumCcp = sumCcp + Ccp(ib3(l),:);
end
sumCcp = repmat(sumCcp,N,1);
Bbar = L*B - sumCcp;

H = [2*((Bbar'*Bbar) + lambda*eye(K)), ones(K,1);
     ones(1,K), 0];
b = [2*Bbar'*sumh; 1.0];
invH = H \ b;

c = invH(1:(end-1),:);

% FINALIZE POSITION
p = Ap \ (ySum - Bsum*c);

% FINALIZE VELOCITY
v = Cv * p;

% MAIN OPTIMIZATION
prob_obj = 0;
for l = 1:L
    for i = 1:N
        obj = R(ib3(l),:)'*y(ib3(i),l) - B(ib3(i),:)*c - p(ib3(l));
        prob_obj = prob_obj + obj' * obj;
    end
end
prob_obj = prob_obj + lambda*(c'*c);
for l = 2:L
    % delta p
    delp = p(ib3(l)) - (p(ib3(l-1)) + v(ib3(l-1))*dt);
    prob_obj = prob_obj + delp'*delp;
    % delta v
    delv = v(ib3(l)) - v(ib3(l-1));
    prob_obj = prob_obj + delv'*delv;
    % delta R
    delR = reshape(R(ib3(l),:) - Rh(ib3(l-1),:),9,1);
    prob_obj = prob_obj + delR'*delR;
    % dR
    deldR = reshape(dR(ib3(l),:) - dR(ib3(l-1),:),9,1);
    prob_obj = prob_obj + deldR'*deldR;
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
        h = [c1; c2];
    end
end

% Rh = R dR constraint
for l = 2:L
    h = [h; reshape(Rh(ib3(l-1),:) - R(ib3(l),:)*dR(ib3(l),:),9,1)];
end

% INEQUALITY
% p in range for just first time (p'*p<=pBoundSq)
% TODO: enforce in range for ALL time steps
pBoundSq = pBound^2;
g_p_first = pBoundSq*L - p(ib3(1))'*p(ib3(1));

% v bound (v'*v<=vBoundSq)
vBoundSq = vBound^2;
g_v = vBoundSq*L - v'*v;

% c bound (0<=c<=1)
cBoundSq = 1.0; % should just be 1
g_c = [cBoundSq - c'*c;c];

% g = [];
% g = [g_c;g_v];
g = [g_p_first; g_v; g_c];

%% Complete problem definition
problem.vars = x;
problem.objective = prob_obj;
problem.equality = h; % equality
problem.inequality = g; % inequsality

%% Relax!
kappa = 1; % relaxation order
[SDP,info] = dense_sdp_relax(problem,kappa);

%% Solve using MOSEK
% load("sdpvars.mat");
prob = convert_sedumi2mosek(SDP.sedumi.At,...
                            SDP.sedumi.b,...
                            SDP.sedumi.c,...
                            SDP.sedumi.K);
addpath(genpath(mosekpath))
[~,res] = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
rmpath(genpath(mosekpath))

% figure; bar(eig(Xopt{1})); % if rank = 1, then relaxation is exact/tight

%% Compare to ground truth
% clip to first eigenvalue
[eigvecs, ~] = eig(Xopt{1});
vecmax = eigvecs(:,end);
% re-normalize so first element is 1
vecmax_normed = vecmax / vecmax(1);
x_est = vecmax_normed(2:end);

% generate a "gt" vector
x_gt = zeros(d,1);

R_cur = problem.R_gt;
for l = 1:L
    x_gt((9*(l-1)+1):(9*(l-1)+9)) = reshape(R_cur,9,1);
    x_gt((9*(L + l-1)+1):(9*(L + l-1)+9)) = reshape(problem.dR_gt,9,1);
    if (l < L)
        x_gt((9*(2*L + l-1)+1):(9*(2*L + l-1)+9)) = reshape(R_cur * problem.dR_gt,9,1);
    end
    R_cur = R_cur * problem.dR_gt;
end

% Technically I should project to SO(3), but this is quick error measure
x_err = norm(x_gt - x_est);

% TODOS:
% 1) project R
% 2) better p estimate
% 3) angular R error

% FOR TESTING
% x_gt = x_gt(2:end);
% c_gt = dmsubs(c,x,x_gt);
% p_gt = dmsubs(p,x,x_gt);
% v_gt = dmsubs(v,x,x_gt);

%% Pack into struct
% raw SDP/MOSEK data
soln.raw.Xopt = Xopt;
soln.raw.yopt = yopt;
soln.raw.Sopt = Sopt;
soln.raw.obj = obj;
soln.raw.relax_info = info;

% convert to c, p, s, v and save
soln.x_est = x_est;
soln.c_est = dmsubs(c,x,x_est);
soln.p_est = dmsubs(p,x,x_est);
soln.v_est = dmsubs(v,x,x_est);

% errors
soln.x_err = x_err;

end

