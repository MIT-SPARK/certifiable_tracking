%% Dense SDP relaxation for certifiable tracking
% Uses exact optimization form.
% 
% Replaced by dense_tracking.m + solvers/solve_full_tracking.m
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all; restoredefaultpath
rng("default")

%% dependencies
spotpath    = '../spotless';
mosekpath   = 'C:/Program Files/Mosek/10.1/toolbox/r2017a';
stridepath  = '../STRIDE';
manoptpath  = '../manopt';
sdpnalpath  = '../../SDPNALv1.0';
path.stridepath = stridepath;
path.mosekpath  = mosekpath;
path.manoptpath = manoptpath;
addpath('../utils')
addpath('./solvers')
addpath(genpath('../spotless')) % Use spotless for defining polynomials
addpath('../SDPRelaxations') % implementations for SDP relaxation

%% Generate random tracking problem
problem.N_VAR = 10;
problem.K = 3;
problem.L = 9;

problem.outlierRatio = 0.0; % TODO: no support for outliers
problem.noiseSigma = 0.01;
problem.intraRadius = 0.2;
problem.translationBound = 20.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

% ground truth for testing
% problem.c_gt = zeros(problem.K,1);
% problem.c_gt(1) = 1.0; % this is actually bad with noise
% problem.p_gt = zeros(3,1);
% problem.v_gt = zeros(3,1);
% problem.R_gt = eye(3);
% problem.dR_gt = eye(3);

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 0.0;


%% Define objective
% constants shorthand
N = problem.N_VAR;
K = problem.K;
L = problem.L;

% optimization vector
d = 9*(3*L - 1) + 3*L + 3*L; % 3L - 1 rotations, 3L rotated positions, 3L positions
x = msspoly('x',d);

% pull out individual variables
r  = x(1:(9*L));
dr = x((9*L + 1):(9*L + 9*L));
rh = x((18*L + 1):(18*L + 9*(L-1)));
p = x((18*L + 9*(L-1) + 1):(27*L - 9 + 3*L));
s = x((30*L - 9 + 1):(30*L - 9 + 3*L));

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

% define auxillary variables
B = problem.B; % 3*N x K matrix of b_i(k)
y = problem.y; % 3*N x L matrix of y_i(t_l)
dt = problem.dt;

% define eliminated variables
% SHAPE (VERIFIED)
sumh = msspoly(zeros(3*N,1));
for l = 1:L
    for i = 1:N
        sumh(ib3(i), :) = sumh(ib3(i), :) + ...
                          R(ib3(l),:)'*y(ib3(i),l) ...
                          - s(ib3(l));
    end
end
H = [2*(L*(B'*B) + lambda*eye(K)), ones(K,1);
     ones(1,K), 0];
b = [2*B'*sumh; 1];
invH = H \ b;

c = invH(1:(end-1),:);

% VELOCITY (VERIFIED)
eye3LL = [zeros(3*(L-1),3), eye(3*(L-1))];
eye3LR = [eye(3*(L-1)), zeros(3*(L-1),3)];
Av = dt*dt*(eye3LR'*eye3LR) + ((eye3LL-eye3LR)'*(eye3LL-eye3LR));

v = Av \ (dt*eye3LR'*(eye3LL-eye3LR) * p);

% MAIN OPTIMIZATION
prob_obj = 0;
for l = 1:L
    for i = 1:N
        obj2 = R(ib3(l),:)'*y(ib3(i),l) - B(ib3(i),:)*c - s(ib3(l));
        prob_obj = prob_obj + obj2' * obj2;
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
    c1 = rot_mat_constraints( R(ib3(l),:));
    c2 = rot_mat_constraints(dR(ib3(l),:));
    if (l < L)
        c3 = rot_mat_constraints(Rh(ib3(l),:));
        h = [h; c1; c2; c3];
    else
        h = [h; c1; c2];
    end
end

% Rh = R dR constraint
for l = 2:L
    h = [h; reshape(Rh(ib3(l-1),:) - R(ib3(l-1),:)*dR(ib3(l-1),:),9,1)];

    % the following constraints do not help
    % h = [h; reshape(Rh(ib3(l-1),:)*dR(ib3(l-1),:)' - R(ib3(l-1),:),9,1)];
    % h = [h; reshape(R(ib3(l-1),:)'*Rh(ib3(l-1),:) - dR(ib3(l-1),:),9,1)];
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
pBound = problem.translationBound;
pBoundSq = pBound^2;
g_p_first = pBoundSq*L - p(ib3(1))'*p(ib3(1));
g_s_first = pBoundSq*L - s(ib3(1))'*s(ib3(1));

% v bound (v'*v<=vBoundSq)
vBound = problem.velocityBound;
vBoundSq = vBound^2;
g_v = vBoundSq*L - v'*v;

% c bound (0<=c<=1)
cBoundSq = 1.0; % should just be 1
g_c = [cBoundSq - c'*c;c];

% g = [];
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

figure; bar(eig(Xopt{1})); % if rank = 1, then relaxation is exact/tight

%% Check answer
% TODO: CLEAN UP
% clip to first eigenvalue
[eigvecs, eigvals] = eig(Xopt{1});
eigmax = eigvals(end,end);
vecmax = eigvecs(:,end);
% re-normalize so first element is 1
vecmax_normed = vecmax / vecmax(1);

% generate a "gt" vector
x_gt = zeros(d,1);

p_cur = problem.p_gt;
R_cur = problem.R_gt;
for l = 1:L
    x_gt((9*(l-1)+1):(9*(l-1)+9)) = reshape(R_cur,9,1);
    x_gt((9*(L + l-1)+1):(9*(L + l-1)+9)) = reshape(problem.dR_gt,9,1);
    if (l < L)
        x_gt((9*(2*L + l-1)+1):(9*(2*L + l-1)+9)) = reshape(R_cur * problem.dR_gt,9,1);
    end

    x_gt((9*(3*L-1) + 3*(l-1) + 1):(9*(3*L-1) + 3*(l-1) + 3)) = p_cur;
    x_gt((9*(3*L-1) + 3*(L + l-1) + 1):(9*(3*L-1) + 3*(L + l-1) + 3)) = R_cur' * p_cur;

    p_cur = p_cur + problem.v_gt * dt;
    R_cur = R_cur * problem.dR_gt;
end
x_gt = [1.0; x_gt];

% Technically I should project to SO(3), but this is quick error measure
x_err = norm(x_gt - vecmax_normed)

% error excluding p
x_gt_no_p = [x_gt(1:(1 + 9*(3*L-1) + 1)); x_gt((1 + 9*(3*L-1) + 3*L + 1):end)];
vecmax_normed_no_p = [vecmax_normed(1:(1 + 9*(3*L-1) + 1)); vecmax_normed((1 + 9*(3*L-1) + 3*L + 1):end)];
x_err_no_p = norm(x_gt_no_p - vecmax_normed_no_p)

% also convert to c, p, v
x_est = vecmax_normed(2:end);
c_est = dmsubs(c,x,x_est);
p_est = dmsubs(p,x,x_est);
s_est = dmsubs(s,x,x_est);
v_est = dmsubs(v,x,x_est);

% FOR TESTING
% x_gt = x_gt(2:end);
% c_gt = dmsubs(c,x,x_gt); % will not match gt with noise
% p_gt = dmsubs(p,x,x_gt);
% v_gt = dmsubs(v,x,x_gt);

% can also test objective->0 in noise-free case
% obj_est = dmsubs(prob_obj,x,x_est)
% obj_gt = dmsubs(prob_obj,x,x_gt)


%% helper function for indexing
function ran = ib3(i)
% index every three elements

ran = (3*(i-1) + 1):(3*(i-1) + 3);

end

%% helper function for rotation matrix constraints
function c = rot_mat_constraints(R)

% R
col1 = R(:,1);
col2 = R(:,2);
col3 = R(:,3);

% cols unit length
unitlen = [1.0 - col1'*col1; ...
           1.0 - col2'*col2; ...
           1.0 - col3'*col3];
% cols orthogonal
orth = [col1'*col2;...
        col2'*col3;...
        col3'*col1];
% cols righthanded
righthand = [cross(col1,col2) - col3;...
             cross(col2,col3) - col1;...
             cross(col3,col1) - col2];
% add to matrix!
c = [unitlen; orth; righthand];

end