%% Testing dense SDP relaxation for certifiable tracking
% Category registration with dense relaxation.
% PACE, but only analytically eliminates c.
% s = R'p included in SDP.
%
% Used mostly for debugging
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all; restoredefaultpath

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
problem.L = 1;

problem.outlierRatio = 0.0; % TODO: no support for outliers
problem.noiseSigma = 0.02;
problem.intraRadius = 0.2;
problem.translationBound = 10.0;
problem.velocityBound = 0.0;
problem.dt = 1.0;

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 1.0;


%% Define objective
% constants shorthand
N = problem.N_VAR;
K = problem.K;

% optimization vector
d = 9 + 3; % 1 rotation, 1 position
x = msspoly('x',d);

% pull out individual variables
% one = x(1);
r  = x(1:(0 + 9));
s = x((9 + 1):(9 + 3));

% convert to useful form
R  = reshape(r ,3,3);

% define auxillary variables
B = problem.B; % 3*N x K matrix of b_i(k)
y = problem.y; % 3*N matrix of y_i

% define eliminated variables
% SHAPE
sumh = msspoly(zeros(3*N,1));
for i = 1:N
    sumh(ib3(i), :) = R'*y(ib3(i)) - s;
end
H = [2*(B'*B + lambda*eye(K)), ones(K,1);
     ones(1,K), 0];
b = [2*B'*sumh; 1.0];
invH = H \ b;

c = invH(1:(end-1),:);

% MAIN OPTIMIZATION
prob_obj = 0;
for i = 1:N
    obj = R'*y(ib3(i)) - B(ib3(i),:)*c - s;
    prob_obj = prob_obj + obj' * obj;
end
prob_obj = prob_obj + lambda*(c'*c);

%% Define constraints
% EQUALITY
h = [];
% SO(3) constraints
c1 = rot_mat_constraints( R );
h = [h; c1];

% INEQUALITY
% s'*s <= sBoundSq
sBound = problem.translationBound;
sBoundSq = sBound^2;
g_s = sBoundSq - s'*s;

% c bound
cBoundSq = 1.0; % should just be 1
g_c = [cBoundSq - c'*c;c];

% g = [];
g = [g_s];
% g = [g_s; g_c];

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
% TODO: CLEANUP
% NOTE: THIS ONLY WORKS IF k = 1!!
% clip to first eigenvalue
[eigvecs, eigvals] = eig(Xopt{1});
eigmax = eigvals(end,end);
vecmax = eigvecs(:,end);
% re-normalize so first element is 1
vecmax_normed = vecmax / vecmax(1);

% Pull out rotation
R_est = reshape(vecmax_normed(2:(1+9)),3,3);
R_est = project2SO3(R_est);
% Pull out position
s_est = vecmax_normed((10+1):(10+3));
% s = R' p
p_est = R_est * s_est;

err_R = norm(R_est - problem.R_gt)
err_p = norm(p_est - problem.p_gt)


% generate a "gt" vector
x_gt = zeros(d,1);

x_gt(1:9) = reshape(problem.R_gt,9,1);
x_gt(10:12) = problem.R_gt'*problem.p_gt;

x_gt = [1.0; x_gt];

% Technically I should project to SO(3), but this is quick error measure
x_err = norm(x_gt - vecmax_normed)

% also convert to c
x_est = vecmax_normed(2:end);
c_est = dmsubs(c,x,x_est);


%% solve using STRIDE
% TODO: this is not working! (no local search method)
addpath(genpath(pgdpath))

pgdopts.pgdStepSize     = 10;
pgdopts.SDPNALpath      = sdpnalpath;
pgdopts.tolADMM         = 10e-5;
pgdopts.phase1          = 1;
pgdopts.rrOpt           = 1:3;
pgdopts.rrFunName       = 'local_search_tracking'; % local search function
pgdopts.rrPar           = info; % need the original POP formulation for local search (TODO: actually?)
pgdopts.maxiterLBFGS    = 1000;
pgdopts.maxiterSGS      = 300;
pgdopts.tolLBFGS        = 1e-12;
pgdopts.tolPGD          = 1e-8;

[outPGD,Xopt,yopt,Sopt]     = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, [], pgdopts);
time_pgd                    = outPGD.totaltime;
% round solutions and check optimality certificate
res = get_performance_bqp(Xopt,yopt,Sopt,SDP,info,pgdpath);

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