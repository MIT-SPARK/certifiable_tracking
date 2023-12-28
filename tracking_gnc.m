%% Dense SDP relaxation for certifiable tracking
%  Version with outlier rejection through GNC
%  TODO: THIS NEEDS TO BE UPDATED / REMOVED
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all; restoredefaultpath
% rng("default")

%% dependencies
% Change paths here
certifiablyrobustperceptionpath = "../CertifiablyRobustPerception";
gncpath = "../GNC-and-ADAPT";
mosekpath   = 'C:/Program Files/Mosek/10.1/toolbox/r2017a';
sdpnalpath  = '../SDPNALv1.0';

% add external paths
spotpath    = certifiablyrobustperceptionpath + '/spotless';
stridepath  = certifiablyrobustperceptionpath + '/STRIDE';
manoptpath  = certifiablyrobustperceptionpath + '/manopt';
addpath(certifiablyrobustperceptionpath + '/utils')
addpath(genpath(spotpath)) % Use spotless for defining polynomials
addpath(certifiablyrobustperceptionpath + '/SDPRelaxations') % implementations for SDP relaxation
addpath(genpath(gncpath))

% add internal paths
addpath('./solvers')
addpath('./visualization')
addpath('./utils')

path.stridepath = stridepath;
path.mosekpath  = mosekpath;
path.manoptpath = manoptpath;

%% Generate random tracking problem
problem.N_VAR = 11; % nr of keypoints
problem.K = 3; % nr of shapes
problem.L = 10; % nr of keyframes in horizon

problem.outlierRatio = 0.1;
problem.noiseSigmaSqrt = 0.01; % [m]
problem.intraRadius = 0.2;
problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity
% problem.velprior = "world";      % constant world frame velocity
% problem.velprior = "grav-world"; % add gravity in z direction

problem.accelerationNoiseBoundSqrt = 0.0;%0.01;
problem.rotationNoiseBound = 0;%pi/32; % rad

problem.N = problem.N_VAR*problem.L; % How many measurements this problem has
problem.outliers = []; % outlier indicies
problem.priors = [];
problem.dof = 3;

% Optional: use a specified velocity trajectory
% problem = make_trajectory(problem);

% add shape, measurements, outliers
problem = gen_random_tracking(problem);
lambda = 0.0;
problem.lambda = lambda;

problem.mosekpath = mosekpath;

%% Solve!
epsilon = chi2inv(0.99, problem.dof)*problem.noiseSigmaSqrt;
[inliers, info] = gnc(problem, @solver_for_gnc, 'NoiseBound', epsilon,'MaxIterations',100,'Debug',false);


%% Check solutions
if isequal(problem.inliers_gt,inliers)
    disp("Correct inliers found after " + string(info.Iterations) + " iterations.");
else
    disp("Inliers not found after " + string(info.Iterations) + " iterations.");
end

% play done sound
load handel.mat
sound(y,2*Fs);
