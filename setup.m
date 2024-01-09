%% Setup for Tracking
% adds all necessary paths

clc; clear; close all; restoredefaultpath

%% Select experiment
experiment = "synthetic";
addpath("experiments/"+experiment);

%% Change paths here
certifiablyrobustperceptionpath = "../CertifiablyRobustPerception";
gncpath = "../GNC-and-ADAPT"; % optional if no outliers
mosekpath   = '../mosek/10.1/toolbox/r2017a';
sdpnalpath  = '../SDPNALv1.0';

%% add external paths
spotpath    = certifiablyrobustperceptionpath + '/spotless';
stridepath  = certifiablyrobustperceptionpath + '/STRIDE';
manoptpath  = certifiablyrobustperceptionpath + '/manopt';
addpath(certifiablyrobustperceptionpath + '/utils')
addpath(genpath(spotpath)) % Use spotless for defining polynomials
addpath(certifiablyrobustperceptionpath + '/SDPRelaxations') % implementations for SDP relaxation
addpath(genpath(gncpath)) % optional if no outliers

addpath(genpath(mosekpath))
addpath(genpath(stridepath))

%% add internal paths
addpath('./outlier_rejection')
addpath('./solvers')
addpath('./utils')
addpath('./visualization')

%% Setup for ROBIN
flag = int32(bitor(2, 8));
py.sys.setdlopenflags(flag);