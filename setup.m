%% Setup for Tracking
% adds all necessary paths

clc; clear; close all; restoredefaultpath

%% Select experiment
% see 'experiments' folder
experiment = "ycbineoat";

%% Change paths here
certifiablyrobustperceptionpath = "../CertifiablyRobustPerception";
gncpath = "../GNC-and-ADAPT"; % optional if no outliers
mosekpath   = '../mosek/10.1/toolbox/r2017a';
sdpnalpath  = '../SDPNALv1.0';

if (experiment == "pascal")
    cadpath = "../datasets/pascal3d";
end

%% add external paths
spotpath    = certifiablyrobustperceptionpath + '/spotless';
stridepath  = certifiablyrobustperceptionpath + '/STRIDE';
manoptpath  = certifiablyrobustperceptionpath + '/manopt';
addpath(certifiablyrobustperceptionpath + '/utils')
addpath(genpath(spotpath)) % Use spotless for defining polynomials
addpath(certifiablyrobustperceptionpath + '/SDPRelaxations') % implementations for SDP relaxation
addpath(genpath(certifiablyrobustperceptionpath + '/CategoryRegistration')) % implementations for SDP relaxation
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

%% Setup experiments
addpath("experiments/"+experiment);

if (experiment == "pascal")
    addpath(cadpath);
end