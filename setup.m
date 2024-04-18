%% Setup for Tracking
% adds all necessary paths

clc; clear; close all; restoredefaultpath

%% Select experiment
% see 'experiments' folder
experiment = "pascal";

%% Change paths here
certifiablyrobustperceptionpath = "../CertifiablyRobustPerception";
gncpath = "../GNC-and-ADAPT"; % optional if no outliers
mosekpath   = '/opt/mosek/10.1/toolbox/r2017a';
sdpnalpath  = '../SDPNALv1.0';
coptpath = '/opt/copt71';
pyvenvpath = '~/research/tracking/trackvenv/bin/python3';

if (experiment == "pascal")
    cadpath = "../datasets/pascal3d";
end
if (experiment == "racecar_offline")
    cadpath = "../datasets/racecar_offline";
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
% addpath("/opt/copt71")
setenv("COPT_LICENSE_DIR", "~/copt")

%% add internal paths
addpath('./outlier_rejection')
addpath('./solvers')
addpath('./utils')
addpath('./visualization')

%% Setup for python
flag = int32(bitor(2, 8));
py.sys.setdlopenflags(flag);
pyenv('Version', pyvenvpath, 'ExecutionMode','OutOfProcess');

%% Setup experiments
addpath("experiments/"+experiment);

if (experiment == "pascal") || (experiment == "racecar_offline")
    addpath(cadpath);
end