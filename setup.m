%% Setup for Tracking
% adds all necessary paths

clc; clear; close all; restoredefaultpath

% Change paths here
certifiablyrobustperceptionpath = "../CertifiablyRobustPerception";
mosekpath   = 'C:/Program Files/Mosek/10.1/toolbox/r2017a';
sdpnalpath  = '../SDPNALv1.0';

% add external paths
spotpath    = certifiablyrobustperceptionpath + '/spotless';
stridepath  = certifiablyrobustperceptionpath + '/STRIDE';
manoptpath  = certifiablyrobustperceptionpath + '/manopt';
addpath(certifiablyrobustperceptionpath + '/utils')
addpath(genpath(spotpath)) % Use spotless for defining polynomials
addpath(certifiablyrobustperceptionpath + '/SDPRelaxations') % implementations for SDP relaxation

addpath(genpath(mosekpath))
addpath(genpath(stridepath))

% add internal paths
addpath('./solvers')
addpath('./visualization')
addpath('./utils')