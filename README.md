# MATLAB Implementation of Certifiable Tracking

## Setup
You can play with a limited set of results without installing dependencies. To run the solver, this repository has the following dependencies:
- [CertifiablyRobustPerception](https://github.com/MIT-SPARK/CertifiablyRobustPerception/tree/master)
- [MOSEK](https://www.mosek.com/)
- [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)
- [GNC](https://github.com/MIT-SPARK/GNC-and-ADAPT) (only for tracking with outliers)

You may set the paths to these dependencies in [setup.m](setup.m).

## Organization
There are two main scripts:
1. [dense_tracking.m](dense_tracking.m) runs outlier-free certifiable tracking. It contains results processing code intended to evaluate the accuracy of our method. This script does not require the GNC dependency.
2. [tracking_gnc.m](tracking_gnc.m) runs certifiable tracking with outliers using GNC.

The subfolders contain helper functions. `solvers` holds functions for our method and comparison methods. `utils` holds a mix of helper functions to generate, process, and solve the SDP. `visualization` holds visualization tools.