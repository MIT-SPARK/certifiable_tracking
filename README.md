# MATLAB Implementation of Certifiable Tracking

## Setup
You can play with a limited set of results without installing dependencies. To run the solver, this repository has the following dependencies:
- [CertifiablyRobustPerception](https://github.com/MIT-SPARK/CertifiablyRobustPerception/tree/master) (you need to run SPOTLESS install script)
- [MOSEK](https://www.mosek.com/)
- [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)
- [GNC](https://github.com/MIT-SPARK/GNC-and-ADAPT) (only for tracking with outliers): TODO: CHANGE TO MY FORK OF GNC
- [ROBIN](https://github.com/MIT-SPARK/ROBIN) (only for robust outlier rejection)

You may set the paths to these dependencies in [setup.m](setup.m).

### ROBIN Setup
ROBIN is a more complicated install. I used Ubuntu 22.04. First, compile ROBIN according to the instructions in the repository. Make sure to turn off the tests in the cmakelists file or it may fail.

The following may not be necessary but worked for me. Run MATLAB from the terminal as follows:
```
LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6 matlab
```
Then, run the following two lines in the MATLAB terminal:
```
flag = int32(bitor(2, 8));
py.sys.setdlopenflags(flag);
```
Now you should be all set!

NOTE: I'm doing something rather strange with ROBIN. Instead of compiling the MATLAB binaries and using those directly, I am using MATLAB's python interface. This is mainly due to the availability of already-written python code and not for any good reason.


## Organization
There are two main scripts:
1. [dense_tracking.m](dense_tracking.m) runs outlier-free certifiable tracking. It contains results processing code intended to evaluate the accuracy of our method. This script does not require the GNC dependency.
2. [tracking_gnc.m](tracking_gnc.m) runs certifiable tracking with outliers using GNC.

The subfolders contain helper functions. `solvers` holds functions for our method and comparison methods. `utils` holds a mix of helper functions to generate, process, and solve the SDP. `visualization` holds visualization tools.