# MATLAB Implementation of PUMO (Pose estimation Under Multiple Observations)

## Setup
You can play with a limited set of results without installing dependencies. To run the solver, this repository has the following dependencies:
- [CertifiablyRobustPerception](https://github.com/MIT-SPARK/CertifiablyRobustPerception/tree/master)
- [MOSEK](https://www.mosek.com/)
- [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)

You may set the paths to these dependencies in [dense_tracking.m](dense_tracking.m).

## Example (no dependencies required)
Sample results are given in [demo_data.mat](data/demo_data.mat). The script [luca_demo.m](luca_demo.m) serves as an example of working with these results.