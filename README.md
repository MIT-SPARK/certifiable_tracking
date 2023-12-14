# MATLAB Implementation of PUMO (Pose estimation Under Multiple Observations)

## Setup
You can play with a limited set of results without installing dependencies. To run the solver, this repository has the following dependencies:
- [CertifiablyRobustPerception](https://github.com/MIT-SPARK/CertifiablyRobustPerception/tree/master)
- [MOSEK](https://www.mosek.com/)
- [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)
- [GNC](https://github.com/MIT-SPARK/GNC-and-ADAPT) (only for tracking with outliers)

You may set the paths to these dependencies in [dense_tracking.m](dense_tracking.m) and [tracking_gnc.m](tracking_gnc.m).

## Example (no dependencies required)
TBD