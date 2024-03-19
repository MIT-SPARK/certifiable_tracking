# Certifiable Tracking (MATLAB Implementation)

## Setup
To run the solver, this repository has the following dependencies:
- [CertifiablyRobustPerception](https://github.com/MIT-SPARK/CertifiablyRobustPerception/tree/master) (you need to run SPOTLESS install script)
- [MOSEK](https://www.mosek.com/)
- [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)

For outlier rejection, we use the COPT solver:
- [COPT](https://www.copt.de/) (we also provide a cvx version that works with MOSEK, it is just slower)

Reproducing results also requires the following datasets:
- [PASCAL3D+](https://cvgl.stanford.edu/projects/pascal3d.html)
- TODO: Link to other datasets

Reproducing ablations requires the following dependencies:
- [UKF-M](https://github.com/CAOR-MINES-ParisTech/ukfm/tree/master) (for PACE+UKF; use `pip install -e .`)
- [ROBIN](https://github.com/MIT-SPARK/ROBIN) (for OURS+ROBIN; be sure turn OFF the unit tests in `cmakelists.txt`)

You may set the paths to these dependencies in [setup.m](setup.m). I recommend cloning this repository and these dependencies in the same parent folder.

### Comments on Python
MATLAB's python integration isn't the easiest to debug. If you run into any python errors, first make sure you have all the dependencies installed.

### Running MATLAB
On my machine (Ubuntu 22.04) I needed to set the `LD_PRELOAD` environment variable. To do this, launch MATLAB with the following line:
```
LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6 matlab
```
This appears to be necessary to run ROBIN.

### Quick Start
Each time you start MATLAB you must run [setup.m](setup.m) to add the necessary paths. The setup file also allows you to select which experiment/dataset you wish to work with. Options are the subfolders within `experiments`.

Suppose you wish to run the experiments that use synthetic data. First run `setup` with `experiments="synthetic"`. Then, from the home directory run [tracking_outlier_free](experiments/synthetic/tracking_outlier_free.m).

## Reproduce Experiments
TODO

## Notes on Organization
General purpose scripts are in the `visualization`, `utils`, `solvers`, and `outlier_rejection` folders. Within the `experiments` folder are scripts specialized to individual experiments. In general, scripts that begin with "tracking" are single-batch test scripts, while scripts that begin with "EXP" are full experiments.