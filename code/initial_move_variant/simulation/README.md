## Description

* `basic_particle` files define the class describing a single individual and implements the movement functions (diffusion and advection).
* `main_BBM.cpp` is used to run all simulations for the Brownian Bug Model. The set of parameters used for nanophytoplankton and microphytoplankton can be chosen at the beginning of the file. In addition to spatial distributions and different indices (pcf, Ripley's K-functions and dominance), the script outputs a summary of the parameters used in the simulation, identified by a number defined at the beginning of the file. Each simulation should be identified by a different number. **Here the model is modified to implement a fixed move by one cell diameter upon birth, in a random direction.**
* `compute_gamma_BBM.cpp` is used to approximate the value of ![\gamma](https://latex.codecogs.com/svg.latex?\gamma) for a given value of _U_ by computing the separation between particles as a function of time (details of this computation are given in Young et al. 2001 and Picoche et al. 2022).

Executable files can be produced with the command `make -f makefile` (for `main_BBM.cpp`, and for `main_PoissonThomas.cpp` with slight modifications) and `make -f makefile_gamma` for `compute_gamma_BBM.cpp`. 

Extraction of relevant values (stretching parameter or pair density) and figures showing the results of these simulations can be found in the folder `../figure`.

In this folder we have additionnally

* `initial_move_test.cpp` which tests the procedure making the move of a cell diameter in a random direction
* `plot_traj_3D.R` which plots it in 3D. 

