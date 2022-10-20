## Description

* `check_convergence.r` is used for Fig. S3 in SI which shows the dominance index without advection for different simulations durations (corresponding file: `convergence_wo_advection.pdf`).
* `compare_advection.r` plots Figs. 3 and 4 in main text, i.e. the dominance index for micro- and nanophytoplankton communities with and without advection (corresponding files: `dominance_diatom_nano_compare_advection*.pdf`). It is also used to plot Fig. 5 in main text (corresponding file: `carac_10species.pdf`).
* `compare_advection_comparison_th.r` plots the same figure as above, but adds a comparison with theoretical values (not shown in main paper for the sake of clarity).
* `dominance_code_*.r` presents examples of spatial distributions, pcf, Ripley's K-function and dominance indices for a Poisson or Thomas point process (corresponding to Figs. S1 and S2 in SI, i.e., `example_*_distribution.pdf`). 
* `dominance_gamma.r` plots Fig. S9 in SI which presents the shift in dominance index due to different ways of computing the value of _U_ (corresponding file: `theoretical_dominance_with_adv.pdf`) 
* `find_gamma.r` computes the value of ![\gamma](https://latex.codecogs.com/svg.latex?\gamma) corresponding to our value of _U_.
* `min_dist_abundances_v2.r` shows the value of the minimum or mean distances to the nearest neighbour as a function of abundances in different plankton communities (Fig. S11 in SI, `dist_abundances_*.pdf`)
* `pcf_bandwidth_BBM.r` shows the impact of the kernel bandwidth value used to compute the pcf for the Brownian Bug Model (Fig. S7 in SI, `bandwidth_BBM.pdf`).
* `plot_K_BBM_micronano.r` plots the comparison between simulation and theory for Ripley's K-function for the Brownian Bug Model (Fig. 2 in main text, `K_micronano.pdf`).
* `quick_view_Thomas.r` plots the comparison between simulation and theory for the pcf and K function for Poisson and Thomas point process (Figs S5 and S6 in SI, `K_PCF_*.pdf`).
* `read_minimum_distances.r` plots Fig. S10 in SI which compares the distance to the nearest neighbour in the BBM and in a uniform spatial distribution (corresponding files: `distrib_distance_*_box_10sp.pdf`).
* `spatial_distribution.r` plots spatial distributions of particles (Fig. 1 in main text and S8 in SI, corresponding files: `spatial_distribution_zoom*.pdf`)
* `theoretical_functions.r` contains the theoretical functions for the pcf and Ripley's K-functions of Thomas and Poisson point processes as well as the 3D Brownian Bug Model.
* `theory_dominance.r` shows the dominance index with and without advection with different simulation durations (Fig. S4 in SI, `theory_dominance.pdf`).
* `utilitary_functions.r` contains a plotting function used in several other scripts (log axes)

It should be noted that most scripts are applied to simulation files that are identified by numbers (most of the time presented through lines similar to `sim_diatom=list(10,11)`). These files are not present in the folder but can be obtained by running `main_BBM.cpp` or `main_PoissonThomas.cpp`.
