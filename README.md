# Local intraspecific aggregation in phytoplankton model communities: spatial scales of occurrence and implications for coexistence

**Authors**: Coralie Picoche, William R. Young and Frédéric Barraquand

This repository contains the mathematical analyses of the three-dimension, multispecies Brownian Bug Model described in the manuscript _Local intraspecific aggregation in phytoplankton model communities: spatial scales of occurrence and implications for coexistence_>, as well as the scripts running corresponding simulations and analyses.

### Organisation

Folders are organised as follows:

* `article` contains all files used to produce the manuscript of the article.
* `code/simulation` contains the code for the simulations of the Brownian bug model. Output files are also stored in this folder.
* `code/figure` contains the codes to produce figures shown in the article, as well as said figures.

### Pre-requisites

Simulations of the Brownian bug model are run in C++ (v.5.4.0, using the GSL Library v2.6) and output files are then treated separately with R (v.3.6.3) to produce figures.

All codes were tested on Ubuntu 16.04 and Ubuntu 20.04. 

### License

Computer codes are released under the GNU General Public License (GPLv3) and data under the Creative Commons Attribution 4.0 International License.

