# Code for the paper *Joint estimation of bipartite networks collections. Application to plant-pollinator networks*

This repository presents the code that generates the results and figures for the colBiSBM method.

The code is written in R.

The scripts are organized as follows:

- `simulations_*.R` : Simulations scripts used to generate the results of the paper. Note that the simulations are run in parallel and need to be run on a cluster. Running them on a local machine may take a long time.
- `application_*.R` : Scripts used to generate the results of the application to plant-pollinator networks. These scripts are also run in parallel and need to be run on a cluster.
- `analyze_*.R` : Scripts used to analyze the results of the simulations and the application. These scripts are used to generate the figures and tables of the paper.
- Other name files : These files are additional scripts that were used throughout the analysis for misceallaneous tasks. They are not used in the paper but may be useful for the reader to understand the code.
