# YFT-lattice-st-CPUE-models
Yellowfin Tuna NOAA spatial simulation experiment - R code to standardize lattice spatio-temporal CPUE data
https://github.com/aaronmberger-nwfsc/Spatial-Assessment-Modeling-Workshop

In this repository we standardize a CPUE from lattice (regular polygons) data via R-INLA. Three different models are used, 1) iid spatial model with time trend, 2) besag spatial model with time trend and 3) besag spatio-temporal interaction model.

Predictions performed via inla.posterior.sample() function.

Code to create figure maps and model effects plotting included.
