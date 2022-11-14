# YFT-lattice-st-CPUE-models
Yellowfin Tuna NOAA spatial simulation experiment - R code to standardize lattice spatio-temporal CPUE data

https://github.com/aaronmberger-nwfsc/Spatial-Assessment-Modeling-Workshop

In this repository we standardize a CPUE from lattice (regular polygons) data via R-INLA. Three different models are applied:

01) iid model with time trend
02) besag spatial model with time trend 
03) besag spatio-temporal interaction model

Spatio-temporal predictions are performed via inla.posterior.sample() function.

R code to create figure maps and model effects plotting included.
