# Data and code for "Mt or not Mt: Temporal variation in detection probability in spatial capture-recapture and occupancy models", Sollmann, R., Peer Community Journal (submitted)
---

DOI: https://doi.org/10.5281/zenodo.8221253

All code necessary to implement the simulation and case studies described in the above references manuscript

## Description of the data and file structure

This repository contains the following code files:

Occu simulation all scenarios.R: R code to generate occupancy data and analyze them under a model that does and one that does not account for temporal variation in detection

Occupancy case study MHB.R: R code to repeat occupancy case study for 10 breeding birds from Switzerland

secr case study Ft Drum bears.R: R code to repeat SCR case study for bears in Ft Drum, NY

secr simulation all scenario random Z.R: R code to generate SCR data and analyze them under a model that does and one that does not account for temporal variation in detection

Summarize occu simulations.R: R code to summarize occupancy simulation output and create main Figure 2 and some appendix figures

Summarize secr sim all scenarios.R: R code to summarize SCR simulation output and create main Figure 1 and some appendix figures


## Sharing/Access information

NA


## Code/Software

Simulations require packages secr and unmarked. Data for case studies comes with R packages oSCR (SCR) and AHMbook (occupancy). 
Please see publication for version information and R scripts for additional packages used in processing data. 
