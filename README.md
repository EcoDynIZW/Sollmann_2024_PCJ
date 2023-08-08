# Data and code for "Mt or not Mt: Temporal variation in detection probability in spatial capture-recapture and occupancy models", Sollmann, R., Peer Community Journal (submitted)
---

DOI (repository): https://doi.org/10.5281/zenodo.8221253

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

Simulations require packages secr (Efford, 2022) and unmarked (Fiske and Chandler, 2011). Datasets for case studies come with R packages oSCR (Sutherland et al. 2018) and AHMbook (Kéry et al., 2022). Please see R scripts for additional packages used in processing data. 


## References

Efford MG (2022) secr: Spatially explicit capture-recapture models. R package version 4.5.6. https://CRAN.R-project.org/package=secr 

Fiske I, Chandler R (2011) Unmarked: an R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software, 43, 1–23. https://doi.org/10.18637/jss.v043.i10 

Kéry M, Royle JA, Meredith M (2022) AHMbook: Functions and Data for the Book “Applied Hierarchical Modeling in Ecology” Vols 1 and 2. R package version 0.2.6. https://CRAN.R-project.org/package=AHMbook

Sutherland C, Royle JA, Linden D (2018) oSCR: Multi-Session Sex-Structured Spatial Capture-Recapture Models. R package version 0.42.0. https://CRAN.R-project.org/package=AHMbook 
