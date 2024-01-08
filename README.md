# Data and code for "Mt or not Mt: Temporal variation in detection probability in spatial capture-recapture and occupancy models", Sollmann, R., Peer Community Journal (2024)
---

DOI (submitted pre-print): https://doi.org/10.1101/2023.08.08.552394

URL (PCJ publication): https://peercommunityjournal.org/articles/10.24072/pcjournal.357/

All code necessary to implement the simulation and case studies described in the above manuscript

## Description of the data and file structure

<ins> **This repository contains the following code files (in folder R)** </ins>:

*Occu simulation all scenarios.R*: R code to generate occupancy data and analyze them under a model that does and one that does not account for temporal variation in detection

*Occupancy case study MHB.R*: R code to repeat occupancy case study for 10 breeding birds from Switzerland

*SCR case study brown tree snakes.R*: R code to repeat SCR case study for brown tree snakes

*SCR case study data prep.R*: R code to format original brown tree snake data for secr analysis; no need to run this script as its output is available in data-raw

*secr simulation all scenario random Z.R*: R code to generate SCR data and analyze them under a model that does and one that does not account for temporal variation in detection

*Summarize occu simulations.R*: R code to summarize occupancy simulation output and create main Figure 2 and some appendix figures

*Summarize secr sim all scenarios.R*: R code to summarize SCR simulation output and create main Figure 1 and some appendix figures

*Helper functions.txt*: R code for two functions, to calculate distances between two sets of coordinates (e2dist) and to generate correlated variables (spcov2)

<ins> **This repository contains the following data files (in folder data-raw)** </ins>:

*Snake_secr_data.rds*: Brown tree snake capture data formatted for secr analysis (original data from Amburgey et al., 2021a; original analysis of the data: Amburgey et al., 2021b) 

## Sharing/Access information

NA


## Code/Software

Simulations require packages secr (Efford, 2022) and unmarked (Fiske and Chandler, 2011). Dataset for occupancy case studies comes with R packages AHMbook (Kéry et al., 2022). Please see R scripts for additional packages used in processing data. 


## References

Efford MG (2022) secr: Spatially explicit capture-recapture models. R package version 4.5.6. https://CRAN.R-project.org/package=secr 

Fiske I, Chandler R (2011) Unmarked: an R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software, 43, 1–23. https://doi.org/10.18637/jss.v043.i10 

Kéry M, Royle JA, Meredith M (2022) AHMbook: Functions and Data for the Book “Applied Hierarchical Modeling in Ecology” Vols 1 and 2. R package version 0.2.6. https://CRAN.R-project.org/package=AHMbook

Amburgey, S.M., Lardner, B., Knox, A.J., Converse, S.J., and A.A. Yackel Adams, 2021a, Brown Treesnake detections on transects using potential attractants of live-mouse lures or fish-spray scent, Guam: U.S. Geological Survey data release, https://doi.org/10.5066/P9G6JHZ3.

Amburgey SM, AA Yackel Adams, B Gardner, B Lardner, AJ Knox, and SJ Converse. 2021b. Tools for increasing visual encounter probabilities for invasive species removal: a case study of brown treesnakes. Neobiota 70:107-122. https://doi.org/10.3897/neobiota.70.71379
