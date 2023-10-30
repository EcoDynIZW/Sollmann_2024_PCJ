################################################################################

## R code to implement SCR case study described in Sollmann, R. (in review). 
## Mt or not Mt: Temporal variation in detection probability in spatial 
## capture-recapture and occupancy models. Peer Community Journal.

## Script uses data already formatted for analysis in secr.
## Raw data obtained from: 
## Amburgey SM, Lardner B, Knox A, Converse SJ, Yackel Adams AA (2021) 
## Brown Treesnake detections on transects using potential attractants of 
## live-mouse lures or fish-spray scent, Guam. U.S. Geological Survey 
## data release. http://doi.org/10.5066/P9G6JHZ3
## Please cite original data source if using the processed data.

## For background and original analysis of the data, see:
## Amburgey SM, AA Yackel Adams, B Gardner, B Lardner, AJ Knox, and SJ Converse. 
## 2021. Tools for increasing visual encounter probabilities for invasive species 
## removal: a case study of brown treesnakes. Neobiota 70:107-122. 
## 10.3897/neobiota.70.71379

rm(list = ls())
library(secr)

## read in data
dat<-readRDS('data-raw/Snake_secr_data.rds')

##make trap object for secr analysis
trp<-read.traps(data=dat$locs, detector='proximity')

## make discrete state space
delta<- 11.874929 #taken from original analysis
mask<-make.mask(trp, buffer=delta, spacing = 5)

# Matrix of effort where locations are active (surveyed = 1)/not active (not surveyed = 0) for all dates
act <- as.matrix(dat$act)
# Number of survey occasions
nocc <- ncol(act)

## add as usage to trap file
usage(trp)<-act

## turn binary scent treatment into 'proportion of active occasions with scent'
## for model Ms
pscent<-apply(dat$stat, 1, sum, na.rm=T)/apply(act, 1, sum)

## set up scent treatment as time varying covariate
## NA=location/date combination not surveyed
## 0 no scent/old scent
## 1 fresh scent

timevaryingcov(trp) <- list(blockt = paste('X',1:nocc, sep=''))
covariates (trp) <- data.frame(dat$stat, pscent=pscent)


## make capthist file
chs<-make.capthist(dat$cap.df, traps=trp, fmt='trapID', noccasions = nocc)


## compress data into binomial counts format (for non-temp models only)
chs.comp<-reduce(chs, outputdetector = 'count', by = 'all')

# #M0: 76.83 seconds
# system.time(
#   (m0<-secr.fit(chs.comp, model=list(D~1, g0~1, sigma~1), mask=mask,
#                 trace=TRUE,details = list(fastproximity = FALSE), binomN = 1))
# )

#Ms, compressed data: 689.42 seconds
system.time(
  (mS<-secr.fit(chs.comp, model=list(D~1, g0~pscent, sigma~1), mask=mask,
                trace=TRUE,details = list(fastproximity = FALSE), binomN = 1))
)

#Ms uncompressed, for AIC comparison with Mst: 7433.71
system.time(
  (mS2<-secr.fit(chs, model=list(D~1, g0~pscent, sigma~1), mask=mask,
                trace=TRUE,details = list(fastproximity = FALSE)))
)

# run model with scent covariate on detection - 2038.68 seconds
## use uncompressed data
## Note: Warning that missing covariate values are set to -1 can be ignored 
##       as these only refer to unsampled occasions, which are taken into
##       account via the usage argument (above)
system.time(
  (mT<-secr.fit(chs, model=list(D~1, g0~blockt, sigma~1), mask=mask,
                trace=TRUE))
)#

## extract AIC
AIC(mS2, mT)

## calculate density per ha
round(predict(mT, realnames=c('D'))[2:5], dig=2)
round(predict(mS2, realnames=c('D'))[2:5], dig=2)


## get detection parameters on natural scale
predict(m0, realnames='g0')
#no scent
plogis(coef(mT)['g0',1])
#scent
plogis(coef(mT)['g0',1]+coef(mT)['g0.blockt',1])

predict(mS2, realnames='sigma')
predict(mT, realnames='sigma')
