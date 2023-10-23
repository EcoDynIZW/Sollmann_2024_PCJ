################################################################################

## R code to implement SCR case study described in Sollmann, R. (in review). 
## Mt or not Mt: Temporal variation in detection probability in spatial 
## capture-recapture and occupancy models. Peer Community Journal.

library(secr)
library(oSCR)

## load data on black bears (see ms for details)
## data comes with the oSCR package

data('beardata')

## extract detector data
tdf<-beardata$trapmat
J <- dim(tdf)[1]

## extract encounter history data
edf<-beardata$edf

## number of sampling occasions
K=8

## format for analysis in sect
X.df<-data.frame(trapID=1:J, x=tdf$V1*1000, y=tdf$V2*1000)
trp<-read.traps(data=X.df, detector='proximity')

## make discrete state space
mask<-make.mask(trp, buffer=6000, spacing = 500)

## set up time covariate (occasion; same for all detectors)
tempcov<-matrix(0:(K-1),J,K, byrow=TRUE)
timevaryingcov(trp) <- list(blockt = paste('X',1:K, sep=''))
covariates (trp) <- data.frame(tempcov)

## format detection data for secr
cap.df<-edf
colnames(cap.df)<-c('session', 'ID', 'occasion', 'trap', 'sex')
chs<-make.capthist(cap.df, traps=trp, fmt='trapID', noccasions = K, covnames='sex')

## run null model 
#7 seconds
system.time(
  (m0<-secr.fit(chs, model=list(D~1, g0~1, sigma~1), mask=mask,
                trace=TRUE,details = list(fastproximity = FALSE)))
)

#run model with quadratic trend on detection - 83 seconds
system.time(
  (mT<-secr.fit(chs, model=list(D~1, g0~blockt+I(blockt^2), sigma~1), mask=mask,
                trace=TRUE))
)#

## calculate density per 100km2
round(predict(mT, realnames=c('D'))[2:5]*10000, dig=2)
round(predict(m0, realnames=c('D'))[2:5]*10000, dig=2)


## extract AIC, get detection parameters on natural scale
AIC(m0, mT)
predict(m0, realnames='g0')
plogis(coef(mT)['g0',1])
plogis(coef(mT)['g0',1]+coef(mT)['g0.blockt',1]*5+
         coef(mT)['g0.I(blockt^2)',1]*25)
predict(m0, realnames='sigma')
predict(mT, realnames='sigma')

