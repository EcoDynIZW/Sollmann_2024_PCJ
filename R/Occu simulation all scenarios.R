################################################################################

## R code to implement occupancy simulation study described in Sollmann, R. 
## (in review). 
## Mt or not Mt: Temporal variation in detection probability in spatial 
## capture-recapture and occupancy models. Peer Community Journal.

## Simulation scenarios
## scenario 1: variable causing variation in p, Z, is identical for all detectors 
## scenario 2: Z varies randomly across traps and occasions
## scenario 3: Z varies randomly across traps and occasions covaries in space
##             with variable causing variation in density
## Each scenario is implemented with low (~3-fold) and high (~6-fold) temporal
## variation in p.
## For each scenario, the data-generating model accounting for temporal variation
## in p0 is compared to an alternative model NOT accounting for temporal variation 
## in detection. 
## Note that for Sc 2 and 3, a null model is also run but is not included in the 
## summary script as the focus is on ignoring temporal, not spatial, variation.

library(unmarked)
source('Helper functions.txt')


##############################################################################
######### constant settings across all scenarios #############################

## number of iterations
niter=250

## detector locations
X<-as.matrix(expand.grid(seq(0, 900, 100), seq(0, 900, 100)))
J<-nrow(X)

## number of sampling occasions
K<-8

##calculate distances between detectors (for spatially correlated predictor, X)
G<-cbind(X[,1], X[,2])
dg <- e2dist(G, G)

## effect of covariate X on occurrence
beta.d <- 1

##intercept occurrence (logit scale)
b0<-0.5

##calculate true occu at X=0
true.occu<-plogis(b0)


#############################################################################
########### Simulation ######################################################

## label for magnitude of change in p over time
mag<-c(6,3)

## set scenario of interest
scenario<-3


##loop over the two magnitude options
for (fold in 1:2){
  
## keep track of estimates, AIC
## here, dens = intercept
AIC.mat<-matrix(NA, niter, 3)
colnames(AIC.mat)<-c('m0', 'mt', 'mt2')
beta.0<-beta.t<-beta.t2<-matrix(NA, niter, 2)
colnames(beta.0)<-colnames(beta.t)<-colnames(beta.t2)<-c('estimate', 'SE')#, 'lower', 'upper')

## track prediction of occu at X=0
occu.0<-occu.t<-occu.t2<-matrix(NA, niter, 4)
colnames(occu.0)<-colnames(occu.t)<-colnames(occu.t2)<-c('estimate', 'SE', 'lower', 'upper')

## track convergence
converged<-matrix(NA, niter, 3)
colnames(converged)<-c('m0', 'mt', 'mt2')

## keep track of data, all important model output
dat.list<-all.out<-list()

set.seed(1234)

## loop over iterations
for (iter in 1:niter){
  
  ## generate spatially correlated covariate
  X.d <- spcov2(dg)
  psi.eff<-plogis(b0+beta.d*X.d)
  
  ## pick detection parameters depending on scenario and magnitude of change in p
  
  if (scenario == 1){
    tempcov<-matrix(0:7,J,K, byrow=TRUE)
    p0<-c(-1, -1.38)[fold] 
    b.tc<-c(-0.3, -0.18)[fold] 
  }

  if (scenario == 2){
    tempcov.mean<-rnorm(J)
    sdt<-c(0.9, 0.5)
    tempcov<-matrix(NA, J,K)
    for (j in 1:J){
      tempcov[j,]<-rnorm(K, tempcov.mean[j],sdt[fold])#0.5
    }
  
    p0<-c(-1.95,-1.8)[fold]
    b.tc<-c(-0.78, -0.9)[fold]
  }

  if (scenario == 3){
 
    tempcov.mean<-rnorm(J, X.d, 0.5)
    sdt<-c(0.9, 0.5)
    tempcov<-matrix(NA, J,K)
    for (j in 1:J){
      tempcov[j,]<-rnorm(K, tempcov.mean[j],sdt[fold])
    }
  
    p0<-c(-1.95,-1.8)[fold]
    b.tc<-c(-0.78, -0.9)[fold]
  }

  ## calculate detection probability
  p.eff<-plogis(p0 + b.tc*tempcov)
  
  ## compile covariate data
  mean.tc=apply(tempcov,1,mean)
  sitecovs<-data.frame(X=X.d, mean.tc=mean.tc)
  obscovs<-list(tempcov=data.frame(tempcov))

  ## generate detection data
  z<-rbinom(J,1,psi.eff)
  obs<-matrix(NA, J,K)
  for (k in 1:K){
    obs[,k]<-rbinom(J,1,p.eff[,k]*z)
  }
  
  ## bundle for unmarked, run models
  umf<-unmarkedFrameOccu(y=obs, siteCovs = sitecovs, obsCovs = obscovs)
  m0<-occu(~1~X, umf)
  mt<-occu(~tempcov~X, umf)
  ## mts makes no sense for scenario 1 but this is so fast that it's easier to
  ## run this than to write the code to skip it
  mt2<-occu(~mean.tc~X, umf)
  
  ##run times are trivial but 2d model is still 3x faster 
  
  ## track AIC
  AIC.mat[iter,]<-c(m0@AIC, mt@AIC, mt2@AIC)
  
  ## track convergence based on SEs, optim and beta>10
  converged[iter,1]<-ifelse(m0@opt$convergence != 0 | 
                              any(is.na(sqrt(diag(solve(m0@opt$hessian))))) |
                              any( abs(coef(m0)) >10), 0, 1)
  converged[iter,2]<-ifelse(mt@opt$convergence != 0| 
                              any(is.na(sqrt(diag(solve(mt@opt$hessian)))))|
                              any( abs(coef(mt)) >10), 0, 1)
  converged[iter,3]<-ifelse(mt2@opt$convergence != 0| 
                              any(is.na(sqrt(diag(solve(mt2@opt$hessian)))))|
                              any( abs(coef(mt2)) >10), 0, 1)
  
  ## track parameter estimates
  beta.0[iter,1:2 ]<-unlist(summary(m0)$state['X',1:2])
  beta.t[iter,1:2 ]<-unlist(summary(mt)$state['X',1:2])
  beta.t2[iter,1:2 ]<-unlist(summary(mt2)$state['X',1:2])
  
  ## get prediction at X=0
  occu.0[iter,]<-unlist(predict(m0, 'state', newdata= data.frame(X=0)))
  occu.t[iter,]<-unlist(predict(mt, 'state', newdata= data.frame(X=0)))
  occu.t2[iter,]<-unlist(predict(mt2, 'state', newdata= data.frame(X=0)))
  
  ## compile all betas and var-covar, just in case
  all.out[[iter]]<-list(vcv=list(m0=solve(m0@opt$hessian), 
                                 mt=solve(mt@opt$hessian), 
                                 mt2=solve(mt2@opt$hessian)),
                        betas=list(m0=coef(m0), mt=coef(mt), mt2=coef(mt2)))
  
  ## compile data
  dat.list[[iter]]<-list(z=z, obs=obs, X.d=X.d, tempcov=tempcov)
  
  ## save everything at each iteration in case R crashes
  ## not strictly necessary
  saveRDS(list(occu.0=occu.0, occu.t=occu.t, 
               occu.t2=occu.t2, 
               beta.0=beta.0, beta.t=beta.t, beta.t2=beta.t2, 
               AIC.mat=AIC.mat, converged=converged),
          paste('out/Results_Scenario_', scenario,'.', mag[fold],'.rds', sep=''))
  saveRDS(all.out, paste('out/Betas_Scenario_', scenario,'.', mag[fold],'.rds', sep=''))
  saveRDS(dat.list, paste('in/Data_Scenario_', scenario,'.', mag[fold],'.rds', sep=''))
  
}
}

## for the curious, quick plots of results
boxplot(beta.t[converged[,2]==1,1],
        beta.t2[converged[,3]==1,1], ylim=c(-1,5))
abline(h=beta.d, col='red')
boxplot(occu.0[converged[,1]==1,1], occu.t[converged[,2]==1,1],
        occu.t2[converged[,3]==1,1])
abline(h=true.occu, col='red')
