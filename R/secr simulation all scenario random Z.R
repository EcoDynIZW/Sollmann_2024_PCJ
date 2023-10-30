################################################################################

## R code to implement SCR simulation study described in Sollmann, R. (in review). 
## Mt or not Mt: Temporal variation in detection probability in spatial 
## capture-recapture and occupancy models. Peer Community Journal.

## Main simulation scenarios
## scenario 1: variable causing variation in p, Z, is identical for all detectors 
## scenario 2: Z varies randomly across traps and occasions
## scenario 3: Z varies randomly across traps and occasions covaries in space
##             with variable causing variation in density
## For each scenario, the data-generating model accounting for temporal variation
## in p0 is compared to an alternative model NOT accounting for temporal variation 
## in detection. 
## Note that for Sc 2 and 3, a null model is also run but is not included in the 
## summary script as the focus is on ignoring temporal, not spatial, variation.

## Addl scenarios - variants of scenario 2
## scenario 4 (2a in ms): varying variability of tempcov across sites
## scenario 5 (2b in ms): different levels/directions of skew in tempcov
## scenario 6 (2c in ms): temporal correlation in variation in p

library(secr)
library(sn)
source('R/Helper functions.txt')

##############################################################################
######### constant settings across all scenarios #############################

## number of iterations
niter=100 

## detector locations, and formatting of detector data for secr
X<-as.matrix(expand.grid(seq(0, 600, 100), seq(0, 700, 100)))
J<-nrow(X)
X.df<-data.frame(trapID=1:J, x=X[,1], y=X[,2])
trp<-read.traps(data=X.df, detector='proximity')

## scale parameter, detection function
sigma<-80

## number of sampling occasions
K<-8

## effect of covariate on density ('beta' in manuscript)
beta.d <- 1

## create discrete state space S and set up for secr
mask<-make.mask(trp, buffer=300, spacing = 50)

xlim<-range(attributes(mask)$boundingbox[,'x'])
ylim<-range(attributes(mask)$boundingbox[,'y'])

## calculate distances between all cells in S for generation of spatially
## correlated covariate
G<-cbind(mask$x, mask$y)
dg <- e2dist(G, G)

## population size in S
N<-60

#############################################################################
########### Simulation ######################################################

## set scenario of interest (1 through 6)
scenario=4


## keep track of estimates, AIC, convergence
dens.0<-dens.t<-dens.t2<-sigma.0<-sigma.t<-sigma.t2<-matrix(NA, niter, 4)
colnames(dens.0)<-colnames(dens.t)<-colnames(dens.t2)<-
  colnames(sigma.0)<-colnames(sigma.t)<-colnames(sigma.t2)<-c('estimate', 'SE', 'lower', 'upper')
AIC.mat<-matrix(NA, niter, 3)
colnames(AIC.mat)<-c('m0', 'mt', 'mt2')
beta.0<-beta.t<-beta.t2<-matrix(NA, niter, 4)
colnames(beta.0)<-colnames(beta.t)<-colnames(beta.t2)<-c('estimate', 'SE', 'lower', 'upper')
converged<-matrix(NA, niter, 3)
colnames(converged)<-c('m0', 'mt', 'mt2')


## keep track of generated data, all important model output
dat.list<-all.out<-list()

## set random seed, same for all scenarios
set.seed(1234)

## loop over iterations
for (iter in 1:niter){
  
  ## generate spatially correlated covariate on density, X
  X.d <- spcov2(dg)
  
  ## add to mask as covariate data
  covariates(mask)<-X.d
  
  ## generate Z, set detection parameters according to scenario
  if (scenario == 1){
    tempcov<-matrix(0:7,J,K, byrow=TRUE)
    p0<--2.1
    b.tc<--0.3
  }
  
  if (scenario == 2){
    tempcov.mean<-rnorm(J)
    tempcov<-matrix(NA, J,K)
    for (j in 1:J){
      tempcov[j,]<-rnorm(K, tempcov.mean[j],0.5)
    }
    
    p0<--3.2
    b.tc<--0.8
  }
  
  if (scenario == 3){
    
    dxg<-e2dist(X, G)
    neigh<-t(apply(dxg, 1, function(x)which(x==min(x))))
    X.j<-apply(matrix(X.d[neigh], J, ncol(neigh)),1, mean) #column order
    
    tempcov.mean<-rnorm(J, X.j, 0.5)
    tempcov<-matrix(NA, J,K)
    for (j in 1:J){
      tempcov[j,]<-rnorm(K, tempcov.mean[j],0.5)
    }
    p0<--3.2
    b.tc<--0.8
  }
  
  if (scenario == 4){
    ##alternative value in sigvec and b.tc - scenario 4 high 
    ##with avg 10-fold variation in p0 per detector
    tempcov.mean<-rnorm(J)
    sigvec<-c(0.2, 0.5, 1, 1.2) #1.5
    tempcov<-matrix(NA, J,K)
    for (j in 1:J){
      tempcov[j,]<-rnorm(K, tempcov.mean[j],sample(sigvec, 1))
    }
    
    p0<--3.2
    b.tc<--0.7 #-0.8
  }
  
  if (scenario == 5){
    
    ##have more or less skewed distributions
    tempcov.mean<-rnorm(J)
    #set varying skew
    alphavec<-c(-10, -2, 0, 2, 10)
    tempcov<-matrix(NA, J,K)
    for (j in 1:J){
      tempcov[j,]<-rsn(K, xi=tempcov.mean[j],omega=1, alpha=sample(alphavec, 1))
    }
    
    p0<--3.2
    b.tc<--0.8
  }
  
  if (scenario == 6){
    ## create temporal 'distance' matrix
    Dt<-matrix(NA, K,K)
    for (i in 1:K){
      if (i ==1) {Dt[i,]<-(i-1):(K-1)} else{
        Dt[i,]<-c((i-1):1, 0:(K-i))}
    }
    
    ##divide D by 3 to get decent amount of correlation (0.66)
    ## equivalent to setting alpha=3
    tempcov<-matrix(NA, J,K)
    for (j in 1:J){
      cov <- spcov2(Dt/3, alpha=1, standardize = F)#
      tempcov[j,]<-cov
    }
    
    p0<--3.2
    b.tc<--0.8
  }
  
  ##compile covariates for secr
  covariates (trp) <- data.frame(tempcov, mean.tc=apply(tempcov,1,mean))
  timevaryingcov(trp) <- list(blockt = paste('X',1:K, sep=''))
  
  ## calculate cell probabilities for activity center placement
  p.cell <- exp(beta.d*X.d)/sum(exp(beta.d*X.d))
  
  ## generate activity centers 
  Sxy<-matrix(NA,N,2) 
  s.multi <- rmultinom(N, 1, p.cell) 
  s.g <- apply(s.multi,2, function(x){which(x==1)})
  Sxy[,1]<-G[s.g,1]
  Sxy[,2]<-G[s.g,2]
  
  ## get distance between activity centers and detectors
  D <- e2dist(Sxy, X)

  ## generate observations
  obs <- array(NA, c(N, J, K))
  for (i in 1:N){
    p.d <- exp(-D[i,]^2/(2*sigma^2))  

    for (k in 1:K){
      p.c<-plogis(p0+b.tc*tempcov[,k])
      p.eff<-p.c*p.d
      obs[i,,k] <- rbinom(J, 1, p.eff)
    }
  }

  ## how many individuals detected?
  n <- sum(apply(obs, 1, sum, na.rm = TRUE) >0)
  
  ## subset detection data to detected individuals only
  seen <- which(apply(obs, 1, sum, na.rm = TRUE)>0)
  Y <- obs[seen,,] 

  #make 'flat' capture history for secr (there is probably a function for that...)
  cap.df<-data.frame(session=NA, ID=NA, occasion=NA, trap=NA)

  for (i in 1:n){
    for (j in 1:J){
     for (k in 1:K){
        if(Y[i,j,k]==0)next
        sub<-data.frame(session=1, ID=i, occasion=k, trap=j)
        cap.df<-rbind(cap.df, sub)
      }
    }
  }

  ## remove first row of NAs
  cap.df<-cap.df[-1,]

  ## convert to secr format
  chs<-make.capthist(cap.df, traps=trp, fmt='trapID', noccasions = K)

  ## run models
  m0<-secr.fit(chs, model=list(D~value, g0~1, sigma~1), mask=mask,
             details = list(fastproximity = FALSE), trace=FALSE)
  mt<-secr.fit(chs, model=list(D~value, g0~blockt, sigma~1), mask=mask,
             trace=FALSE)

  ## calculate N for S
  N0<-region.N(m0)
  Nt<-region.N(mt)

  ## save N output
  dens.t[iter,]<-unlist(Nt['R.N',1:4])
  dens.0[iter,]<-unlist(N0['R.N',1:4])

  ## save sigma
  out.mt<-predict(mt, realnames=c('sigma'))
  sigma.t[iter,]<-unlist(out.mt['sigma',2:5])
  sigma.0[iter,]<-unlist(summary(m0)$predicted['sigma',2:5])

  ## save coefficient density covariate, beta
  beta.0[iter,]<-unlist(coef(m0)['D.value',1:4])
  beta.t[iter,]<-unlist(coef(mt)['D.value',1:4])
  
  ## save AICc
  AIC.mat[iter,]<-c(AIC(m0)[,'AICc'], AIC(mt)[,'AICc'],NA)

  ## flag non-converged and non-SE models
  converged[iter,1]<-ifelse(m0$fit$code !=1 | any(is.na(coef(m0)[,'SE.beta'])),
                          0,1)
  converged[iter,2]<-ifelse(mt$fit$code !=1 | any(is.na(coef(mt)[,'SE.beta'])),
                          0,1)
  
  ## save varcov and raw betas for scenario 1 (no mt2)
  if(scenario ==1){
  all.out[[iter]]<-list(vcv=list(m0=m0$beta.vcv, mt=mt$beta.vcv),
                        betas=list(m0=coef(m0), mt=coef(mt)))
  
  ## saving and over-writing at each iteration avoids loss if the loop is
  ## interrupted, but is not strictly necessary
  saveRDS(all.out, paste('out/Betas_Scenario_', scenario,'.rds', sep=''))
  }
  
  #run Ms for all but scenario 1, save output as above
  if(scenario >1){
    mt2<-secr.fit(chs, model=list(D~value, g0~mean.tc, sigma~1),mask=mask,
              details = list(fastproximity = FALSE), trace=FALSE)

    Nt2<-region.N(mt2)
    dens.t2[iter,]<-unlist(Nt2['R.N',1:4])
    sigma.t2[iter,]<-unlist(summary(mt2)$predicted['sigma',2:5])
    beta.t2[iter,]<-unlist(coef(mt2)['D.value',1:4])
    AIC.mat[iter,3]<-AIC(mt2)[,'AICc']
    converged[iter,3]<-ifelse(mt2$fit$code !=1 | any(is.na(coef(mt2)[,'SE.beta'])),
                          0,1)

    all.out[[iter]]<-list(vcv=list(m0=m0$beta.vcv, mt=mt$beta.vcv, mt2=mt2$beta.vcv),
                          betas=list(m0=coef(m0), mt=coef(mt), mt2=coef(mt2)))
    saveRDS(all.out, paste('out/Betas_Scenario_', scenario,'.rds', sep=''))
  }
  
  #save data
  dat.list[[iter]]<-list(Sxy=Sxy, obs=obs, X.d=X.d, tempcov=tempcov)
  saveRDS(dat.list, paste('in/Data_Scenario_', scenario,'.rds', sep=''))
  
  #save estimates, AIC, convergence results
  saveRDS(list(dens.0=dens.0, dens.t=dens.t, 
             dens.t2=dens.t2, sigma.0=sigma.0, 
             sigma.t=sigma.t, sigma.t2=sigma.t2,
             beta.0=beta.0, beta.t=beta.t, beta.t2=beta.t2, 
             AIC.mat=AIC.mat, converged=converged),
        paste('out/Results_Scenario_', scenario,'.rds', sep=''))
  
  ## keep track of where in the iteration loop you are
  print(iter)
}#end iter loop


##quick check for the curious - estimates of beta from 3 models
boxplot(beta.0[,1], beta.t[,1],beta.t2[,1])
abline(h=beta.d, col='red')
