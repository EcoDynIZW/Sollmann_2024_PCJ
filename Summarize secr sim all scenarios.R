################################################################################

## R code to summarize output from the SCR simulation study described in 
## Sollmann, R. (in review). 
## Mt or not Mt: Temporal variation in detection probability in spatial 
## capture-recapture and occupancy models. Peer Community Journal.

## DISCLAIMER: because scenario 2a with 10-fold variation was an afterthought,
## there is some awkwardness in the code in combining results from all scenarios


library(ggplot2)
library(writexl)
library(sn)
library(plotrix)
library(gridExtra)
library(cowplot)
library(ggpubr)


##############################################################################
######### constant settings across all scenarios #############################

niter=100
X<-as.matrix(expand.grid(seq(0, 600, 100), seq(0, 700, 100)))
J<-nrow(X)
sigma<-80
K<-8

##effect of covariate on density
beta.d <- 1
N<-60

###############################################################################
###############################################################################

## make data frame with estimates of all parameters under all scenarios
out.all<-data.frame(matrix(NA, nrow=0, ncol=7))
colnames(out.all)<-c('estimate', 'SE', 'Parameter', 'Model', 'Bias', 'CV', 'Scenario')
nconverged<-NULL


for (i in 1:6){
  sub<-readRDS(paste('out/Results_Scenario_',i, '.rds', sep=''))
  
  ## get parameter estimates only
  sub2<-lapply(sub[1:9], function(x)x[,1:2])
  xx<-as.data.frame(do.call(rbind, sub2))
  
  ## look at convergence of ALL models for given iteration
  ## for scenarios >1, disregard M0
  cmat<-matrix(sub$converged, 100, 3)
  if(i>1){
    cmat[,1]<-NA
  }
  cvec<-apply(cmat, 1, function(x){any(x==2, na.rm=TRUE)})
  
  ## repeat 3x for all 3 models, then again 3 times for all parameters
  nconverged[i]<-sum(!cvec)
  cvec2<-rep(rep(cvec, 3), 3)
  
  ## if cvec==TRUE, not converged, ie, remove all params for all models
  ## for that iteration
  xx[cvec2,]<-NA
  
  ## create labels for parameter and model
  xx$Parameter<-c(rep('N', 300), rep('sigma', 300), rep('beta', 300))
  xx$Model<-rep(rep(c('M0', 'Mt', 'Mt2'), each=100),3)
  
  ## calculate bias and CV
  xx$Bias<-(xx$estimate- c(rep(N, 300),rep(sigma, 300), rep(beta.d, 300)))/
    c(rep(N, 300),rep(sigma, 300), rep(beta.d, 300))
  xx$CV<-(xx$SE/xx$estimate)*100
  
  ##create label for scenario
  xx$Scenario<-rep(i, 900)

  ##for scenario 1, remove mt2 (not run), for others remove m0 (not of interest)
  if(i==1){
    xx<-xx[xx$Model != 'Mt2',]
  } else {
    xx<-xx[xx$Model != 'M0',]
    
  }
  
  out.all<-rbind(out.all, xx)

}


## to add in scenario 4high, rename all scenarios 5 and higher
out.all$Scenario[out.all$Scenario>4]<-out.all$Scenario[out.all$Scenario>4]+1
nconverged<-c(nconverged[1:4], NA, nconverged[5:6])


## add scenario 4high, call it scenario 5
## code as in loop above
sub<-readRDS('out/Results_Scenario_4.high.rds')
sub2<-lapply(sub[1:9], function(x)x[,1:2])
xx<-as.data.frame(do.call(rbind, sub2))
cmat<-matrix(sub$converged, 100, 3)
cmat[,1]<-NA
cvec<-apply(cmat, 1, function(x){any(x==2, na.rm=TRUE)})
nconverged[5]<-sum(!cvec)
cvec2<-rep(rep(cvec, 3), 3)
xx[cvec2,]<-NA
xx$Parameter<-c(rep('N', 300), rep('sigma', 300), rep('beta', 300))
xx$Model<-rep(rep(c('M0', 'Mt', 'Mt2'), each=100),3)
xx$Bias<-(xx$estimate- c(rep(N, 300),rep(sigma, 300), rep(beta.d, 300)))/
  c(rep(N, 300),rep(sigma, 300), rep(beta.d, 300))
xx$CV<-(xx$SE/xx$estimate)*100
xx$Scenario<-rep(5, 900)
xx<-xx[xx$Model != 'M0',]
out.all<-rbind(out.all, xx)


################################################################################
###### make summary table ######################################################

## calculate mean and median (median used in ms)
mean.bias<-aggregate(cbind(Bias, CV)~Parameter + Scenario+Model, data=out.all, FUN=mean)
median.bias<-aggregate(cbind(Bias, CV)~Parameter + Scenario+Model, data=out.all, FUN=median)

all.summ<-cbind(mean.bias, median.bias$Bias, median.bias$CV)


## get RMSE, divide rmse by truth to get relative rmse

truth<-rep(c(1,60,80), 14) #beta, N, sigma
all.summ$RMSE<-NA
for (jj in 1:nrow(all.summ)){
  
  ## get data pertaining to unique combo of parm, sc, model
  sub<-out.all[out.all$Parameter == all.summ$Parameter[jj] &
                 out.all$Scenario == all.summ$Scenario[jj] &
                 out.all$Model == all.summ$Model[jj],]
  ## calculate absolute bias
  abs.bias<-sub$Bias*truth[jj]
  
  ## calculate RMSE
  all.summ$RMSE[jj]<-sqrt(sum(abs.bias^2, na.rm=TRUE)/sum(!is.na(sub$Bias)))/truth[jj]
}


## order by scenario first, then parameter, then model
all.summ.df<-all.summ[with(all.summ, order(Scenario, Parameter)),]
colnames(all.summ.df)<-c('Parameter', 'Scenario', 'Model', 'Bias.mean','CV.mean',
                         'Bias.median','CV.median', 'RMSE')

## round to 2 digits
for (jj in 4:8){
  all.summ.df[,jj]<-round(all.summ.df[,jj], dig=2)
}

##change some labels in table to match ms

##models: M0, Ms, Mt or Mst 
##For scenarios >1
##Mt is data generating model, which should be labelled Mst (space and time)
##Mt2 is model ignoring temporal variation and should be labelled Ms
all.summ.df$Model[all.summ.df$Scenario>1 & all.summ.df$Model == 'Mt']<- 'Mst'
all.summ.df$Model[all.summ.df$Scenario>1 & all.summ.df$Model == 'Mt2']<- 'Ms'

##scenarios: 4 should be 2a, 5 2a.high, 6 and 7 2b and c
all.summ.df$Scenario<-as.character(all.summ.df$Scenario)
all.summ.df$Scenario[all.summ.df$Scenario == '4']<- '2a'
all.summ.df$Scenario[all.summ.df$Scenario == '5']<- '2a.high'
all.summ.df$Scenario[all.summ.df$Scenario == '6']<- '2b'
all.summ.df$Scenario[all.summ.df$Scenario == '7']<- '2c'
all.summ.df<-all.summ.df[with(all.summ.df, order(Scenario, Parameter)),]

## finally, switch order for scenario 1 so that data generating model always 
## comes first
all.summ.df<-all.summ.df[c(2,1,4,3,6,5, 7:nrow(all.summ.df)),]

## save table
write_xlsx(all.summ.df, 'Summary results all scenarios.xlsx')



################################################################################
######## make AICc table #######################################################

prop<-matrix(NA, 6,2)
colnames(prop)<-c('Data generating', 'Alternative')

for (i in 1:6){
  sub<-readRDS(paste('out/Results_Scenario_',i, '.rds', sep=''))
  ## get AIC only
  sub2<-sub$AIC.mat

  ## get converged, disregard m0 from convergence and AIC
  cmat<-matrix(sub$converged, 100, 3)
  if(i>1){
    cmat[,1]<-NA
    sub2[,'m0']<-NA
  }
  cvec<-apply(cmat, 1, function(x){any(x==2, na.rm=TRUE)})
  sub2[cvec,]<-NA
  
  ## calculate deltaAIC
  if(i==1){
    d1<-sub2[,'mt']-sub2[,'m0']
    #if negative, mt is better; if more than 2 units, unequivocally better
    #if positive, m0 is better, if >2 etc etc
  } else {
    d1<-sub2[,'mt']-sub2[,'mt2']
  }
  
  ## get proportion top model 
  prop[i,]<-c( sum(d1< (-2), na.rm=TRUE)/nconverged[i],
               sum(d1> (2), na.rm=TRUE) /nconverged[i])

}

## make room for 2a.high
prop<-rbind(prop[1:4,], c(NA, NA), prop[5:6,])

## add in 2a.high
sub<-readRDS('out/Results_Scenario_4.high.rds')

sub2<-sub$AIC.mat
cmat<-matrix(sub$converged, 100, 3)
cmat[,1]<-NA
sub2[,'m0']<-NA
cvec<-apply(cmat, 1, function(x){any(x==2, na.rm=TRUE)})
sub2[cvec,]<-NA
d1<-sub2[,'mt']-sub2[,'mt2']
prop[5,]<-c( sum(d1< (-2), na.rm=TRUE)/nconverged[5],
             sum(d1> (2), na.rm=TRUE) /nconverged[5])


## round, add number of converged iterations, rename scenarios and save
prop<-data.frame(round(prop, dig=2))
prop$N.converged<-nconverged
prop$Scenario<-c(as.character(1:3),'2a', '2a.high',paste(2, c('b', 'c'), sep='') )
write_xlsx(prop, 'AIC.All.Scenarios.xlsx')



################################################################################
###### Make Figure 1: plot of bias in estimates ################################

##DISCLAIMER: This is most likely not the cleanest or most efficient way to make
##            and combine this series of plots. Code provided only for transpa-
##            rency, not as template for similar figures


## add correct model names
out.all$Model[out.all$Model == 'Mt' & out.all$Scenario >1]<-'Mst'
out.all$Model[out.all$Model == 'Mt2' & out.all$Scenario >1]<-'Ms'

## add correct scenario names
sc.name<-prop$Scenario
out.all$Scenario.name<-sc.name[out.all$Scenario]

## add indicator of data generating vs alternative model
out.all$ModelCode<-out.all$Model

out.all$Model[out.all$Scenario == 1 & out.all$ModelCode == 'Mt']<-'D.g.'
out.all$Model[out.all$Scenario == 1 & out.all$ModelCode == 'M0']<-'Alt.'

out.all$Model[out.all$Scenario > 1 & out.all$ModelCode == 'Mst']<-'D.g.'
out.all$Model[out.all$Scenario > 1 & out.all$ModelCode == 'Ms']<-'Alt.'

## change level order 
out.all$Model<-as.factor(out.all$Model)
out.all$Model<-relevel(out.all$Model, 'D.g.')


## Strategy: per parameter, plot bias for 3 main scenarios and 4 addl scenarios
##           then, arrange on one page

## for horizontal line at 0
v_line <- data.frame(
  yintercept = c(0,0,0,0)
)

##set colors for data generating and alternative models
colr<-c('grey80', 'grey50')


#############################################################################
### plot N ##################################################################

## main scenarios
bias.Nmain<-out.all[out.all$Parameter == 'N' & out.all$Scenario<=3, ]

bNm<-ggplot(bias.Nmain, aes(x = Scenario.name, y = Bias*100, fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + 
  labs(x = NULL, y = "Relative bias N") +
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-55,100)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = 'none')#+ 


## addl scenarios
bias.Nsub<-out.all[out.all$Parameter == 'N' & out.all$Scenario>3, ]

bNs<-ggplot(bias.Nsub, aes(x = Scenario.name, y = Bias*100, fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + 
  labs(x = NULL, y = NULL) +
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-55,100)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),        
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position = 'none')


################################################################################
###### plot beta ###############################################################

## main scenarios
bias.bmain<-out.all[out.all$Parameter == 'beta' & out.all$Scenario<=3, ]

bbm<-ggplot(bias.bmain, aes(x = Scenario.name, y = Bias*100, fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + 
  labs(x = NULL, y = expression(paste("Relative bias  ", beta, sep=''))) +
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-250,200)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = 'none')#+ 

## addl. scenarios
bias.bsub<-out.all[out.all$Parameter == 'beta' & out.all$Scenario>3, ]

bbs<-ggplot(bias.bsub, aes(x = Scenario.name, y = Bias*100, fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + 
  labs(x = NULL, y = NULL) +
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-250,200)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position = 'none')


################################################################################
###### plot sigma ##############################################################

## main scenarios
bias.smain<-out.all[out.all$Parameter == 'sigma' & out.all$Scenario<=3, ]

bsm<-ggplot(bias.smain, aes(x = Scenario.name, y = Bias*100, fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  labs(x = 'Scenario', y = expression(paste("Relative bias  ", sigma, sep=''))) +
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-50,50)) +
  theme_bw()+
  theme(legend.position = 'none')#+ 

## addl. scenarios
bias.ssub<-out.all[out.all$Parameter == 'sigma' & out.all$Scenario>3, ]

bss<-ggplot(bias.ssub, aes(x = Scenario.name, y = Bias*100, fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  labs(x = 'Scenario', y = NULL) +
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-50,50)) +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position = 'none')


## combine all on a single page

## align plots first by column
plotfull.1<-plot_grid(bNm, bbm, bsm, align = 'v', ncol=1 , 
                     rel_heights = c(1,1,1.2))
plotfull.2<-plot_grid(bNs, bbs, bss, align = 'v', ncol=1 , 
                      rel_heights = c(1,1,1.2))
## then align both columns
plotfull<-plot_grid(plotfull.1, plotfull.2, ncol=2, nrow=1,
                    align='h', rel_widths=c(1, 0.9))

# ## Save in appropriate size to fit on one page with caption
# pdf('Figure 1 Bias SCR.pdf', paper='letter', width = 6.5, height = 9)
# plotfull
# dev.off()

## Save in appropriate size to fit on page
jpeg('Manuscript/Figure 1 Bias SCR.jpg', width=10.4, height = 12.33, units='cm',
     res=600)
plotfull
dev.off()

################################################################################
### repeat for CV, for Appendix ################################################

#############################################################################
### plot N ##################################################################

cNm<-ggplot(bias.Nmain, aes(x = Scenario.name, y = CV, fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  labs(x = NULL, y = "CV(N)") +
  coord_cartesian(ylim=c(15,60)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = 'none')#+ 


cNs<-ggplot(bias.Nsub, aes(x = Scenario.name, y = CV, fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim=c(15,60)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),        
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position = 'none')


#############################################################################
### plot beta #################################################################

cbm<-ggplot(bias.bmain, aes(x = Scenario.name, y = abs(CV), fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  labs(x = NULL, y = expression(paste("CV(", beta,')',sep=''))) +
  coord_cartesian(ylim=c(15,300)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = 'none')#+ 


cbs<-ggplot(bias.bsub, aes(x = Scenario.name, y = abs(CV), fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim=c(15,300)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),        
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position = 'none')


#############################################################################
### plot sigma ##################################################################

csm<-ggplot(bias.smain, aes(x = Scenario.name, y = abs(CV), fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  labs(x = 'Scenario', y = expression(paste("CV(", sigma,')', sep=''))) +
  coord_cartesian(ylim=c(8,35)) +
  theme_bw()+
  theme(legend.position = 'none')#+ 


css<-ggplot(bias.ssub, aes(x = Scenario.name, y = abs(CV), fill=Model)) +
  geom_boxplot(color = "black") + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  labs(x = 'Scenario', y = NULL) +
  coord_cartesian(ylim=c(8,35)) +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position = 'none')


###combine on a single page
plotfull.1b<-plot_grid(cNm, cbm, csm, align = 'v', ncol=1 , 
                      rel_heights = c(1,1,1.2))
plotfull.2b<-plot_grid(cNs, cbs, css, align = 'v', ncol=1 , 
                      rel_heights = c(1,1,1.2))

plotfullb<-plot_grid(plotfull.1b, plotfull.2b, ncol=2, nrow=1,
                    align='h', rel_widths=c(1, 0.89))

##suppl figure: jpeg, with room for caption
jpeg('Figure S1 CV SCR.jpg', width=15, height = 18, units='cm',
     res=600)
plotfullb
dev.off()


################################################################################
##### visualize distribution of Z at location for 2, 2a and 2b #################

## not for manuscript

# zz<-seq(-3,3,0.01)
# 
# ##symmetric, same variance
# plot(zz, dnorm(zz, 0, 0.5), type='l', xlab='', ylab='', lwd=2, col=coll[1],
#      axes=F)
# box()
# abline(v=0, lwd=2)
# ##symmetric, different variances
# sigvec<-c(0.2, 0.5, 1, 1.5)
# coll<-topo.colors(n=6)
# plot(zz, dnorm(zz, 0, sigvec[1]), type='l', xlab='', ylab='', lwd=2, col=coll[1],
#      axes=F)
# for (i in 2:length(sigvec)){
#   points(zz, dnorm(zz, 0, sigvec[i]), type='l', lwd=2, col=coll[i])
# }
# box()
# abline(v=0, lwd=2)
# 
# ##asymmetric, same variances
# alphavec<-c(-10, -2, 0, 2, 10)
# #rsn(K, xi=tempcov.mean[j],omega=0.5, alpha=sample(alphavec, 1))
# plot(zz, dsn(zz, xi=0, omega=1, alphavec[1]), type='l', xlab='', 
#      ylab='', lwd=2,col=coll[1], axes=FALSE)
# for (i in 2:length(alphavec)){
#   points(zz, dsn(zz, xi=0, omega=1, alphavec[i]), type='l', lwd=2, col=coll[i])
# }
# box()
# abline(v=0, lwd=2)
