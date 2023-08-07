################################################################################

## R code to summarize output from the SCR simulation study described in 
## Sollmann, R. (in review). 
## Mt or not Mt: Temporal variation in detection probability in spatial 
## capture-recapture and occupancy models. Peer Community Journal.

library(ggplot2)
library(writexl)
library(gridExtra)
library(cowplot)
library(ggpubr)

##############################################################################
######### constant settings across all scenarios #############################

niter=250
mag<-c(6,3)

X<-as.matrix(expand.grid(seq(0, 900, 100), seq(0, 900, 100)))
J<-nrow(X)

K<-8
beta.d <- 1
b0<-0.5
true.occu<-plogis(b0)

## create list of all scenarios
scen<-expand.grid((1:3), mag)

## read in results, data (for occupancy predictions)
res<-list()
dat<-list()
for (i in 1:nrow(scen)){
  res[[i]]<-readRDS(paste('out/Results_Scenario_',scen[i,1],'.', scen[i,2], '.rds', sep=''))
  dat[[i]]<-readRDS(paste('in/Data_Scenario_',scen[i,1],'.', scen[i,2], '.rds', sep=''))
}


## make data frame with estimates of all parameters under all scenarios
out.all<-data.frame(matrix(NA, nrow=0, ncol=7))
colnames(out.all)<-c('estimate', 'SE', 'Parameter', 'Model', 'Bias', 'CV', 'Scenario')
nconverged<-NULL

## matrices for difference in predictions
mean.d<-max.d<-nmax.d<-matrix(NA, nrow(scen), niter)

for (i in 1:nrow(scen)){
  sub<-res[[i]]
  
  ## extract estimates with SEs
  sub2<-lapply(sub[1:6], function(x)x[,1:2])
  xx<-as.data.frame(do.call(rbind, sub2))
  
  ##look at convergence of ALL models for given iteration
  ##but for scenarios >1, disregard M0; for 1, disregard mt2
  cmat<-sub$converged
  if(scen[i,1]>1){
    cmat[,1]<-NA
  } else {
    cmat[,3]<-NA
  }
  
  cvec<-apply(cmat, 1, function(x){any(x==0, na.rm=TRUE)})
  nconverged[i]<-sum(!cvec)
  
  ##if cvec==TRUE, not converged, ie, remove all params for all models
  ##for that iteration
  ##repeat convergence info 3x for the 3 models, then that twice for the 2 parameters
  cvec2<-rep(rep(cvec, 3), 2)
  xx[cvec2,1:2]<-NA
  
  ## make label for parameter and model
  xx$Parameter<-c(rep('occu', 3*niter), rep('beta', 3*niter))
  xx$Model<-rep(rep(c('M0', 'Mt', 'Mt2'), each=niter),2)
  
  ## calculate relative bias, CV
  xx$Bias<-(xx$estimate- c(rep(true.occu, 3*niter), rep(beta.d, 3*niter)))/
            c(rep(true.occu, 3*niter), rep(beta.d, 3*niter))
  xx$CV<-(xx$SE/xx$estimate)*100
  
  ## make label for scenario and magnitude
  xx$Scenario<-rep(scen[i,1], length(xx$CV))
  xx$Change<-scen[i,2]
  
  ## for sc 1 remove mt2, for others remove m0
  if (scen[i,1] == 1){
    xx<-xx[xx$Model!='Mt2',]
  }else {
    xx<-xx[xx$Model != 'M0',]
    
  }

  out.all<-rbind(out.all, xx)
  
  ## make occu predictions 
  ## first, get covariate info for all iterations, remove non-converged
  x.l<-lapply(dat[[i]], function(x)x$X.d)
  allx<-do.call(rbind,x.l)
  
  ## next, get intercept, coefficients, data generating (dg) and alternative (a)
  ## model
  if (scen[i,1] == 1){
    b0.dg<-log((xx$estimate[xx$Parameter == 'occu' & xx$Model == 'Mt'])/
      (1-xx$estimate[xx$Parameter == 'occu' & xx$Model == 'Mt']))
    b0.a<-log((xx$estimate[xx$Parameter == 'occu' & xx$Model == 'M0'])/
      (1-xx$estimate[xx$Parameter == 'occu' & xx$Model == 'M0']))
    b.dg<-xx$estimate[xx$Parameter == 'beta' & xx$Model == 'Mt']
    b.a<-xx$estimate[xx$Parameter == 'beta' & xx$Model == 'M0']
  }else {
    b0.dg<-log((xx$estimate[xx$Parameter == 'occu' & xx$Model == 'Mt'])/
      (1-xx$estimate[xx$Parameter == 'occu' & xx$Model == 'Mt']))
    b0.a<-log((xx$estimate[xx$Parameter == 'occu' & xx$Model == 'Mt2'])/
      (1-xx$estimate[xx$Parameter == 'occu' & xx$Model == 'Mt2']))
    b.dg<-xx$estimate[xx$Parameter == 'beta' & xx$Model == 'Mt']
    b.a<-xx$estimate[xx$Parameter == 'beta' & xx$Model == 'Mt2']
  }
  
  ## next calculate occu predictions, mean and max difference for each sim
  ## (only max used in ms) and number of locations with difference >0.1
  psi.dg<-plogis(matrix(b0.dg, niter, J) + matrix(b.dg, niter, J)*allx)
  psi.a<-plogis(matrix(b0.a, niter, J) + matrix(b.a, niter, J)*allx)
  
  mean.d[i,]<-apply(psi.a-psi.dg, 1, mean, na.rm=TRUE)
  max.d[i,]<-apply(psi.a-psi.dg, 1, function(x){
    unique(x[abs(x)==max(abs(x), na.rm=TRUE)])})
  nmax.d[i,]<-apply(psi.a-psi.dg, 1, function(x){
  sum(abs(x)>0.1)  })
}
## warnings can be ignored!


#############################################################################
### make summary table ######################################################

## calculate mean and median bias (median used in ms)
mean.bias<-aggregate(cbind(Bias, CV)~Parameter + Scenario+Model+Change, data=out.all, FUN=mean)
median.bias<-aggregate(cbind(Bias, CV)~Parameter + Scenario+Model+Change, data=out.all, FUN=median)

all.summ<-cbind(mean.bias, median.bias$Bias, median.bias$CV)


## get RMSE
all.summ$RMSE<-NA
truth<-rep(c(beta.d,true.occu), 12) 

for (jj in 1:24){
  ## get data pertaining to unique combo of parm, sc, model, fold
  sub<-out.all[out.all$Parameter == all.summ$Parameter[jj] &
               out.all$Scenario == all.summ$Scenario[jj] &
              out.all$Model == all.summ$Model[jj] &
                out.all$Change == all.summ$Change[jj] ,]
  ## calculate absolute bias
  abs.bias<-sub$Bias*truth[jj]
  ## calculate RMSE
  all.summ$RMSE[jj]<-sqrt(sum(abs.bias^2, na.rm=TRUE)/sum(!is.na(sub$Bias)))/truth[jj]
}

## order by scenario first, then parameter, then model
all.summ.df<-all.summ[with(all.summ, order(Change, Scenario, Parameter)),]
colnames(all.summ.df)<-c('Parameter', 'Scenario', 'Model','Change', 'Bias.mean',
                         'CV.mean','Bias.median', 'CV.median', 'RMSE')

## round
for (jj in 5:9){
  all.summ.df[,jj]<-round(all.summ.df[,jj], dig=2)
}

##change some labels in table to match ms

##models: M0, Ms, Mt or Mst 
##For scenarios >1
##Mt is data generating model, which should be labelled Mst (space and time)
##Mt2 is model ignoring temporal variation and should be labelled Ms
all.summ.df$Model[all.summ.df$Scenario>1 & all.summ.df$Model == 'Mt']<- 'Mst'
all.summ.df$Model[all.summ.df$Scenario>1 & all.summ.df$Model == 'Mt2']<- 'Ms'

##finally, switch order for scenario 1 so that data generating model always 
##comes first
all.summ.df<-all.summ.df[c(2,1,4,3, 5:12, 14,13,16,15, 17:nrow(all.summ.df)),]

##calculate median diff (max and avg) per scenario; mean number of locations
## with diff >0.1
all.summ.df$mean.diff<-NA
all.summ.df$max.diff<-NA
all.summ.df$n.max.diff<-NA
all.summ.df$mean.diff[seq(1, 24, by=4)]<-round(apply(mean.d, 1, median, na.rm=TRUE)
                                               [c(4:6, 1:3)], dig=2)
all.summ.df$max.diff[seq(1, 24, by=4)]<-round(apply(max.d, 1, median, na.rm=TRUE)
                                              [c(4:6, 1:3)], dig=2)
all.summ.df$n.max.diff[seq(1, 24, by=4)]<-round(apply(nmax.d, 1, mean, na.rm=TRUE)
                                              [c(4:6, 1:3)], dig=2)
#write_xlsx(all.summ.df, 'BiasOccu.xlsx')


############################################################################
#### make AICc table #######################################################

prop<-matrix(NA, 6,2)
colnames(prop)<-c('Data generating', 'Alternative')

for (i in 1:nrow(scen)){
  
  ## get results
  sub<-readRDS(paste('out/Results_Scenario_',scen[i,1],'.', scen[i,2], '.rds', sep=''))

  ## get AIC only
  sub2<-sub$AIC.mat
  
  ## get converged, disregard m0/mt2 where appropriate
  cmat<-sub$converged
  if(scen[i,1]>1){
    cmat[,1]<-NA
    sub2[,'m0']<-NA
  } else {
    cmat[,3]<-NA
    sub2[,'mt2']<-NA
  }
  
  cvec<-apply(cmat, 1, function(x){any(x==0, na.rm=TRUE)})
  sub2[cvec,]<-NA
  
  ## calculate deltaAIC
  if(scen[i,1]==1){
    d1<-sub2[,'mt']-sub2[,'m0']
    #if negative, mt is better; if more than 2 units, unequivocally better
    #if positive, m0 is better, if >2 etc etc
  } else {
    d1<-sub2[,'mt']-sub2[,'mt2']
  }
  
  ## et proportion top model 
  prop[i,]<-c( sum(d1< (-2), na.rm=TRUE)/nconverged[i],
               sum(d1> (2), na.rm=TRUE) /nconverged[i])
  
}

## combine with number of converged iterations
prop<-data.frame(round(prop, dig=2))
prop$N.converged<-nconverged
prop$Scenario<-c(paste(1:3, ' (', rep(c('3-fold', '6-fold'), each=3),')', sep='') )
#write_xlsx(prop, 'AIC.All.Scenarios.xlsx')


################################################################################
###### Make Figure 2: plot of bias in estimates ################################

##DISCLAIMER: This is most likely not the cleanest or most efficient way to make
##            and combine this series of plots. Code provided only for transpa-
##            rency, not as template for similar figures

## add correct model names
out.all$Model[out.all$Model == 'Mt' & out.all$Scenario >1]<-'Mst'
out.all$Model[out.all$Model == 'Mt2' & out.all$Scenario >1]<-'Ms'


## add indicator of data generating vs alternative model
out.all$ModelCode<-out.all$Model

out.all$Model[out.all$Scenario == 1 & out.all$ModelCode == 'Mt']<-'D.g.'
out.all$Model[out.all$Scenario == 1 & out.all$ModelCode == 'M0']<-'Alt.'

out.all$Model[out.all$Scenario > 1 & out.all$ModelCode == 'Mst']<-'D.g.'
out.all$Model[out.all$Scenario > 1 & out.all$ModelCode == 'Ms']<-'Alt.'

## change level order 
out.all$Model<-as.factor(out.all$Model)
out.all$Model<-relevel(out.all$Model, 'D.g.')

out.all$Scenario<-as.factor(out.all$Scenario)


## Strategy: per parameter, plot bias for 3 main scenarios and 2 magnitudes
##           then, arrange on one page

v_line <- data.frame(
  yintercept = c(0,0,0,0)
)

## set colors for data generating and alternative models
colr<-c('grey80', 'grey50')

## set width for median line and outline of boxplots
lwd.m<-1


###############################################################################
#### plot occu ################################################################

## manitude 3

bias.Occu.low<-out.all[out.all$Parameter == 'occu' & out.all$Change == 3, ]

pol<-ggplot(bias.Occu.low, aes(x = Scenario, y = Bias*100, fill=Model,
                               color=Model)) +
  geom_boxplot(fatten = lwd.m) + #color = "black"
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + 
  scale_color_manual(values = c(`D.g.` =colr[2], 
                                `Alt.` = "black"))+
  labs(x = NULL, y = expression(paste("Relative bias  ", bar(psi), sep=''))) +
  scale_x_discrete(breaks=c("1", "2", "3"),
                     labels=c("1 (low)", "2 (low)", "3 (low)"))+
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-50,60)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = 'none')

## manitude 6
bias.Occu.high<-out.all[out.all$Parameter == 'occu' & out.all$Change == 6, ]

poh<-ggplot(bias.Occu.high, aes(x = Scenario, y = Bias*100, fill=Model,
                                color=Model)) +
  geom_boxplot(fatten = lwd.m) + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + 
  scale_color_manual(values = c(`D.g.` =colr[2], 
                                `Alt.` = "black"))+
  labs(x = NULL, y = NULL) +
  scale_x_discrete(breaks=c("1", "2", "3"),
                   labels=c("1 (high)", "2 (high)", "3 (high)"))+
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-50,60)) +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = 'none')


#### beta #####################################################################

## manitude 3
bias.b.low<-out.all[out.all$Parameter == 'beta' & out.all$Change == 3, ]

pbl<-ggplot(bias.b.low, aes(x = Scenario, y = Bias*100, fill=Model,
                            color=Model)) +
  geom_boxplot(fatten = lwd.m) + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  scale_color_manual(values = c(`D.g.` =colr[2], 
                                `Alt.` = "black"))+
  labs(x = NULL, y = expression(paste("Relative bias ", beta, sep=''))) +
  scale_x_discrete(breaks=c("1", "2", "3"),
                   labels=c("1 (low)", "2 (low)", "3 (low)"))+
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-120,400)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = 'none')

## manitude 6
bias.b.high<-out.all[out.all$Parameter == 'beta' & out.all$Change == 6, ]

pbh<-ggplot(bias.b.high, aes(x = Scenario, y = Bias*100, fill=Model,
                             color=Model)) +
  geom_boxplot(fatten = lwd.m) + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  scale_color_manual(values = c(`D.g.` =colr[2], 
                                `Alt.` = "black"))+
  labs(x = NULL, y = NULL) +
  scale_x_discrete(breaks=c("1", "2", "3"),
                   labels=c("1 (high)", "2 (high)", "3 (high)"))+
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-120,400)) +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = 'none')


### max difference ###########################################################
scena<-as.factor(rep(scen[,1], 250))
magn<-as.factor(rep(c('high', 'high', 'high', 'low', 'low', 'low'), 250))
magn<-relevel(magn, 'low')

diff.df<-data.frame(mean.d=c(mean.d),
                    Scenario=scena,
                    Magnitude=magn)
## manitude 3
ml<-ggplot(diff.df[diff.df$Magnitude == 'low',], 
           aes(x = Scenario, y = mean.d)) +
  geom_boxplot(fatten = lwd.m) + 
  labs(x = 'Scenario', 
       y = expression(paste("Difference ", widehat(psi), sep=''))) +
  scale_x_discrete(breaks=c("1", "2", "3"),
                   labels=c("1 (low)", "2 (low)", "3 (low)"))+
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-0.1,0.25)) +
  theme_bw()

## manitude 6
mh<-ggplot(diff.df[diff.df$Magnitude == 'high',], 
           aes(x = Scenario, y = mean.d)) +
  geom_boxplot(fatten = lwd.m) + 
  labs(x = 'Scenario', y = NULL) +
  scale_x_discrete(breaks=c("1", "2", "3"),
                   labels=c("1 (high)", "2 (high)", "3 (high)")) +
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', linewidth=0.5)+
  coord_cartesian(ylim=c(-0.1,0.25)) +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank())


## align, first by column, then both columns
pfb1<-plot_grid(pol, pbl, ml, align = 'v', ncol=1 , 
                       rel_heights = c(1,1,1.2))
pfb2<-plot_grid(poh, pbh, mh, align = 'v', ncol=1 , 
                       rel_heights = c(1,1,1.2))

pfb<-plot_grid(pfb1, pfb2, ncol=2, nrow=1,
                     align='h', rel_widths=c(1, 0.78))

##figure: letter size with 1-in margin, pdf or eps
# pdf('Figure 2 Bias Occu.pdf', paper='letter', width = 6.5, height = 9)
# pfb
# dev.off()

## Save in appropriate size to fit on page
jpeg('Figure 2 Bias Occu.jpg', width=10.4, height = 12.33, units='cm',
     res=600)
pfb
dev.off()

################################################################################
### repeat for CV, for Appendix ################################################

### occu ###################################################################
col<-ggplot(bias.Occu.low, aes(x = Scenario, y = abs(CV), fill=Model,
                               color=Model)) +
  geom_boxplot(fatten = lwd.m) + #color = "black"
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + 
  scale_color_manual(values = c(`D.g.` =colr[2], 
                                `Alt.` = "black"))+
  labs(x = NULL, y = expression(paste("CV(", bar(psi),')', sep=''))) +
  coord_cartesian(ylim=c(0,35)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none')


coh<-ggplot(bias.Occu.high, aes(x = Scenario, y = abs(CV), fill=Model,
                                color=Model)) +
  geom_boxplot(fatten = lwd.m) + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + 
  scale_color_manual(values = c(`D.g.` =colr[2], 
                                `Alt.` = "black"))+
  labs(x = NULL, y = NULL) +
  coord_cartesian(ylim=c(0,35)) +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = 'none')


#### beta #####################################################################

cbl<-ggplot(bias.b.low, aes(x = Scenario, y = abs(CV), fill=Model,
                            color=Model)) +
  geom_boxplot(fatten = lwd.m) + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  scale_color_manual(values = c(`D.g.` =colr[2], 
                                `Alt.` = "black"))+
  labs(x = 'Scenario', y = expression(paste("CV(", beta,')', sep=''))) +
  scale_x_discrete(breaks=c("1", "2", "3"),
                   labels=c("1 (low)", "2 (low)", "3 (low)")) +
  coord_cartesian(ylim=c(0,400)) +
  theme_bw()+
  theme(legend.position = 'none')


cbh<-ggplot(bias.b.high, aes(x = Scenario, y = abs(CV), fill=Model,
                             color=Model)) +
  geom_boxplot(fatten = lwd.m) + 
  scale_fill_manual(values = c(`D.g.` = colr[1],
                               `Alt.` = colr[2])) + # alpha = .3,
  scale_color_manual(values = c(`D.g.` =colr[2], 
                                `Alt.` = "black"))+
  labs(x = 'Scenario', y = NULL) +
  scale_x_discrete(breaks=c("1", "2", "3"),
                   labels=c("1 (high)", "2 (high)", "3 (high)")) +
  coord_cartesian(ylim=c(0,400)) +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = 'none')


### combine on single page
pfc1<-plot_grid(col, cbl, align = 'v', ncol=1 , 
                rel_heights = c(1,1.1))
pfc2<-plot_grid(coh, cbh, align = 'v', ncol=1 , 
                rel_heights = c(1,1.1))

pfc<-plot_grid(pfc1, pfc2, ncol=2, nrow=1,
               align='h', rel_widths=c(1, 0.89))


##suppl figure: jpeg
jpeg('Figure S2 CV Occu.jpg', width=6.5, height = 6.5, units='in',
    res=600)
pfc
dev.off()

