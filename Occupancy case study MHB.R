################################################################################

## R code to implement occupancy case study described in Sollmann, R. (in review). 
## Mt or not Mt: Temporal variation in detection probability in spatial 
## capture-recapture and occupancy models. Peer Community Journal.

library(AHMbook)
library(unmarked)
library(writexl)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(lemon)

################################################################################
## load common bird data Switzerland
data("MHB2014")

## extract and scale site covariatess
covs.raw<-MHB2014$sites[,-1]
covs<-data.frame(scale(covs.raw))

## extract and scale detection covariates
date<-MHB2014$date
date<-(date-mean(date, na.rm=TRUE))/sd(date, na.rm=TRUE)
dur<-MHB2014$dur
dur<-(dur-mean(dur, na.rm=TRUE))/sd(dur, na.rm=TRUE)

## calculate location-average detection covariates
covs$avg.date<-scale(apply(date, 1, mean, na.rm=TRUE))
covs$avg.dur<-scale(apply(dur, 1, mean, na.rm=TRUE))

## extract observations, covert to binary
obs<-MHB2014$counts
obs[obs>1]<-1

################################################################################
## Step 1: identify species suitable for analysis ##############################

## Select species with 50+ detections
n<-apply(obs,3,sum, na.rm=TRUE)
inn<-which(n>50)
length(inn)

## Based on target model, identify species with sufficient temporal variation 
## in p 
## Here: ratio of max:min p is (on average across sampling locations) >2
## Also check on p-value of detection covariates

pval<-dif<-NULL

for (ii in 1:length(inn)){
y<-obs[,,inn[ii]]

umf<-unmarkedFrameOccu(y=y, siteCovs = covs,
                       obsCovs = list(date=date, dur=dur))
m.full<-occu(~date+ I(date^2)+dur~elev + I(elev^2)+forest, data=umf)

log <- capture.output({
  res <- summary(m.full@estimates)$det;
})
log <- capture.output({
  res2 <- summary(m.full@estimates)$state;
})

## filter out extreme or singular model fits
if(any(is.na(res))|any(res[,c(1,2)]>10)|
   any(is.na(res2))|any(res2[,c(1,2)]>10)) next

pval[ii]<-min(res$`P(>|z|)`[2:3])

p.pred<-predict(m.full, 'det')

pmat<-matrix(p.pred[,1], nrow(date)-1, ncol(date), byrow=TRUE)
dp<-apply(pmat, 1, function(x){
  max(x)/min(x)
})
dif[ii]<-mean(dp, na.rm=TRUE)

}

## check how many species have a max:min p ratio >2
which(dif>2)
##ten with dif>2

###############################################################################
### For selected species, fit models ##########################################

keep<-which(dif>2)

## create objects to store model results
delta<-NULL #dAIC Ms to Mst

## occupancy predictions at sampling locations
psi.pred<-array(NA, c(length(keep),nrow(covs)-1,4, 2))
dimnames(psi.pred)<-list(names(inn[keep]),
                         NULL,
                         c('Predicted', 'SE', 'lower', 'upper'),
                         c('M.st', 'M.s'))

## coefficients, occupancy
coef.mat.psi<-array(NA, c(length(keep),4,6, 2)) #coefficients, psi, p
dimnames(coef.mat.psi)<-list(names(inn[keep]),
                            c('(Intercept)', 'elev','I(elev^2)','forest'),
                            c('Estimate', 'SE', 'z', 'p', 'lower', 'upper'),
                            c('M.st', 'M.s'))

## coefficients, p
coef.mat.p<-array(NA, c(length(keep),4,6, 2))
dimnames(coef.mat.p)<-list(names(inn[keep]),
                           c('(Intercept)', 'date', 
                             'I(date^2)','dur'),
                           c('Estimate', 'SE', 'z', 'p', 'lower', 'upper'),
                           c('M.st', 'M.s'))

## for response curves (occupancy and detection)
## create ranges of covariates and compile in data frames
pfor<-seq(min(covs$forest), max(covs$forest), length.out=100)
pele<-seq(min(covs$elev), max(covs$elev), length.out=100)
pdate<-seq(min(date, na.rm=TRUE), max(date, na.rm=TRUE), length.out=100)
pdur<-seq(min(dur, na.rm=TRUE), max(dur, na.rm=TRUE), length.out=100)
pav.date<-seq(min(covs$avg.date, na.rm=TRUE), max(covs$avg.date, na.rm=TRUE), length.out=100)
pav.dur<-seq(min(covs$avg.dur, na.rm=TRUE), max(covs$avg.dur, na.rm=TRUE), length.out=100)

newdat.for<-data.frame(forest=pfor, elev=rep(0, length(pfor)))
newdat.ele<-data.frame(elev=pele, forest=rep(0, length(pfor)))

newdat.date<-data.frame(date=pdate, dur=rep(0, length(pdate)))
newdat.dur<-data.frame(dur=pdur, date=rep(0, length(pdate)))

newdat.avgdate<-data.frame(avg.date=pav.date, avg.dur=rep(0, length(pav.date)))
newdat.avgdur<-data.frame(avg.dur=pav.dur, avg.date=rep(0, length(pav.date)))

## objects to hold response curve calculations
psi.response.for<-psi.response.ele<-array(NA, c(length(keep), 
                                                nrow(newdat.ele), 4, 2))
dimnames(psi.response.for)<-dimnames(psi.response.ele)<-list(c(names(inn[keep])),
                             NULL,
                             c('Predicted', 'SE', 'lower', 'upper'),
                             c('M.st', 'M.s'))

p.response.date<-p.response.dur<-array(NA, c(length(keep), 
                                              nrow(newdat.date), 4, 2))
dimnames(p.response.date)<-dimnames(p.response.dur)<-list(c(names(inn[keep])),
                                                             NULL,
                                                             c('Predicted', 'SE', 'lower', 'upper'),
                                                             c('M.st', 'M.s'))
## track number of sites with detections
nsites<-NULL

for (jj in 1:length(keep)){

  ## get detections
  y<-obs[,,inn[keep[jj]]]

  ## track number of sites with detections
  nsites[jj]<-sum(apply(y,1,sum, na.rm=TRUE)>0)

  ## compile everything in unmarked frame
  umf<-unmarkedFrameOccu(y=y, siteCovs = covs,
                       obsCovs = list(date=date, dur=dur))
  
  ## fit models with (m.full) and without (m.s) temporal variation in p
  m.full<-occu(~date+I(date^2)+dur~elev + I(elev^2)+forest, data=umf)
  m.s<-occu(~avg.date+I(avg.date^2)+avg.dur~elev + I(elev^2)+forest, data=umf)

  ## get delta AIC
  delta[jj]<-m.s@AIC-m.full@AIC

  ## make predictions at sampled locations
  psi.pred[jj,,,1]<-as.matrix(predict(m.full, 'state'))
  psi.pred[jj,,,2]<-as.matrix(predict(m.s, 'state'))

  ## extract coefficient estimates
  log <- capture.output({
    res.f <- summary(m.full@estimates)$state;#[,4];
  })
  log <- capture.output({
    res.s <- summary(m.s@estimates)$state;#[,4];
  })
  log <- capture.output({
    res.f.p <- summary(m.full@estimates)$det;#[,4];
  })
  log <- capture.output({
    res.s.p <- summary(m.s@estimates)$det;#[,4];
  })

  coef.mat.psi[jj,,1:4,1]<-as.matrix(res.f)
  coef.mat.psi[jj,,1:4,2]<-as.matrix(res.s)
  coef.mat.psi[jj,,5:6,1]<-as.matrix(confint(m.full, type='state'))
  coef.mat.psi[jj,,5:6,2]<-as.matrix(confint(m.s, type='state'))

  coef.mat.p[jj,,1:4,1]<-as.matrix(res.f.p)
  coef.mat.p[jj,,1:4,2]<-as.matrix(res.s.p)
  coef.mat.p[jj,,5:6,1]<-as.matrix(confint(m.full, type='det'))
  coef.mat.p[jj,,5:6,2]<-as.matrix(confint(m.s, type='det'))

  ## make predictions for response curve, occupancy
  psi.response.for[jj,,,1]<-as.matrix(predict(m.full, type='state', 
                                        newdata=newdat.for))
  psi.response.for[jj,,,2]<-as.matrix(predict(m.s, type='state', 
                                        newdata=newdat.for))
  psi.response.ele[jj,,,1]<-as.matrix(predict(m.full, type='state', 
                                            newdata=newdat.ele))
  psi.response.ele[jj,,,2]<-as.matrix(predict(m.s, type='state', 
                                            newdata=newdat.ele))

  ## make predictions for response curve, p
  p.response.date[jj,,,1]<-as.matrix(predict(m.full, type='det', 
                                            newdata=newdat.date))
  p.response.date[jj,,,2]<-as.matrix(predict(m.s, type='det', 
                                            newdata=newdat.avgdate))
  p.response.dur[jj,,,1]<-as.matrix(predict(m.full, type='det', 
                                            newdata=newdat.dur))
  p.response.dur[jj,,,2]<-as.matrix(predict(m.s, type='det', 
                                            newdata=newdat.avgdur))

} #end species loop


## compile parameter estimates in tables

diff.list<-p.list<-list()

for (ii in 1:length(keep)){
  
  ## calculate differences in coefficients for psi
  dd<-round(((coef.mat.psi[ii,,1,2]-coef.mat.psi[ii,,1,1])/
            coef.mat.psi[ii,,1,1]), dig=2)

  ## compile occupancy estimates
  diff.list[[ii]]<-data.frame(Parameter=names(coef.mat.psi[ii,,1,1]),
                            `Estimate(M.st)`=round(coef.mat.psi[ii,,1,1],dig=2),
                            `SE(M.st)`=round(coef.mat.psi[ii,,2,1],dig=2),
                            `lower(M.st)`=round(coef.mat.psi[ii,,5,1],dig=2),
                            `upper(M.st)`=round(coef.mat.psi[ii,,6,1],dig=2),
                            `Estimate(M.s)`=round(coef.mat.psi[ii,,1,2],dig=2),
                            `SE(M.s)`=round(coef.mat.psi[ii,,2,2],dig=2),
                            `lower(M.s)`=round(coef.mat.psi[ii,,5,2],dig=2),
                            `upper(M.s)`=round(coef.mat.psi[ii,,6,2],dig=2),
                            Rel.diff=dd,
                            Coverage=coef.mat.psi[ii,,1,2]>coef.mat.psi[ii,,5,1]&
                              coef.mat.psi[ii,,1,2]<coef.mat.psi[ii,,6,1], 
                            check.names = FALSE)

  ## compile detection estimates
  p.list[[ii]]<-data.frame(Parameter=names(coef.mat.p[ii,,1,1]),
                            `Estimate(M.st)`=round(coef.mat.p[ii,,1,1],dig=2),
                            `SE(M.st)`=round(coef.mat.p[ii,,2,1],dig=2),
                            `lower(M.st)`=round(coef.mat.p[ii,,5,1],dig=2),
                            `upper(M.st)`=round(coef.mat.p[ii,,6,1],dig=2),
                            `Estimate(M.s)`=round(coef.mat.p[ii,,1,2],dig=2),
                            `SE(M.s)`=round(coef.mat.p[ii,,2,2],dig=2),
                            `lower(M.s)`=round(coef.mat.p[ii,,5,2],dig=2),
                            `upper(M.s)`=round(coef.mat.p[ii,,6,2],dig=2),
                            check.names = FALSE)
}

## add species names
names(diff.list)<-names(p.list)<-names(inn[keep])

## write out tables
write_xlsx(diff.list, 'OccuExampleDiffCoef.xlsx')
write_xlsx(p.list, 'OccuExampleCoefP.xlsx')


##also write out as one long table for Appendix

full.frame<-data.frame(Species=rep(names(inn[keep]),each=nrow(diff.list[[1]])), 
                       do.call(rbind, diff.list))
write_xlsx(full.frame, 'OccuCoefAppendix.xlsx')

full.frame.p<-data.frame(Species=rep(names(inn[keep]),each=nrow(p.list[[1]])), 
                       do.call(rbind, p.list))
write_xlsx(full.frame.p, 'PCoefAppendix.xlsx')

##get mean and max difference in predictions

max.d<-round(apply(psi.pred, 1, function(x){
  xx<-x[,1,2]-x[,1,1]
  unique(xx[abs(xx)==max(abs(xx))])
}), dig=2)

avg.d<-round(apply(psi.pred, 1, function(x){
 mean(x[,1,2]-x[,1,1])
}), dig=2)


#double check that predictions under M.s are inside confint of M.st
apply(psi.pred, 1, function(x){
  sum(x[,1,2]>x[,4,1])
})


##make table with species specific results (for 2 App tables)
spec.res<-data.frame(Species=names(inn[keep]),
                     `p.min:p.max`=round(na.omit(dif[dif>2]), dig=2),
                     D.mean=avg.d,
                     D.max=max.d,
                     dAIC=round(delta, dig=2),
                     n=n[inn[keep]],
                     psi=round(nsites/(nrow(covs)-1), dig=2))
write_xlsx(spec.res, 'OccuSpeciesResults.xlsx')


################################################################################
###### make plots of predictions ###############################################

##ignore warnings, some arrows have 0 length because predictions are 0

## compile all necessary output in one data frame

for (ii in 1:length(keep)){
  if (ii ==1){
    df<-data.frame(pred.st=psi.pred[ii,,1,1],
                   pred.s=psi.pred[ii,,1,2],
                   lower = psi.pred[ii,,3,1], 
                   upper = psi.pred[ii,,4,1],
                   spec=names(inn[keep])[ii])
  }else{
    df2<-data.frame(pred.st=psi.pred[ii,,1,1],
                    pred.s=psi.pred[ii,,1,2],
                    lower = psi.pred[ii,,3,1], 
                    upper = psi.pred[ii,,4,1],
                    spec=names(inn[keep])[ii])
    df<-rbind(df, df2)
  }
  
}

## make panel plot of predictions under Ms against predictions under Mst,
## with 95% CI for Mst; if line of equivalence goes through CI, difference in
## predictions is not statistically significant

p.pred<-ggplot(df, aes(y=pred.s, x=pred.st))+
  geom_errorbar(aes(xmin=lower, xmax=upper), width=0.05, color='lightgrey')+
  labs(x = expression(paste(psi, ' under ', 'M'[st]), sep=''), 
       y = expression(paste(psi, ' under ', 'M'[s]), sep='')) +
  geom_abline(intercept=0, slope=1, color = 'red')+
  geom_point(size=0.5)+
  facet_wrap(~spec, nrow=4, ncol=, scales="free")+
  theme_bw()

### Make Appendix 3 Figure S1
jpeg('PredictionCorrelationsFigS1.jpg', width=6.5, height = 8, units='in',
     res=600)
p.pred
dev.off()



################################################################################
###### make plots of response curves for psi ###################################

## set dimensions, colors
Tt<-dim(psi.response.ele)[2]
col.for<-c('red', 'lightcoral')
col.ele<-c('black', 'lightgrey')

##compile data for all species
for(ii in 1:length(keep)){
  
  if (ii ==1){
    resp.df<-data.frame(resp.ele=psi.response.ele[ii,,1,1],
                        resp.ele.a=psi.response.ele[ii,,1,2],
                        resp.ele.l=psi.response.ele[ii,,3,1],
                        resp.ele.h=psi.response.ele[ii,,4,1],
                        ele=pele,
                        resp.for=psi.response.for[ii,,1,1],
                        resp.for.a=psi.response.for[ii,,1,2],
                        resp.for.l=psi.response.for[ii,,3,1],
                        resp.for.h=psi.response.for[ii,,4,1],
                        forr=pfor,
                        spec = names(inn[keep])[ii])
  } else {
    resp.df2<-data.frame(resp.ele=psi.response.ele[ii,,1,1],
                         resp.ele.a=psi.response.ele[ii,,1,2],
                         resp.ele.l=psi.response.ele[ii,,3,1],
                         resp.ele.h=psi.response.ele[ii,,4,1],
                         ele=pele,
                         resp.for=psi.response.for[ii,,1,1],
                         resp.for.a=psi.response.for[ii,,1,2],
                         resp.for.l=psi.response.for[ii,,3,1],
                         resp.for.h=psi.response.for[ii,,4,1],
                         forr=pfor,
                         spec = names(inn[keep])[ii])
    resp.df<-rbind(resp.df, resp.df2)
  }
}

## plot response curves under Mst with 95%CI, and add response curves under Ms
## for both predictors

##use facet.wrap
pp.x<-ggplot(resp.df, aes(y=resp.ele, x=ele))+
  geom_line(aes(color = 'Elevation', linetype = 'Mst'))+
  geom_line(aes(y=resp.for, x=forr,color='%Forest', linetype = 'Mst'))+
  coord_cartesian(ylim = c(0, 1))+
  labs(x = 'Predictor (scaled)',
       y = 'Occupancy') +
  theme_bw()+
  geom_ribbon(aes(ymin=resp.ele.l, ymax=resp.ele.h),
              alpha=0.3, color=col.ele[2], fill=col.ele[2])+
  geom_ribbon(aes(ymin=resp.for.l, ymax=resp.for.h),
              alpha=0.3, color=col.for[2], fill=col.for[2])+
  geom_line(aes(y=resp.for.a, x=forr, linetype='Ms'), color=col.for[1])+
  geom_line(aes(y=resp.ele.a, x=ele, linetype='Ms'), color=col.ele[1])+
  theme(axis.title = element_text(size = 9))+
  scale_color_manual(name='Predictor',
                     breaks=c('Elevation', '%Forest'),
                     values=c('Elevation'='black', '%Forest'="red"))+
  scale_linetype_manual(name='Model',
                        values = c('Mst'='solid', 
                                   'Ms'='dashed'))+
  theme(legend.direction = "vertical", legend.box = "horizontal",
        axis.title = element_text(size = 9), legend.position = 'right')+
  facet_wrap(~spec, nrow=4, ncol=3)

##find empty panels
pnls <- cowplot::plot_to_gtable(pp.x) %>% gtable::gtable_filter("panel") %>%
  with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))

###write out with repositioned legend - full page figure
pdf('Figure 3 Response curves.pdf', paper='letter', width=6.5, height = 9,
    onefile=FALSE)
lemon::reposition_legend( pp.x, "center", panel=names(pnls) )
dev.off()

## write out with repositioned legend, adjust size to fit page


## Save in appropriate size to fit on page
jpeg('Figure 3 Response curves.jpg', width=15, height = 18, units='cm',
     res=600)
lemon::reposition_legend( pp.x, "center", panel=names(pnls) )
dev.off()


################################################################################
###### Repeat for p for Appendix ###############################################

Tt2<-dim(p.response.dur)[2]

for(ii in 1:length(keep)){
  
  ##make data frame with all necessary info
  if(ii == 1){
  resp.df<-data.frame(resp.date=p.response.date[ii,,1,1],
                      resp.date.a=p.response.date[ii,,1,2],
                      resp.date.l=p.response.date[ii,,3,1],
                      resp.date.h=p.response.date[ii,,4,1],
                      date=pdate,
                      resp.dur=p.response.dur[ii,,1,1],
                      resp.dur.a=p.response.dur[ii,,1,2],
                      resp.dur.l=p.response.dur[ii,,3,1],
                      resp.dur.h=p.response.dur[ii,,4,1],
                      dur=pdur,
                      spec = names(inn[keep])[ii])
  }else {
    resp.df2<-data.frame(resp.date=p.response.date[ii,,1,1],
                        resp.date.a=p.response.date[ii,,1,2],
                        resp.date.l=p.response.date[ii,,3,1],
                        resp.date.h=p.response.date[ii,,4,1],
                        date=pdate,
                        resp.dur=p.response.dur[ii,,1,1],
                        resp.dur.a=p.response.dur[ii,,1,2],
                        resp.dur.l=p.response.dur[ii,,3,1],
                        resp.dur.h=p.response.dur[ii,,4,1],
                        dur=pdur,
                        spec = names(inn[keep])[ii])
    resp.df<-rbind(resp.df, resp.df2)
  }
}

    ppp<-ggplot(resp.df, aes(y=resp.date, x=date))+
      geom_line(aes(color = 'Date', linetype = 'Mst'))+
      geom_line(aes(y=resp.dur, x=dur,color='Duration', linetype = 'Mst'))+
      coord_cartesian(ylim = c(0, 1))+
      labs(x = 'Predictor (scaled)', 
           y = 'Detection') +
      theme_bw()+
      geom_ribbon(aes(ymin=resp.date.l, ymax=resp.date.h),
                  alpha=0.3, color=col.ele[2], fill=col.ele[2])+
      geom_ribbon(aes(ymin=resp.dur.l, ymax=resp.dur.h),
                  alpha=0.3, color=col.for[2], fill=col.for[2])+
      geom_line(aes(y=resp.dur.a, x=dur, linetype='Ms'), color=col.for[1])+
      geom_line(aes(y=resp.date.a, x=date, linetype='Ms'), color=col.ele[1])+
      scale_color_manual(name='Predictor',
                         breaks=c('Date', 'Duration'),
                         values=c('Date'='black', 'Duration'="red"))+
      scale_linetype_manual(name='Model',
                            values = c('Mst'='solid', 
                                       'Ms'='dashed'))+
      theme(legend.direction = "vertical", legend.box = "horizontal",
            axis.title = element_text(size = 9))+
      facet_wrap(~spec, nrow=4, ncol=3)

    ##find empty panels
    pnls2 <- cowplot::plot_to_gtable(ppp) %>% gtable::gtable_filter("panel") %>%
      with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))
    
jpeg('Figure A3S1 Detection curves.jpg', width=6.5, height = 7.5, units='in',
     res=600)
lemon::reposition_legend( ppp, "center", panel=names(pnls2) )

dev.off()
