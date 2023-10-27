################################################################################

## R code to format data for SCR case study described in Sollmann, R. (in review). 
## Mt or not Mt: Temporal variation in detection probability in spatial 
## capture-recapture and occupancy models. Peer Community Journal.

## Files CapturesScent.csv and SurveysScent.csv (loaded in the code below)
## can be freely downloaded from 
## https://www.sciencebase.gov/catalog/item/6181a72bd34e9f2789e44877
## Place in data-raw folder of R Project or adjust path when reading in files

## Data citation:
## Amburgey SM, Lardner B, Knox A, Converse SJ, Yackel Adams AA (2021) 
## Brown Treesnake detections on transects using potential attractants of 
## live-mouse lures or fish-spray scent, Guam. U.S. Geological Survey 
## data release. http://doi.org/10.5066/P9G6JHZ3
## Please cite original data source if using the processed data.

## The script builds on code and makes use of a custom data prep function obtained from
## github.com/Quantitative-Conservation-Lab/Amburgey_etal_2021_NeoBiota/tree/main


rm(list = ls())


## source function to format data from github
## note: this will load a series of R packages that may have to be installed
## for the function to work

source('https://raw.githubusercontent.com/Quantitative-Conservation-Lab/Amburgey_etal_2021_NeoBiota/main/Data/PrepDataScent.R')


# Read in capture data and survey data and keep only relevant columns
# See metadata in ScienceBase repository in README
caps <- read_csv("data-raw/CapturesScent.csv")[,c("EFFORTID","Date","PITTAG","SVL","TOTAL","WEIGHT","SEX","TRANSECT","LOCATION")]
survs <- read_csv("data-raw/SurveysScent.csv")
# Restrict to 2-month period to meet assumptions of population closure
caps <- caps %>%
  filter(!grepl('Oct', Date)) %>%
  filter(!grepl('Jan', Date))
survs <- survs %>%
  select(-contains("Oct")) %>%
  select(-contains("Jan"))

# Format captures and survey effort for traditional SCR analysis
dat <- PrepDat(caps,survs)

# Create trapping grid of Closed Population (CP) dimensions (5 ha, 50,000 m2)
locs <- secr::make.grid(nx = 13, ny = 27, spacex = 16, spacey = 8)

## Number of locations (lure/no lure points)
J <- nrow(locs)

## format for use in secr
locs<-data.frame(trapID=1:J, x=locs[,1], y=locs[,2])


# Status of surveys (1 = not active, 2 = active but no scent, 
##                   3 = active and fresh scent, 4 = active and old scent)
stat <- as.matrix(dat$scent)
colnames(stat)<-NULL

## Change to no scent/old scent vs fresh scent
## (based on original analysis that suggests differences in detection between
## these two categories but not no scent vs old scent)

## set inactive status to NA; inactive site-occasion combinations are not used 
## in the analysis (indicated to secr via 'usage' argument)
stat[stat==1]<-NA
stat[stat==3]<-1
stat[stat %in% c(2,4)]<-0


## format detection data for secr
cap.df<-data.frame(session=NA, 
                   ID=NA, 
                   occasion=NA,
                   trap=NA)

##individual maximum #captures in trap and occasion ==1 --> binary data

for ( i in 1:dim(dat$snks)[1]){
  
  ## subset to data for individual i
  sub<-dat$snks[i,,]
  
  ## for each trap, identify occasions with captures
  test<-apply(sub, 1, function(x)which(x>0))
  test<-test[which(sapply(test, length)>0)]
  
  ##get traps
  trpid<-pmatch(rep(names(test), sapply(test, length)), rownames(sub), duplicates.ok = TRUE)
  
  ##get occasions
  occid<-unlist(test)
  
  ##compile in data frame
  new.df<-data.frame(session = rep(1, length(trpid)),
                     ID= rep(i, length(trpid)),
                     occasion = occid,
                     trap = trpid)
  cap.df<-rbind(cap.df, new.df)
  
}
## remove starter row from cap.df
cap.df<-cap.df[-1, ]


##write out data for secr analysis
dat.list<-list(locs=locs, act=as.matrix(dat$act),
               stat=stat, cap.df=cap.df)
write_rds(dat.list, 'data-raw/Snake_secr_data.rds')
