## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(readxl,dplyr,tidyverse,optmatch,lpSolve,
               randomForest,survival,ggsurvfit,survminer,MatchIt,boot)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions")
source("ilp_matching.R")
source("alignment_psm_matching_c.R")
source("exact_matching.R")
source("boot_matching.R")
source("risk_ratio.R")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
pop=readRDS("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\paper\\main_analysis\\with duration\\time_to_event\\population_rr_0.5.xlsx")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
trt=pop %>% filter(trt==1) 
untrt=pop %>% filter(trt==0) 


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(03111963)
dftrial=rbind(trt[sample(1:nrow(trt),20),],
            untrt[sample(1:nrow(untrt),50),])


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=matching(data=dftrial,treatment=trt,ageofonset=ageonset,sex=sex,
                 duration=duration,time=time,method='optimal.ilp',
                 caliper=0.5,step=0.125,maxcaliper=0.5,ratio=2,
                 replacement=F,history=F,warning=T)
# extract dataset
matcheddf=matched$dataset
# see optimal pairs
#matched$optimalres


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(time_trial,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,
              treatment="trt",estimand="ATT",landmark=730)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       full_data=dftrial,
       estimand="ATT",
       statistic=boot_matching,
       estimate="surv_rr",
       landmark=730,
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(time_trial,event)~trt+calendar,
       duration="duration",time="time",
       method="exact",
       search=list(exact="optimal.ilp"))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0[1]
lower=bate-1.96*sd(B$t[,1])
upper=bate+1.96*sd(B$t[,2])
# converged proportion?
sum(B$t[,2])/length(B$t[,2])
# check all
B$t


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=psm_matching(data=dftrial,
                     formula=trt~ageonset+sex,
                     ratio=2,
                     replacement=F,
                     distance='glm',
                     method='optimal.ilp',
                     #tune=TRUE,
                     #caliper=0.1,
                     duration=duration,
                     time=time)
# extract dataset
matcheddf=matched$dataset


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(time_trial,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,treatment="trt",estimand="ATT",landmark=730)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       full_data=dftrial,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand="ATT",
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(time_trial,event)~trt+calendar,
       duration="duration",time="time",
       method="psm",search=list(psm="greedy"),
       psm_arguments=list(#caliper=0.1,
                          distance="glm",
                          seed=123,tune=FALSE))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0[1]
lower=bate-1.96*sd(B$t[,1])
upper=bate+1.96*sd(B$t[,1])
# converged proportion?
sum(B$t[,2])/length(B$t[,2])


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=psm_matching(data=dftrial,
                     formula=trt~ageonset+sex,
                     ratio=2,
                     replacement=F,
                     distance='random forest',
                     method='greedy',
                     tune=FALSE,
                     caliper=NULL,
                     duration=duration,
                     time=time)
# extract dataset
matcheddf=matched$dataset


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(time_trial,event)~trt+calendar,x=T,model=T,
              data=matcheddf,robust=T,cluster=cluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(data=matcheddf,model=model,
              treatment="trt",estimand="ATT",landmark=730)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       full_data=dftrial,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand="ATT",
       search=list(psm="greedy"),
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(time_trial,event)~trt+calendar,
       duration="duration",time="time",
       psm_arguments=list(tune=FALSE,
                          seed=181299,
                          distance="random forest"),
       method="psm")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
rrest=B$t0[1]
lower=B$t0-1.96*sd(B$t[,1])
upper=B$t0+1.96*sd(B$t[,1])
# converged proportion?
sum(B$t[,2])/length(B$t[,2])


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=psm_matching(data=dftrial,
                     formula=trt~ageonset+sex,
                     ratio=2,
                     replacement=F,
                     distance='mahalanobis',
                     method='optimal.ilp',
                     #tune=TRUE,
                     #caliper=0.1,
                     duration=duration,
                     time=time)
# extract dataset
matcheddf=matched$dataset


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(time_trial,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(data=matcheddf,model=model,
              treatment="trt",estimand="ATT",landmark=730)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       #indices=1:nrow(dftrial),
       full_data=dftrial,
       #replacement=T,
       statistic=boot_matching,
       estimand="ATT",
       ratio=2,
       landmark=730,
       estimate="surv_rr",
       search=list(psm="optimal.ilp"),
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(time_trial,event)~trt+calendar,
       duration="duration",time="time",
       psm_arguments=list(tune=FALSE,distance="mahalanobis",caliper=NULL),
       method="psm")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0[1]
lower=bate-1.96*sd(B$t[,1])
upper=bate+1.96*sd(B$t[,1])
# converged proportion?
sum(B$t[,2])/length(B$t[,2])


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\main_analysis\\with duration\\time_to_event\\risk\\single_methods\\single_methods.qmd")

