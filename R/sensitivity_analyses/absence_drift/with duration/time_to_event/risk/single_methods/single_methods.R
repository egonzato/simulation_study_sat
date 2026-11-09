## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(readxl,dplyr,tidyverse,optmatch,lpSolve,
               randomForest,survival,ggsurvfit,survminer,MatchIt,boot)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions")
source("ilp_matching.R")
source("alignment_psm_matching.R")
source("exact_matching.R")
source("boot_matching.R")
source("risk_ratio.R")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pop=read_xlsx("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_dur.xlsx")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
trt=pop %>% filter(trt==1)
untrt=pop %>% filter(trt==0)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(03111963)
dftrial=rbind(trt[sample(1:nrow(trt),15),],
            untrt[sample(1:nrow(untrt),150),])


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=matching(data=dftrial,treatment=trt,ageofonset=ageonset,sex=sex,
                 duration=duration,time=time,method='greedy',
                 caliper=0.5,step=0.125,maxcaliper=0.5,ratio=2,
                 replacement=F,history=F,warning=T)
# extract dataset
matcheddf=matched$dataset
# see optimal pairs
#matched$optimalres


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(time_trial,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,treatment="trt",estimand="ATT",landmark=365)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       full_data=dftrial,
       estimand="ATT",
       statistic=boot_matching,
       estimate="surv_rr",
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(time_trial,event)~trt+calendar,
       duration="duration",time="time",
       method="exact")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rrest=B$t0
lower=B$t0-1.96*sd(B$t)
upper=B$t0+1.96*sd(B$t)
# difference
rr-rrest


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=psm_matching(data=dftrial,
                     formula=trt~ageonset+sex,
                     ratio=2,
                     replacement=F,
                     distance='glm',
                     method='greedy',
                     #tune=TRUE,
                     caliper=0.1,
                     duration=duration,
                     time=time)
# extract dataset
matcheddf=matched$dataset


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(time_trial,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,treatment="trt",estimand="ATT",landmark=365)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       full_data=dftrial,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand="ATT",
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(time_trial,event)~trt+calendar,
       duration="duration",time="time",
       method="psm",
       psm_arguments=list(caliper=0.1,
                          distance="glm",
                          seed=123,tune=FALSE))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rrest=B$t0
lower=B$t0-1.96*sd(B$t)
upper=B$t0+1.96*sd(B$t)
# difference
rr-rrest


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(time_trial,event)~trt+calendar,x=T,model=T,
              data=matcheddf,robust=T,cluster=cluster)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(data=matcheddf,model=model,
              treatment="trt",estimand="ATT",landmark=365)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rrest=B$t0
lower=B$t0-1.96*sd(B$t)
upper=B$t0+1.96*sd(B$t)
# difference
rr-rrest


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=psm_matching(data=dftrial,
                     formula=trt~ageonset+sex,
                     ratio=2,
                     replacement=F,
                     distance='mahalanobis',
                     method='greedy',
                     #tune=TRUE,
                     #caliper=0.1,
                     duration=duration,
                     time=time)
# extract dataset
matcheddf=matched$dataset


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(time_trial,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(data=matcheddf,model=model,
              treatment="trt",estimand="ATT",landmark=365)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       #indices=1:nrow(dftrial),
       full_data=dftrial,
       #replacement=T,
       statistic=boot_matching,
       estimand="ATT",
       #ratio=2,
       estimate="surv_rr",
       search=list(psm="greedy"),
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(time_trial,event)~trt+calendar,
       duration="duration",time="time",
       psm_arguments=list(tune=FALSE,distance="mahalanobis"),
       method="psm")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rrest=B$t0
lower=B$t0-1.96*sd(B$t)
upper=B$t0+1.96*sd(B$t)
# difference
rr-rrest


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\with duration\\time_to_event\\risk\\single_methods\\single_methods.qmd")

