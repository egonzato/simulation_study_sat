## ----------------------------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(readxl,dplyr,tidyverse,optmatch,lpSolve,randomForest,survival,ggsurvfit,survminer,MatchIt,boot)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
pop=read_xlsx("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_nodur.xlsx")


## ----------------------------------------------------------------------------------------------------------------------------------------------------
trt=pop %>% filter(trt==1)
untrt=pop %>% filter(trt==0)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(181299)
dftrial=rbind(trt[sample(1:nrow(trt),33),],
            untrt[sample(1:nrow(untrt),123),]) %>% 
  mutate(event=ifelse(timetotal>365,0,1),
         timetotal=ifelse(timetotal>365,365,timetotal),
         sex=ifelse(sex=="female",0,1))


## ----------------------------------------------------------------------------------------------------------------------------------------------------
setwd("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions")
source("ilp_matching.R")
source("exact_matching.R")
source("boot_matching.R")
source("iptw.R")
source("boot_iptw.R")
source("gcomp.R")
source("risk_ratio.R")


## ----------------------------------------------------------------------------------------------------------------------------------------------------
matched=matching(data=dftrial,treatment=trt,ageofonset=ageonset,sex=sex,
                 duration=NULL,time=timetotal,method='greedy',
                 caliper=0.5,step=0.1,maxcaliper=0.5,ratio=2,
                 replacement=FALSE,history=F,warning=TRUE)
# extract dataset
matcheddf=matched$dataset


## ----------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,treatment="trt",estimand='ATT',landmark=365)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       full_data=dftrial,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand='ATT',
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(timetotal,event)~trt+calendar,
       duration=NULL,time="timetotal",
       method="exact")


## ----------------------------------------------------------------------------------------------------------------------------------------------------
est=B$t0
lower=B$t0-1.96*sd(B$t)
upper=B$t0+1.96*sd(B$t)
rr-est


## ----------------------------------------------------------------------------------------------------------------------------------------------------
matched=matchit(trt~ageonset+sex,
                ratio=2,
                #method='optimal',
                caliper=0.05,
                #exact=~ageonset,
                dftrial)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
matcheddf=match.data(matched) %>% mutate(cluster=subclass)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,treatment="trt",estimand='ATT',landmark=365)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       full_data=dftrial,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand='ATT',
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(timetotal,event)~trt+calendar,
       psm_arguments=list(caliper=0.05,distance="glm"),
       duration=NULL,time="timetotal",
       method="psm")


## ----------------------------------------------------------------------------------------------------------------------------------------------------
est=B$t0
lower=B$t0-1.96*sd(B$t)
upper=B$t0+1.96*sd(B$t)
est-rr


## ----------------------------------------------------------------------------------------------------------------------------------------------------
matched=matchit(trt~ageonset+sex,
                ratio=2,
                distance='mahalanobis',
                #method='optimal',
                caliper=NULL,
                dftrial)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
matcheddf=match.data(matched) %>% mutate(cluster=subclass)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,treatment="trt",estimand='ATT',landmark=365)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       full_data=dftrial,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand='ATT',
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(timetotal,event)~trt+calendar,
       psm_arguments=list(caliper=NULL,distance="mahalanobis"),
       duration=NULL,time="timetotal",
       method="psm")


## ----------------------------------------------------------------------------------------------------------------------------------------------------
est=B$t0
lower=B$t0-1.96*sd(B$t)
upper=B$t0+1.96*sd(B$t)
rr-est


## ----------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(181299)
matched=matchit(trt~ageonset+sex,
                ratio=2,
                distance='randomforest',
                #method='optimal',
                caliper=NULL,
                dftrial)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
matcheddf=match.data(matched) %>% mutate(cluster=subclass)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,treatment="trt",estimand='ATT',landmark=365)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       full_data=dftrial,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand='ATT',
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(timetotal,event)~trt+calendar,
       psm_arguments=list(caliper=NULL,distance="randomforest"),
       duration=NULL,time="timetotal",
       method="psm")


## ----------------------------------------------------------------------------------------------------------------------------------------------------
est=B$t0
lower=B$t0-1.96*sd(B$t)
upper=B$t0+1.96*sd(B$t)
rr-est


## ----------------------------------------------------------------------------------------------------------------------------------------------------
weights=iptw(data=dftrial,type='stabilised',estimand='ATT',
             formula=trt~ageonset+sex)
# dataset
dftrialw=weights$data


## ----------------------------------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,weights=w,data=dftrialw)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(dftrialw,model=model,treatment="trt",estimand='ATT',landmark=365)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
B=boot(data=dftrial,
       statistic=boot_iptw,
       landmark=365,
       estimate="surv_rr",
       R=5,
       trt_model=list(formula=trt~ageonset+sex,
                      type='stabilised',
                      estimand='ATT',
                      truncate=NULL,
                      distance='glm',
                      tune=FALSE),
       outcome_model=list(formula=Surv(timetotal,event)~trt+calendar))

# see results
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
rr-bate


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------
# rr=gcomp(data=dftrial,
#             estimand=c("ATT"),
#             treatment="trt",
#             estimate="surv_rr",
#             formula=Surv(timetotal,event)~trt+calendar+ageonset+sex,
#             landmark=365)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
B=boot(data=dftrial,
       statistic=gcomp,
       R=5,
       estimand="ATT",
       treatment="trt",
       estimate="surv_rr",
       formula=Surv(timetotal,event)~trt+calendar+ageonset+sex,
       landmark=365)
# see results
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
rr-bate


## ----------------------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\without duration\\time_to_event\\risk\\single_methods\\single_methods.qmd")

