## ------------------------------------------------------------------------------------------------------------------------
pacman::p_load(readxl,dplyr,tidyverse,optmatch,lpSolve,randomForest,survival,ggsurvfit,survminer,MatchIt,boot)


## ------------------------------------------------------------------------------------------------------------------------
pop=readRDS("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\paper\\main_analysis\\without duration\\time_to_event\\population_rr_0.5.xlsx")


## ------------------------------------------------------------------------------------------------------------------------
trt=pop %>% filter(trt==1)
untrt=pop %>% filter(trt==0)


## ------------------------------------------------------------------------------------------------------------------------
set.seed(11112000)
dftrial=rbind(trt[sample(1:nrow(trt),20),],
            untrt[sample(1:nrow(untrt),40),]) %>% 
  mutate(event=ifelse(timetotal>730,0,1),
         timetotal=ifelse(timetotal>730,730,timetotal),
         sex=ifelse(sex=="female",0,1))


## ------------------------------------------------------------------------------------------------------------------------
setwd("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions")
source("ilp_matching_full.R")
source("exact_matching.R")
source("boot_matching.R")
source("alignment_psm_matching_optimal.R")
source("iptw.R")
source("boot_iptw.R")
source("gcomp.R")
source("risk_ratio.R")


## ------------------------------------------------------------------------------------------------------------------------
matched=matching(data=dftrial,treatment=trt,ageofonset=ageonset,sex=sex,
                 duration=NULL,time=timetotal,method='optimal',
                 caliper=NULL,step=NULL,maxcaliper=NULL,
                 #caliper=0.5,step=0.1,maxcaliper=0.5,
                 ratio=2,
                 replacement=FALSE,history=F,warning=TRUE)
# extract dataset
matcheddf=matched$dataset
# number of treated?
table(matcheddf$trt); table(table(matcheddf$cluster))


## ------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,
              treatment="trt",estimand='ATT',landmark=730)


## ------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       strata=dftrial[["trt"]],
       full_data=dftrial,
       statistic=boot_matching,
       estimate="surv_rr",
       landmark=730,
       exact_arguments=list(caliper=NULL,
                            step=NULL,
                            maxcaliper=NULL),
       estimand='ATT',
       ratio=2,
       search=list(exact="optimal"),
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(timetotal,event)~trt+calendar,
       duration=NULL,time="timetotal",
       method="exact")


## ------------------------------------------------------------------------------------------------------------------------
est=B$t0[1]
lower=B$t0-1.96*sd(B$t[,1])
upper=B$t0+1.96*sd(B$t[,1])
# see runs
B$t


## ----eval=FALSE----------------------------------------------------------------------------------------------------------
# set.seed(151515)
# # set iterations
# B=5
# # collect estimates
# rrs=c(); n=c(); ntrt=c(); nuntrt=c()
# # loop through
# for (i in 1:B){
#   # resample ids
#   idx_trt=sample(nrow(dftrial[dftrial$trt==1,]),
#                  size=nrow(dftrial[dftrial$trt==1,]),replace=TRUE)
#   idx_untrt=sample(nrow(dftrial[dftrial$trt==0,]),
#                  size=nrow(dftrial[dftrial$trt==0,]),replace=TRUE)
#   # put together
#   dfboot_trt=dftrial[idx_trt,]; dfboot_untrt=dftrial[idx_untrt,]
#   # bind
#   dfboot=rbind(dfboot_trt,dfboot_untrt)
#   # redefine ids
#   dfboot=dfboot %>% mutate(id=row_number())
#   # run matching
#   matched=matching(data=dfboot,treatment=trt,ageofonset=ageonset,sex=sex,
#                    duration=NULL,time=timetotal,method='optimal',
#                    caliper=NULL,step=NULL,maxcaliper=NULL,
#                    #caliper=0.5,step=0.1,maxcaliper=0.5,
#                    ratio=1,
#                    replacement=FALSE,history=F,warning=TRUE)
#   # extract dataset
#   matcheddf=matched$dataset
#   # see numbers
#   n=c(n,nrow(matcheddf)); ntrt=c(ntrt,table(matcheddf$trt)[2]);nuntrt=c(ntrt,table(matcheddf$trt)[1])
#   # model
#   model=coxph(Surv(time_trial,event)~trt+calendar,x=T,
#               data=matcheddf,robust=T,cluster=cluster)
#   # caclulate risk ratio
#   rr=risk_ratio(matcheddf,model=model,
#                 treatment="trt",estimand="ATT",landmark=730)
#   # append
#   rrs=c(rrs,rr)
# }


## ------------------------------------------------------------------------------------------------------------------------
matched=psm_matching(data=dftrial,
                     formula=trt~ageonset+sex,
                     ratio=2,
                     replacement=F,
                     distance='glm',
                     method='optimal',
                     #tune=TRUE,
                     caliper=NULL,
                     duration=NULL,
                     time=timetotal)
# extract dataset
matcheddf=matched$dataset


## ------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,
              treatment="trt",estimand='ATT',landmark=730)


## ------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       strata=dftrial[["trt"]],
       full_data=dftrial,
       landmark=730,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand='ATT',
       ratio=2,
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(timetotal,event)~trt+calendar,
       psm_arguments=list(caliper=NULL, 
                          distance="glm", 
                          seed=123, 
                          tune=FALSE),
       duration=NULL,
       time="timetotal",
       method="psm",
       search=list(psm="optimal"))


## ------------------------------------------------------------------------------------------------------------------------
est=B$t0[1]
lower=B$t0-1.96*sd(B$t[,1])
upper=B$t0+1.96*sd(B$t[,1])
# see B
B$t


## ------------------------------------------------------------------------------------------------------------------------
# parameters
betas=c(); B=5; nuntrt=c(); ntrt=c(); n=c(); rrs=c()
# loop
for (b in 1:B){
  # resample ids
  ## treated
  dfboot_trt=dftrial %>%
    filter(trt == 1) %>%
    slice_sample(n = 20, replace = TRUE)
  ## untreated
  dfboot_untrt=dftrial %>%
    filter(trt == 0) %>%
    slice_sample(n = 50, replace = TRUE)
  # bind together
  dfboot=bind_rows(dfboot_trt, dfboot_untrt) %>%
    mutate(id = row_number())
  # redefine ids
  dfboot=dfboot %>% mutate(id=row_number())
  # run matching
    ## own function
  matched=psm_matching(data=dfboot,
                       formula=trt~ageonset+sex,
                       ratio=2,
                       replacement=F,
                       distance='glm',
                       method='optimal',
                       duration=NULL,
                       time=timetotal)
  ## extract dataset
  matcheddfp=matched$dataset
  # model
  model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddfp,robust=T,cluster=cluster)
  ## append quantities
  n=c(n,nrow(matcheddfp)); nuntrt=c(nuntrt,table(matcheddfp$trt)[1]); 
  ntrt=c(ntrt,table(matcheddfp$trt)[2]);   betas=c(betas,coef(model)[1])
  rrs=c(rrs,
        risk_ratio(matcheddfp,model=model,
                   treatment="trt",estimand='ATT',landmark=730))
}


## ------------------------------------------------------------------------------------------------------------------------
matched=psm_matching(data=dftrial,
                     formula=trt~ageonset+sex,
                     ratio=2,
                     replacement=F,
                     distance='mahalanobis',
                     method='optimal',
                     #tune=TRUE,
                     #caliper=0.1,
                     duration=NULL,
                     time=timetotal)
# extract dataset
matcheddf=matched$dataset


## ------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,
              treatment="trt",estimand='ATT',landmark=730)


## ------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       strata=dftrial[["trt"]],
       full_data=dftrial,
       landmark=730,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand='ATT',
       ratio=2,
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(timetotal,event)~trt+calendar,
       psm_arguments=list(caliper=NULL, distance="mahalanobis", 
                          seed=123, tune=FALSE),
       duration=NULL,time="timetotal",
       method="psm",search=list(psm="optimal"))


## ------------------------------------------------------------------------------------------------------------------------
est=B$t0[1]
lower=B$t0-1.96*sd(B$t[,1])
upper=B$t0+1.96*sd(B$t[,1])
# see B
B$t


## ------------------------------------------------------------------------------------------------------------------------
# parameters
betas=c(); B=5; nuntrt=c(); ntrt=c(); n=c(); rrs=c()
# loop
for (b in 1:B){
  # resample ids
  ## treated
  dfboot_trt=dftrial %>%
    filter(trt == 1) %>%
    slice_sample(n = 20, replace = TRUE)
  ## untreated
  dfboot_untrt=dftrial %>%
    filter(trt == 0) %>%
    slice_sample(n = 50, replace = TRUE)
  # bind together
  dfboot=bind_rows(dfboot_trt, dfboot_untrt) %>%
    mutate(id = row_number())
  # redefine ids
  dfboot=dfboot %>% mutate(id=row_number())
  # run matching
    ## own function
  matched=psm_matching(data=dfboot,
                       formula=trt~ageonset+sex,
                       ratio=2,
                       replacement=F,
                       distance='mahalanobis',
                       method='optimal',
                       duration=NULL,
                       time=timetotal)
  ## extract dataset
  matcheddfp=matched$dataset
  # model
  model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddfp,robust=T,cluster=cluster)
  ## append quantities
  n=c(n,nrow(matcheddfp)); nuntrt=c(nuntrt,table(matcheddfp$trt)[1]); 
  ntrt=c(ntrt,table(matcheddfp$trt)[2]);   betas=c(betas,coef(model)[1])
  rrs=c(rrs,
        risk_ratio(matcheddfp,model=model,
                   treatment="trt",estimand='ATT',landmark=730))
}


## ------------------------------------------------------------------------------------------------------------------------
set.seed(181299)
matched=psm_matching(data=dftrial,
                     formula=trt~ageonset+sex,
                     ratio=2,
                     replacement=F,
                     distance='random forest',
                     method='optimal',
                     tune=TRUE,
                     #caliper=0.1,
                     duration=NULL,
                     time=timetotal)


## ------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,x=T,
              data=matcheddf,robust=T,cluster=cluster)


## ------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(matcheddf,model=model,
              treatment="trt",estimand='ATT',landmark=730)


## ------------------------------------------------------------------------------------------------------------------------
B=boot(data=dftrial,
       R=5,
       strata=dftrial[["trt"]],
       full_data=dftrial,
       landmark=730,
       statistic=boot_matching,
       estimate="surv_rr",
       estimand='ATT',
       ratio=2,
       treatment_model=trt~ageonset+sex,
       outcome_model=Surv(time_trial,event)~trt+calendar,
       psm_arguments=list(caliper=NULL, distance="random forest", 
                          seed=181299, tune=FALSE),
       duration=NULL,time="timetotal",
       method="psm")


## ------------------------------------------------------------------------------------------------------------------------
est=B$t0
lower=B$t0-1.96*sd(B$t)
upper=B$t0+1.96*sd(B$t)
B$t


## ------------------------------------------------------------------------------------------------------------------------
weights=iptw(data=dftrial,type='stabilised',estimand='ATT',
             formula=trt~ageonset+sex)
# dataset
dftrialw=weights$data


## ------------------------------------------------------------------------------------------------------------------------
model=coxph(Surv(timetotal,event)~trt+calendar,weights=w,data=dftrialw)


## ------------------------------------------------------------------------------------------------------------------------
rr=risk_ratio(dftrialw,model=model,treatment="trt",estimand='ATT',landmark=730)


## ------------------------------------------------------------------------------------------------------------------------
set.seed(123)
B=boot(data=dftrial,
       strata=dftrial[["trt"]],
       statistic=boot_iptw,
       landmark=730,
       estimate="surv_rr",
       R=100,
       trt_model=list(formula=trt~ageonset+sex,
                      type='stabilised',
                      estimand='ATT',
                      truncate=NULL,
                      distance='glm',
                      tune=FALSE),
       outcome_model=list(formula=Surv(timetotal,event)~trt+calendar))

# see results
bate=B$t0[1]
lower=bate-1.96*sd(B$t[,1])
upper=bate+1.96*sd(B$t[,1])
# see number of convergence
mean(B$t[,2])


## ----eval=FALSE----------------------------------------------------------------------------------------------------------
# rr=gcomp(data=dftrial,
#             estimand=c("ATT"),
#             treatment="trt",
#             estimate="surv_rr",
#             formula=Surv(timetotal,event)~trt+calendar+ageonset+sex,
#             landmark=730)


## ------------------------------------------------------------------------------------------------------------------------
set.seed(123)
B=boot(data=dftrial,
       statistic=gcomp,
       R=5,
       estimand="ATT",
       treatment="trt",
       estimate="surv_rr",
       formula=Surv(timetotal,event)~trt+calendar+ageonset+sex,
       landmark=730)
# see results
bate=B$t0[1]
lower=bate-1.96*sd(B$t[,1])
upper=bate+1.96*sd(B$t[,1])
mean(B$t[,2])


## ------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\main_analysis\\without duration\\time_to_event\\risk\\single_methods\\single_methods.qmd")

