## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr,tidyr,tidyverse,ggplot2,lpSolve,
               optmatch,ipw,lme4,lmerTest,emmeans,boot,randomForest)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions")
source("ilp_matching.R")
source("boot_matching.R")
source("exact_matching_interpolation.R")
source("alignment_psm_matching_interpolation.R")


## ----warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------
df=readxl::read_xlsx("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_hui_like_dur.xlsx")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# set seed for reproducibility
set.seed(11112000)
# sample untreated
untrtid=unique(df$id[which(df$trt_overall==0 & df$dateonset<df$trialdate)])
# sample treated
trtid=unique(df$id[which(df$trt_overall==1 & df$dateonset<df$trialdate)])
# bind together
ids=c(sample(untrtid,99),sample(trtid,33))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
trialdf=df %>% 
  filter(miss==0) %>% 
  filter(id %in% ids) 


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
trialdf=trialdf %>% 
  mutate(across(c(date,dateonset,starttrial,lastvisit,trialdate),as.Date)) %>% 
  mutate(time=time/365.25) 


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
treateddf=trialdf %>% 
  filter(trt_overall==1)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=matching_interpolation(data=trialdf,identifier=id,y=hui_linear_trt,
                 treatment=trt_overall,ageofonset=aoo,
                 sex=sex,duration=duration,time=timen,
                 ratio=2,caliper=0.5,step=0.1,maxcaliper=0.5,method='greedy')


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(scipen=999)
trial=matched$interpolated %>% 
  filter(time_aligned<=2) %>%
  lmer(y~(1|cluster)+(1|cluster:id)+trt*time_aligned+yoo,.)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimand=emmeans(trial, ~ trt | time_aligned,
                              at = list(time_aligned=2))
# see contrast
ate=summary(contrast(estimand, method="trt.vs.ctrl", ref=1))$estimate[1]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pick baseline
trialdfb=trialdf %>% group_by(id) %>% slice(1) %>% ungroup()
# run bootstrap
B=boot(data=trialdfb,
       full_data=trialdf,
       R=5,
       landmark=2,
       duration="duration",time="timen",
       statistic=boot_matching,
       treatment_model=trt_overall~aoo+sex,
       outcome_model=hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*time_aligned+yoo,
       method="exact",
       estimate="ate",
       exact_arguments=list(caliper=0.5,
                            step=0.1,
                            maxcaliper=0.5),
       longitudinal_arguments=list(identifier="id",
                                   y="hui_linear_trt",
                                   time="time_aligned"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
# difference?
bate-ate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=psm_matching_interpolation(data=trialdf,identifier=id,outcome=hui_linear_trt,
                     formula=trt_overall~aoo+sex,duration=duration,time=timen,
                     ratio=2,replacement=F,tune=F,
                     distance='glm',method='greedy',caliper=0.1)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(scipen=999)
trial=matched$interpolated %>% 
  filter(time_aligned<=2) %>%
  lmer(y~(1|cluster)+(1|cluster:id)+trt*time_aligned+yoo,.)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimand=emmeans(trial, ~ trt | time_aligned,
                              at = list(time_aligned=2))
# see contrast
ate=summary(contrast(estimand, method="trt.vs.ctrl", ref=1))$estimate[1]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pick baseline
trialdfb=trialdf %>% group_by(id) %>% slice(1) %>% ungroup()
# run bootstrap
B=boot(data=trialdfb,
       full_data=trialdf,
       R=5,
       landmark=2,
       duration="duration",time="timen",
       statistic=boot_matching,
       treatment_model=trt_overall~aoo+sex,
       outcome_model=hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*time_aligned+yoo,
       method="psm",
       estimate="ate",
       psm_arguments=list(tune=FALSE,
                          caliper=0.1,
                          distance="glm"),
       longitudinal_arguments=list(identifier="id",
                                   y="hui_linear_trt",
                                   time="time_aligned"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
# difference?
bate-ate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=psm_matching_interpolation(data=trialdf,identifier=id,outcome=hui_linear_trt,
                     formula=trt_overall~aoo+sex,duration=duration,time=timen,
                     ratio=2,replacement=F,tune=F,seed=181299,
                     distance='random forest',method='greedy',caliper=0.1)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(scipen=999)
trial=matched$interpolated %>% 
  filter(time_aligned<=2) %>%
  lmer(y~(1|cluster)+(1|cluster:id)+trt*time_aligned+yoo,.)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimand=emmeans(trial, ~ trt | time_aligned,
                              at = list(time_aligned=2),
                 cov.reduce=list(yoo=mean))
# see contrast
ate=summary(contrast(estimand, method="trt.vs.ctrl", ref=1))$estimate[1]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pick baseline
trialdfb=trialdf %>% group_by(id) %>% slice(1) %>% ungroup()
# run bootstrap
B=boot(data=trialdfb,
       full_data=trialdf,
       R=5,
       landmark=2,
       duration="duration",time="timen",
       statistic=boot_matching,
       treatment_model=trt_overall~aoo+sex,
       outcome_model=hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*time_aligned+yoo,
       method="psm",
       estimate="ate",
       psm_arguments=list(tune=FALSE,
                          seed=181299,
                          caliper=0.1,
                          distance="random forest"),
       longitudinal_arguments=list(identifier="id",
                                   y="hui_linear_trt",
                                   time="time_aligned"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
# difference?
bate-ate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=psm_matching_interpolation(data=trialdf,identifier=id,outcome=hui_linear_trt,
                     formula=trt_overall~aoo+sex,duration=duration,time=timen,
                     ratio=2,replacement=F,tune=F,
                     distance='mahalanobis',method='greedy',caliper=0.1)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(scipen=999)
trial=matched$interpolated %>% 
  filter(time_aligned<=2) %>%
  lmer(y~(1|cluster)+(1|cluster:id)+trt*time_aligned+yoo,.)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimand=emmeans(trial, ~ trt | time_aligned,
                              at = list(time_aligned=2),
                 cov.reduce=list(yoo=mean))
# see contrast
ate=summary(contrast(estimand, method="trt.vs.ctrl", ref=1))$estimate[1]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pick baseline
trialdfb=trialdf %>% group_by(id) %>% slice(1) %>% ungroup()
# run bootstrap
B=boot(data=trialdfb,
       full_data=trialdf,
       R=5,
       landmark=2,
       duration="duration",time="timen",
       statistic=boot_matching,
       treatment_model=trt_overall~aoo+sex,
       outcome_model=hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*time_aligned+yoo,
       method="psm",
       estimate="ate",
       psm_arguments=list(tune=FALSE,
                          seed=181299,
                          caliper=0.1,
                          distance="mahalanobis"),
       longitudinal_arguments=list(identifier="id",
                                   y="hui_linear_trt",
                                   time="time_aligned"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
# difference?
bate-ate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\main_analysis\\with duration\\continuous\\single_methods\\single_methods.qmd")

