## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr,tidyr,tidyverse,ggplot2,lpSolve,MatchIt,
               optmatch,ipw,lme4,lmerTest,emmeans,boot,randomForest)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions")
source("ilp_matching.R")
source("exact_matching.R")
source("boot_matching.R")
source("gcomp.R")
source("iptw.R")
source("boot_iptw.R")


## ----warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------
df=readxl::read_xlsx("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_hui_like_without_dur.xlsx")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# set seed for reproducibility
set.seed(11112000)
# sample untreated
untrtid=unique(df$id[which(df$trt_overall==0)])
# sample treated
trtid=unique(df$id[which(df$trt_overall==1)])
# bind together
ids=c(sample(untrtid,99),sample(trtid,28))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
trialdf=df %>% 
  filter(miss==0) %>% 
  filter(id %in% ids) %>% 
  group_by(id) %>% 
  mutate(id=cur_group_id())


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
trialdf=trialdf %>% 
  mutate(across(c(date,dateonset,starttrial,lastvisit,trialdate),as.Date)) %>% 
  mutate(time=time/365.25)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
baseline=trialdf %>% 
  group_by(id) %>% 
  filter(row_number()==1) %>% 
  ungroup()


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=matching(data=baseline,treatment=trt_overall,ageofonset=aoo,sex=sex,
                 duration=NULL,time=timen,method='greedy',
                 caliper=0.5,step=0.1,maxcaliper=0.5,ratio=2,
                 replacement=F,history=F,warning=T)
# inspect
matcheddf=matched$dataset
# get data out and keep only cluster and id
matches=matched$dataset %>% select(id,cluster)
# merge with initial dataset
trialdf_matched=merge(trialdf,matches,by='id') %>% 
  arrange(id,timen)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# model
trial=trialdf_matched %>% 
  lmer(hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*timen+yoo,.)
# marginalize
estimand=emmeans(trial, ~ trt | timen,
                              at = list(timen=2),
                 cov.reduce=list(yoo=mean))
# see contrast
ate=summary(pairs(estimand))$estimate[1]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# select random cluster
clusters=sample(unique(trialdf_matched$cluster),5)
# table
table(trialdf_matched$trt_overall)
# plot cluster
trialdf_matched %>% 
  filter(cluster %in% clusters) %>% 
  filter(timen<=2) %>% 
  ggplot(.,aes(x=timen,y=hui_linear_trt,col=factor(trt_overall),group = id))+
  geom_point()+
  geom_line()+
  guides(col='none')+
  facet_wrap(~cluster)+
  theme_bw()


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pick baseline
trialdfb=trialdf %>% group_by(id) %>% slice(1) %>% ungroup()
# run bootstrap
B=boot(data=trialdfb,
       full_data=trialdf,
       R=5,
       landmark=2,
       duration=NULL,time="timen",
       statistic=boot_matching,
       treatment_model=trt_overall~aoo+sex,
       outcome_model=hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*timen+yoo,
       method="exact",
       estimate="ate",
       exact_arguments=list(caliper=0.5,
                            step=0.1,
                            maxcaliper=0.5),
       longitudinal_arguments=list(identifier="id",
                                   y="hui_linear_trt",
                                   time="timen"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
# differences
ate-bate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=matchit(trt_overall~aoo+sex,
                ratio=2,
                #method='optimal',
                caliper=0.1,
                #exact=~ageonset,
                baseline)
# get data out and keep only cluster and id
matches=match.data(matched) %>% select(id,subclass)
# merge with initial dataset
trialdf_matched=merge(trialdf,matches,by='id') %>% 
  mutate(cluster=subclass) %>% 
  select(-subclass)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# select random cluster
clusters=sample(unique(trialdf_matched$cluster),5)
# table
table(trialdf_matched$trt_overall)
# plot cluster
trialdf_matched %>% 
  filter(cluster %in% clusters) %>% 
  filter(timen<=2) %>% 
  ggplot(.,aes(x=timen,y=hui_linear_trt,col=factor(trt_overall),group = id))+
  geom_point()+
  geom_line()+
  guides(col='none')+
  facet_wrap(~cluster)+
  theme_bw()


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# model
trial=trialdf_matched %>% 
  lmer(hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*timen+yoo,.)
# marginalize
estimand=emmeans(trial, ~ trt | timen,
                              at = list(timen=2),
                 cov.reduce=list(yoo=mean))
# see contrast
ate=summary(pairs(estimand))$estimate[1]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pick baseline
trialdfb=trialdf %>% group_by(id) %>% slice(1) %>% ungroup()
# run bootstrap
B=boot(data=trialdfb,
       full_data=trialdf,
       R=5,
       landmark=2,
       duration=NULL,time="timen",
       statistic=boot_matching,
       treatment_model=trt_overall~aoo+sex,
       outcome_model=hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*timen+yoo,
       method="psm",
       estimate="ate",
       psm_arguments=list(caliper=0.1,
                          distance="glm",
                          seed=123,
                          tune=FALSE),
       longitudinal_arguments=list(identifier="id",
                                   y="hui_linear_trt",
                                   time="timen"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
# differences
ate-bate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=matchit(trt_overall~aoo+sex,
                ratio=2,
                distance='randomforest',
                #method='optimal',
                caliper=NULL,
                #exact=~ageonset,
                baseline)
# get data out and keep only cluster and id
matches=match.data(matched) %>% select(id,subclass)
# merge with initial dataset
trialdf_matched=merge(trialdf,matches,by='id') %>% 
  mutate(cluster=subclass) %>% 
  select(-subclass)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# model
trial=trialdf_matched %>% 
  lmer(hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*timen+yoo,.)
# marginalize
estimand=emmeans(trial, ~ trt | timen,
                              at = list(timen=2),
                 cov.reduce=list(yoo=mean))
# see contrast
ate=summary(pairs(estimand))$estimate[1]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pick baseline
trialdfb=trialdf %>% group_by(id) %>% slice(1) %>% ungroup()
# run bootstrap
B=boot(data=trialdfb,
       full_data=trialdf,
       R=5,
       landmark=2,
       duration=NULL,time="timen",
       statistic=boot_matching,
       treatment_model=trt_overall~aoo+sex,
       outcome_model=hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*timen+yoo,
       method="psm",
       estimate="ate",
       psm_arguments=list(caliper=NULL,
                          distance="randomforest"),
       longitudinal_arguments=list(identifier="id",
                                   y="hui_linear_trt",
                                   time="timen"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
ate-bate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
matched=matchit(trt_overall~aoo+sex,
                ratio=2,
                distance='mahalanobis',
                #method='optimal',
                #exact=~ageonset,
                baseline)
# get data out and keep only cluster and id
matches=match.data(matched) %>% select(id,subclass)
# merge with initial dataset
trialdf_matched=merge(trialdf,matches,by='id') %>% 
  mutate(cluster=subclass)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# model
trial=trialdf_matched %>% 
  lmer(hui_linear_trt~(1|subclass)+(1|subclass:id)+trt*timen+yoo,.)
# marginalize
estimand=emmeans(trial, ~ trt | timen,
                              at = list(timen=2),
                 cov.reduce=list(yoo=mean))
# see contrast
ate=summary(pairs(estimand))$estimate[1]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pick baseline
trialdfb=trialdf %>% group_by(id) %>% slice(1) %>% ungroup()
# run bootstrap
B=boot(data=trialdfb,
       full_data=trialdf,
       R=5,
       landmark=2,
       duration=NULL,time="timen",
       statistic=boot_matching,
       treatment_model=trt_overall~aoo+sex,
       outcome_model=hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*timen+yoo,
       method="psm",
       estimate="ate",
       psm_arguments=list(caliper=NULL,
                          distance="mahalanobis"),
       longitudinal_arguments=list(identifier="id",
                                   y="hui_linear_trt",
                                   time="timen"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
ate-bate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
weights=iptw(data=baseline,estimand="ATT",
             formula=trt_overall~aoo+sex)
# extract weights
dataw=weights$data %>% 
  select(id,w)
# merge dataset
trialdfw=merge(trialdf,dataw,by='id')


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(scipen=999)
fit=trialdfw %>% 
  #filter(time_aligned<=2) %>% 
  lmer(hui_linear_trt~(1|id)+trt_overall*timen+yoo,weights = w,data=.)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
estimand=emmeans(fit, ~ trt_overall | timen,
                              at = list(timen=2),
                 cov.reduce=list(yoo=mean))
# see contrast
ate=summary(pairs(estimand))$estimate[1]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
B=boot(data=trialdf,
       statistic=boot_iptw,
       landmark=2,
       time="timen",
       estimate="ate",
       R=5,
       trt_model=list(formula=trt_overall~aoo+sex,
                             type='stabilised',
                             estimand='ATT',
                             truncate=NULL,
                             distance='glm',
                             tune=FALSE),
       outcome_model=list(formula=hui_linear_trt~(1|id)+trt_overall*timen+yoo))

# see results
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
ate-bate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ate=gcomp(trialdf,
          estimate="ate",
          estimand="ATT",
          treatment="trt_overall",
          time="timen",
          formula=hui_linear_trt~(1|id)+trt_overall*timen+aoo+sex+yoo,
          landmark=2)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
B=boot(trialdf,
       statistic=gcomp,
       estimand="ATT",
       estimate="ate",
       treatment="trt_overall",
       time="timen",
       R=5,
       formula=hui_linear_trt~(1|id)+trt_overall*timen+aoo+sex+yoo,
       landmark=2)
# see results
bate=B$t0
lower=bate-1.96*sd(B$t)
upper=bate+1.96*sd(B$t)
# difference
ate-bate


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\without duration\\continuous\\single_methods\\single_methods.qmd")

