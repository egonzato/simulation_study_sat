## ----------------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(dplyr,tidyr,tidyverse,ggplot2,lpSolve,smd,ESS,
               optmatch,ipw,lme4,lmerTest,emmeans,boot,randomForest,boot,
               foreach,doParallel,doRNG,future,progressr,doFuture)


## ----warning=FALSE-----------------------------------------------------------------------------------------------------------------------
df=readxl::read_xlsx("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\population_hui_like_without_dur.xlsx")


## ----------------------------------------------------------------------------------------------------------------------------------------
# get ids out of function
trtid=df %>% filter(trt_overall==1) %>% select(id) %>% unique()
untrtid=df %>% filter(trt_overall==0) %>% select(id) %>% unique()


## ----warning=FALSE-----------------------------------------------------------------------------------------------------------------------
setwd("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions")
source("iptw.R")
source("boot_iptw.R")
source("gcomp.R")
source("dre.R")
source("exact_matching.R")


## ----------------------------------------------------------------------------------------------------------------------------------------
# grid of combinations
grid=expand.grid(ntrt=c(28),nctrl=c(90,150))
n_reps=3
tasks=grid[rep(1:nrow(grid),each=n_reps),]
tasks$rep=rep(1:n_reps,times=nrow(grid))
total_steps=nrow(tasks)
# register components to run in parallel
registerDoFuture()
plan(multisession, workers = parallel::detectCores()-1)
handlers("txtprogressbar")
registerDoRNG(181299)
# run loop
with_progress({
  # define number of steps
  p=progressor(steps = total_steps)
  # define results object
  results=foreach(t=1:nrow(tasks), .combine = rbind,
                     .packages = c("survival","MatchIt","boot","dplyr","smd","WeighIt",
                                   "optmatch","lpSolve","randomForest"),
                     .export = c("trt","untrt","exact_matching_interpolation","alignment_psm_matching_interpolation",
                            "bootiptw","gcomp","dre","psm_matching")) %dopar% {
      # define number of treated and untreated
      ntrti=tasks$ntrt[t]
      nctrli=tasks$nctrl[t]
      # print to see where we are
      p(sprintf("ntrt=%d, nctrl=%d, i=%d",ntrti,nctrli,tasks$rep[t]))
        # sample ids
          ids=c(sample(untrtid$id,90),sample(trtid$id,28))
          # retain patients sampled
          trialdf=df %>% 
            filter(miss==0) %>% 
            filter(id %in% ids) %>% 
          group_by(id) %>% 
          mutate(id=cur_group_id())
          # get baseline
          baseline=trialdf %>% group_by(id) %>% slice(1) %>% ungroup()
          ## calculate metrics
          smdagebi=smd(baseline$aoo,baseline$trt)$estimate
          smdsexbi=smd(baseline$sex,baseline$trt)$estimate
          # exact matching
          ## match
          matched=matching(data=baseline,treatment=trt_overall,ageofonset=aoo,sex=NULL,
                         duration=NULL,time=timen,method='greedy',
                         caliper=0.5,step=0.1,maxcaliper=0.5,ratio=2,
                         replacement=F,history=F,warning=F)
          # inspect
          matcheddf=matched$dataset
          # get data out and keep only cluster and id
          matches=matched$dataset %>% select(id,cluster)
          # merge with initial dataset
          trialdf_matched=merge(trialdf,matches,by='id') %>% 
            arrange(id,timen) %>% 
            filter(timen<=2)
          ## get baseline matched
          baselinem=trialdf_matched %>% group_by(id) %>% slice(1)
          ## run model
          trial=lmer(hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*timen+yoo,trialdf_matched)
          ## marginalize based on the model
          estimand=emmeans(trial, ~ trt | timen,
                                      at = list(timen=2))
          ## append quantities
          res_exact=data.frame(beta=summary(pairs(estimand))$estimate,
                             se=summary(pairs(estimand))$SE, 
                             method="exact",
                             smdageb=smdagebi,
                             smdsexb=smdsexbi,
                             smdagea=smd(baselinem$aoo,baselinem$trt)$estimate,
                             smdsexa=smd(baselinem$sex,baselinem$trt)$estimate,
                             ntrt=table(baselinem$trt_overall)[2],
                             nctrl=table(baselinem$trt_overall)[1]) %>% 
            mutate(ess=ntrt+nctrl)
          # propensity score matching
          ## match
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
            select(-subclass) %>% 
            filter(timen<=2)
          ## run model
          trial=lmer(hui_linear_trt~(1|cluster)+(1|cluster:id)+trt*timen+yoo,trialdf_matched)
          ## marginalize based on the model
          estimand=emmeans(trial, ~ trt | timen,
                                      at = list(timen=2))
          ## append quantities
          res_psm=data.frame(beta=summary(pairs(estimand))$estimate,
                             se=summary(pairs(estimand))$SE, 
                             method="psm",
                             smdageb=smdagebi,
                             smdsexb=smdsexbi,
                             smdagea=smd(baselinem$aoo,baselinem$trt)$estimate,
                             smdsexa=smd(baselinem$sex,baselinem$trt)$estimate,
                             ntrt=table(baselinem$trt_overall)[2],
                             nctrl=table(baselinem$trt_overall)[1]) %>% 
            mutate(ess=ntrt+nctrl)
          # iptw
          ## point estimate and standard error
          # boot_res=boot(data=dftrial,
          #             statistic=bootiptw,
          #             #parallel="snow", ncpus=parallel::detectCores()-1,
          #             R=500,
          #             trt_model=list(formula=trt~ageonset+sex,
          #                            type='stabilised',
          #                            estimand='ATT',
          #                            truncate=NULL,
          #                            distance='glm',
          #                            tune=FALSE),
          #             outcome_model=list(formula=Surv(timetotal,event)~trt+calendar,
          #                                treatment='trt'))
          # ## append quantities
          # res_iptw=data.frame(beta=boot_res$t0, se=sd(boot_res$t), method="iptw")
          # append all quantities
        rbind(res_exact, res_psm)  
        out=rbind(res_exact, res_psm)
        out$ntrt_pool=ntrti
        out$nctrl_pool=nctrli
        out
                     }
})


## ----------------------------------------------------------------------------------------------------------------------------------------
results=results %>% 
  mutate(truth=0.1,
         lower=beta-1.96*se,
         upper=beta+1.96*se)


## ----------------------------------------------------------------------------------------------------------------------------------------
writexl::write_xlsx(results,
                    "C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\paper\\without duration\\continuous\\results_01.xlsx")


## ----------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\without duration\\continuous\\analysis\\analysis.qmd")

