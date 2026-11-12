## ---------------------------------------------------------------------------------------------------------------------
pacman::p_load(readxl,dplyr,tidyverse,optmatch,lpSolve,smd,WeightIt,lme4,
               emmeans,
               randomForest,survival,ggsurvfit,survminer,MatchIt,boot,
               foreach,doParallel,doRNG,future,progressr,doFuture)


## ---------------------------------------------------------------------------------------------------------------------
# missingness mechanism: "mcar" or "mar" (populations created in dgm/missingness.qmd)
mechanism="mcar"
df=readRDS(paste0("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\paper\\sensitivity_analysis\\missing_data\\with duration\\population_diff_0.15_",mechanism,".rds")) %>% 
  select(id,timen,trt,trt_overall,duration,aoo,sex,hui_linear_trt,yoo)
# n rows
nrow(df)


## ----warning=F--------------------------------------------------------------------------------------------------------
function_dirs=c("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\wrapper\\continuous\\with duration",
                "C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\matching",
                "C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\functions\\regression")
for (i in 1:length(function_dirs)){
  files=list.files(function_dirs[i], pattern = "\\.R$", full.names = TRUE)
sapply(files, source)
}


## ----warning=F--------------------------------------------------------------------------------------------------------
# grid of combinations
grid=expand.grid(ntrt=c(20,30,40),nctrl=c(50,100,150))
#grid=expand.grid(ntrt=c(2000),nctrl=c(5000))
ratio=c(1,2,3)
matching_methods=c("duration_matching",
                   "exact_greedy",
                   "exact_optimal",
                   "psm_p_greedy_noc",
                   "psm_p_greedy_caliper_01",
                   "psm_p_greedy_caliper_005",
                   "psm_m_greedy_noc",
                   "psm_m_greedy_caliper_01",
                   "psm_m_greedy_caliper_005",
                   "psm_p_optimal","psm_m_optimal")
# set number of repetitions for iterations and bootstrap
n_reps=50; b_reps=20
tasks=grid[rep(1:nrow(grid),each=n_reps),]
tasks$rep=rep(1:n_reps,times=nrow(grid))
total_steps=nrow(tasks)
# register components to run in parallel
options(future.globals.maxSize = 2 * 1024^3)
registerDoFuture()
plan(multisession, workers = parallel::detectCores()-1)
handlers("progress")
registerDoRNG(181299)
# time 
start.time=proc.time()
# run loop
with_progress({
  # define number of steps
  p=progressor(steps = total_steps)
  # define results object
  results=foreach(t=1:nrow(tasks), .combine = rbind,
                     .packages = c("survival","MatchIt","boot",
                                   "dplyr","smd","WeightIt",
                                   "optmatch","lpSolve","randomForest",
                                   "lme4","emmeans"),
                     .export = c("df",
                                 "matching","psm_matching","risk_ratio",
                                 "trial_sampler","bootstrap_sampler",
                                 "analyse_call",
                                 "matching_call","regression_call",
                                 "run_matching_methods","run_regression_methods",
                                 "run_all_matching","run_all_regression")) %dopar% {
      # define number of treated, untreated and ratio
      ntrti=tasks$ntrt[t]
      nctrli=tasks$nctrl[t]
      # print to see where we are
      p(sprintf("ntrt=%d, nctrl=%d, i=%d",ntrti,nctrli,tasks$rep[t]))
      # sample
      dftrial=trial_sampler(df,id,trt_overall,ntrti,nctrli) %>% 
        group_by(id) %>% 
        mutate(id=cur_group_id()) %>% 
        ungroup()
      # baseline
      baseline=dftrial %>% 
        group_by(id) %>% 
        filter(row_number()==1) %>% 
        ungroup()
      # smd before
      smdageb=smd(baseline$aoo,baseline$trt_overall)$estimate
      smdsexb=smd(baseline$sex,baseline$trt_overall)$estimate
      # available ratios=ratio
      available=ratio[ntrti*ratio<=nctrli]
      # run for point estimates
      ## matching
      point_results=bind_rows(
        lapply(available,function(ratioi){
          run_all_matching(data=dftrial,
                           ratio=ratioi,
                           id=id,
                           methods=matching_methods,
                           landmark=2)
        })
      )
      # bootstrap
      ## initialize list
      bootstrap_results_matching=vector("list",b_reps)
      ## loop
      for (b in seq_len(b_reps)){
        # perform bootstrap
        boot_data=bootstrap_sampler(dftrial,id,trt_overall,ntrti,nctrli)
        # run all method on the same bootstrapped dataset
        ## matching
        bootstrap_results_matching[[b]]=bind_rows(
          lapply(available,function(ratioi){
            run_all_matching(data=boot_data,
                             ratio=ratioi,
                             id=id,
                             methods=matching_methods,
                             landmark=2)
          })
        ) %>% mutate(boot=b)
      }
      # bind_together results from bootstrap,for both methods
      bootstrap_results=bind_rows(bootstrap_results_matching)
      ## no bootstrap (b_reps=0): empty table so that se is NA
      if (nrow(bootstrap_results)==0) bootstrap_results=data.frame(method=character(),ratio=numeric(),result=numeric(),converged=logical())
      # obtain standard error
      se_results=bootstrap_results %>%
        group_by(method,ratio) %>%
        summarise(se=sd(result,na.rm=TRUE),
                  mean_converged=mean(converged),
                  .groups = "drop")
      # obtain standard error
        point_results %>%
          left_join(se_results,
                    by=c("method","ratio")) %>%
          mutate(sim=t,
                 ntrt=ntrti,
                 nctrl=nctrli,
                 smdageb=smdageb,
                 smdsexb=smdsexb)
        }
})
# time? 
elapsed=proc.time()-start.time
elapsed; elapsed["elapsed"]


## ---------------------------------------------------------------------------------------------------------------------
saveRDS(results,paste0("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Data\\paper\\sensitivity_analysis\\missing_data\\with duration\\results_diff_0.15_",mechanism,".rds"))


## ---------------------------------------------------------------------------------------------------------------------
knitr::purl("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project 5 - Simulation study in SATs\\Code\\paper\\R\\sensitivity_analyses\\continuous\\missing_data\\with_duration\\analysis\\analysis_diff_0.15.qmd")

