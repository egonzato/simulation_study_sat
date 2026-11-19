## =============================================================================================
## Anchored continuous branch: analysis for the HPC (HELIOS), also runs on the laptop
##
## Same analysis as analysis_diff_0.15.qmd. What differs:
##  - paths are built from one project root (SIM_ROOT), so the script runs on Windows and Linux
##  - settings come from the command line, so the same script serves every run
##  - the cores are the ones Slurm gives to the job, not all the cores of the node
##  - with a Slurm job array, every array job runs one part of the grid and saves its own file;
##    combine_results_hpc.R puts the parts together
##
## Usage
##   laptop (RStudio): set the defaults below and source the file; with no Slurm array it runs the
##     whole grid, with the same seed as analysis_diff_0.15.qmd
##   terminal / HPC:   Rscript analysis_diff_0.15_hpc.R <effect> <n_reps> <b_reps> <chunk_size>
##     effect     : effect size in the file names of the populations ("0", "0.05", "0.15", "0.3")
##     n_reps     : repetitions per sample-size scenario (grid cell)
##     b_reps     : bootstrap samples per repetition
##     chunk_size : repetitions per array job (default n_reps: one array job per grid cell)
##   Slurm array:      job k runs grid cell ceiling(k / n_chunks), chunk ((k-1) %% n_chunks) + 1,
##                     with n_chunks = ceiling(n_reps / chunk_size); --array=1-(9 * n_chunks)
## =============================================================================================

## ---- packages (install them once on the HPC, see the Slurm script) --------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(emmeans)
  library(smd)
  library(optmatch)
  library(foreach)
  library(doRNG)
  library(future)
  library(doFuture)
  library(parallelly)
})

## ---- settings ---------------------------------------------------------------------------------
args=commandArgs(trailingOnly=TRUE)
effect    =if (length(args)>=1) args[1] else "0.15"
n_reps    =if (length(args)>=2) as.integer(args[2]) else 10
b_reps    =if (length(args)>=3) as.integer(args[3]) else 1
chunk_size=if (length(args)>=4) as.integer(args[4]) else n_reps
# project root: on the HPC set in the Slurm script (export SIM_ROOT=...), on the laptop the OneDrive folder
root=Sys.getenv("SIM_ROOT",unset="C:/Users/P095206/OneDrive - Amsterdam UMC/Shared material with Elia/PhD/Project 5 - Simulation study in SATs")
code_dir=file.path(root,"Code","paper","R")
data_dir=file.path(root,"Data","paper","main_analysis","without duration","continuous")
out_dir =file.path(data_dir,"hpc")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
# base seed, the same as analysis_diff_0.15.qmd
seed=11022004

## ---- which part of the grid this job runs -----------------------------------------------------
grid=expand.grid(ntrt=c(20,30,40),nctrl=c(50,100,150))
n_chunks=ceiling(n_reps/chunk_size)
array_id=as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",unset="0"))
if (array_id==0){
  # no job array: the whole grid in one go (laptop), same seed as the original analysis
  cells=seq_len(nrow(grid)); reps=seq_len(n_reps); job_seed=seed
  out_file=file.path(out_dir,paste0("results_diff_",effect,"_all.rds"))
} else {
  if (array_id>nrow(grid)*n_chunks) stop("SLURM_ARRAY_TASK_ID larger than 9 x number of chunks")
  cell=ceiling(array_id/n_chunks)
  chunk=(array_id-1)%%n_chunks+1
  cells=cell
  reps=((chunk-1)*chunk_size+1):min(chunk*chunk_size,n_reps)
  # one seed per grid cell and chunk: the same part always gives the same results
  job_seed=seed+1000*cell+chunk
  out_file=file.path(out_dir,sprintf("results_diff_%s_cell%d_chunk%03d.rds",effect,cell,chunk))
}
cat(sprintf("effect %s | cells %s | repetitions %d-%d | bootstrap %d | seed %d\n",
            effect,paste(cells,collapse=","),min(reps),max(reps),b_reps,job_seed))

## ---- data and functions -----------------------------------------------------------------------
df=readRDS(file.path(data_dir,paste0("population_diff_",effect,".rds"))) %>%
  filter(timen<=2) %>%
  select(id,timen,trt_overall,trt,aoo,sex,hui_trt,yoo) %>%
  group_by(id,trt) %>%
  filter(!(trt==1 & !(row_number() %in% c(1,2,3,5,7,9)))) %>%
  ungroup()
for (d in c(file.path(code_dir,"functions","wrapper","continuous","without duration"),
            file.path(code_dir,"functions","matching"),
            file.path(code_dir,"functions","regression"))){
  for (f in list.files(d,pattern="\\.R$",full.names=TRUE)) source(f)
}

## ---- methods ----------------------------------------------------------------------------------
ratio=c(1,2,3)
regression_methods=c("raw",
                     "iptw, unstabilised",
                     "iptw, stabilised",
                     "iptw, stabilised, truncated 0.01",
                     "iptw, stabilised, truncated 0.1",
                     "gcomp")
matching_methods=c("exact_greedy",
                   "exact_optimal",
                   "psm_p_greedy_noc",
                   "psm_p_greedy_caliper_01",
                   "psm_p_greedy_caliper_005",
                   "psm_m_greedy_noc",
                   "psm_m_greedy_caliper_01","psm_m_greedy_caliper_005",
                   "psm_p_optimal",
                   "psm_m_optimal")

## ---- tasks --------------------------------------------------------------------------------------
# the same numbering of repetitions (sim) as in the full grid, whichever part this job runs
tasks=expand.grid(rep=reps,cell=cells) %>%
  mutate(ntrt=grid$ntrt[cell],nctrl=grid$nctrl[cell],sim=(cell-1)*n_reps+rep)

## ---- parallel backend: the cores of this job ----------------------------------------------------
# availableCores() follows SLURM_CPUS_PER_TASK on the HPC; on the laptop it leaves one core free
workers=if (array_id==0) max(1,availableCores()-1) else availableCores()
options(future.globals.maxSize=2*1024^3)
registerDoFuture()
if (.Platform$OS.type=="unix") plan(multicore,workers=workers) else plan(multisession,workers=workers)
registerDoRNG(job_seed)
cat(sprintf("workers: %d\n",workers))

## ---- run ----------------------------------------------------------------------------------------
start.time=proc.time()
results=foreach(t=seq_len(nrow(tasks)),.combine=rbind,
                .packages=c("dplyr","smd","optmatch","lme4","emmeans")) %dopar% {
  ntrti=tasks$ntrt[t]
  nctrli=tasks$nctrl[t]
  # sample
  dftrial=trial_sampler(df,id,trt,ntrti,nctrli) %>%
    group_by(id) %>%
    mutate(id=cur_group_id())
  # baseline
  baseline=dftrial %>%
    group_by(id) %>%
    filter(row_number()==1) %>%
    ungroup()
  # smd before
  smdageb=smd(baseline$aoo,baseline$trt)$estimate
  smdsexb=smd(baseline$sex,baseline$trt)$estimate
  # available ratios
  available=ratio[ntrti*ratio<=nctrli]
  # point estimates
  point_regression=bind_rows(run_all_regression(data=dftrial,methods=regression_methods,landmark=2))
  point_matching=bind_rows(lapply(available,function(ratioi){
    run_all_matching(data=dftrial,ratio=ratioi,methods=matching_methods,landmark=2)
  }))
  point_results=bind_rows(point_regression,point_matching)
  # bootstrap
  bootstrap_results_regression=vector("list",b_reps)
  bootstrap_results_matching=vector("list",b_reps)
  for (b in seq_len(b_reps)){
    boot_data=bootstrap_sampler(dftrial,id,trt,ntrti,nctrli)
    bootstrap_results_regression[[b]]=bind_rows(
      run_all_regression(data=boot_data,id=id,methods=regression_methods,landmark=2))
    bootstrap_results_matching[[b]]=bind_rows(lapply(available,function(ratioi){
      run_all_matching(data=boot_data,ratio=ratioi,id=id,methods=matching_methods,landmark=2)
    })) %>% mutate(boot=b)
  }
  bootstrap_results=bind_rows(bind_rows(bootstrap_results_regression),bind_rows(bootstrap_results_matching))
  # no bootstrap (b_reps=0): empty table so that se is NA
  if (nrow(bootstrap_results)==0) bootstrap_results=data.frame(method=character(),ratio=numeric(),result=numeric(),converged=logical())
  se_results=bootstrap_results %>%
    group_by(method,ratio) %>%
    summarise(se=sd(result,na.rm=TRUE),
              mean_converged=mean(converged),
              .groups="drop")
  point_results %>%
    left_join(se_results,by=c("method","ratio")) %>%
    mutate(sim=tasks$sim[t],
           ntrt=ntrti,
           nctrl=nctrli,
           smdageb=smdageb,
           smdsexb=smdsexb)
}
elapsed=proc.time()-start.time
cat(sprintf("done: %d repetitions in %.1f minutes\n",nrow(tasks),elapsed["elapsed"]/60))

## ---- save ---------------------------------------------------------------------------------------
saveRDS(results,out_file)
cat("saved:",out_file,"\n")
