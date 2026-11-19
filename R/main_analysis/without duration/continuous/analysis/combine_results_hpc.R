## =============================================================================================
## Combine the results of the array jobs of analysis_diff_0.15_hpc.R into one file
##
## Usage: Rscript combine_results_hpc.R <effect> <n_reps>
##   checks that every grid cell has n_reps repetitions and saves results_diff_<effect>_hpc.rds,
##   which report.qmd can read like the results of the laptop runs
## =============================================================================================
suppressPackageStartupMessages(library(dplyr))
args=commandArgs(trailingOnly=TRUE)
effect=if (length(args)>=1) args[1] else "0.15"
n_reps=if (length(args)>=2) as.integer(args[2]) else NA
root=Sys.getenv("SIM_ROOT",unset="C:/Users/P095206/OneDrive - Amsterdam UMC/Shared material with Elia/PhD/Project 5 - Simulation study in SATs")
data_dir=file.path(root,"Data","paper","main_analysis","without duration","continuous")
out_dir=file.path(data_dir,"hpc")

files=list.files(out_dir,pattern=paste0("^results_diff_",gsub(".","\\.",effect,fixed=TRUE),"_cell[0-9]+_chunk[0-9]+\\.rds$"),full.names=TRUE)
if (length(files)==0) stop("no results files found in ",out_dir)
results=bind_rows(lapply(files,readRDS))
cat(sprintf("%d files, %d rows\n",length(files),nrow(results)))

# repetitions per sample-size scenario: every scenario should have n_reps distinct repetitions
check=results %>% group_by(ntrt,nctrl) %>% summarise(repetitions=n_distinct(sim),.groups="drop")
print(check)
if (!is.na(n_reps) && any(check$repetitions!=n_reps)) warning("some scenarios do not have ",n_reps," repetitions: check for failed array jobs")
if (any(duplicated(results[,c("sim","method","ratio")]))) warning("duplicated repetitions: the same part was run twice")

out_file=file.path(data_dir,paste0("results_diff_",effect,"_hpc.rds"))
saveRDS(results,out_file)
cat("saved:",out_file,"\n")
