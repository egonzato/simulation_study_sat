#!/bin/bash
## =============================================================================================
## Slurm job array for analysis_diff_0.15_hpc.R (anchored continuous branch)
##
## Submit from a login node:   sbatch analysis_diff_0.15_hpc.slurm
## Check the jobs:             squeue -u $USER
## After all jobs finished:    Rscript combine_results_hpc.R 0.15 50
##
## To confirm with HELIOS support: the R module name, partition, time and memory limits,
## and where the project folder (SIM_ROOT) should live (home, T1 BeeGFS or T0 scratch).
## =============================================================================================
#SBATCH --job-name=anch_cont_0.15
## one array job per grid cell and chunk: 9 cells x ceiling(N_REPS / CHUNK_SIZE) chunks
#SBATCH --array=1-9
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
## the logs folder must exist before submitting (mkdir -p logs)
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
##SBATCH --partition=...         # to confirm with HELIOS support

## ---- run settings: change these, and --array above to 9 x number of chunks ------------------
EFFECT=0.15        # "0", "0.05", "0.15" or "0.3"
N_REPS=50          # repetitions per grid cell
B_REPS=20          # bootstrap samples per repetition
CHUNK_SIZE=50      # repetitions per array job (N_REPS here: one job per grid cell)

## ---- environment ------------------------------------------------------------------------------
module load R                  # exact module name to confirm with HELIOS support
# folder that contains Code/paper/R and Data/paper (the same layout as on OneDrive)
export SIM_ROOT=/path/to/project
# R packages installed once in a personal library, e.g. ~/R/library:
# Rscript -e 'install.packages(c("dplyr","lme4","emmeans","smd","optmatch","foreach","doRNG","future","doFuture","parallelly"), lib="~/R/library")'
export R_LIBS_USER=~/R/library

SCRIPT="$SIM_ROOT/Code/paper/R/main_analysis/without duration/continuous/analysis/analysis_diff_0.15_hpc.R"
srun Rscript "$SCRIPT" "$EFFECT" "$N_REPS" "$B_REPS" "$CHUNK_SIZE"
