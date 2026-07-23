#!/bin/bash -l



#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-36
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --job-name=tpc
#SBATCH --output=logs/tpc-%a.out
#SBATCH --error=logs/tpc-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL


# IMPORTANT -- SET THIS FOR THE MAX ARRAY NUMBER ABOVE:
# (I'm doing it this way--instead of using existing SLURM_ARRAY_TASK_MAX object--
# so that I can restart a job with just one array index--
# e.g., `sbatch --array=X big-ace-test.sh`--if one job fails for a weird reason)
export MAX_ARRAY_INDEX=36



# I first moved this script and the preamble over to bioHPC using the following:
#
# cd ~/GitHub/Stanford/TPCdesign/_testing
# scp big-ace-test.R big-ace-test.sh \
#   lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/tpc/testing/
#
# This was then run on BioHPC in a non-interactive batch job started with the following:
#
# cd /home2/lan68/tpc/testing/
# mkdir -p logs
# sbatch big-ace-test.sh
#
# Then, when the jobs are done (assuming you're back on your desktop in
# directory `~/GitHub/Stanford/TPCdesign/_testing`):
#
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/tpc/testing/big-ace-test*.csv \
#     ./interm-data/ace-test
#
#
# To test:
# cd /home2/lan68/tpc/testing/
# srun -N 1 -n 1 -c 4 --mem=20G --time=2:00:00 --job-name="ace" --pty bash -i
# export OMP_NUM_THREADS=1
# R --vanilla
#
#


# I don't need OpenMP multithreading:
export OMP_NUM_THREADS=1




Rscript - << EOF

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(tidyverse)
    library(nls.multstart)
    library(rTPC)
    library(TPCdesign)
    library(acebayes)
    library(lhs)
    library(future.apply)
})

n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_threads)) stop("n_threads cannot be NA")

options("readr.num_threads" = n_threads)
plan(multisession, workers = n_threads)



# Contains functions one_combo_fits and one_test_fit
source("big-ace-test.R")



input_df <- crossing(n_temps = 5:10,
                   n_reps = 3:10,
                   b = c(0.2, 0.5, 1, 2),
                   obs_cv = 0.2 * c(0.5, 1, 2),
                   # Things that don't vary:
                   ctmin = 5,
                   ctmax = 40,
                   a = 1)


# --------------*
# Read inputs from job manager:

curr_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")

max_idx <- suppressWarnings(as.integer(Sys.getenv("MAX_ARRAY_INDEX")))
if (is.na(max_idx)) stop("MAX_ARRAY_INDEX must be an integer")
if (curr_idx > max_idx) stop("SLURM_ARRAY_TASK_ID cannot exceed MAX_ARRAY_INDEX")

# Get seed for this set of simulations:
seed <- c(1708876859,  771729491,  106836587, 1666747557, 1835919150,
          1650012022, 1550552152, 1104661569, 1339491109,  641105046,
          1451392045,  751711886, 2022559069,  258770169, 1085677670,
           774475743, 2074678884,  775811129, 1073723116,  652825051,
           455347514,  416484533, 1596904740, 2004673084,  159899192,
           757246500,  621837632, 1433482434, 1710243287, 1988468718,
          1625207889,  805328827,  796557836, 1470116913, 2111356371,
           712313943)[curr_idx]

# Verify that indices align with 'input_df':
stopifnot(nrow(input_df) %% max_idx == 0L)

# number of rows to simulate per job:
n_rows_pj <- nrow(input_df) %/% max_idx

# Start and stop for rows to simulate this job:
curr_start <- (curr_idx - 1L) * n_rows_pj + 1L
curr_stop <- curr_idx * n_rows_pj


indices <- curr_start:curr_stop





# --------------*
# Main simulator function:


all_combo_fits <- function(indices, input_df, .seed) {

    out <- future_lapply(indices,
                         one_combo_fits,
                         input_df = input_df,
                         future.seed = .seed,
                         future.globals = c("one_combo_fits",
                                            "one_test_fit"),
                         future.packages = c("tidyverse",
                                             "nls.multstart",
                                             "TPCdesign",
                                             "acebayes",
                                             "lhs")) |>
        list_rbind()

    return(out)
}


# Took < 90 min per job, where each job processed 16 rows using 4 threads

cat("Starting simulations...\n")
t0 <- Sys.time()
out_df <- all_combo_fits(indices, input_df, seed)
t1 <- Sys.time()
difftime(t1, t0, units = "min")

# --------------*


cat("Writing output...\n")


write_csv(out_df, sprintf("big-ace-test-%02i.csv", curr_idx))



cat("\nFINISHED!!\n\n")



EOF
