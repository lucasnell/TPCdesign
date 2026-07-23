#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-40
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=06:00:00
#SBATCH --job-name=ace-tune
#SBATCH --output=logs/ace-tune-%a.out
#SBATCH --error=logs/ace-tune-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL


# IMPORTANT -- SET THIS FOR THE MAX ARRAY NUMBER ABOVE:
# (I'm doing it this way--instead of using existing SLURM_ARRAY_TASK_MAX object--
# so that I can restart a job with just one array index--
# e.g., `sbatch --array=X ace-tuning.sh`--if one job fails for a weird reason)
export MAX_ARRAY_INDEX=40



# I first moved this script and the preamble over to bioHPC using the following:
#
# cd ~/GitHub/Stanford/TPCdesign/_testing
# scp ace-tuning.* lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/tpc/tuning/
#
# This was then run on BioHPC in a non-interactive batch job started with the following:
#
# cd /home2/lan68/tpc/tuning/
# mkdir -p logs
# sbatch ace-tuning.sh
#
# Then, when the jobs are done (assuming you're back on your desktop in
# directory `~/GitHub/Stanford/TPCdesign/_testing`):
#
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/tpc/tuning/ace-tuning*.csv \
#     ./interm-data/ace-tuning/
#
#
# To test:
# cd /home2/lan68/tpc/tuning
# srun -N 1 -n 1 -c 6 --mem=20G --time=2:00:00 --job-name="ace-tune" --pty bash -i
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
    library(lhs)
    library(future.apply)
    library(progressr)
})

n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_threads)) stop("n_threads cannot be NA")

options("readr.num_threads" = n_threads)
plan(multisession, workers = n_threads)
handlers(handler_cli(format = paste("{cli::pb_bar} {cli::pb_percent} |",
                                        "{cli::pb_elapsed} | ETA: {cli::pb_eta}"),
                         clear = FALSE))


# Contains functions one_combo_fits and one_test_fit
source("ace-tuning.R")


set.seed(1981443946)
lhs_df <- maximinLHS(n = 1000, k = 9) |>
    as.data.frame() |>
    set_names(c("n_temps", "b", "ctmin_eps", "ctmax_eps",
                "logb_eps", "min_sep", "n_filler", "n_draws", "n_starts")) |>
    as_tibble() |>
    mutate(n_temps = n_temps |> qinteger(5, 10) |> as.integer(),
           b = c(0.2, 0.5, 1)[qinteger(b, 1, 3)],
           ctmin_eps = 2 * 2.5 * (2 * ctmin_eps - 1), # [-5, +5]
           ctmax_eps = 2 * 1.5 * (2 * ctmax_eps - 1), # [-3, +3]
           logb_eps = 2 * 0.26 * (2 * logb_eps - 1), # [-0.52, +0.52]
           min_sep = min_sep * (2 - 0.5) + 0.5,
           n_filler = n_filler |> qinteger(0, 3) |> as.integer(),
           n_draws = n_draws |> qinteger(50, 500) |> as.integer(),
           n_starts = n_starts |> qinteger(5, 10) |> as.integer()) |>
    # Variables that don't change:
    mutate(n_reps = 6L,
           obs_cv = 0.2,
           ctmin = 5,
           ctmax = 40,
           a = 1,
           # to keep track of specific combinations:
           combo = 1:n()) |>
    select(combo, everything()) |>
    # Now there's one row per combo + rep:
    crossing(tibble(rep = 1:100))





# --------------*
# Read inputs from job manager:

curr_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")

max_idx <- suppressWarnings(as.integer(Sys.getenv("MAX_ARRAY_INDEX")))
if (is.na(max_idx)) stop("MAX_ARRAY_INDEX must be an integer")
if (curr_idx > max_idx) stop("SLURM_ARRAY_TASK_ID cannot exceed MAX_ARRAY_INDEX")

# Get seed for this set of simulations:
seed <- c(1066802577, 1373261426,   51277529,  210376342,  920159360,
          1146973941, 1007152868,  534147735,  970437820,  142493451,
          1128770799,  777660544,  356311217, 1548414489,  253378587,
          1600791213, 1043286887, 1242410279,  230438103,  513892003,
           682731035, 1822575486, 1401602719, 1417269504, 1315811922,
          1302120418, 2083447359, 1747511432,  586618815, 1999373684,
          1496113079, 1761135690,  658172835, 1460341010, 1948846207,
          1669298813,  477760098, 1696519555,  848624250, 111549496)[curr_idx]

# Verify that indices align with 'lhs_df':
stopifnot(nrow(lhs_df) %% max_idx == 0L)

# number of rows to simulate per job:
n_rows_pj <- nrow(lhs_df) %/% max_idx

# Start and stop for rows to simulate this job:
curr_start <- (curr_idx - 1L) * n_rows_pj + 1L
curr_stop <- curr_idx * n_rows_pj


lhs_df <- lhs_df[curr_start:curr_stop,]





# --------------*
# Main simulator function:


all_combo_rep_rmse <- function(lhs_df, .seed) {

    n_rows <- nrow(lhs_df)

    with_progress({
            p <- progressor(n_rows)
            out <- future_lapply(1:n_rows,
                                 one_combo_fits,
                                 input_df = lhs_df,
                                 prog = p,
                                 future.seed = .seed,
                                 future.globals = c("one_combo_fits",
                                                    "one_fit_rmse"),
                                 future.packages = c("tidyverse",
                                                     "nls.multstart",
                                                     "rTPC",
                                                     "TPCdesign",
                                                     "acebayes",
                                                     "lhs")) |>
                list_c()
    })

    return(out)
}




# Took  min per job, where each job processed 2500 rows using 6 threads

cat("Starting simulations...\n")
t0 <- Sys.time()
lhs_df[["rmse"]] <- all_combo_rep_rmse(lhs_df, seed)
t1 <- Sys.time()
difftime(t1, t0, units = "min")

# --------------*


cat("Writing output...\n")


write_csv(lhs_df, sprintf("ace-tuning-%02i.csv", curr_idx))



cat("\nFINISHED!!\n\n")



EOF
