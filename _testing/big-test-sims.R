
source("_scripts/00-preamble.R")

library(lhs)



# Number of fits per combination of parameters:
n_test_fits <- 100L






# Coefficients from one fit for one combo of parameters
one_test_fit <- function(i, temps, n_reps, obs_cv, ctmin, ctmax, a, b) {


    obs <- sim_gamma_data(temps, n_reps, obs_cv, ctmin, ctmax, a, b,
                          scale_tpc = FALSE)
    fit <- nls_multstart(
        formula     = y ~ a * temp * (temp - ctmin) *
            pmax(ctmax - temp, 0)^exp(b),
        data        = obs,
        start_lower = c(a = 0,  ctmin = 0,  ctmax = 30, b = log(0.01)),
        start_upper = c(a = 2,  ctmin = 15, ctmax = 50, b = log(3)),
        iter        = 50,
        supp_errors = "Y",
        control = list(maxfev = 5e3, maxiter = 1e3),
        lhstype = "improved")
    if (is.null(fit)) {
        fitted <- tibble(ctmin = NA_real_, ctmax = NA_real_,
                         a = NA_real_, b = NA_real_, Topt = NA_real_,
                         converged = FALSE, rmse = NA_real_)
    } else {
        fitted <- as_tibble(as.list(coef(fit))[c("ctmin", "ctmax", "a", "b")])
        fitted[["b"]] <- exp(fitted[["b"]])
        fitted[["Topt"]] <- briere2_tpc_Topt(fitted[["ctmin"]], fitted[["ctmax"]],
                                             fitted[["b"]])
        fitted[["converged"]] <- fit$convInfo$isConv
        # Add RMSE
        test_temps <- seq(ctmin, ctmax, length.out = 101)
        obs_y <- briere2_tpc(test_temps, fitted[["ctmin"]], fitted[["ctmax"]],
                             fitted[["a"]], fitted[["b"]], FALSE)
        tru_y <- briere2_tpc(test_temps, ctmin, ctmax, a, b, FALSE)
        fitted[["rmse"]] <- sqrt(mean((obs_y - tru_y)^2))
    }

    fitted[["rep"]] <- i
    fitted[["type"]] <- "fitted"

    real <- tibble(ctmin = .env$ctmin, ctmax = .env$ctmax,
                   a = .env$a, b = .env$b)
    real[["Topt"]] <- briere2_tpc_Topt(ctmin, ctmax, b)
    real[["converged"]] <- NA
    real[["rep"]] <- i
    real[["type"]] <- "real"
    real[["rmse"]] <- NA_real_

    return(bind_rows(fitted, real))
}


# Do all fits for a single combination of parameters
one_combo_fits <- function(j, input_df, prog) {

    n_temps <- input_df[["n_temps"]][[j]]
    n_reps <- input_df[["n_reps"]][[j]]
    z0 <- input_df[["z0"]][[j]]
    z1 <- input_df[["z1"]][[j]]
    obs_cv <- input_df[["obs_cv"]][[j]]
    ctmin <- input_df[["ctmin"]][[j]]
    ctmax <- input_df[["ctmax"]][[j]]
    a <- input_df[["a"]][[j]]
    b <- input_df[["b"]][[j]]
    temp_buffer <- input_df[["temp_buffer"]][[j]]

    temps <- deriv_design_temps(n_temps, z0, z1, ctmin, ctmax, a, b,
                                temp_buffer)

    coef_df <- map(1:n_test_fits, one_test_fit, temps = temps, n_reps = n_reps,
                    obs_cv = obs_cv, ctmin = ctmin, ctmax = ctmax, a = a, b = b) |>
        list_rbind() |>
        mutate(n_temps = .env$n_temps,
               n_reps = .env$n_reps,
               z0 = .env$z0,
               z1 = .env$z1,
               obs_cv = .env$obs_cv,
               temp_buffer = .env$temp_buffer) |>
        mutate(combo = j) |>
        select(combo, rep, type, n_temps:temp_buffer, everything())

    if (!isTRUE(is.null(prog))) prog()

    return(coef_df)

}


# Do all fits for all combinations
all_combo_fits <- function(input_df) {

    suppressPackageStartupMessages(library(future.apply))
    suppressPackageStartupMessages(library(progressr))
    # handlers("progress")
    handlers(handler_cli(format = "{cli::pb_bar} {cli::pb_percent} | {cli::pb_elapsed} | ETA: {cli::pb_eta}",
                         clear = FALSE))
    with(plan(multisession, workers = options()[["mc.cores"]], gc = TRUE),
         local = TRUE)

    n_rows  <- nrow(input_df)
    # stopifnot(n_rows <= nrow(input_df))

    with_progress({
        p <- progressor(n_rows)
        out <- future_lapply(1:n_rows,
                             one_combo_fits,
                             prog = p,
                             input_df = input_df,
                             future.seed = 1634404298,
                             future.globals = c("one_combo_fits",
                                                "one_test_fit",
                                                "n_test_fits"),
                             future.packages = c("tidyverse",
                                                 "nls.multstart",
                                                 "TPCdesign")) |>
            list_rbind()
    })
    return(out)
}





# =============================================================================*
# =============================================================================*
# Simulations
# =============================================================================*
# =============================================================================*


#
# Factors that I'll vary by simulation:
#
par_ranges <- list(n_temps = c(5, 10),
                   n_reps = c(3, 10),
                   z1 = c(0, 1),
                   b = c(0.2, 2),
                   obs_cv = 0.2 * c(0.5, 2),
                   temp_buffer = c(1, 10))

# These need to be treated differently:
int_pars <- c("n_temps", "n_reps", "temp_buffer")

# Create LHS data frame, then transform variables to proper scales
set.seed(191621615)
lhs_df <- improvedLHS(n = 1e3, k = length(par_ranges)) |>
    (\(x) {colnames(x) <- names(par_ranges); return(x)})() |>
    as_tibble()
for (n in names(par_ranges)) {
    p1 <- min(par_ranges[[n]])
    p2 <- max(par_ranges[[n]])
    if (! n %in% int_pars) {
        lhs_df[[n]] <- p1 + lhs_df[[n]] * (p2 - p1)
    } else {
        lhs_df[[n]] <- qinteger(lhs_df[[n]], p1, p2)
    }
}; rm(n, p1, p2)

# Now add variables that won't vary:
lhs_df[["z0"]] <- 0
lhs_df[["ctmin"]] <- 5
lhs_df[["ctmax"]] <- 40
lhs_df[["a"]] <- 1



if (file.exists("_testing/test-fits.csv")) {
    # Takes ~40 min
    fit_df <- all_combo_fits(lhs_df)
    write_csv(fit_df, "_testing/test-fits.csv")
} else {
    fit_df <- read_csv("_testing/test-fits.csv", col_types = "iiciidddiddddd")
}


#
# Lots of noise here, but we can conclude a couple of things:
#
# 1. z1 = 1 produces lowest rmse
# 2. temp_buffer = 6 because of the following:
#     a. higher values produce low levels of non-convergence but are more
#        likely to have divergent (especially high RMSE) fits
#     b. 6 gives good balance of relatively high convergence and low divergence
#     c. Note: this parameter has little overall effect on rmse
#

fit_df |>
    filter(type == "fitted") |>
    group_by(combo, n_temps, n_reps, z1, obs_cv, temp_buffer) |>
    summarize(rmse = mean(log10(rmse)), .groups = "drop") |>
    ggplot(aes(z1, (rmse))) +
    geom_jitter(aes(), shape = 1) +
    stat_smooth(method = "gam", formula = y ~ s(x)) +
    scale_color_viridis_c()

fit_df |>
    filter(type == "fitted") |>
    mutate(temp_buffer = factor(temp_buffer)) |>
    group_by(combo, n_temps, n_reps, z1, obs_cv, temp_buffer) |>
    summarize(rmse = sd((rmse)) / mean(rmse), .groups = "drop") |>
    ggplot(aes(temp_buffer, rmse)) +
    geom_jitter(aes(), shape = 1) +
    stat_summary(fun = "mean", geom = "point", color = "red", size = 3) +
    scale_color_viridis_c()

fit_df |>
    filter(type == "fitted") |>
    group_by(combo, n_temps, n_reps, z1, obs_cv, temp_buffer) |>
    summarize(p_conv = mean(converged), .groups = "drop") |>
    mutate(temp_buffer = factor(temp_buffer)) |>
    ggplot(aes(temp_buffer, p_conv)) +
    geom_jitter(aes(), shape = 1) +
    stat_summary(fun = "mean", geom = "point", color = "red", size = 3) +
    scale_color_viridis_c()


