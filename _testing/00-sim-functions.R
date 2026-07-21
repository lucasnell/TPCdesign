
source("_testing/ace.R")

# Number of fits per combination of parameters:
n_test_fits <- 100L






# Coefficients from one fit for one combo of parameters
one_test_fit <- function(i, temps, n_reps, obs_cv, ctmin, ctmax, a, b,
                         .scale_tpc = FALSE) {


    obs <- sim_gamma_data(temps, n_reps, obs_cv, ctmin, ctmax, a, b,
                          scale_tpc = .scale_tpc)
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
                             fitted[["a"]], fitted[["b"]], .scale_tpc)
        tru_y <- briere2_tpc(test_temps, ctmin, ctmax, a, b, .scale_tpc)
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
one_combo_fits <- function(j, input_df, prog, .scale_tpc = FALSE) {

    # n_temps <- 6L; n_reps <- 5L
    # z0 <- 0; z1 <- 1; obs_cv <- 0.2; ctmin <- 10; ctmax <- 50
    # a <- 1; b <- 0.2; temp_buffer <- 6L
    # curve_off <- TRUE
    # rm(n_temps, n_reps, z0, z1, obs_cv, ctmin, ctmax, a, b, temp_buffer, curve_off)
    # rm(temps, coef_df)

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
    curve_off <- input_df[["curve_off"]][[j]]


    if (curve_off) {
        # genus-level relatives was 2.50°C for Tmin, 1.50°C for Tmax,
        # and 0.26 for log(b)
        temps <- map(1:n_test_fits, \(i) {
            .ctmin <- ctmin + sample(c(-1,1),1) * 2.50
            .ctmax <- ctmax + sample(c(-1,1),1) * 1.50
            .b <- exp(log(b) + sample(c(-1,1),1) * 0.26)
            deriv_design_temps(n_temps, z0, z1, .ctmin, .ctmax, a, .b,
                               temp_buffer)
        })
    } else {
        temps <- deriv_design_temps(n_temps, z0, z1, ctmin, ctmax, a, b,
                                    temp_buffer)
        temps <- rep(list(temps), n_test_fits)
    }


    coef_df <- map2(1:n_test_fits, temps,
                    \(i, temps) {
                        one_test_fit(i = i, temps = temps,
                                     n_reps = n_reps, obs_cv = obs_cv,
                                     ctmin = ctmin, ctmax = ctmax, a = a, b = b,
                                     .scale_tpc = .scale_tpc)
                    }) |>
        list_rbind() |>
        mutate(n_temps = .env$n_temps,
               n_reps = .env$n_reps,
               z0 = .env$z0,
               z1 = .env$z1,
               obs_cv = .env$obs_cv,
               temp_buffer = .env$temp_buffer,
               curve_off = .env$curve_off) |>
        mutate(combo = j) |>
        select(combo, rep, type, n_temps:curve_off, everything())

    if (!isTRUE(is.null(prog))) prog()

    return(coef_df)

}


# Do all fits for all combinations
all_combo_fits <- function(input_df, .seed, .scale_tpc = FALSE) {

    suppressPackageStartupMessages(library(future.apply))
    suppressPackageStartupMessages(library(progressr))
    # handlers("progress")
    handlers(handler_cli(format = paste("{cli::pb_bar} {cli::pb_percent} |",
                                        "{cli::pb_elapsed} | ETA: {cli::pb_eta}"),
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
                             .scale_tpc = .scale_tpc,
                             future.seed = .seed,
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


