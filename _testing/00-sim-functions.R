
source("_testing/ace.R")








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
one_combo_fits <- function(j, input_df, prog,
                           .scale_tpc = FALSE,
                           .ctmin_err = 2.50,
                           .ctmax_err = 1.50,
                           .logb_err = 0.26,
                           .n_test_fits = 50L) {

    # j = 1L; n_temps = 6L; n_reps = 5L; obs_cv = 0.2; ctmin = 10; ctmax = 50
    # a = 1; b = 0.2
    # .scale_tpc = FALSE; .ctmin_err = 2.50; .ctmax_err = 1.50
    # .logb_err = 0.26; .n_test_fits = 100L
    #
    # rm(n_temps, n_reps, obs_cv, ctmin, ctmax, a, b)
    # rm(.scale_tpc, .ctmin_err, .ctmax_err, .logb_err, .n_test_fits)
    # rm(temps, coef_df)

    n_temps <- input_df[["n_temps"]][[j]]
    n_reps <- input_df[["n_reps"]][[j]]
    obs_cv <- input_df[["obs_cv"]][[j]]
    ctmin <- input_df[["ctmin"]][[j]]
    ctmax <- input_df[["ctmax"]][[j]]
    a <- input_df[["a"]][[j]]
    b <- input_df[["b"]][[j]]


    # genus-level relatives was 2.50°C for Tmin, 1.50°C for Tmax,
    # and 0.26 for log(b)
    temps <- map(1:.n_test_fits, \(i) {
        .ctmin <- ctmin + sample(c(-1,1),1) * .ctmin_err
        .ctmax <- ctmax + sample(c(-1,1),1) * .ctmax_err
        .b <- exp(log(b) + sample(c(-1,1),1) * .logb_err)
        dtemps <- ace_design_temps(n_temps, .ctmin, .ctmax, a, .b)
        etemps <- seq(.ctmin-5, .ctmax+5, length.out = n_temps+2L) |>
            head(-1) |> tail(-1) |> round(2)
        return(list(design = dtemps, even = etemps))
    })



    coef_df <- map2(1:.n_test_fits, temps,
                    \(i, temp) {
                        imap(temp, \(tmp, m) {
                            one_test_fit(i = i, temps = tmp,
                                         n_reps = n_reps, obs_cv = obs_cv,
                                         ctmin = ctmin, ctmax = ctmax, a = a, b = b,
                                         .scale_tpc = .scale_tpc) |>
                                mutate(method = m)
                        }) |>
                            list_rbind() |>
                            select(rep, method, type, everything()) |>
                            filter(!(method == "design" & type == "real")) |>
                            mutate(method = ifelse(type == "real", NA, method))
                    }) |>
        list_rbind() |>
        mutate(n_temps = .env$n_temps,
               n_reps = .env$n_reps,
               obs_cv = .env$obs_cv) |>
        mutate(combo = j) |>
        select(combo, rep, method, type, n_temps:obs_cv, everything())

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
                                                "ace_design_temps",
                                                "gap_filler"),
                             future.packages = c("tidyverse",
                                                 "nls.multstart",
                                                 "TPCdesign",
                                                 "acebayes",
                                                 "lhs")) |>
            list_rbind()
    })
    return(out)
}


