#'
#' Preamble for big ACE test
#' Meant to be run on cluster
#'
#'





# Coefficients from one fit for one combo of parameters
one_test_fit <- function(i, temps, n_reps, obs_cv, ctmin, ctmax, a, b) {


    obs <- sim_gamma_data(temps, n_reps, obs_cv, ctmin, ctmax, a, b)

    starts_lo <- c(a = log(1e-6),  ctmin = -5,  ctmax = 30, b = log(0.01))
    starts_up <- c(a = log(2),  ctmin = 15, ctmax = 50, b = log(3))

    fit <- nls_multstart(
        formula     = y ~ exp(a) * temp * pmax(temp - ctmin, 0) *
            pmax(ctmax - temp, 0)^exp(b),
        data        = obs,
        start_lower = starts_lo,
        start_upper = starts_up,
        supp_errors = "Y",
        control = list(maxfev = 5e3, maxiter = 1e3),
        iter        = 500)

    if (is.null(fit)) {
        fitted <- tibble(ctmin = NA_real_, ctmax = NA_real_,
                         a = NA_real_, b = NA_real_, Topt = NA_real_,
                         converged = FALSE, rmse = NA_real_)
    } else {
        fitted <- as_tibble(as.list(coef(fit))[c("ctmin", "ctmax", "a", "b")])
        for (x in c("a", "b")) fitted[[x]] <- exp(fitted[[x]])
        fitted[["Topt"]] <- briere2_tpc_Topt(fitted[["ctmin"]], fitted[["ctmax"]],
                                             fitted[["b"]])
        fitted[["converged"]] <- fit$convInfo$isConv
        # Add RMSE
        test_temps <- seq(ctmin, ctmax, length.out = 101)
        obs_y <- briere2_tpc(test_temps, fitted[["ctmin"]], fitted[["ctmax"]],
                             fitted[["a"]], fitted[["b"]])
        tru_y <- briere2_tpc(test_temps, ctmin, ctmax, a, b)
        fitted[["rmse"]] <- sqrt(mean((obs_y - tru_y)^2))
    }

    fitted[["rep"]] <- i
    return(fitted)
}




# Do all fits for a single combination of parameters
one_combo_fits <- function(j, input_df,
                           .ctmin_err = 2.50,
                           .ctmax_err = 1.50,
                           .logb_err = 0.26,
                           .n_test_fits = 100L,
                           .parent_prog = function(...) NULL,
                           .parallel = FALSE,
                           .parallel_seed = TRUE) {

    # j = 7L; .ctmin_err = 2.50; .ctmax_err = 1.50; .logb_err = 0.26
    # .n_test_fits = 30L
    # rm(j, input_df, .ctmin_err, .ctmax_err, .logb_err, .n_test_fits)
    # rm(ace_args, x, .n_par_fits, par_lhs_list, temps, Topt, fit_df)

    # Create argument list for extra parameters for ace_design_temps:
    ace_args <- slice(input_df, j) |> as.list()
    for (x in c("combo", "rep", "rmse")) ace_args[[x]] <- NULL
    for (x in c("n_temps", "n_reps", "obs_cv", "ctmin", "ctmax", "a", "b")) {
        assign(x, ace_args[[x]])
        # These are needed for ace_design_temps (and are not modified):
        if (!x %in% c("n_temps", "a")) ace_args[[x]] <- NULL
    }

    Topt <- briere2_tpc_Topt(ctmin, ctmax, b) # used in output

    if (!identical(as.integer(.n_test_fits) %% 10L, 0L)) {
        stop(".n_test_fits must be a multiple of 10")
    }
    .n_par_fits <- as.integer(.n_test_fits) %/% 10L
    par_lhs_list <- lhs::optimumLHS(n = .n_par_fits, k = 3) |>
        as.data.frame() |>
        set_names(c("ctmin_eps", "ctmax_eps", "logb_eps")) |>
        mutate(ctmin_eps = 2 * 2.5 * (2 * ctmin_eps - 1),    # [-5, +5]
               ctmax_eps = 2 * 1.5 * (2 * ctmax_eps - 1),    # [-3, +3]
               logb_eps = 2 * 0.26 * (2 * logb_eps - 1)) |>  # [-0.52, +0.52]
        asplit(1) |>
        map(as.list)


    one_temp_fun <- function(i) {
        eps  <- par_lhs_list[[(i-1L) %/% 10L + 1L]]
        args <- ace_args
        args[["ctmin"]] <- ctmin + eps[["ctmin_eps"]]
        args[["ctmax"]] <- ctmax + eps[["ctmax_eps"]]
        args[["b"]] <- exp(log(b) + eps[["logb_eps"]])
        dtemps <- do.call(ace_design_temps, args)
        utemps <- seq(args$ctmin-5, args$ctmax+5, length.out = n_temps+2L) |>
            head(-1) |> tail(-1) |> round(2)

        dfit <- one_test_fit(i = i, temps = dtemps,
                             n_reps = n_reps, obs_cv = obs_cv,
                             ctmin = ctmin, ctmax = ctmax,
                             a = a, b = b) |>
            mutate(method = "design")
        ufit <- one_test_fit(i = i, temps = utemps,
                             n_reps = n_reps, obs_cv = obs_cv,
                             ctmin = ctmin, ctmax = ctmax,
                             a = a, b = b) |>
            mutate(method = "uniform")

        out <- bind_rows(dfit, ufit) |>
            mutate(ctmin_eps = eps[["ctmin_eps"]],
                   ctmax_eps = eps[["ctmax_eps"]],
                   logb_eps = eps[["logb_eps"]]) |>
            select(rep, method, everything())

        return(out)
    }

    if (.parallel) {
        fit_df <- future_lapply(1:.n_test_fits, one_temp_fun,
                               future.seed = .parallel_seed,
                               future.globals = c("par_lhs_list", "ace_args",
                                                  "ctmin", "ctmax", "a", "b",
                                                  "n_reps", "obs_cv", "n_temps",
                                                  "one_test_fit"),
                               future.packages = c("TPCdesign", "dplyr",
                                                   "nls.multstart", "tibble"))
    } else {
        fit_df <- map(1:.n_test_fits, one_temp_fun)
    }

    fit_df <- fit_df |>
        list_rbind() |>
        add_row(rep = NA, method  = "real",
                ctmin = .env$ctmin, ctmax = .env$ctmax,
                a = .env$a, b = .env$b, Topt = .env$Topt) |>
        mutate(n_temps = .env$n_temps,
               n_reps = .env$n_reps,
               obs_cv = .env$obs_cv) |>
        mutate(combo = j) |>
        select(combo, rep, method, n_temps:obs_cv, everything())

    .parent_prog()


    return(fit_df)

}





