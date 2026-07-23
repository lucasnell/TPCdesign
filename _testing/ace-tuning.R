#'
#' Preamble for ACE tuning
#' Meant to be run on cluster
#'
#'

# RMSE for one fit for one combo of parameters
one_fit_rmse <- function(n_temps, n_reps, obs_cv, ctmin, ctmax, a, b,
                         ctmin_eps, ctmax_eps, logb_eps,
                         ...) {

    # Design parameters (can be inaccurate)
    dsn_ctmin <- ctmin + ctmin_eps
    dsn_ctmax <- ctmax + ctmax_eps
    dsn_b <- exp(log(b) + logb_eps)
    temps <- ace_design_temps(n_temps, dsn_ctmin, dsn_ctmax, a, dsn_b, ...)

    obs <- sim_gamma_data(temps, n_reps, obs_cv, ctmin, ctmax, a, b)

    starts_lo <- c(a = 0,  ctmin = 0,  ctmax = 30, b = 0.01)
    starts_up <- c(a = 2,  ctmin = 15, ctmax = 50, b = 3)

    # Rougly based on rTPC:::briere2_1999.lower[upper]_lims
    lims <- list(lo = c(ctmin = 0, ctmax = 5, a = 0, b = 0),
                 up = c(ctmin = 40, ctmax = 400, a = 10, b = 30))

    fit <- nls_multstart(
        formula = y ~ briere2_tpc(temp, ctmin, ctmax, a, b),
        data        = obs,
        start_lower = starts_lo,
        start_upper = starts_up,
        lower = lims$lo, upper = lims$up,
        supp_errors = "Y",
        control = list(maxfev = 5e3, maxiter = 1e3),
        lhstype = "improved",
        iter        = 500)

    if (is.null(fit)) {
        rmse <- NA_real_
    } else {
        fitted <- as.list(coef(fit))[c("ctmin", "ctmax", "a", "b")]
        test_temps <- seq(ctmin, ctmax, length.out = 101)
        obs_y <- briere2_tpc(test_temps, fitted[["ctmin"]], fitted[["ctmax"]],
                             fitted[["a"]], fitted[["b"]])
        tru_y <- briere2_tpc(test_temps, ctmin, ctmax, a, b)
        rmse <- sqrt(mean((obs_y - tru_y)^2))
    }

    return(rmse)
}





# Do one fits for a single combination of parameters
one_combo_fits <- function(j, input_df, prog) {

    args <- slice(input_df, j) |> select(-combo, -rep) |> as.list()

    rmse <- do.call(one_fit_rmse, args)

    if (!isTRUE(is.null(prog))) prog()

    return(rmse)

}





