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

    if (is.null(fit)) return(NA_real_)

    fitted <- as.list(coef(fit))[c("ctmin", "ctmax", "a", "b")]
    for (x in c("a", "b")) fitted[[x]] <- exp(fitted[[x]])
    test_temps <- seq(ctmin, ctmax, length.out = 101)
    obs_y <- briere2_tpc(test_temps, fitted[["ctmin"]], fitted[["ctmax"]],
                         fitted[["a"]], fitted[["b"]])
    tru_y <- briere2_tpc(test_temps, ctmin, ctmax, a, b)
    rmse <- sqrt(mean((obs_y - tru_y)^2))

    return(rmse)
}





# Do one fits for a single combination of parameters
one_combo_fits <- function(j, input_df, prog) {

    args <- slice(input_df, j) |> as.list()
    for (x in c("combo", "rep", "rmse")) args[[x]] <- NULL

    rmse <- do.call(one_fit_rmse, args)

    if (!isTRUE(is.null(prog))) prog()

    return(rmse)

}





