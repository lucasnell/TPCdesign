
# =============================================================================*
# =============================================================================*
# Fit one TPC ----
# =============================================================================*
# =============================================================================*


#' Fit parameters for an individual TPC
#'
#' @noRd
#'
fit_tpc <- function(obs,
                    lower_bounds,
                    upper_bounds,
                    iter,
                    optim_pkg,
                    optim_fun,
                    optim_args,
                    lhs_mapper,
                    scale_tpc,
                    verbose) {

    lhs_mat <- lhs_mapper(iter, 4L)
    for (i in 1:4) {
        lhs_mat[,i] <- lower_bounds[i] + (upper_bounds[i] - lower_bounds[i]) *
            lhs_mat[,i]
    }

    other_args <- add_args(list(fn = rmse_objective, y = obs$y, temp = obs$temp,
                                scale = scale_tpc),
                           optim_args, naughty_pars = c("par", "x0"))

    # Switch names for nloptr package:
    name_par0 <- ifelse(optim_pkg == "nloptr", "x0", "par")
    get_converged <- function(op) {
        if (optim_pkg == "nloptr") return(op$convergence > 0)
        return(op$convergence == 0)
    }

    fits <- lapply(1:nrow(lhs_mat),
                   \(i) {
                       args <- other_args
                       args[[name_par0]] <- lhs_mat[i,]
                       op <- do.call(optim_fun, args)
                       return(op)
                   })
    values <- sapply(fits, \(x) x$value)
    convergences <- sapply(fits, get_converged)
    best_fit_idx <- which(values == min(values[convergences], na.rm = TRUE) &
                              convergences)
    if (length(best_fit_idx) == 0 && sum(convergences) == 0) {
        if (verbose) warning("no fits converged")
        best_fit_idx <- which(values == min(values, na.rm = TRUE))
        if (length(best_fit_idx) == 0) stop("no non-NA fits")
    }
    best_fits <- fits[best_fit_idx]

    best_pars <- lapply(best_fits, \(x) x$par) |>
        do.call(what = rbind) |>
        colMeans(na.rm = TRUE)
    best_pars[3:4] <- exp(best_pars[3:4])
    names(best_pars) <- c("ctmin", "ctmax", "a", "b")

    return(best_pars)

}



# =============================================================================*
# =============================================================================*
# Fit temperatures ----
# =============================================================================*
# =============================================================================*




#' Objective function for parameters related to sampling temperatures
#'
#'
#'
#' @noRd
#'
temp_objective <- function(params,
                           n_reps,
                           obs_cv,
                           temp_min,
                           temp_max,
                           ctmin,
                           ctmax,
                           a,
                           b,
                           output,
                           iter,
                           optim_pkg,
                           optim_fun,
                           optim_args,
                           lhs_mapper,
                           n_grid_temps,
                           scale_tpc,
                           verbose) {


    design_temps <- make_temps(params, temp_min, temp_max)

    obs <- sim_gamma_data(design_temps, n_reps, obs_cv, ctmin, ctmax, a, b,
                          scale_tpc)

    lb  <- c(temp_min, 0.75 * temp_min + 0.25 * temp_max, -10, -10)
    ub <- c(0.25 * temp_min + 0.75 * temp_max, temp_max, 10, 10)

    tpc_pars <- fit_tpc(obs, lower_bounds = lb, upper_bounds = ub,
                        iter = iter,
                        optim_pkg = optim_pkg,
                        optim_fun = optim_fun,
                        optim_args = optim_args,
                        lhs_mapper = lhs_mapper,
                        verbose = verbose,
                        scale_tpc = scale_tpc)
    # This only happens when coefs couldn't be fit
    if (is.null(tpc_pars)) return(1e10)


    if (output == "rmse") {

        temp_grid <- seq(temp_min, temp_max, length.out = n_grid_temps)

        true_grid <- briere2_tpc(temp = temp_grid, a = a, ctmin = ctmin,
                                 ctmax = ctmax, b = b,
                                 scale = scale_tpc)

        pred_grid <- briere2_tpc(temp = temp_grid,
                                 a = tpc_pars[["a"]],
                                 ctmin = tpc_pars[["ctmin"]],
                                 ctmax = tpc_pars[["ctmax"]],
                                 b = tpc_pars[["b"]],
                                 scale = scale_tpc)

        out <- sqrt(mean((true_grid - pred_grid)^2))
        if (is.infinite(out) || is.na(out)) out <- 1e10

    } else {

        out <- abs(eval(parse(text = output)) - tpc_pars[[output]])

    }


    return(out)

}






#' Design experimental temperatures
#'
#' @param n_temps Single integer for the number of temperatures to sample from.
#' @param n_reps Single integer for the number of experimental reps per temperature.
#' @param obs_cv Single numeric for the observation error coefficient of variation.
#'     Must be >= 0.
#'     Observation error is simulated using a gamma distribution, except for
#'     when performance values are below zero, in which case it uses
#'     a normal approximation.
#' @param temp_min Single numeric for the minimum possible temperature allowed.
#' @param temp_max Single numeric for the maximum possible temperature allowed.
#' @param ctmin Single numeric for parameter `ctmin`.
#' @param ctmax Single numeric for parameter `ctmax`.
#' @param a Single numeric for parameter `a`.
#' @param b Single numeric for parameter `b`.
#' @param tpc_measure Single string for the measure to use to fit temperatures.
#'     Options are `"rmse"`, `"ctmin"`, `"ctmax"`, `"b"`, or `"a"`, and
#'     capitalization is ignored.
#'     Option `"rmse"` returns the root mean square deviation between the
#'     between "observed" (simulated) and expected performance values for a
#'     given set of parameters.
#'     The other options are absolute differences between true and fitted
#'     values for individual TPC parameter values.
#'     Defaults to `"rmse"`.
#' @param tpc_optim_fun A function to use for optimization in TPC fitting.
#'     Options are `stats::optim` or functions from the `nloptr` package.
#'     Defaults to `stats::optim`.
#' @param tpc_optim_args A named list of additional arguments to pass
#'     to `tpc_optim_fun`.
#'     Defaults `list()`.
#' @param tpc_iter Single integer specifying the
#'     "number of partitions (simulations or design points or rows)"
#'     to use in Latin Hypercube Sampling for TPC fitting.
#' @param tpc_lhstype Single string specifying the type of Latin Hypercube
#'     Sampling to use (via the `lhs` package) for TPC fitting. Options are
#'     `"random"` (uses [lhs::randomLHS])
#'     `"improved"` (uses [lhs::improvedLHS])
#'     `"maximin"` (uses [lhs::maximinLHS])
#'     `"genetic"` (uses [lhs::geneticLHS]).
#'     Defaults to `"improved"`.
#' @param tpc_n_grid_temps Single interger indicating the number of
#'     temperatures to use to compute RMSE in TPC fitting. Must be `>= 10`.
#'     Defaults to `101L`.
#' @param scale_tpc Single logical for whether to scale performance values
#'     to make them have a max value of 1.
#'     Defaults to `TRUE`.
#' @param lower_bounds Numeric vectors of length `n_temps+1L` indicating
#'     the lower bounds in starting conditions used in
#'     `estimatePMR::winnowing_optim`.
#'     Values must be > 0 and < 1.
#'     If not provided, it defaults to `rep(1e-9, n_temps+1L)`.
#' @param upper_bounds Numeric vectors of length `n_temps+1L` indicating
#'     the upper bounds in starting conditions used in
#'     `estimatePMR::winnowing_optim`.
#'     Values must be > 0 and < 1.
#'     If not provided, it defaults to `rep(1 - 1e-9, n_temps+1L)`.
#' @param verbose Single logical for whether to warn of poor fits.
#'     Useful for debugging, but results in potentially many warnings.
#'     Defaults to `FALSE`.
#' @param ... Other arguments to pass to [estimatePMR::winnowing_optim()].
#'     Arguments `fn` and `fn_args` are not allowed to be passed here because
#'     they are handled by this function explicitly.
#'     If nothing is provided here, `winnowing_optim` is run with
#'     `n_bevals = 100L`,
#'     `n_boxes = 10L`,
#'     `n_outputs = 1L`,
#'     `controls = list(list(maxit = 1000, reltol = 1e-08))`,
#'     and
#'     `optimizers = c(nloptr::newuoa)`.
#'
#' @returns Typically the output will be a single object of the type returned
#'     by `optim_fun`. The exception is if the `n_outputs` is provided
#'     (and passed to `winnowing_optim()`), and the last item in `n_outputs`
#'     is greater than 1.
#'     In that case, it will be a list of items returned from `optim_fun`,
#'     and the list's length will be `tail(n_outputs, 1L)`.
#'     Regardless of length, each item in the output has an extra field added
#'     to it (`"temps"`) indicating the temperatures from this fit.
#'
#' @export
#'
design_temps <- function(n_temps,
                         n_reps,
                         obs_cv,
                         temp_min,
                         temp_max,
                         ctmin,
                         ctmax,
                         a,
                         b,
                         tpc_measure = "rmse",
                         tpc_optim_fun = stats::optim,
                         tpc_optim_args = list(),
                         tpc_iter = 10L,
                         tpc_lhstype = "improved",
                         tpc_n_grid_temps = 101L,
                         scale_tpc = TRUE,
                         lower_bounds = NULL,
                         upper_bounds = NULL,
                         verbose = FALSE,
                         ...) {

    single_integer(n_temps, "n_temps", .min = 2L)
    single_integer(n_reps, "n_reps", .min = 1L)
    single_number(obs_cv, "obs_cv", .min = .Machine$double.eps)
    single_number(temp_min, "temp_min")
    single_number(temp_max, "temp_max", .min = temp_min)
    if (temp_max == temp_min) stop("temp_min cannot equal temp_max")
    single_number(ctmin, "ctmin", .min = temp_min)
    single_number(ctmax, "ctmax", .min = ctmin, .max = temp_max)
    if (ctmax == ctmin) stop("ctmin cannot equal ctmax")
    single_number(a, "a", .min = .Machine$double.eps)
    single_number(b, "b", .min = .Machine$double.eps)
    single_string(tpc_measure, "tpc_measure")
    is_type(tpc_optim_fun, "tpc_optim_fun", "function")
    is_type(tpc_optim_args, "tpc_optim_args", "list")
    single_integer(tpc_iter, "tpc_iter", .min = 1L)
    single_string(tpc_lhstype, "tpc_lhstype")
    single_integer(tpc_n_grid_temps, "tpc_n_grid_temps", .min = 10L)
    single_logical(scale_tpc, "scale_tpc")
    single_logical(verbose, "verbose")

    if (is.null(lower_bounds)) {
        lower_bounds <- rep(logit(1e-9), n_temps+1L)
    } else {
        is_type(lower_bounds, "lower_bounds", is.numeric,
                len_min = n_temps+1, len_max = n_temps + 1,
                .min = .Machine$double.eps, .max = 1 - .Machine$double.eps)
    }
    if (is.null(upper_bounds)) {
        upper_bounds <- rep(logit(1-1e-9), n_temps+1L)
    } else {
        is_type(upper_bounds, "upper_bounds", is.numeric,
                len_min = n_temps+1, len_max = n_temps + 1,
                .min = .Machine$double.eps, .max = 1 - .Machine$double.eps)
    }

    tpc_measure <- match.arg(tolower(tpc_measure),
                             c("rmse", "ctmin", "ctmax", "b", "a"))
    tpc_lhstype <- match.arg(tolower(tpc_lhstype),
                             c("random", "improved", "maximin", "genetic"))

    # From R package nls.multistarts
    # permalink: https://github.com/padpadpadpad/nls.multstart/blob/4842e857d5c5cd182c82c30db6c702db4a4d5da2/R/nls_multstart.R#L204-L209
    lhs_mapper <- switch(tpc_lhstype,
                         random = lhs::randomLHS,
                         improved = lhs::improvedLHS,
                         maximin = lhs::maximinLHS,
                         genetic = lhs::geneticLHS)
    if (is.null(lhs_mapper)) {
        err_msg("tpc_lhstype", "\"random\", \"improved\", \"maximin\", or",
                "\"genetic\"")
    }

    tpc_optim_pkg <- packageName(environment(tpc_optim_fun))
    if (! tpc_optim_pkg %in% c("stats", "nloptr")) {
        err_msg("tpc_optim_fun", "from stats or nloptr packages")
    }
    if (tpc_optim_pkg == "stats" && !identical(tpc_optim_fun, stats::optim)) {
        err_msg("tpc_optim_fun", "`optim` if it's from the stats package")
    }


    args <- list(fn = temp_objective,
                 lower_bounds = lower_bounds,
                 upper_bounds = upper_bounds,
                 fn_args = list(n_reps = n_reps,
                                obs_cv = obs_cv,
                                temp_min = temp_min,
                                temp_max = temp_max,
                                ctmin = ctmin,
                                ctmax = ctmax,
                                a = a,
                                b = b,
                                output = tpc_measure,
                                iter = tpc_iter,
                                optim_pkg = tpc_optim_pkg,
                                optim_fun = tpc_optim_fun,
                                optim_args = tpc_optim_args,
                                lhs_mapper = lhs_mapper,
                                n_grid_temps = tpc_n_grid_temps,
                                scale_tpc = scale_tpc,
                                verbose = verbose),
                 n_bevals = 100L,
                 n_boxes = 10L,
                 n_outputs = 1L,
                 controls = list(list(maxit = 1000, reltol = 1e-08)),
                 optimizers = c(nloptr::newuoa))

    args <- add_args(args, list(...), c("fn", "fn_args"))

    op <- do.call(estimatePMR::winnowing_optim, args)

    for (i in 1:length(op)) {
        op[[i]][["temps"]] <- make_temps(op[[i]][["par"]], temp_min, temp_max)
    }

    if (length(op) == 1L) op <- op[[1]]

    return(op)

}
