
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
                    temp_min,
                    temp_max,
                    optim_args = list(),
                    nls = FALSE,
                    nls_args = list(),
                    verbose = FALSE) {


    if (! inherits(obs, "data.frame")) stop("obs must be a data.frame")
    if (! "y" %in% colnames(obs)) stop("obs doesn't contain the column \"y\"")
    if (! "temp" %in% colnames(obs)) stop("obs doesn't contain the column \"temp\"")

    if (! is.list(optim_args)) stop("optim_args must be a list")
    if (! is.list(nls_args)) stop("nls_args must be a list")

    if (length(verbose) != 1 || !inherits(verbose, "logical"))
        stop("verbose must be a single logical")

    if (nls) {

        if (!requireNamespace("nls.multstart", quietly = TRUE)) {
            stop("Package \"nls.multstart\" must be installed to use nls.")
        }

        args <- list(formula = y ~ a * temp * (temp - CTmin) * (CTmax - temp)^b,
                     data = obs,
                     # start_lower = c(a = 0, CTmin = 0,  CTmax = 30, b = 0.01),
                     # start_upper = c(a = 2, CTmin = 15, CTmax = 50, b = 3),
                     iter        = 50,
                     supp_errors = "Y",
                     control = list(maxfev = 5e3, maxiter = 1e3),
                     lhstype = "improved")
        if (length(nls_args) > 0) {
            if (is.null(names(nls_args)) || any(names(nls_args) == "")) {
                stop("nls_args must only contain named elements")
            }
            for (n in names(nls_args)) args[[n]] <- nls_args[[n]]
        }
        fit <- do.call(nls.multstart::nls_multstart, args)

        if (is.null(fit)) {
            if (verbose) warning("nls_multistart returned NULL")
            # This will result in objective function returning very large
            # output if fit fails
            return(NULL)
        }

        if (!fit$convInfo$isConv) {
            if (verbose) warning("fit didn't converge")
            # This will result in objective function returning very large
            # output if fit fails
            return(NULL)
        }

        coefs <- coef(fit)
        names(coefs) <- tolower(names(coef))

        return(coef(fit)[c("ctmin", "ctmax", "a", "b")])

    }

    args <- list(fn = rmse_objective,
                 lower_bounds = c(temp_min, 0.75 * temp_min + 0.25 * temp_max,
                                  -10, -10),
                 upper_bounds = c(0.25 * temp_min + 0.75 * temp_max, temp_max,
                                  10, 10),
                 fn_args = list(y = obs$y, temp = obs$temp),
                 n_bevals = 100L,
                 n_boxes = 100L,
                 n_outputs = 1L,
                 controls = list(list(maxit = 1000, reltol = 1e-08)),
                 optimizers = c(stats::optim))
    if (length(optim_args) > 0) {
        if (is.null(names(optim_args)) || any(names(optim_args) == "")) {
            stop("optim_args must only contain named elements")
        }
        for (n in names(optim_args)) args[[n]] <- optim_args[[n]]
    }

    optim_fun <- args[["optimizers"]][[length(args[["optimizers"]])]]
    if (! inherits(optim_fun, "function")) {
        stop("optim_args$optimizers, if provided, must be a function")
    }
    optim_pkg <- packageName(environment(optim_fun))
    if (! optim_pkg %in% c("stats", "nloptr")) {
        stop("optim_args$optimizers, if provided, can only be from stats or nloptr packages")
    }
    if (optim_pkg == "stats" && !identical(optim_fun, stats::optim)) {
        stop("optim is the only stats function that optim_args$optimizers can be")
    }

    op <- do.call(estimatePMR::winnowing_optim, args)[[1]]

    if (optim_pkg == "nloptr") {
        # This will result in objective function returning very large
        # output if fit fails
        if (op$convergence < 0) return(NULL)
    } else {
        # This will result in objective function returning very large
        # output if fit fails
        if (op$convergence != 0) return(NULL)
    }


    coefs <- op$par
    coefs[3:4] <- exp(coefs[3:4])
    names(coefs) <- c("ctmin", "ctmax", "a", "b")
    return(coefs)

}



# =============================================================================*
# =============================================================================*
# Fit temperatures ----
# =============================================================================*
# =============================================================================*




#' Objective function for parameters related to sampling temperatures
#'
#' @noRd
#'
temp_objective <- function(params,
                           n_reps,
                           obs_cv,
                           temp_min,
                           temp_max,
                           CTmin,
                           CTmax,
                           a,
                           b,
                           output = "RMSE",
                           optim_args = list(),
                           nls = FALSE,
                           nls_args = list(),
                           n_grid_temps = 101L,
                           verbose = FALSE) {

    is_type(params, "params", "numeric", len_min = 2L)
    single_integer(n_reps, "n_reps", .min = 1L)
    single_number(obs_cv, "obs_cv", .min = .Machine$double.eps)
    single_number(temp_min, "temp_min")
    single_number(temp_max, "temp_max", .min = temp_min)
    if (temp_max == temp_min) stop("temp_min cannot equal temp_max")
    single_number(CTmin, "CTmin", .min = temp_min)
    single_number(CTmax, "CTmax", .min = CTmin, .max = temp_max)
    if (CTmax == CTmin) stop("CTmin cannot equal CTmax")
    single_number(a, "a", .min = .Machine$double.eps)
    single_number(b, "b", .min = .Machine$double.eps)
    single_string(output, "output")
    is_type(optim_args, "optim_args", "list")
    single_logical(nls, "nls")
    is_type(nls_args, "nls_args", "list")
    single_logical(verbose, "verbose")


    output <- match.arg(tolower(output), c("rmse", "ctmin", "ctmax", "b", "a"))

    design_temps <- make_temps(params, temp_min, temp_max)

    obs <- sim_gamma_data(design_temps, n_reps, obs_cv, CTmin, CTmax, a, b)

    tpc_pars <- fit_tpc(obs, temp_min, temp_max, optim_args, nls, nls_args,
                        verbose)
    # This only happens when coefs couldn't be fit
    if (is.null(tpc_pars)) return(1e10)


    if (output == "rmse") {

        temp_grid <- seq(temp_min, temp_max, length.out = n_grid_temps)

        true_grid <- briere2_tpc(temp = temp_grid, a = a, CTmin = CTmin,
                                 CTmax = CTmax, b = b)

        pred_grid <- briere2_tpc(temp = temp_grid,
                                 a = tpc_pars[["a"]],
                                 CTmin = tpc_pars[["ctmin"]],
                                 CTmax = tpc_pars[["ctmax"]],
                                 b = tpc_pars[["b"]],
                                 scale = FALSE)

        out <- sqrt(mean((true_grid - pred_grid)^2))
        if (is.infinite(out)) out <- 1e10

    } else {

        out <- tpc_pars[[output]]

    }


    return(out)

}





