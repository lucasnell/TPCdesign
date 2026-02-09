
#' Fit parameters for an individual TPC
#'
#' @noRd
#'
fit_tpc <- function(obs,
                    params0,
                    optim_args = list(),
                    optim_fun = stats::optim,
                    nls = FALSE,
                    nls_args = list(),
                    verbose = FALSE) {


    if (! inherits(obs, "data.frame")) stop("obs must be a data.frame")
    if (! "y" %in% colnames(obs)) stop("obs doesn't contain the column \"y\"")
    if (! "temp" %in% colnames(obs)) stop("obs doesn't contain the column \"temp\"")

    if (!is.numeric(params0) || !is.vector(params0)) stop("params0 must be numeric vector")

    if (! inherits(optim_fun, "function")) stop("optim_fun must be a function")
    optim_pkg <- packageName(environment(optim_fun))
    if (! optim_pkg %in% c("stats", "nloptr")) {
        stop("optim_fun can only be from stats or nloptr packages")
    }
    if (optim_pkg == "stats") {
        if (gsub(".*::", "", deparse(substitute(optim_fun))) != "optim") {
            stop("optim is the only stats function that optim_fun can be")
        }
    }
    if (! is.list(optim_args)) stop("optim_args must be a list")
    if (! is.list(nls_args)) stop("nls_args must be a list")

    if (length(verbose) != 1 || !inherits(verbose, "logical"))
        stop("verbose must be a single logical")

    if (nls) {

        if (!requireNamespace("nls.multstart", quietly = TRUE)) {
            stop("Package \"nls.multstart\" must be installed to use nls.")
        }

        args <- list(formula = y ~ a * temp * (temp - CTmin) *
                         pmax(CTmax - temp, 0)^b,
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

        return(coef(fit)[c("CTmin", "CTmax", "a", "b")])

    }

    args <- list(par = params0,
                 fn = rmse_objective,
                 y = obs$y, temp = obs$temp)
    if (length(optim_args) > 0) {
        if (is.null(names(optim_args)) || any(names(optim_args) == "")) {
            stop("optim_args must only contain named elements")
        }
        for (n in names(optim_args)) args[[n]] <- optim_args[[n]]
    }
    # Switch names for nloptr package:
    if (optim_pkg == "nloptr") {
        args[["x0"]] <- params0
        args[["par"]] <- NULL
    }
    op <- do.call(optim_fun, args)

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
    coefs[1:3] <- exp(coefs[1:3])
    names(coefs) <- c("CTmin", "CTmax", "a", "b")
    return(coefs)

}
