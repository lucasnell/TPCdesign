


# Version of acebayes::pace that allows multithreading using future.apply
# Some options that we don't use are removed.
tpc_pace <- function (utility, start.d, B, lower, upper, limits) {

    # non-default args to acebayes::pace that won't change
    N2 = 0
    deterministic = TRUE

    # args that are defaults
    Q = 20
    N1 = 20
    binary = FALSE
    n.assess = 20

    ptm <- proc.time()[3]
    C <- length(start.d)
    # if (missing(B) && identical(deterministic, FALSE)) {
    #     B <- c(20000, 1000)
    # }
    # if (missing(B) && identical(deterministic, TRUE)) {
    #     B <- NULL
    # }
    innerace <- function(d) {
        out <- acebayes::ace(utility = utility, start.d = d, B = B, Q = Q,
                             N1 = N1, N2 = N2, lower = lower, upper = upper, limits = limits,
                             binary = binary, deterministic = deterministic, progress = FALSE)
        list(phase2.d = out$phase2.d, phase1.trace = out$phase1.trace,
             phase2.trace = out$phase2.trace)
    }

    do_parallel <- requireNamespace("future.apply", quietly = TRUE) &&
        requireNamespace("future", quietly = TRUE) &&
        !methods::is(future::plan(), "sequential")

    if (do_parallel) {
        globals <- c("utility", "B", "Q", "N1", "N2", "lower", "upper", "limits",
                     "binary", "deterministic")
        rout <- future.apply::future_lapply(start.d, innerace,
                                            future.seed = TRUE,
                                            future.globals = globals,
                                            future.packages = c("TPCdesign", "acebayes"))
    } else {
        rout <- lapply(start.d, innerace)
    }

    # if (.Platform$OS.type == "unix") {
    #     rout <- mclapply(start.d, innerace, mc.cores = mc.cores)
    # }
    # else {
    #     if (mc.cores == 1) {
    #         rout <- lapply(start.d, innerace)
    #     }
    #     if (mc.cores > 1) {
    #         warning("mc.cores > 1 not currently supported under a non-Unix OS. Proceeding with mc.cores = 1 \n")
    #         rout <- lapply(start.d, innerace)
    #     }
    # }
    routd <- list()
    routphase1 <- list()
    routphase2 <- list()
    for (i in 1:C) {
        routd[[i]] <- rout[[i]]$phase2.d
        routphase1[[i]] <- rout[[i]]$phase1.trace
        routphase2[[i]] <- rout[[i]]$phase2.trace
    }
    # if (!deterministic) {
    #     inner <- function(d) {
    #         evals <- rep(0, n.assess)
    #         for (i in 1:n.assess) {
    #             evals[i] <- mean(utility(d = d, B = B[1]))
    #         }
    #         evals
    #     }
    #     if (.Platform$OS.type == "unix") {
    #         fout <- mclapply(routd, inner, mc.cores = mc.cores)
    #     }
    #     else {
    #         if (mc.cores == 1) {
    #             fout <- lapply(routd, inner)
    #         }
    #         if (mc.cores > 1) {
    #             fout <- lapply(routd, inner)
    #         }
    #     }
    #     mout <- lapply(fout, mean)
    #     besti <- which.max(mout)
    # }
    # if (deterministic) {
    inner <- function(d) {
        utility(d = d, B = B)
    }
    if (do_parallel) {
        globals <- c("utility", "B")
        fout <- future.apply::future_lapply(routd, inner,
                                            future.seed = TRUE,
                                            future.globals = globals,
                                            future.packages = c("TPCdesign", "acebayes"))
    } else {
        fout <- lapply(routd, inner)
    }
    # if (.Platform$OS.type == "unix") {
    #     fout <- mclapply(routd, inner, mc.cores = mc.cores)
    # }
    # else {
    #     if (mc.cores == 1) {
    #         fout <- lapply(routd, inner)
    #     }
    #     if (mc.cores > 1) {
    #         fout <- lapply(routd, inner)
    #     }
    # }
    besti <- which.max(fout)
    # }
    ptm <- proc.time()[3] - ptm
    phase1.trace <- NULL
    phase2.trace <- NULL
    if (N1 > 0) {
        phase1.trace <- routphase1[[besti]]
    }
    if (N2 > 0) {
        phase2.trace <- routphase2[[besti]]
    }
    output <- list(d = routd[[besti]], phase1.trace = phase1.trace,
                   phase2.trace = phase2.trace, eval = fout[[besti]], utility = utility,
                   start.d = start.d, final.d = routd, besti = besti, B = B,
                   Q = Q, N1 = N1, N2 = N2, glm = FALSE, nlm = FALSE, criterion = "NA",
                   prior = "NA", time = ptm, binary = binary, deterministic = deterministic)
    class(output) <- "pace"
    output
}



# Fills in gaps to hedge against curve mismatch
gap_filler <- function(n_filler, opt_temps, min_temp, max_temp, digits) {

    if (n_filler <= 0) return(opt_temps)

    points <- sort(opt_temps)

    # Include the domain boundaries as candidate gap edges
    points <- c(min_temp, points, max_temp)

    for (i in 1:n_filler) {
        gaps <- diff(points)
        largest_gap_idx <- which.max(gaps)
        new_point <- mean(points[largest_gap_idx:(largest_gap_idx + 1)])
        points <- sort(c(points, new_point))
    }

    # Drop the boundary markers if they weren't part of the original optimal set
    final <- setdiff(points, c(min_temp, max_temp))
    # unless lower/upper were themselves chosen as optimal points, add them back
    if (min_temp %in% opt_temps) final <- c(min_temp, final)
    if (max_temp %in% opt_temps) final <- c(final, max_temp)

    final <- sort(round(final, digits))
    return(final)
}




#' Design experimental temperatures
#'
#'
#'
#' @param n_temps Single integer for the number of temperatures to sample from.
#'     Must be >= 2.
#' @param ctmin Single numeric for parameter `ctmin`.
#' @param ctmax Single numeric for parameter `ctmax`. Must be > `ctmin`.
#' @param a Single numeric for parameter `a`. Must be > 0.
#' @param b Single numeric for parameter `b`. Must be > 0.
#' @param ctmin_err Single numeric for error for `ctmin`.
#'     This means priors for `ctmin` will be generated from a uniform
#'     distribution with minimum `ctmin - ctmin_err` and
#'     maximum `ctmin + ctmin_err`. Must be `>= 0`.
#'     Defaults to `2.50`.
#' @param ctmax_err Single numeric for error for `ctmax`.
#'     This means priors for `ctmax` will be generated from a uniform
#'     distribution with minimum `ctmax - ctmax_err` and
#'     maximum `ctmax + ctmax_err`. Must be `>= 0`.
#'     Defaults to `1.50`.
#' @param logb_err Single numeric for error for `log(b)`.
#'     This means priors for `log(b)` will be generated from a uniform
#'     distribution with minimum `log(b) - logb_err` and
#'     maximum `log(b) + logb_err`. Must be `>= 0`.
#'     Defaults to `0.26`.
#' @param min_sep Single numeric specifying the minimum separation between
#'     optimized temperatures. The step where equally-spaced temperatures are
#'     added (only happens when `n_filler >= 1`) is not affected by this
#'     argument, and unless the number of points is very high relative
#'     to the difference in maximum (`ctmax + ctmax_err`) and
#'     minimum (`ctmin + ctmin_err`) temperatures surveyed, having a minimum
#'     separation shouldn't change this step anyway. Defaults to `1`.
#' @param n_filler Single integer specifying the number of temperatures
#'     that are equally spaced versus optimized.
#'     Equally spacing points can be a good bet hedging strategy if you're
#'     quite unsure of the TPC parameters (`ctmin`, `ctmax`, `a`, `b`)
#'     you're using. Defaults to `1L`.
#' @param n_draws Single integer specifying the number of prior draws to
#'     use. Increase if your certainty is low (i.e., `*_err` parameters are
#'     high). Defaults to `250L`.
#' @param n_starts Single integer specifying the number of starts
#'     Defaults to `7L`.
#' @param digits Single integer specifying the digits to round temperatures to.
#'     Defaults to `2L`.
#'
#'
#'
#' @returns A vector of temperatures at which to sample.
#'
#' @examples
#' design_temps(n_temps = 10, ctmin = 5, ctmax = 40, a = 1, b = 0.2)
#'
#' # To run in parallel:
#' # library(future.apply)
#' # plan(multisession, workers = 4L)
#' # design_temps(n_temps = 10, ctmin = 5, ctmax = 40, a = 1, b = 0.2)
#'
#' @export
#'
design_temps <- function(n_temps, ctmin, ctmax, a, b,
                         ctmin_err = 2.50,
                         ctmax_err = 1.50,
                         logb_err = 0.26,
                         min_sep = 1,
                         n_filler = 1L,
                         n_draws = 250L,
                         n_starts = 7L,
                         digits = 2L) {

    # Check args:
    single_integer(n_temps, "n_temps", .min = 2L)
    single_number(ctmin, "ctmin")
    single_number(ctmax, "ctmax", .min = ctmin + .Machine$double.eps)
    single_number(a, "a", .min = .Machine$double.eps)
    single_number(b, "b", .min = .Machine$double.eps)
    single_number(ctmin_err, "ctmin_err", .min = 0)
    single_number(ctmax_err, "ctmax_err", .min = 0)
    single_number(logb_err, "logb_err", .min = 0)
    single_number(min_sep, "min_sep", .min = 0)
    single_integer(n_filler, "n_filler", .min = 0L, .max = n_temps)
    single_integer(n_draws, "n_draws", .min = 1L)
    single_integer(n_starts, "n_starts", .min = 1L)
    single_integer(digits, "digits", .min = 0L)

    n_optimal <- n_temps - n_filler

    prior_ctmin <- ctmin + c(-1,1) * ctmin_err
    prior_ctmax <- ctmax + c(-1,1) * ctmax_err
    prior_lb <- log(b) + c(-1,1) * logb_err

    theta_draws <- cbind(ctmin = stats::runif(n_draws, prior_ctmin[1], prior_ctmin[2]),
                         ctmax = stats::runif(n_draws, prior_ctmax[1], prior_ctmax[2]),
                         b = exp(stats::runif(n_draws, prior_lb[1], prior_lb[2])),
                         a = rep(a, n_draws))


    if ((n_optimal - 1L) * min_sep > (ctmax - ctmin)) {
        stop("min_sep is too large for given ctmin and ctmax temperatures and the ",
             "number of desired optimized temps")
    }

    start_temp_list <- lapply(1:n_starts, function(i) {
        # We want to be more conservative with starting values (compared to
        # min_temp and max_temp below)
        start_min <- ctmin + min_sep
        start_max <- ctmax - min_sep
        d <- lhs::randomLHS(n = n_optimal, k = 1) * (start_max - start_min) +
            start_min
        d <- d[order(d),,drop=FALSE]
        colnames(d) <- "temp"
        return(d)
    })

    min_temp <- min(prior_ctmin)
    max_temp <- max(prior_ctmax)

    # Limit temps to have a minimum separation:
    limits_minsep <- function(d, i, j) {
        grid <- seq(from = min_temp, to = max_temp, length.out = 2000)
        other_points <- as.vector(d)[-i]
        for (s in other_points) {
            grid <- grid[(grid < (s - min_sep)) | (grid > (s + min_sep))]
        }
        return(grid)
    }

    design_robust <- tpc_pace(utility = utility_briere2,
                                    start.d = start_temp_list,
                                    B = rep(list(theta = theta_draws), 2),
                                    lower = min_temp,
                                    upper = max_temp,
                                    limits = limits_minsep)

    opt_temps <- sort(round(design_robust$d, digits))

    final_temps <- gap_filler(n_filler, opt_temps, min_temp, max_temp, digits)

    return(final_temps)

}

