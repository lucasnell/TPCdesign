




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
    if (max_temp %in% opt_temps) final <- c(max_temp, upper)

    final <- sort(round(final, digits))
    return(final)
}




#' Design experimental temperatures using ACE method
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
#'     maximum `ctmin + ctmin_err`. Defaults to `2.50`.
#' @param ctmax_err Single numeric for error for `ctmax`.
#'     This means priors for `ctmax` will be generated from a uniform
#'     distribution with minimum `ctmax - ctmax_err` and
#'     maximum `ctmax + ctmax_err`. Defaults to `1.50`.
#' @param logb_err Single numeric for error for `log(b)`.
#'     This means priors for `log(b)` will be generated from a uniform
#'     distribution with minimum `log(b) - logb_err` and
#'     maximum `log(b) + logb_err`. Defaults to `0.26`.
#' @param min_sep Single numeric specifying the minimum separation between
#'     optimized temperatures. The step where equally-spaced temperatures are
#'     added (only happens when `n_filler >= 1`) is not affected by this
#'     argument, and unless the number of points is very high relative
#'     to the difference in maximum (`ctmax + ctmax_err`) and
#'     minimum (`ctmin + ctmin_err`) temperatures surveyed, having a minimum
#'     separation shouldn't change this step anyway. Defaults to `0.5`.
#' @param n_filler Single integer specifying the number of temperatures
#'     that are equally spaced versus optimized.
#'     Equally spacing points can be a good bet hedging strategy if you're
#'     quite unsure of the TPC parameters (`ctmin`, `ctmax`, `a`, `b`)
#'     you're using. Defaults to `0`.
#' @param n_draws Single integer specifying the number of prior draws to
#'     use. Increase if your certainty is low (i.e., `*_err` parameters are
#'     high). Defaults to `100L`.
#' @param n_starts Single integer specifying the number of starts
#'     Defaults to `5L`.
#' @param digits Single integer specifying the digits to round temperatures to.
#'     Defaults to `2L`.
#' @param n_threads Single integer for the number of threads to use.
#'     Multithreading is done across starts, so there is no need for `n_threads`
#'     to exceed `n_starts`. Defaults to `1L`.
#'
#'
#'
#' @returns A vector of temperatures at which to sample.
#'
#' @export
#'
ace_design_temps <- function(n_temps, ctmin, ctmax, a, b,
                             ctmin_err = 2.50,
                             ctmax_err = 1.50,
                             logb_err = 0.26,
                             min_sep = 0.5,
                             n_filler = 0,
                             n_draws = 100L,
                             n_starts = 5L,
                             digits = 2L,
                             n_threads = 1L) {

    # Check args:
    single_integer(n_temps, "n_temps", .min = 2L)
    single_number(n_temps, "n_temps", .min = 2L)
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
    single_integer(n_threads, "n_threads", .min = 1L)


    n_optimal <- n_temps - n_filler

    # genus-level relatives was 2.50°C for Tmin, 1.50°C for Tmax,
    # and 0.26 for log(b)
    prior_ctmin <- ctmin + c(-1,1) * ctmin_err
    prior_ctmax <- ctmax + c(-1,1) * ctmax_err
    prior_lb <- log(b) + c(-1,1) * logb_err


    theta_draws <- cbind(ctmin = runif(n_draws, prior_ctmin[1], prior_ctmin[2]),
                         ctmax = runif(n_draws, prior_ctmax[1], prior_ctmax[2]),
                         b = exp(runif(n_draws, prior_lb[1], prior_lb[2])),
                         a = rep(a, n_draws))


    min_temp <- min(prior_ctmin)
    max_temp <- max(prior_ctmax)

    if ((n_optimal - 1L) * min_sep > (max_temp - min_temp)) {
        stop("min_sep is too large for given min and max temperatures and the ",
             "number of desired optimized temps")
    }

    start_temp_list <- lapply(1:n_starts, function(i) {
        d <- lhs::randomLHS(n = n_optimal, k = 1) * (max_temp - min_temp) + min_temp
        colnames(d) <- "temp"
        return(d)
    })

    # Limit temps to have a minimum separation:
    limits_minsep <- function(d, i, j) {
        grid <- seq(from = min_temp, to = max_temp, length.out = 2000)
        other_points <- as.vector(d)[-i]
        for (s in other_points) {
            grid <- grid[(grid < (s - min_sep)) | (grid > (s + min_sep))]
        }
        return(grid)
    }

    design_robust <- acebayes::pace(utility = utility_briere2D,
                                    start.d = start_temp_list,
                                    B = rep(list(theta = theta_draws), 2),
                                    lower = min_temp,
                                    upper = max_temp,
                                    limits = limits_minsep,
                                    N2 = 0L,
                                    deterministic = TRUE,
                                    mc.cores = n_threads)

    opt_temps <- sort(round(design_robust$d, digits))

    final_temps <- gap_filler(n_filler, opt_temps, min_temp, max_temp, digits)

    return(final_temps)

}

