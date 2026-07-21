



#' Design experimental temperatures
#'
#'
#'
#' @param n_temps Single integer for the number of temperatures to sample from.
#'     Must be >= 2.
#' @param z0 Single numeric for the degree to which sampled temperatures
#'     are derivative based versus uniform at the lowest possible
#'     sampled temperature. Values of `0` indicate uniform, and values of `1`
#'     indicate totally derivative based.
#'     Must be >= 0.
#' @param z1 Single numeric for the degree to which sampled temperatures
#'     are derivative based versus uniform at the highest possible
#'     sampled temperature. Values of `0` indicate uniform, and values of `1`
#'     indicate totally derivative based.
#'     Must be >= 0.
#' @param ctmin Single numeric for parameter `ctmin`.
#' @param ctmax Single numeric for parameter `ctmax`. Must be > `ctmin`.
#' @param a Single numeric for parameter `a`. Must be > 0.
#' @param b Single numeric for parameter `b`. Must be > 0.
#' @param temp_buffer Single numeric for the "buffer" around `ctmin` and `ctmax`
#'     to have for the possible sampling temperatures.
#'     More specifically, the minimum possible sampling temperature is
#'     `ctmin - temp_buffer`, and the maximum is `ctmax + temp_buffer`
#'     (both are rounded to have the same number of decimal places as
#'     `temp_by`).
#'     Must be >= 0.
#'     Must not have more decimal places than those in `temp_by`.
#'     Defaults to `5`.
#' @param temp_by Single numeric for the increment to generate possible
#'     sampling temperatures by.
#'     Specifically, possible temperatures will be
#'     `seq(ctmin - temp_buffer, ctmax + temp_buffer, temp_by)`.
#'     This number should reflect the precision of the equipment used to
#'     generate temperatures.
#'     Must be > 1e-6. Defaults to `0.1`.
#'
#' @returns A vector of temperatures at which to sample.
#'
#' @export
#'
deriv_design_temps <- function(n_temps, z0, z1, ctmin, ctmax, a, b,
                               temp_buffer = 5, temp_by = 0.1) {

    single_integer(n_temps, "n_temps", .min = 2L)
    single_number(z0, "z0", .min = 0)
    single_number(z1, "z1", .min = 0)
    single_number(ctmin, "ctmin")
    single_number(ctmax, "ctmax", .min = ctmin + .Machine$double.eps)
    single_number(a, "a", .min = .Machine$double.eps)
    single_number(b, "b", .min = .Machine$double.eps)
    single_number(temp_buffer, "temp_buffer", .min = 0)
    single_number(temp_by, "temp_by", .min = 1e-6)

    temp_digits <- decimalplaces(temp_by)
    if (decimalplaces(temp_buffer) > temp_digits)
        stop("`temp_buffer` must not have more decimal places than those in `temp_by`.")

    temp_min <- round(ctmin - temp_buffer, temp_digits)
    temp_max <- round(ctmax + temp_buffer, temp_digits)

    all_temps <- seq(temp_min, temp_max, temp_by) |>
        round(temp_digits)

    n_total_temps <- length(all_temps)

    abs_derivs <- abs(briere2_tpc_deriv(all_temps, ctmin, ctmax, a, b))
    deriv_wts <- abs_derivs / sum(abs_derivs)

    unif_wts <- rep(1 / n_total_temps, n_total_temps)

    z <- seq(z0, z1, length.out = n_total_temps)

    wts <- z * deriv_wts + (1 - z) * unif_wts
    wts <- wts / sum(wts)
    cs_wts <- cumsum(wts)

    temp_idx <- sapply(1:n_temps, \(x) which(cs_wts > x/(n_temps + 1))[1])

    return(all_temps[temp_idx])
}

