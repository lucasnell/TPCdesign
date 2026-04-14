

#' Thermal optimum for Brière-2 thermal performance curve (TPC)
#'
#' Note that parameter `a` is not necessary for this calculation.
#'
#'
#' @param ctmin Numeric vector for parameter `ctmin`.
#' @param ctmax Numeric vector for parameter `ctmax`.
#' @param b Numeric vector for parameter `b`.
#'
#' @returns A numeric vector of optimum temperatures.
#'
#' @export
#'
briere2_tpc_Topt <- function(ctmin, ctmax, b) {

    is_type(ctmin, "ctmin", is.numeric, len_min = 1)
    is_type(ctmax, "ctmax", is.numeric, len_min = 1)
    is_type(b, "b", is.numeric, len_min = 1, .min = .Machine$double.eps)

    n <- list(ctmin, ctmax, b) |>
        sapply(length) |>
        max()

    if (! length(ctmin) %in% c(1L, n)) stop("`ctmin` must be of length 1 or ", n)
    if (! length(ctmax) %in% c(1L, n)) stop("`ctmax` must be of length 1 or ", n)
    if (! length(b) %in% c(1L, n)) stop("`b` must be of length 1 or ", n)

    (2 * ctmax + ctmin + b * ctmin + sqrt(-4 * (2 + b) * ctmax * ctmin +
                                          (-2 * ctmax - ctmin - b * ctmin)^2)) /
        (2 * (2 + b))

    # (2 ctmax + ctmin + b ctmin +
    #     Sqrt[-4 (2 + b) ctmax ctmin + (-2 ctmax - ctmin - b ctmin)^2]) /
    # (2 (2 + b))
    #
}

