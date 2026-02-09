# For avoiding warnings for comparing nonsensible inputs
comparable <- function(x) {
    return(!is.null(x) && (is.vector(x) || is.array(x)) && !any(is.na(x)))
}

# Check for a single, whole number, perhaps in range
single_integer <- function(x, par, .min = NULL, .max = NULL) {
    all_good <- TRUE
    if (!comparable(x)) {
        all_good <- FALSE
    } else {
        all_good <- is.numeric(x) && length(x) == 1 && x %% 1 == 0
        if (all_good && !is.null(.min)) all_good <- x >= .min
        if (all_good && !is.null(.max)) all_good <- x <= .max
    }
    if (!all_good) {
        err_suff <- "a single integer"
        if (!is.null(.min)) err_suff <- paste0(err_suff, ", min = ", .min)
        if (!is.null(.max)) err_suff <- paste0(err_suff, ", max = ", .max)
        err_msg(par, err_suff)
    }
    invisible(NULL)
}
# Check for a single number, perhaps in range
single_number <- function(x, par, .min = NULL, .max = NULL) {
    all_good <- TRUE
    if (!comparable(x)) {
        all_good <- FALSE
    } else {
        all_good <- is.numeric(x) && length(x) == 1
        if (all_good && !is.null(.min)) all_good <- x >= .min
        if (all_good && !is.null(.max)) all_good <- x <= .max
    }
    if (!all_good) {
        err_suff <- "a single number"
        if (!is.null(.min)) err_suff <- paste0(err_suff, ", min = ", .min)
        if (!is.null(.max)) err_suff <- paste0(err_suff, ", max = ", .max)
        err_msg(par, err_suff)
    }
    invisible(NULL)
}
single_string <- function(x, par) {
    all_good <- comparable(x) && is.character(x) && length(x) == 1
    if (!all_good) err_msg(par, "a single string")
    invisible(NULL)
}
single_logical <- function(x, par) {
    all_good <- comparable(x) && is.logical(x) && length(x) == 1
    if (!all_good) err_msg(par, "a single logical")
    invisible(NULL)
}
# general type checker
is_type <- function(x, par, type, len_min = NULL, len_max = NULL, .min = NULL, .max = NULL) {
    all_good <- TRUE
    if (!comparable(x)) {
        all_good <- FALSE
    } else {
        if (is.character(type)) all_good <- inherits(x, type)
        if (is.function(type)) all_good <- type(x)
        if (all_good && !is.null(len_min)) all_good <- length(x) >= len_min
        if (all_good && !is.null(len_max)) all_good <- length(x) <= len_max
        if (all_good && !is.null(.min)) all_good <- all(x >= .min)
        if (all_good && !is.null(.max)) all_good <- all(x <= .max)
    }
    if (!all_good) {
        err_suff <- paste("an object of type", type)
        if (!is.null(L)) err_suff <- paste(err_suff, "with possible length(s) =",
                                           paste(L, collapse = ", "))
        if (!is.null(.min)) err_suff <- paste0(err_suff, "; min value(s) = ", .min)
        if (!is.null(.max)) err_suff <- paste0(err_suff, "; max value(s) = ", .max)
        err_msg(par, err_suff)
    }
    invisible(NULL)
}

# Standard way to show error messages (also to make input-checking less verbose).
#
#
err_msg <- function(par, ...) {
    stop(sprintf("\nArgument `%s` must be %s.",
                 par, paste(...)), call. = FALSE)
}

