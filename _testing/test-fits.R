
library(tidyverse)
library(nls.multstart)
library(estimatePMR)
# library(distionary)
# library(distplyr)
library(TPCdesign)





make_temps <- TPCdesign:::make_temps
make_temps2 <- TPCdesign:::make_temps2





fit_coefs <- function(obs, verbose, nls = TRUE) {

    if (nls) {

        # 4) Fit model to these n_temps points
        # First try nls_multstart which is fast but often has issues converging
        fit <- nls_multstart(
            formula     = y ~ a * temp * (temp - CTmin) * pmax(CTmax - temp, 0)^b,
            data        = obs,
            start_lower = c(a = 0, CTmin = 0,  CTmax = 30, b = 0.01),
            start_upper = c(a = 2, CTmin = 15, CTmax = 50, b = 3),
            iter        = 50,
            # supp_errors = "Y",
            control = list(maxfev = 5e3, maxiter = 1e3),
            lhstype = "improved")

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

    fn <- function(params, obs) {
        y <- briere2_tpc(obs$temp,
                         CTmin = exp(params[1]),
                         CTmax = exp(params[2]),
                         a = exp(params[3]),
                         b = params[4])
        return(sqrt(mean((obs$y - y)^2)))
    }

    op <- optim(par = c(log(7.5), log(40), log(1), 0.15), fn = fn, obs = obs)
    coefs <- op$par
    coefs[1:3] <- exp(coefs[1:3])
    names(coefs) <- c("CTmin", "CTmax", "a", "b")
    return(coefs)


    # # If that doesn't work, use winnowing_optim which is slower but more reliable
    # if (!fit$convInfo$isConv) {
    #     fit <- winnowing_optim(\(p, obs) {
    #         y <- briere2_tpc(obs$temp,
    #                          a = exp(p[[1]]),
    #                          CTmin = exp(p[[2]]),
    #                          CTmax = exp(p[[3]]),
    #                          b = p[[4]],
    #                          scale = FALSE)
    #         return(sqrt(mean((obs$y - y)^2)))
    #     },
    #     lower_bounds = c(a = log(1e-9),
    #                      CTmin = log(1),
    #                      CTmax = log(30),
    #                      b = 0.01),
    #     upper_bounds = c(a = log(2),
    #                      CTmin = log(15),
    #                      CTmax = log(50),
    #                      b = 3),
    #     fn_args = list(obs = obs),
    #     n_bevals = 100L,
    #     n_boxes = 1000L,
    #     n_outputs = 1L,
    #     controls = list(list(maxit = 1000, reltol = 1e-08)),
    #     optimizers = c(optim))
    #
    #     # Extract fitted coefficients
    #     coefs <- fit[[1]][["par"]]
    #     names(coefs) <- c("a", "CTmin", "CTmax", "b")
    #     for (i in 1:3) coefs[[i]] <- exp(coefs[[i]])
    #
    #     if (fit[[1]][["convergence"]] != 0L) {
    #         if (verbose) warning("fit didn't converge")
    #         # This will result in objective function returning very large
    #         # output if fit fails
    #         return(NULL)
    #     }
    #
    # } else {
    #
    #     # Extract fitted coefficients
    #     coefs <- coef(fit)
    #
    # }
    #
    # return(coefs)
}





# ============================================================================*
# ============================================================================*
# ============================================================================*


scenario_bs <- c(
    severe_skew = 0.2,
    mild_skew   = 0.5,
    near_sym    = 2
)

# combos of curve parameters (not naming param to distinguish from fitted parameters)
set.seed(567834)
curve_combos <- tibble(skew = names(scenario_bs)) |>
    mutate(a = 1,
           b = map_dbl(skew, \(x) scenario_bs[[x]])) |>
    mutate(CTmin = map(1:n(), \(i) rnorm(100, mean = 5, sd = 2)),
           CTmax = map(1:n(), \(i) rnorm(100, mean = 40, sd = 2))) |>
    unnest(c(CTmin, CTmax))




temp_fn <- function(params,
                    n_temps,
                    n_reps,
                    obs_cv,
                    temp_min,
                    temp_max,
                    CTmin,
                    CTmax,
                    b,
                    a,
                    nls = FALSE,
                    output = "RMSE",
                    verbose = FALSE) {

    # params = c(1, -3, 0, 0)
    # n_temps = 6; n_reps = 5; obs_cv = 0.2; temp_min = 5; temp_max = 40
    # CTmin = curve_combos$CTmin[[1]]; CTmax = curve_combos$CTmax[[1]]
    # b = curve_combos$b[[1]]; a = curve_combos$a[[1]]
    # nls = FALSE; output = "RMSE"; verbose = FALSE
    # rm(params, n_temps, n_reps, obs_cv, temp_min, temp_max, CTmin, CTmax, b, a)
    # rm(nls, output, verbose)

    output <- match.arg(tolower(output), c("rmse", "ctmin", "ctmax", "b", "a"))

    # design_temps <- make_temps(params, temp_min, temp_max, n_temps)
    design_temps <- make_temps2(params, temp_min, temp_max)

    # 2) Compute "true" unscaled performance at design points (with a = 1)
    true_vals <- briere2_tpc(temp = design_temps, a = a, CTmin = CTmin,
                             CTmax = CTmax, b = b)

    # 3) Simulate gamma data at these design points
    obs <- map2_dfr(design_temps, true_vals, \(temp, val) {
        tibble(temp = rep(temp, n_reps),
               y = rgamma(n_reps, shape = 1 / (obs_cv^2), scale = val * (obs_cv^2)))
    })

    coefs <- fit_coefs(obs, verbose, nls)
    # This only happens when coefs couldn't be fit
    if (is.null(coefs)) return(1e10)



    if (output == "rmse") {
        # 5) Evaluate RMSE across a dense temperature grid =====#
        # Construct a dense grid
        Tgrid <- seq(temp_min, temp_max, length.out = 101)

        # Calculate "True" curve across that grid, using the same param set & scaling
        #    We'll recompute unscaled with a=1, then apply the same scale_factor
        true_grid <- briere2_tpc(temp = Tgrid, a = a, CTmin = CTmin,
                                 CTmax = CTmax, b = b)

        # 6) Fitted curve across Tgrid
        #    We must pass each T value to predict(), but we need a new data frame
        #    with columns that match the formula's variables (esp. T).
        #    Also note we need a column for 'a', 'CTmin', 'CTmax','b' if they're in the formula.
        #    But nls won't require them if they're coefficients. We'll do partial approach:

        # a, CTmin, CTmax
        # We'll compute predicted curve by replicating the same Brière expression:
        pred_grid <- briere2_tpc(temp = Tgrid,
                                 a = coefs[["a"]],
                                 CTmin = coefs[["CTmin"]],
                                 CTmax = coefs[["CTmax"]],
                                 b = coefs[["b"]],
                                 scale = FALSE)

        out <- sqrt(mean((true_grid - pred_grid)^2))
        if (is.infinite(out)) out <- 1e10

    } else if (output == "ctmin") {

        out <- abs(CTmin - coefs[["CTmin"]])

    } else if (output == "ctmax") {

        out <- abs(CTmax - coefs[["CTmax"]])

    } else if (output == "b") {

        out <- abs(b - coefs[["b"]])

    } else {

        out <- abs(a - coefs[["a"]])

    }

    return(out)

}






# f <- function() {
#     temp_fn(params = c(0, 1, 1, 1),
#                     n_temps = 6,
#                     n_reps = 5,
#                     obs_cv = 0.2,
#                     temp_min = 0,
#                     temp_max = 50,
#                     CTmin = curve_combos$CTmin[[1]],
#                     CTmax = curve_combos$CTmax[[1]],
#                     b = curve_combos$b[[1]],
#                     a = curve_combos$a[[1]])
# }
#
# # bench::mark(f(), min_iterations = 10)


n_temps = 6L
n_reps = 5L
obs_cv = 0.2
temp_min = 0
temp_max = 50
CTmin = curve_combos$CTmin[[1]]
CTmax = curve_combos$CTmax[[1]]
b = curve_combos$b[[1]]
a = curve_combos$a[[1]]



t0 <- Sys.time()
# op <- optim(par = c(1, -3, 0, 0), fn = temp_fn,
op <- optim(par = rep(0, n_temps+1L), fn = temp_fn,
             n_temps = n_temps,
             n_reps = n_reps,
             obs_cv = obs_cv,
             temp_min = temp_min,
             temp_max = temp_max,
             CTmin = CTmin,
             CTmax = CTmax,
             b = b,
             a = a)
t1 <- Sys.time()
t1 - t0
# Time difference of 1.963169 secs


t0 <- Sys.time()
# op2 <- optim(par = c(1, -3, 0, 0), fn = temp_fn,
op2 <- optim(par = rep(0, n_temps+1L), fn = temp_fn,
             n_temps = n_temps,
             n_reps = n_reps,
             obs_cv = obs_cv,
             temp_min = temp_min,
             temp_max = temp_max,
             CTmin = CTmin,
             CTmax = CTmax,
             b = b,
             a = a,
             nls = FALSE)
t1 <- Sys.time()
t1 - t0
# Time difference of 0.04184008 secs



t0 <- Sys.time()
op3 <- winnowing_optim(temp_fn,
                       # lower_bounds = c(logit(c(1e-6, 1e-6)), rep(log(1e-6), 2)),
                       # upper_bounds = c(logit(1 - c(1e-6, 1e-6)), rep(log(200), 2)),
                       lower_bounds = rep(logit(1e-9), n_temps+1L),
                       upper_bounds = rep(logit(1-1e-9), n_temps+1L),
                       fn_args = list(n_temps = n_temps,
                                      n_reps = n_reps,
                                      obs_cv = obs_cv,
                                      temp_min = temp_min,
                                      temp_max = temp_max,
                                      CTmin = CTmin,
                                      CTmax = CTmax,
                                      b = b,
                                      a = a),
                       n_bevals = 100L,
                       n_boxes = 100L,
                       n_outputs = 1L,
                       controls = list(list(maxit = 1000, reltol = 1e-08)),
                       optimizers = c(optim))[[1]]
t1 <- Sys.time()
t1 - t0
# Time difference of 4.336942 mins
# Didn't converge!


t0 <- Sys.time()
op4 <- winnowing_optim(temp_fn,
                       # lower_bounds = c(logit(c(1e-6, 1e-6)), rep(log(1e-6), 2)),
                       # upper_bounds = c(logit(1 - c(1e-6, 1e-6)), rep(log(200), 2)),
                       lower_bounds = rep(logit(1e-9), n_temps+1L),
                       upper_bounds = rep(logit(1-1e-9), n_temps+1L),
                       fn_args = list(n_temps = n_temps,
                                      n_reps = n_reps,
                                      obs_cv = obs_cv,
                                      temp_min = temp_min,
                                      temp_max = temp_max,
                                      CTmin = CTmin,
                                      CTmax = CTmax,
                                      b = b,
                                      a = a,
                                      nls = FALSE),
                       n_bevals = 100L,
                       n_boxes = 100L,
                       n_outputs = 1L,
                       controls = list(list(maxit = 1000, reltol = 1e-08)),
                       optimizers = c(optim))[[1]]
t1 <- Sys.time()
t1 - t0
# Time difference of 3.944689 mins
# Also did not converge!



fits <- map(list("op" = op, "op2" = op2, "op3" = op3, "op4" = op4, "even" = NA),
            \(op) {
                if (isTRUE(is.na(op))) {
                    temps <- seq(temp_min, temp_max, length.out = n_temps)
                } else {
                    temps <- make_temps2(op$par, temp_min, temp_max)
                }
                true_vals <- briere2_tpc(temp = temps, CTmin = CTmin,
                                         CTmax = CTmax, b = b, a = a)
                obs <- map2_dfr(temps, true_vals, \(temp, val) {
                    tibble(temp = rep(temp, n_reps),
                           y = rgamma(n_reps, shape = 1 / (obs_cv^2),
                                      scale = val * (obs_cv^2)))
                })
                fit <- nls_multstart(
                    formula     = y ~ a * temp * (temp - CTmin) *
                        pmax(CTmax - temp, 0)^b,
                    data        = obs,
                    start_lower = c(a = 0,  CTmin = 0,  CTmax = 30, b = 0.01),
                    start_upper = c(a = 2,  CTmin = 15, CTmax = 50, b = 3),
                    iter        = 50,
                    supp_errors = "Y",
                    control = list(maxfev = 5e3, maxiter = 1e3),
                    lhstype = "improved")
                return(fit)
            })

coefs <- map(fits, \(x) coef(x)[c("CTmin", "CTmax", "a", "b")])


Tgrid <- seq(temp_min, temp_max, length.out = 101)

# Calculate "True" curve across that grid, using the same param set & scaling
#    We'll recompute unscaled with a=1, then apply the same scale_factor
true_grid <- briere2_tpc(temp = Tgrid, CTmin = CTmin, CTmax = CTmax,
                         b = b, a = a)
# We'll compute predicted curve by replicating the same Brière expression:
pred_grids <- imap(coefs, \(cfs, type) {
    y <- briere2_tpc(temp = Tgrid, a = cfs[["a"]], CTmin = cfs[["CTmin"]],
                     CTmax = cfs[["CTmax"]], b = cfs[["b"]])
    tibble(type = .env$type, temp = Tgrid, y = .env$y)
    }) |>
    list_rbind() |>
    bind_rows(tibble(type = "true", temp = Tgrid, y = true_grid)) |>
    mutate(type = factor(type, levels = c(names(coefs), "true")))

pred_grids |>
    ggplot(aes(temp, y, color = type, linetype = type, linewidth = type)) +
    geom_line() +
    scale_linewidth_manual(values = c(rep(1, 5), 1.5)) +
    scale_linetype_manual(values = c(rep("solid", 5), "22")) +
    scale_color_manual(values = c(viridisLite::plasma(5, begin = 0.2), "black"))


cbind(true = c(CTmin, CTmax, a, b), do.call(cbind, coefs))


