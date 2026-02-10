


library(tidyverse)
library(nls.multstart)
library(estimatePMR)
# library(distionary)
# library(distplyr)
library(TPCdesign)




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
    mutate(ctmin = map(1:n(), \(i) rnorm(100, mean = 5, sd = 2)),
           ctmax = map(1:n(), \(i) rnorm(100, mean = 40, sd = 2))) |>
    unnest(c(ctmin, ctmax))


n_temps = 6L
n_reps = 5L
obs_cv = 0.2
temp_min = 0
temp_max = 50
ctmin = curve_combos$ctmin[[1]]
ctmax = curve_combos$ctmax[[1]]
a = curve_combos$a[[1]]
b = curve_combos$b[[1]]
# rm(ctmin, ctmax, a, b)




# Takes ~10 sec
op <- design_temps(n_temps = n_temps,
                   n_reps = n_reps,
                   obs_cv = obs_cv,
                   temp_min = temp_min,
                   temp_max = temp_max,
                   ctmin = ctmin,
                   ctmax = ctmax,
                   a = a,
                   b = b,
                   scale_tpc = TRUE)
op


# ==================================================*
# how well do these temps work compared to evenly distributed?
# ==================================================*


fits <- map(list(designed = op, even = NA),
            \(op) {
                if (isTRUE(is.na(op))) {
                    temps <- seq(temp_min, temp_max, length.out = n_temps+2L)[-1]
                    temps <- head(temps, -1L)
                } else {
                    temps <- op$temps
                }
                obs <- sim_gamma_data(temps, n_reps, obs_cv, ctmin, ctmax, a, b,
                                      scale_tpc = TRUE)
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

coefs <- map(fits, \(x) {
    z <- coef(x)[c("CTmin", "CTmax", "a", "b")]
    names(z) <- tolower(names(z))
    return(z)
})

Tgrid <- seq(temp_min, temp_max, length.out = 101)

true_grid <- briere2_tpc(temp = Tgrid, ctmin = ctmin, ctmax = ctmax, b = b, a = a,
                         scale = TRUE)

# We'll compute predicted curve by replicating the same Brière expression:
pred_grids <- imap(coefs, \(cfs, type) {
    y <- briere2_tpc(temp = Tgrid,
                     ctmin = cfs[["ctmin"]],
                     ctmax = cfs[["ctmax"]],
                     a = cfs[["a"]],
                     b = cfs[["b"]],
                     scale = TRUE)
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
    scale_color_manual(values = c(viridisLite::plasma(length(coefs),
                                                      begin = 0.2), "black")) +
    coord_cartesian(ylim = c(0, NA))


cbind(true = c(ctmin, ctmax, a, b), do.call(cbind, coefs))




# ============================================================================*
# ============================================================================*
# ============================================================================*





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
    design_temps <- make_temps(params, temp_min, temp_max)

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
                    temps <- make_temps(op$par, temp_min, temp_max)
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


