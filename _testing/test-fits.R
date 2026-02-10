


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


# --------------------*
# overall curves:

pred_grids |>
    ggplot(aes(temp, y, color = type, linetype = type, linewidth = type)) +
    geom_line() +
    scale_linewidth_manual(values = c(rep(1, 5), 1.5)) +
    scale_linetype_manual(values = c(rep("solid", 5), "22")) +
    scale_color_manual(values = c(viridisLite::plasma(length(coefs),
                                                      begin = 0.2), "black")) +
    coord_cartesian(ylim = c(0, NA))

# --------------------*
# individual parameters:

tibble(param = names(coefs[[1]]),
       true = c(ctmin, ctmax, a, b),
       designed = coefs$designed,
       even = coefs$even)



