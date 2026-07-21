
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



# # combos of curve parameters (not naming param to distinguish from fitted parameters)
# set.seed(567834)
# curve_combos <- tibble(skew = names(scenario_bs)) |>
#     mutate(a = 1,
#            b = map_dbl(skew, \(x) scenario_bs[[x]])) |>
#     mutate(ctmin = map(1:n(), \(i) rnorm(100, mean = 5, sd = 2)),
#            ctmax = map(1:n(), \(i) rnorm(100, mean = 40, sd = 2))) |>
#     unnest(c(ctmin, ctmax))


# SCENARIOS WHERE FITTING OFTEN FAILS
#
# CTmin = 5, CTmax = 40, a = 1, b=0.5, CV = 0.2, 6 temperature treatments, 4 replicate observations per temperature
# CTmin = 5, CTmax = 40, a = 1, b=0.5, CV = 0.2, 5 temperature treatments, 3 replicate observations per temperature
# CTmin = 5, CTmax = 40, a = 1, b=2, CV = 0.2, 5 temperature treatments, 4 replicate observations per temperature


ctmin = 5
ctmax = 40
a = 1
b = 0.5 # also 2
obs_cv = 0.2
n_temps = 5L
n_reps = 3L


temp_min = ctmin - 5
temp_max = ctmax + 5
# ctmin = curve_combos$ctmin[[1]]
# ctmax = curve_combos$ctmax[[1]]
# a = curve_combos$a[[1]]
# b = curve_combos$b[[1]]
# rm(ctmin, ctmax, a, b)

# -----------------------*
# Derivative approach
# -----------------------*

deriv_temp_maker <- function(q, ctmin, ctmax, a, b, temp_buffer = 5) {
    stopifnot(length(temp_buffer) == 1 && is.numeric(temp_buffer) && temp_buffer >= 0)
    temp_min <- ctmin - temp_buffer
    temp_max <- ctmax + temp_buffer
    all_temps <- temp_min:temp_max
    abs_derivs <- abs(briere2_tpc_deriv(all_temps, ctmin, ctmax, a, b))
    wts <- q * abs_derivs / sum(abs_derivs) +
        (1 - q) * rep(1 / length(abs_derivs), length(abs_derivs))
    cs_wts <- cumsum(wts)
    temp_idx <- sapply(1:n_temps, \(x) which(cs_wts > x/(n_temps + 1))[1])
    return(all_temps[temp_idx])
}


deriv_temp_maker(q = 0.2, ctmin, ctmax, a, b, temp_buffer = 5)
deriv_temp_maker(q = 0.1, ctmin, ctmax, a, b, temp_buffer = 5)
bench::mark(deriv = deriv_temp_maker(q = 0, ctmin, ctmax, a, b, temp_buffer = 5),
            optim = design_temps(n_temps = n_temps,
                                 n_reps = n_reps,
                                 obs_cv = obs_cv,
                                 temp_min = ctmin - 5,
                                 temp_max = ctmax + 5,
                                 ctmin = ctmin,
                                 ctmax = ctmax,
                                 a = a,
                                 b = b,
                                 scale_tpc = TRUE),
            check = FALSE, memory = FALSE)


# -----------------------*
# Optimization approach
# -----------------------*
# Takes ~3 sec
op <- design_temps(n_temps = n_temps,
                   n_reps = n_reps,
                   obs_cv = obs_cv,
                   temp_min = ctmin - 5,
                   temp_max = ctmax + 5,
                   ctmin = ctmin,
                   ctmax = ctmax,
                   a = a,
                   b = b,
                   scale_tpc = TRUE)
op

# ==================================================*
# how well do these temps work compared to evenly distributed?
# ==================================================*


# Takes ~1.5 min
fits <- list(designed = op$temps,
             q4 = deriv_temp_maker(q = 0.4, ctmin, ctmax, a, b, temp_buffer = 5),
             q2 = deriv_temp_maker(q = 0.2, ctmin, ctmax, a, b, temp_buffer = 5),
             q1 = deriv_temp_maker(q = 0.1, ctmin, ctmax, a, b, temp_buffer = 5),
             q0 = deriv_temp_maker(q = 0, ctmin, ctmax, a, b, temp_buffer = 5),
             even = head(seq(temp_min, temp_max, length.out = n_temps+2L)[-1], -1)) |>
    map(\(temps) {
        temps <- round(temps, 2)
        lapply(1:100, \(i) {
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
    })



coefs <- map(fits, \(x) {
    z <- t(sapply(x, coef))[, c("CTmin", "CTmax", "a", "b")]
    colnames(z) <- tolower(colnames(z))
    return(z)
})

Tgrid <- seq(temp_min, temp_max, length.out = 101)

true_grid <- briere2_tpc(temp = Tgrid, ctmin = ctmin, ctmax = ctmax, b = b, a = a,
                         scale = TRUE)

# We'll compute predicted curve by replicating the same Brière expression:
pred_grids <- imap(coefs, \(cfs, type) {
    y <- briere2_tpc(temp = Tgrid,
                     ctmin = median(cfs[,"ctmin"]),
                     ctmax = median(cfs[,"ctmax"]),
                     a = median(cfs[,"a"]),
                     b = median(cfs[,"b"]),
                     scale = TRUE)
    tibble(type = .env$type, temp = Tgrid, y = .env$y)
}) |>
    list_rbind() |>
    bind_rows(tibble(type = "true", temp = Tgrid, y = true_grid)) |>
    mutate(type = factor(type, levels = c(names(coefs), "true")))


# --------------------*
# overall curves:

pred_grids |>
    # filter(!grepl("q", paste(type))) |>
    (\(x) {
        n_coefs <<- length(unique(x$type)) - 1L
        return(x)
    })() |>
    ggplot(aes(temp, y, color = type, linetype = type, linewidth = type)) +
    geom_line() +
    scale_linewidth_manual(values = c(rep(1, n_coefs), 1.5)) +
    scale_linetype_manual(values = c(rep("solid", n_coefs), "22")) +
    scale_color_manual(values = c(viridisLite::plasma(n_coefs,
                                                      begin = 0.2), "black")) +
    # coord_cartesian(ylim = c(0, NA))
    coord_cartesian(xlim = c(39, 40), ylim = c(0, NA))

# --------------------*
# individual parameters:




imap(coefs, \(x, n) {
    as_tibble(x) |> mutate(type = n)
}) |>
    list_rbind() |>
    mutate(type = factor(type, levels = c(names(coefs), "true"))) |>
    pivot_longer(ctmin:b, names_to = "param") |>
    mutate(param = factor(param, levels = colnames(coefs[[1]]))) |>
    ggplot(aes(value, after_stat(density), color = type)) +
    geom_freqpoly(bins = 20) +
    geom_vline(data = tibble(param = factor(colnames(coefs[[1]]),
                                            levels = colnames(coefs[[1]])),
                             true = c(ctmin, ctmax, a, b)),
               aes(xintercept = true), color = "black", linetype = "22",
               linewidth = 1.5) +
    facet_wrap(~ param, scales = "free") +
    scale_color_brewer(palette = "Dark2")


