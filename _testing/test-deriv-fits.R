
library(tidyverse)
library(nls.multstart)
library(TPCdesign)



set_theme(theme_classic() +
              theme(strip.background = element_blank()))




# -----------------------*
# Derivative approach
# -----------------------*

deriv_temp_maker <- function(n_temps, z0, z1, ctmin, ctmax, a, b, temp_buffer = 5,
                             temp_by = 0.1) {
    # n_temps = 5; z0 = 1; z1 = 0; ctmin = 5; ctmax = 40; a = 1; b = 0.2
    # temp_buffer = 5; temp_by = 0.1
    # rm(n_temps, z0, z1, ctmin, ctmax, a, b, temp_buffer, temp_by)
    # rm(z_by_T, deriv_wts, unif_wts, abs_derivs, all_temps, temp_digits)
    # rm(temp_max, temp_min, z, wts, cs_wts, temp_idx)
    stopifnot(length(temp_buffer) == 1 && is.numeric(temp_buffer) && temp_buffer >= 0)
    temp_min <- ctmin - temp_buffer
    temp_max <- ctmax + temp_buffer
    temp_digits <- pmax(0, -1 * floor(log10(c(temp_min, temp_max, temp_by) %% 1)))
    all_temps <- seq(temp_min, temp_max, temp_by) |>
        round(max(temp_digits))
    abs_derivs <- abs(briere2_tpc_deriv(all_temps, ctmin, ctmax, a, b))
    deriv_wts <- abs_derivs / sum(abs_derivs)
    unif_wts <- rep(1 / length(abs_derivs), length(abs_derivs))
    z <- seq(z0, z1, length.out = length(abs_derivs))
    wts <- z * deriv_wts + (1 - z) * unif_wts
    wts <- wts / sum(wts)
    cs_wts <- cumsum(wts)
    temp_idx <- sapply(1:n_temps, \(x) which(cs_wts > x/(n_temps + 1))[1])
    return(all_temps[temp_idx])
}


# Totally derivative-based:
deriv_temp_maker(n_temps = 5, z0 = 1, z1 = 1, ctmin = 5, ctmax = 40, a = 1, b = 0.2)
# Totally uniform:
deriv_temp_maker(n_temps = 5, z0 = 0, z1 = 0, ctmin = 5, ctmax = 40, a = 1, b = 0.2)
# Linear increase:
deriv_temp_maker(n_temps = 5, z0 = 0, z1 = 1, ctmin = 5, ctmax = 40, a = 1, b = 0.2)
# Totally increase, starting at 0.5:
deriv_temp_maker(n_temps = 5, z0 = 0.5, z1 = 1, ctmin = 5, ctmax = 40, a = 1, b = 0.2)




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
b = 2
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

# ==================================================*
# how well do these temps work compared to evenly distributed?
# ==================================================*


# Takes ~4 min
set.seed(2039799578)
fits <- c(0, 0.5, 1) |>
    round(1) |>
    set_names() |>
    map(\(z) {
        temps <- deriv_temp_maker(n_temps, 0, z, ctmin, ctmax, a, b, temp_buffer = 5)
        lapply(1:100, \(i) {
            obs <- sim_gamma_data(temps, n_reps, obs_cv, ctmin, ctmax, a, b,
                                  scale_tpc = FALSE)
            fit <- nls_multstart(
                formula     = y ~ a * temp * (temp - CTmin) *
                    pmax(CTmax - temp, 0)^exp(b),
                data        = obs,
                start_lower = c(a = 0,  CTmin = 0,  CTmax = 30, b = log(0.01)),
                start_upper = c(a = 2,  CTmin = 15, CTmax = 50, b = log(3)),
                iter        = 50,
                supp_errors = "Y",
                control = list(maxfev = 5e3, maxiter = 1e3),
                lhstype = "improved")
            return(fit)
        })
    }, .progress = TRUE)



coefs <- map(fits, \(x) {
    z <- t(sapply(x, coef))[, c("CTmin", "CTmax", "a", "b")]
    colnames(z) <- tolower(colnames(z))
    z[,"b"] <- exp(z[,"b"])
    return(z)
})

Tgrid <- seq(temp_min, temp_max, length.out = 101)

true_grid <- briere2_tpc(temp = Tgrid, ctmin = ctmin, ctmax = ctmax, b = b, a = a,
                         scale = FALSE)

# We'll compute predicted curve by replicating the same Brière expression:
pred_grids <- imap(coefs, \(cfs, z) {
    # cfs = coefs[[1]]; z = names(coefs)[[1]]
    # rm(cfs, z, y, nT, nR, k, i)
    nT <- length(Tgrid)
    nR <- nrow(cfs)
    y <- numeric(nT * nR)
    k <- 1
    for (i in 1:nR) {
        y[k:(k+nT-1L)] <- briere2_tpc(temp = Tgrid,
                                      ctmin = cfs[i,"ctmin"],
                                      ctmax = cfs[i,"ctmax"],
                                      a = cfs[i,"a"],
                                      b = cfs[i,"b"],
                                      scale = FALSE)

        k <- k + nT
    }
    tibble(z = as.numeric(.env$z),
           temp = rep(Tgrid, nR),
           rep = rep(1:nR, each = nT),
           y = .env$y)
}) |>
    list_rbind() |>
    mutate(z_fct = sprintf("z = %.1f", z) |>
               (\(x) ifelse(x == "z = 0.0", "uniform", x))() |>
               factor())



# --------------------*
# overall curves:

pred_grids |>
    mutate(id = interaction(z_fct, rep, drop = TRUE)) |>
    ggplot(aes(temp, y)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(aes(group = id), alpha = 0.1, color = "dodgerblue") +
    geom_line(data = pred_grids |>
                  group_by(z, z_fct, temp) |>
                  summarize(y = mean(y), .groups = "drop"),
              color = "dodgerblue", linewidth = 1) +
    geom_line(data = tibble(temp = Tgrid, y = true_grid),
              linewidth = 1.5, linetype = "22") +
    # coord_cartesian(ylim = c(0, NA)) +
    # coord_cartesian(xlim = c(39, 40), ylim = c(0, NA))
    facet_wrap(~ z_fct)



# --------------------*
# individual parameters:




imap(coefs, \(x, n) {
    as_tibble(x) |> mutate(z = as.numeric(n))
}) |>
    list_rbind() |>
    # filter(z %in% c(0, 0.5, 1)) |>
    mutate(z_fct = sprintf("z = %.1f", z) |>
               (\(x) ifelse(x == "z = 0.0", "uniform", x))() |>
               factor()) |>
    pivot_longer(ctmin:b, names_to = "param") |>
    mutate(param = factor(param, levels = colnames(coefs[[1]]))) |>
    ggplot(aes(z, value)) +
    geom_violin(aes(color = z_fct)) +
    geom_jitter(aes(color = z_fct), shape = 1, alpha = 0.2, width = 0.02) +
    geom_hline(data = tibble(param = factor(colnames(coefs[[1]]),
                                            levels = colnames(coefs[[1]])),
                             y = c(ctmin, ctmax, a, b)),
               aes(yintercept = y), color = "black", linetype = "22",
               linewidth = 1.5) +
    stat_summary(fun = "mean", geom = "point") +
    # stat_smooth(formula = y ~ s(x), se = FALSE, method = mgcv::gam) +
    facet_wrap(~ param, scales = "free") +
    labs(x = "z parameter (0: uniform, 1: derivative)") +
    # scale_color_brewer(palette = "Dark2")
    scale_color_viridis_d(end = 0.9, guide = "none")


