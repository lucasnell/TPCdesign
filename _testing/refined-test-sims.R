

source("_testing/00-preamble.R")
source("_testing/00-sim-functions.R")


# Output files:
csv_files <- list(perfect = "_testing/refined-test-fits.csv",
                  inaccurate = "_testing/refined-off-test-fits.csv")



# =============================================================================*
# =============================================================================*
# Simulations ----
# =============================================================================*
# =============================================================================*



# # ------------------------------------------------*
# # ... perfect parameters ----
# # ------------------------------------------------*
#
# if (!file.exists(csv_files$perfect)) {
#
#     # Takes ~50 min
#     fit_df <- crossing(n_temps = 5:10,
#                        n_reps = 3:10,
#                        b = c(0.2, 0.5, 1, 2),
#                        obs_cv = 0.2 * c(0.5, 1, 2),
#                        # To compare to even temps (z0 = 0):
#                        z1 = c(0, 1),
#                        # These curves are completely accurate:
#                        curve_off = FALSE,
#                        # Things that don't vary:
#                        temp_buffer = 6L,
#                        z0 = 0,
#                        ctmin = 5,
#                        ctmax = 40,
#                        a = 1) |>
#         all_combo_fits(.seed = 1634404298) |>
#         select(-curve_off)
#     write_csv(fit_df, csv_files$perfect)
#
# } else {
#     fit_df <- read_csv(csv_files$perfect, col_types = "iiciidddidddddld")
# }
#
#
#
# fit_df |>
#     # make `b` all the real version:
#     group_by(combo, rep) |>
#     mutate(b = b[type == "real"]) |>
#     ungroup() |>
#     filter(type == "fitted") |>
#     select(combo, n_temps, n_reps, b, obs_cv, z1, rmse) |>
#     group_by(combo, n_temps, n_reps, b, obs_cv, z1) |>
#     summarize(rmse = mean(rmse), .groups = "drop") |>
#     group_by(n_temps, n_reps, b, obs_cv) |>
#     # Values < 0 (if using `-`) or < 1 (if using `\`) are good!
#     summarize(diff_rmse = rmse[z1 == 1] / rmse[z1 == 0], .groups = "drop") |>
#     getElement("diff_rmse") |>
#     (\(x) {print(mean(x < 1)); return(x)})() |>
#     # (\(x) sign(x) * log10(abs(x)))() |>
#     hist(xlab = "Derivative-based RMSE / uniform RMSE", main = NULL)
#
#
#
#
#
# fit_df |>
#     # make `b` all the real version:
#     group_by(combo, rep) |>
#     mutate(b = b[type == "real"]) |>
#     ungroup() |>
#     filter(type == "fitted") |>
#     select(combo, n_temps, n_reps, b, obs_cv, z1, rmse) |>
#     group_by(combo, n_temps, n_reps, b, obs_cv, z1) |>
#     summarize(rmse = mean(rmse), .groups = "drop") |>
#     group_by(n_temps, n_reps, b, obs_cv) |>
#     # Values < 0 (if using `-`) or < 1 (if using `\`) are good!
#     summarize(diff_rmse = rmse[z1 == 1] / rmse[z1 == 0], .groups = "drop") |>
#     mutate(across(n_temps:n_reps, factor),
#            b = factor(b, labels = sprintf("<i>b</i> = %.1f", sort(unique(b)))),
#            obs_cv = factor(obs_cv, labels = sprintf("<i>CV</i> = %.1f", sort(unique(obs_cv))))) |>
#     ggplot(aes(n_temps, n_reps, fill = log2(diff_rmse))) +
#     geom_raster() +
#     facet_grid(b ~ obs_cv) +
#     # scale_fill_scico("log<sub>2</sub>(RMSE<sub>deriv.</sub> / RMSE<sub>unif.</sub>)",
#     scale_fill_scico(expression(log[2](frac(RMSE[deriv], RMSE[unif]))),
#                      palette = "vik", midpoint = 0) +
#     labs(x = "Number of temperature treatments",
#          y = "Replicates per temperature") +
#     theme(# legend.title = element_markdown(),
#           strip.text.y = element_markdown(angle = 0),
#           strip.text.x = element_markdown(),
#           axis.ticks = element_blank(),
#           axis.line = element_blank())
#
#
#
# fit_df |>
#     # make `b` all the real version:
#     group_by(combo, rep) |>
#     mutate(b = b[type == "real"]) |>
#     ungroup() |>
#     filter(type == "fitted") |>
#     select(combo, n_temps, n_reps, b, obs_cv, z1, rmse) |>
#     group_by(combo, n_temps, n_reps, b, obs_cv, z1) |>
#     summarize(rmse = mean(log10(rmse)), .groups = "drop") |>
#     filter(obs_cv == 0.2) |>
#     mutate(across(n_temps:n_reps, factor),
#            b = factor(b, labels = sprintf("<i>b</i> = %.1f", sort(unique(b)))),
#            obs_cv = factor(obs_cv, labels = sprintf("<i>CV</i> = %.1f", sort(unique(obs_cv)))),
#            z1 = factor(as.integer(z1), levels = 0:1,
#                        labels = c("uniform", "deriv."))) |>
#     ggplot(aes(n_temps, n_reps, fill = rmse)) +
#     geom_raster() +
#     facet_grid(b ~ z1) +
#     scale_fill_viridis_c(option = "magma") +
#     theme(strip.text.y = element_markdown(family = "serif", angle = 0),
#           strip.text.x = element_markdown(family = "serif"))






# ------------------------------------------------*
# ... inaccurate parameters ----
# ------------------------------------------------*


# Takes ~3.5 min per row
set.seed(75359468)
off_fit_small_df <- tibble(n_temps = 5, n_reps = 10, b = 0.2, obs_cv = 0.2,
                           ctmin = 5, ctmax = 50, a = 1) |>
    one_combo_fits(j = 1L, prog = NULL)

off_fit_small_df |>
    filter(!is.na(method)) |>
    group_by(method) |>
    summarize(rmse = mean(rmse))
# # A tibble: 2 × 2
#   method  rmse
#   <chr>  <dbl>
# 1 design  242.
# 2 even    476.


set.seed(1977579122)
off_fit_small_df2 <- map(1:2, \(j) {
    one_combo_fits(j = j,
                   input_df = crossing(n_temps = 9, n_reps = 9, b = c(0.2, 2),
                                       obs_cv = 0.2, ctmin = 5, ctmax = 50,
                                       a = 1),
                   prog = NULL)
    }) |>
    list_rbind()

off_fit_small_df2 |>
    filter(!is.na(method)) |>
    group_by(combo, method) |>
    summarize(rmse = mean(rmse), .groups = "drop")
# # A tibble: 4 × 3
#   combo method  rmse
#   <int> <chr>  <dbl>
# 1     1 design  162.
# 2     1 even    366.
# 3     2 design 8187.
# 4     2 even   9137.



off_fit_small_df2 |>
    filter(!is.na(method), combo == 2) |>
    mutate(diff = ifelse(method == "design", abs(rmse - 9387), abs(rmse - 7965))) |>
    group_by(rep) |>
    summarize(diff = mean(diff)) |>
    arrange(diff)

off_fit_small_df2 |>
    filter(!is.na(method), combo == 2, rep == 22)


curve(briere2_tpc(x, 5, 50, 1, 2, FALSE), 0, 55, ylab = NA)
curve(briere2_tpc(x, 4.90, 49.7, 0.872, 2.04, FALSE), add = TRUE, col = "red")
curve(briere2_tpc(x, 4.3, 49.0, 1.65, 1.84, FALSE), add = TRUE, col = "red", lty = 2)


temps <- ace_design_temps(15L, 5, 50, 1, 2)

curve(briere2_tpc(x, 5, 50, 1, 2, FALSE), 0, 55, ylab = NA)
abline(v = temps, lty = 2, col = "red")


y = "ctmax"


off_fit_small_df2 |>
    filter(!is.na(method), combo == 2) |>
    ggplot(aes(.data[[y]])) +
    geom_histogram(bins = 25, fill = "dodgerblue", alpha = 0.5) +
    geom_vline(data = off_fit_small_df2 |>
                   filter(!is.na(method), combo == 2) |>
                   group_by(method) |>
                   summarize(across(all_of(y), mean)),
               aes(xintercept = .data[[y]]),
               color = "dodgerblue") +
    geom_vline(data = off_fit_small_df2 |> filter(is.na(method), combo == 2) |>
                   slice(1) |> select(-method),
               aes(xintercept = .data[[y]]), linetype = "22", color = "red") +
    facet_wrap(~ method, ncol = 1)









# FULL SIMULATIONS ----

if (!file.exists(csv_files$inaccurate)) {

    # Takes ~1 hr
    off_fit_df <- crossing(n_temps = 5:10,
                           n_reps = 3:10,
                           b = c(0.2, 0.5, 1, 2),
                           obs_cv = 0.2 * c(0.5, 1, 2),
                           # Things that don't vary:
                           ctmin = 5,
                           ctmax = 40,
                           a = 1) |>
        all_combo_fits(.seed = 2140098495)
    write_csv(off_fit_df, csv_files$inaccurate)

} else {

    off_fit_df <- read_csv(csv_files$inaccurate, col_types = "iiciidddidddddld")

}





off_fit_df |>
    # make `b` all the real version:
    group_by(combo, rep) |>
    mutate(b = b[type == "real"]) |>
    ungroup() |>
    filter(type == "fitted") |>
    select(combo, n_temps, n_reps, b, obs_cv, z1, converged) |>
    group_by(combo, n_temps, n_reps, b, obs_cv, z1) |>
    summarize(conv = mean(converged), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values > 0 (if using `-`) or > 1 (if using `\`) are good!
    summarize(diff_conv = conv[z1 == 1] / conv[z1 == 0], .groups = "drop") |>
    getElement("diff_conv") |>
    (\(x) {print(mean(x > 1)); return(x)})() |>
    hist(xlab = "Derivative-based converged / uniform converged", main = NULL)


off_fit_df |>
    # make `b` all the real version:
    group_by(combo, rep) |>
    mutate(b = b[type == "real"]) |>
    ungroup() |>
    filter(type == "fitted") |>
    select(combo, n_temps, n_reps, b, obs_cv, z1, rmse) |>
    group_by(combo, n_temps, n_reps, b, obs_cv, z1) |>
    summarize(rmse = mean(rmse), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values < 0 (if using `-`) or < 1 (if using `\`) are good!
    summarize(diff_rmse = rmse[z1 == 1] / rmse[z1 == 0], .groups = "drop") |>
    getElement("diff_rmse") |>
    (\(x) {print(mean(x < 1)); return(x)})() |>
    # (\(x) sign(x) * log10(abs(x)))() |>
    hist(xlab = "Derivative-based RMSE / uniform RMSE", main = NULL)



off_fit_df |>
    # make `b` all the real version:
    group_by(combo, rep) |>
    mutate(b = b[type == "real"]) |>
    ungroup() |>
    filter(type == "fitted") |>
    select(combo, n_temps, n_reps, b, obs_cv, z1, rmse) |>
    group_by(combo, n_temps, n_reps, b, obs_cv, z1) |>
    summarize(rmse = mean(rmse), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values < 0 (if using `-`) or < 1 (if using `\`) are good!
    summarize(diff_rmse = rmse[z1 == 1] / rmse[z1 == 0], .groups = "drop") |>
    mutate(across(n_temps:n_reps, factor),
           b = factor(b, labels = sprintf("<i>b</i> = %.1f", sort(unique(b)))),
           obs_cv = factor(obs_cv, labels = sprintf("<i>CV</i> = %.1f", sort(unique(obs_cv))))) |>
    ggplot(aes(n_temps, n_reps, fill = log2(diff_rmse))) +
    geom_raster() +
    facet_grid(b ~ obs_cv) +
    # scale_fill_scico("log<sub>2</sub>(RMSE<sub>deriv.</sub> / RMSE<sub>unif.</sub>)",
    scale_fill_scico(expression(log[2](frac(RMSE[deriv], RMSE[unif]))),
                     palette = "vik", midpoint = 0) +
    labs(x = "Number of temperature treatments",
         y = "Replicates per temperature") +
    theme(# legend.title = element_markdown(),
        strip.text.y = element_markdown(angle = 0),
        strip.text.x = element_markdown(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

off_fit_df |>
    # make `b` all the real version:
    group_by(combo, rep) |>
    mutate(b = b[type == "real"]) |>
    ungroup() |>
    filter(type == "fitted") |>
    filter(b == 0.2) |>
    filter(obs_cv == 0.2) |>
    select(combo, n_temps, n_reps, b, obs_cv, z1, rmse) |>
    group_by(combo, n_temps, n_reps, b, obs_cv, z1) |>
    summarize(rmse = mean((rmse)), .groups = "drop") |>
    mutate(across(n_temps:n_reps, factor),
           b = factor(b, labels = sprintf("<i>b</i> = %.1f", sort(unique(b)))),
           obs_cv = factor(obs_cv, labels = sprintf("<i>CV</i> = %.1f", sort(unique(obs_cv)))),
           z1 = factor(as.integer(z1), levels = 0:1,
                       labels = c("uniform", "deriv."))) |>
    ggplot(aes(n_temps, n_reps, fill = rmse)) +
    geom_raster() +
    facet_grid(b ~ z1) +
    scale_fill_viridis_c(option = "magma") +
    labs(x = "Number of temperature treatments",
         y = "Replicates per temperature") +
    theme(strip.text.y = element_markdown(family = "serif", angle = 0),
          strip.text.x = element_markdown(family = "serif"))





# Exact true values:
tru <- list(ctmin = 5,
            ctmax = 40,
            a = 1,
            b = 0.2,
            n_temps = 8L,
            n_reps = 7L,
            obs_cv = 0.2)


combos <- off_fit_df |>
    filter(type == "real", n_temps == tru$n_temps, n_reps == tru$n_reps,
           b == tru$b, obs_cv == tru$obs_cv) |>
    getElement("combo") |> unique()

off_fit_df |>
    filter(combo %in% combos) |>
    filter(type == "fitted") |>
    filter(rmse == max(rmse))


off_fit_df |>
    filter(combo %in% combos) |>
    pivot_longer(ctmin:Topt, names_to = "parameter") |>
    group_by(z1, combo, rep, parameter) |>
    summarize(value = value[type == "fitted"] - value[type == "real"],
              .groups = "drop") |>
    filter(parameter == "Topt") |>
    filter(abs(value) == max(abs(value)))

off_fit_df |>
    filter(combo %in% combos) |>
    select(-n_temps, -n_reps, -b, -obs_cv)
off_fit_df |>
    filter(combo == 676, rep == 14) |>
    select(-n_temps, -n_reps, -b, -obs_cv)




# For RMSE
test_temps <- seq(tru$ctmin, tru$ctmax, length.out = 101)
tru_y <- briere2_tpc(test_temps, tru$ctmin, tru$ctmax, tru$a, tru$b, FALSE)


topt_split_design_temps <- function(n_temps, ctmin, ctmax, b) {
    # Trying out innacurate Topt, split down middle.
    topt <- briere2_tpc_Topt(ctmin, ctmax, b)
    n_above <- n_temps %/% 2L
    n_below <- n_temps - n_above + 1L
    temps <- c(seq(ctmin, topt, length.out = n_below),
               seq(topt, ctmax, length.out = n_above)) |>
        unique() |>
        sort()
    stopifnot(length(temps) == n_temps)
    return(temps)
}

# # Versions used for temps (not necessarily true):
# tmp <- tru
# tmp$ctmin <- tru$ctmin + 2.50
# tmp$ctmax <- tru$ctmax + 1.50
# tmp$b <- exp(log(tru$b) + 0.26)

if (!file.exists("_testing/rmse-test.rds")) {
    # Takes ~ 15 min
    set.seed(680555590)
    rmse_df <- crossing(z1 = c(0, 1, -1),
                        ctmin = tru$ctmin + (-1:1) * 2.50,
                        ctmax = tru$ctmax + (-1:1) * 1.50,
                        b = exp(log(tru$b) + (-1:1) * 0.26),
                        seed = sample.int(2^31-1, 100)) |>
        mutate(rmse = pmap_dbl(across(everything()), \(z1, ctmin, ctmax, b, seed) {
            set.seed(seed)
            if (z1 < 0) {
                temps <- topt_split_design_temps(tru$n_temps, ctmin, ctmax, b)
            } else {
                temps <- deriv_design_temps(tru$n_temps, z0 = 0, z1 = z1, ctmin, ctmax,
                                            tru$a, b, 6L)
            }
            obs <- sim_gamma_data(temps, n_reps = tru$n_reps, obs_cv = tru$obs_cv,
                                  ctmin = tru$ctmin, ctmax = tru$ctmax, a = tru$a,
                                  b = tru$b, scale_tpc = FALSE)
            fit <- nls_multstart(
                formula     = y ~ a * temp * (temp - ctmin) *
                    pmax(ctmax - temp, 0)^exp(b),
                data        = obs,
                start_lower = c(a = 0,  ctmin = 0,  ctmax = 30, b = log(0.01)),
                start_upper = c(a = 2,  ctmin = 15, ctmax = 50, b = log(3)),
                iter        = 50,
                supp_errors = "Y",
                control = list(maxfev = 5e3, maxiter = 1e3),
                lhstype = "improved")
            fitted <- as.list(coef(fit))[c("ctmin", "ctmax", "a", "b")]
            fitted[["b"]] <- exp(fitted[["b"]])
            obs_y <- briere2_tpc(test_temps, fitted[["ctmin"]], fitted[["ctmax"]],
                                 fitted[["a"]], fitted[["b"]], FALSE)
            # RMSE
            return(sqrt(mean((obs_y - tru_y)^2)))
        }, .progress = TRUE))
    write_rds(rmse_df, "_testing/rmse-test.rds")

} else {

    rmse_df <- read_rds("_testing/rmse-test.rds")

}




rmse_df |>
    filter(ctmin == tru$ctmin) |>
    group_by(z1, b, ctmax) |>
    summarize(rmse = mean(rmse), .groups = "drop") |>
    mutate(b = factor(b, levels = sort(unique(b)),
                      labels = sprintf("%.3f", sort(unique(b)))),
           z1 = factor(z1, levels = -1:1,
                       labels = c("mid", "unif", "deriv"))) |>
    ggplot(aes(b, ctmax)) +
    geom_raster(aes(fill = rmse)) +
    scale_fill_viridis_c(option = "magma") +
    facet_wrap(~ z1)

rmse_df |>
    filter(ctmax == tru$ctmax) |>
    group_by(z1, b, ctmin) |>
    summarize(rmse = mean(rmse), .groups = "drop") |>
    mutate(b = factor(b, levels = sort(unique(b)),
                      labels = sprintf("%.3f", sort(unique(b)))),
           z1 = factor(z1, levels = -1:1,
                       labels = c("mid", "unif", "deriv"))) |>
    ggplot(aes(b, ctmin)) +
    geom_raster(aes(fill = rmse)) +
    scale_fill_viridis_c(option = "magma") +
    facet_wrap(~ z1)




ex_temps <- list(
    s = topt_split_design_temps(ctmin = tru$ctmin, ctmax = tru$ctmax + 1.5,
                                n_temps = tru$n_temps, b = tru$b),
    d = deriv_design_temps(n_temps = tru$n_temps, z0 = 0, z1 = 1,
                           ctmin = tru$ctmin, ctmax = tru$ctmax + 1.5,
                           a = tru$a, b = tru$b, temp_buffer = 6),
    u = deriv_design_temps(n_temps = tru$n_temps, z0 = 0, z1 = 0,
                           ctmin = tru$ctmin, ctmax = tru$ctmax + 1.5,
                           a = tru$a, b = tru$b, temp_buffer = 6)
)


ex_temps |>
    imap(\(x, i) {
        col <- list(s = "dodgerblue", d = "firebrick", u = "gray70")[[i]]
        tibble(temp = seq(tru$ctmin-5, tru$ctmax+5, length.out = 101)) |>
            mutate(y = briere2_tpc(temp, tru$ctmin, tru$ctmax, tru$a, tru$b, FALSE)) |>
            ggplot(aes(temp, y)) +
            geom_line() +
            geom_vline(xintercept = x, linetype = "22", color = col)
    }) |>
    wrap_plots(ncol = 1)




curve(briere2_tpc(x, tru$ctmin, tru$ctmax, tru$a, tru$b,
                  TRUE), 5, 45, ylab = NA)
curve(briere2_tpc(x, tru$ctmin, tru$ctmax + 1.50, tru$a, tru$b,
                  TRUE), add = TRUE, col = "red")
curve(briere2_tpc(x, tru$ctmin, tru$ctmax, tru$a, exp(log(tru$b) - 0.26),
                  TRUE), add = TRUE, lty = 2)
curve(briere2_tpc(x, tru$ctmin, tru$ctmax + 1.50, tru$a, exp(log(tru$b) - 0.26),
                  TRUE), add = TRUE, lty = 2, col = "red")


dxdt <- function(temp, ctmin, ctmax, a, b) {
    z <- abs(briere2_tpc_deriv(temp, ctmin, ctmax, a, b))
    return(z / max(z))
}

curve(dxdt(x, tru$ctmin, tru$ctmax, tru$a, tru$b), 5, 45, ylab = NA)
curve(dxdt(x, tru$ctmin, tru$ctmax + 1.50, tru$a, tru$b), add = TRUE, col = "red")
curve(dxdt(x, tru$ctmin, tru$ctmax, tru$a, exp(log(tru$b) - 0.26)),
      add = TRUE, lty = 2)
curve(dxdt(x, tru$ctmin, tru$ctmax + 1.50, tru$a, exp(log(tru$b) - 0.26)),
      add = TRUE, lty = 2, col = "red")



rmse_df |>
  filter(z1 == 1, ctmin == tru$ctmin, ctmax == max(ctmax)) |>
  group_by(b) |>
  summarize(rmse = mean(rmse), .groups = "drop")

TPCdesign::briere2_tpc_Topt(tru$ctmin, tru$ctmax, tru$b)
TPCdesign::briere2_tpc_Topt(tru$ctmin, tru$ctmax + 1.50, tru$b)

# Trying out innacurate Topt, split down middle.
topt <- TPCdesign::briere2_tpc_Topt(tru$ctmin, tru$ctmax + 1.50, tru$b)
tru$n_temps

temps <- c(seq(tru$ctmin, topt, length.out = tru$n_temps - (tru$n_temps %/% 2L) + 1L),
  seq(topt, tru$ctmax + 1.50, length.out = tru$n_temps %/% 2L)) |>
  unique() |>
  sort()

# LEFT OFF FOR REAL ----
# TRY OPTIMIZATION ON ABOVE TEMPS

obs <- sim_gamma_data(temps, tru$n_reps, tru$obs_cv, tru$ctmin, tru$ctmax,
                      tru$a, tru$b, scale_tpc = FALSE)
fit <- nls_multstart(
    formula     = y ~ a * temp * (temp - ctmin) *
        pmax(ctmax - temp, 0)^exp(b),
    data        = obs,
    start_lower = c(a = 0,  ctmin = 0,  ctmax = 30, b = log(0.01)),
    start_upper = c(a = 2,  ctmin = 15, ctmax = 50, b = log(3)),
    iter        = 50,
    supp_errors = "Y",
    control = list(maxfev = 5e3, maxiter = 1e3),
    lhstype = "improved")
fit
fitted <- as.list(coef(fit))[c("ctmin", "ctmax", "a", "b")]
fitted[["b"]] <- exp(fitted[["b"]])
# RMSE
obs_y <- briere2_tpc(test_temps, fitted[["ctmin"]], fitted[["ctmax"]],
                     fitted[["a"]], fitted[["b"]], FALSE)
sqrt(mean((obs_y - tru_y)^2))





# LEFT OFF -----
# Trying to figure out why in the example above, with a higher (inaccurate) ctmax,
# lower (incaccurate) b gives better rmse
# --> real result or stochastic outcome?



off_fit_df |>
    filter(combo %in% combos) |>
    pivot_longer(ctmin:Topt, names_to = "parameter") |>
    group_by(z1, combo, rep, parameter) |>
    summarize(value = value[type == "fitted"] - value[type == "real"],
              .groups = "drop") |>
    mutate(z1 = factor(as.integer(z1), levels = 0:1,
                       labels = c("uniform", "deriv.")))|>
    ggplot(aes(value)) +
    geom_freqpoly(aes(color = z1), bins = 20, linewidth = 1) +
    scale_color_manual(values = c(uniform = "gray60", `deriv.` = "gold")) +
    facet_wrap(~ parameter, scales = "free")


input_df <- off_fit_df |>
    filter(combo == 682, rep == 8, type == "real") |>
    select(n_temps, n_reps, z0, z1, obs_cv, ctmin, ctmax, a, b, temp_buffer) |>
    mutate(curve_off = TRUE)

# Takes ~13 sec each
fits0 <- one_combo_fits(1, input_df |> mutate(z1 = 0), prog = NULL)
fits <- one_combo_fits(1, input_df, prog = NULL)
fits20 <- one_combo_fits(1, input_df |> mutate(z1 = 0), prog = NULL, .scale_tpc = TRUE)
fits2 <- one_combo_fits(1, input_df, prog = NULL, .scale_tpc = TRUE)

list(fits0, fits, fits20, fits2) |>
    map_dbl(\(x) filter(x, type == "fitted") |> getElement("rmse") |> mean())






# =============================================================================*
# =============================================================================*
# Phylogeny ----
# =============================================================================*
# =============================================================================*



# https://github.com/HuckleyLab/ThermalStress/blob/master/data/tpcs.csv
read_csv("~/Library/CloudStorage/Box-Box/TPCdesign/tpcs.csv",
         col_types = "cccdddcddcc") |>
    filter(!is.na(species), !is.na(genus)) |>
    select(genus, species) |>
    mutate(species = species |>
               str_remove_all("\\ .*") |>
               # remove everything before and including first _:
               str_remove("^.*?_") |>
               # remove everything after and including first _:
               str_remove("_.*")) |>
    distinct(genus, species) |>
    arrange(genus, species) |>
    pmap_chr(\(genus, species) paste(genus, species)) |>
    write_lines("~/Library/CloudStorage/Box-Box/TPCdesign/tpc-species.txt")



