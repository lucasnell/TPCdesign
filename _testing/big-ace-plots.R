

source("_testing/00-preamble.R")



ace_df <- list.files("_testing/interm-data/ace-test/opt_temp_p-0.7", "big-ace-test.*.csv", full.names = TRUE) |>
    map(\(x) {
        read_csv(x, col_types = "iiciiddddddld")
    }) |>
    list_rbind()

if (length(unique(ace_df$combo)) < 560) {
    # correct combo vector:
    ace_df <- ace_df|>
        mutate(combo = rep(1:576, each = 300))
}


# ============================================================================*
# ============================================================================*
# Plots ----
# ============================================================================*
# ============================================================================*


ace_df |>
    # make `b` all the real version:
    group_by(combo, rep) |>
    mutate(b = b[method == "real"]) |>
    ungroup() |>
    filter(method != "real") |>
    select(method, combo, n_temps, n_reps, b, obs_cv, converged) |>
    group_by(method, combo, n_temps, n_reps, b, obs_cv) |>
    summarize(conv = mean(converged), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values > 1 are good!
    summarize(diff_conv = conv[method == "design"] / conv[method == "even"], .groups = "drop") |>
    getElement("diff_conv") |>
    (\(x) {print(mean(x > 1)); return(x)})() |>
    hist(xlab = "Optimized-based converged / uniform converged", main = NULL)


ace_df |>
    # make `b` all the real version:
    group_by(combo, rep) |>
    mutate(b = b[method == "real"]) |>
    ungroup() |>
    filter(method != "real") |>
    select(method, combo, n_temps, n_reps, b, obs_cv, rmse) |>
    group_by(method, combo, n_temps, n_reps, b, obs_cv) |>
    summarize(rmse = mean(rmse, na.rm = TRUE), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values < 1 are good!
    summarize(diff_rmse = rmse[method == "design"] / rmse[method == "even"], .groups = "drop") |>
    getElement("diff_rmse") |>
    (\(x) {print(mean(x < 1)); return(x)})() |>
    # (\(x) sign(x) * log10(abs(x)))() |>
    hist(xlab = "Optimized-based RMSE / uniform RMSE", main = NULL)



rel_p <- ace_df |>
    # make `b` all the real version:
    group_by(combo, rep) |>
    mutate(b = b[method == "real"]) |>
    ungroup() |>
    filter(method != "real") |>
    select(method, combo, n_temps, n_reps, b, obs_cv, rmse) |>
    group_by(method, combo, n_temps, n_reps, b, obs_cv) |>
    summarize(rmse = mean(rmse, na.rm = TRUE), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values < 1 are good!
    summarize(diff_rmse = rmse[method == "design"] / rmse[method == "even"], .groups = "drop") |>
    mutate(across(n_temps:n_reps, factor),
           b = factor(b, labels = sprintf("<i>b</i> = %.1f", sort(unique(b)))),
           obs_cv = factor(obs_cv, labels = sprintf("<i>CV</i> = %.1f", sort(unique(obs_cv))))) |>
    ggplot(aes(n_temps, n_reps, fill = log2(diff_rmse))) +
    geom_raster() +
    facet_grid(b ~ obs_cv) +
    # scale_fill_scico("log<sub>2</sub>(RMSE<sub>deriv.</sub> / RMSE<sub>unif.</sub>)",
    scale_fill_scico(expression(log[2](frac(RMSE[optim], RMSE[unif]))),
                     palette = "vik", midpoint = 0) +
    labs(x = "Number of temperature treatments",
         y = "Replicates per temperature") +
    theme(# legend.title = element_markdown(),
        strip.text.y = element_markdown(angle = 0),
        strip.text.x = element_markdown(),
        axis.ticks = element_blank(),
        axis.line = element_blank())



val_p <- ace_df |>
    # make `b` all the real version:
    group_by(combo, rep) |>
    mutate(b = b[method == "real"]) |>
    ungroup() |>
    filter(method != "real") |>
    filter(b == 1) |>
    filter(obs_cv == 0.1) |>
    select(method, combo, n_temps, n_reps, b, obs_cv, rmse) |>
    group_by(method, combo, n_temps, n_reps, b, obs_cv) |>
    summarize(rmse = mean((rmse), na.rm = TRUE), .groups = "drop") |>
    mutate(across(n_temps:n_reps, factor),
           b = factor(b, labels = sprintf("<i>b</i> = %.1f", sort(unique(b)))),
           obs_cv = factor(obs_cv, labels = sprintf("<i>CV</i> = %.1f", sort(unique(obs_cv)))),
           method = factor(method, levels = c("even", "design"),
                           labels = c("uniform", "optimized"))) |>
    ggplot(aes(n_temps, n_reps, fill = rmse)) +
    geom_raster() +
    facet_grid(b ~ method) +
    scale_fill_viridis_c(option = "magma") +
    labs(x = "Number of temperature treatments",
         y = "Replicates per temperature",
         subtitle = "*CV* = 0.1") +
    theme(strip.text.y = element_markdown(family = "serif", angle = 0),
          strip.text.x = element_markdown(family = "serif"),
          plot.subtitle = element_markdown())


rel_val_p <- rel_p + val_p + plot_layout(ncol = 1, heights = c(1, 0.8))
# rel_val_p

ggsave("_testing/_plots/rel-val-optim1.pdf", rel_val_p, width = 6, height = 6)


# ============================================================================*
# ============================================================================*
# Testing problem areas ----
# ============================================================================*
# ============================================================================*

source("_testing/big-ace-test.R")




# Takes ~1 min
x0 <- one_combo_fits(1L, tibble(n_temps = 8L, n_reps = 3L, b = 1, obs_cv = 0.1,
                                 ctmin = 5, ctmax = 40, a = 1),
                    .opt_temp_p = 1, n_threads = 5L)
x0 |>
    filter(method != "real") |>
    group_by(method) |>
    summarize(rmse = mean(rmse)) |>
    (\(xx) {print(xx); return(xx)})() |>
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "even"]))
# # A tibble: 2 × 2
#   method  rmse
#   <chr>  <dbl>
# 1 design  372.
# 2 even    262.
# # A tibble: 1 × 1
#     rel
#   <dbl>
# 1 0.507


# Takes ~1 min
x2 <- one_combo_fits(1L, tibble(n_temps = 8L, n_reps = 3L, b = 1, obs_cv = 0.1,
                                 ctmin = 5, ctmax = 40, a = 1),
                    .opt_temp_p = 1, min_sep = 2, n_threads = 5L)
x2 |>
    filter(method != "real") |>
    group_by(method) |>
    summarize(rmse = mean(rmse)) |>
    (\(xx) {print(xx); return(xx)})() |>
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "even"]))
# # A tibble: 2 × 2
#   method  rmse
#   <chr>  <dbl>
# 1 design  181.
# 2 even    228.
# # A tibble: 1 × 1
#      rel
#    <dbl>
# 1 -0.328



y0 <- one_combo_fits(1L, tibble(n_temps = 8L, n_reps = 3L, b = 0.2, obs_cv = 0.1,
                                ctmin = 5, ctmax = 40, a = 1),
                    .opt_temp_p = 1, n_threads = 5L)
y0 |>
    filter(method != "real") |>
    group_by(method) |>
    summarize(rmse = mean(rmse)) |>
    (\(xx) {print(xx); return(xx)})() |>
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "even"]))
# # A tibble: 2 × 2
#   method  rmse
#   <chr>  <dbl>
# 1 design  98.5
# 2 even   233.
# # A tibble: 1 × 1
#     rel
#   <dbl>
# 1 -1.24



y2 <- one_combo_fits(1L, crossing(n_temps = 5L, n_reps = 10L, b = 1, obs_cv = 0.1,
                                 ctmin = 5, ctmax = 40, a = 1),
                    .opt_temp_p = 1, min_sep = 2, n_threads = 5L)
y2 |>
    filter(method != "real") |>
    group_by(method) |>
    summarize(rmse = mean(rmse)) |>
    (\(xx) {print(xx); return(xx)})() |>
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "even"]))
# # A tibble: 2 × 2
#   method  rmse
#   <chr>  <dbl>
# 1 design  143.
# 2 even    180.
# # A tibble: 1 × 1
#      rel
#    <dbl>
# 1 -0.334



one_combo_fits(1L, crossing(n_temps = 10L, n_reps = 3L, b = 0.2, obs_cv = 0.4,
                            ctmin = 5, ctmax = 40, a = 1),
               .opt_temp_p = 0.7, n_threads = 5L) |>
    filter(method != "real") |>
    group_by(method) |>
    summarize(rmse = mean(rmse)) |>
    (\(xx) {print(xx); return(xx)})() |>
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "even"]))


one_combo_fits(1L, crossing(n_temps = 10L, n_reps = 3L, b = 0.2, obs_cv = 0.4,
                            ctmin = 5, ctmax = 40, a = 1),
               .opt_temp_p = 1, n_threads = 5L) |>
    filter(method != "real") |>
    group_by(method) |>
    summarize(rmse = mean(rmse)) |>
    (\(xx) {print(xx); return(xx)})() |>
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "even"]))
