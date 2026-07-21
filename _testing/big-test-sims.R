
source("_scripts/00-preamble.R")
source("_scripts/00-sim-functions.R")

library(lhs)






# =============================================================================*
# =============================================================================*
# Simulations ----
# =============================================================================*
# =============================================================================*


#
# Factors that I'll vary by simulation:
#
par_ranges <- list(n_temps = c(5, 10),
                   n_reps = c(3, 10),
                   z1 = c(0, 1),
                   b = c(0.2, 2),
                   obs_cv = 0.2 * c(0.5, 2),
                   temp_buffer = c(1, 10))

# Create LHS data frame, then transform variables to proper scales
set.seed(191621615)
lhs_df <- improvedLHS(n = 1e3, k = length(par_ranges)) |>
    (\(x) {colnames(x) <- names(par_ranges); return(x)})() |>
    as_tibble()
for (n in names(par_ranges)) {
    p1 <- min(par_ranges[[n]])
    p2 <- max(par_ranges[[n]])
    if (! n %in% c("n_temps", "n_reps", "temp_buffer")) {
        lhs_df[[n]] <- p1 + lhs_df[[n]] * (p2 - p1)
    } else {
        # These need to be treated differently:
        lhs_df[[n]] <- qinteger(lhs_df[[n]], p1, p2)
    }
}; rm(n, p1, p2)

# Now add variables that won't vary:
lhs_df[["z0"]] <- 0
lhs_df[["ctmin"]] <- 5
lhs_df[["ctmax"]] <- 40
lhs_df[["a"]] <- 1



if (!file.exists("_testing/test-fits.csv")) {

    # Takes ~40 min
    fit_df <- all_combo_fits(lhs_df)
    write_csv(fit_df, "_testing/test-fits.csv")

} else {
    fit_df <- read_csv("_testing/test-fits.csv", col_types = "iiciidddiddddd")
}


#
# Lots of noise here, but we can conclude a couple of things:
#
# 1. z1 = 1 produces lowest rmse
# 2. temp_buffer = 6 because of the following:
#     a. higher values produce low levels of non-convergence but are more
#        likely to have divergent (especially high RMSE) fits
#     b. 6 gives good balance of relatively high convergence and low divergence
#     c. Note: this parameter has little overall effect on rmse
#

fit_df |>
    filter(type == "fitted") |>
    group_by(combo, n_temps, n_reps, z1, obs_cv, temp_buffer) |>
    summarize(rmse = mean(log10(rmse)), .groups = "drop") |>
    ggplot(aes(z1, (rmse))) +
    geom_jitter(aes(), shape = 1) +
    stat_smooth(method = "gam", formula = y ~ s(x)) +
    scale_color_viridis_c()

fit_df |>
    filter(type == "fitted") |>
    mutate(temp_buffer = factor(temp_buffer)) |>
    group_by(combo, n_temps, n_reps, z1, obs_cv, temp_buffer) |>
    summarize(rmse = sd((rmse)) / mean(rmse), .groups = "drop") |>
    ggplot(aes(temp_buffer, rmse)) +
    geom_jitter(aes(), shape = 1) +
    stat_summary(fun = "mean", geom = "point", color = "red", size = 3) +
    scale_color_viridis_c()

fit_df |>
    filter(type == "fitted") |>
    group_by(combo, n_temps, n_reps, z1, obs_cv, temp_buffer) |>
    summarize(p_conv = mean(converged), .groups = "drop") |>
    mutate(temp_buffer = factor(temp_buffer)) |>
    ggplot(aes(temp_buffer, p_conv)) +
    geom_jitter(aes(), shape = 1) +
    stat_summary(fun = "mean", geom = "point", color = "red", size = 3) +
    scale_color_viridis_c()



