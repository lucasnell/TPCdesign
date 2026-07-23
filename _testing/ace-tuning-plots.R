

source("_testing/00-preamble.R")

library(randomForest)
library(iml)


out_files <- list(mtry = "rf-mtry.rds") |>
    map(\(x) paste0("_testing/interm-data/randomForest/", x))

ace_df <- list.files("_testing/interm-data/ace-tuning", "ace-tuning.*.csv", full.names = TRUE) |>
    read_csv(col_types = "iidddddiiiiddddid")


# Columns that I varied across sims:
# ace_df |>
#     select(-combo, -rep, -rmse) |>
#     select(where(\(x) n_distinct(x) > 1)) |>
#     colnames()
# # [1] "n_temps"   "b"         "ctmin_eps" "ctmax_eps" "logb_eps"  "min_sep"
# # [7] "n_filler"  "n_draws"   "n_starts"


ace_summ_df <- ace_df |>
    group_by(across(n_temps:n_starts)) |>
    summarize(sd_rmse = sd(rmse, na.rm = TRUE),
              log_rmse = mean(log10(rmse), na.rm = TRUE),
              rmse = mean(rmse, na.rm = TRUE),
              .groups = "drop")


# Takeaways with n_temps = 5, b = 0.2:
# - n_filler == 1 is very helpful to prevent awful (rmse > 1e6) fits, but
#   more makes the fits slightly worse
# - *_eps don't seem to matter
# - n_starts doesn't do much
# - min_sep = 2 is generally better
#
#


# Overall (only filter is n_filler >= 1) results:
# - n_draws = 250 seems to be mildly best
# - n_starts has little effect
# - n_filler = 1 works well for all
#
#

ace_summ_df |>
    # filter(b == 0.2, n_temps == 5) |>
    filter(n_filler >= 1) |>
    ggplot(aes(min_sep, rmse)) +
    geom_point(color = "gray60", size = 2, shape = 1) +
    stat_smooth(formula = y ~ s(x, bs = "cs"), method = "gam", se = TRUE, linewidth = 1.5) +
    # stat_summary(fun = "mean", geom = "point", size = 4, shape = 5, stroke = 1) +
    facet_wrap(~ interaction(b, n_temps, sep = " - "), scales = "free_y", ncol = 3) +
    scale_color_viridis_c(option = "plasma", begin = 0.2, end = 0.95) +
    theme(panel.spacing = unit(0, "lines"),
          strip.text = element_text())



