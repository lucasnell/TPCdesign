

source("_testing/00-preamble.R")



ace_df <- list.files("_testing/interm-data/ace-test", "big-ace-test.*.csv", full.names = TRUE) |>
    map(\(x) {
        read_csv(x, col_types = "iiciiddddddld")
    }) |>
    list_rbind()



# ============================================================================*
# ============================================================================*
# Plots ----
# ============================================================================*
# ============================================================================*


ace_df |>
    # make `b` all the real version:
    group_by(combo) |>
    mutate(b = b[method == "real"]) |>
    ungroup() |>
    filter(method != "real") |>
    select(method, combo, n_temps, n_reps, b, obs_cv, converged) |>
    group_by(method, combo, n_temps, n_reps, b, obs_cv) |>
    summarize(conv = mean(converged), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values > 1 are good!
    summarize(diff_conv = conv[method == "design"] / conv[method == "uniform"], .groups = "drop") |>
    getElement("diff_conv") |>
    (\(x) {print(mean(x > 1)); return(x)})() |>
    hist(xlab = "Optimized-based converged / uniform converged", main = NULL)


ace_df |>
    # make `b` all the real version:
    group_by(combo) |>
    mutate(b = b[method == "real"]) |>
    ungroup() |>
    filter(method != "real") |>
    select(method, combo, n_temps, n_reps, b, obs_cv, rmse) |>
    group_by(method, combo, n_temps, n_reps, b, obs_cv) |>
    summarize(rmse = mean(rmse, na.rm = TRUE), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values < 1 are good!
    summarize(diff_rmse = rmse[method == "design"] / rmse[method == "uniform"], .groups = "drop") |>
    getElement("diff_rmse") |>
    (\(x) {print(mean(x < 1)); return(x)})() |>
    # (\(x) sign(x) * log10(abs(x)))() |>
    hist(xlab = "Optimized-based RMSE / uniform RMSE", main = NULL)



ace_df |>
    # make `b` all the real version:
    group_by(combo) |>
    mutate(b = b[method == "real"]) |>
    ungroup() |>
    filter(method != "real") |>
    select(method, combo, n_temps, n_reps, b, obs_cv, rmse) |>
    group_by(method, combo, n_temps, n_reps, b, obs_cv) |>
    summarize(rmse = mean(rmse, na.rm = TRUE), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values < 1 are good!
    summarize(diff_rmse = rmse[method == "design"] / rmse[method == "uniform"],
              .groups = "drop") |>
    group_by(b) |>
    summarize(prob_better = mean(diff_rmse < 1))



rel_p <- ace_df |>
    # make `b` all the real version:
    group_by(combo) |>
    mutate(b = b[method == "real"]) |>
    ungroup() |>
    filter(method != "real") |>
    select(method, combo, n_temps, n_reps, b, obs_cv, rmse) |>
    group_by(method, combo, n_temps, n_reps, b, obs_cv) |>
    summarize(rmse = mean(rmse, na.rm = TRUE), .groups = "drop") |>
    group_by(n_temps, n_reps, b, obs_cv) |>
    # Values < 1 are good!
    summarize(diff_rmse = rmse[method == "design"] / rmse[method == "uniform"], .groups = "drop") |>
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
    group_by(combo) |>
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
           method = factor(method, levels = c("uniform", "design"),
                           labels = c("uniform", "optimized"))) |>
    ggplot(aes(n_temps, n_reps, fill = rmse)) +
    geom_raster() +
    facet_grid(b ~ method) +
    scale_fill_viridis_c(option = "magma", begin = 0.2) +
    labs(x = "Number of temperature treatments",
         y = "Replicates per temperature",
         subtitle = "*CV* = 0.1") +
    theme(strip.text.y = element_markdown(family = "serif", angle = 0),
          strip.text.x = element_markdown(family = "serif"),
          plot.subtitle = element_markdown())


rel_val_p <- rel_p + val_p + plot_layout(ncol = 1, heights = c(1, 0.8))
# rel_val_p

ggsave("_testing/_plots/rel-val-optim1.pdf", rel_val_p, width = 6, height = 6)
ggsave("_testing/_plots/rel-ace-optim.pdf", rel_p, width = 6, height = 6)


# ============================================================================*
# ============================================================================*
# Testing problem areas ----
# ============================================================================*
# ============================================================================*

source("_testing/big-ace-test.R")

library(future.apply)
library(progressr)

plan(multisession, workers = 4L)
handlers(handler_cli(format = paste("{cli::pb_bar} {cli::pb_percent} |",
                                    "{cli::pb_elapsed} | ETA: {cli::pb_eta}"),
                     clear = FALSE))





foo <- function(i, prog) {

    # .ctmin_err = 2.50; .ctmax_err = 1.50; .logb_err = 0.26
    .ctmin_err = 0; .ctmax_err = 0; .logb_err = 0

    args <- list(n_temps = 5L, b = 1,
                 # n_reps = 3L, obs_cv = 0.1,
                 ctmin = 5, ctmax = 40, a = 1)
    args[["ctmin"]] <- runif(1, args$ctmin - 2 * .ctmin_err,
                             args$ctmin + 2 * .ctmin_err)
    args[["ctmax"]] <- runif(1, args$ctmax - 2 * .ctmax_err,
                             args$ctmax + 2 * .ctmax_err)
    args[["b"]] <- exp(runif(1, log(args$b) - 2 * .logb_err,
                             log(args$b) + 2 * .logb_err))

    prog()

    return(do.call(ace_design_temps, args))
}


goo <- function(.seed) {

    .n_test_fits = 100L

    with_progress({
        p <- progressor(.n_test_fits)
        out <- future_lapply(1:.n_test_fits,
                             foo,
                             prog = p,
                             future.seed = .seed,
                             future.packages = c("tidyverse",
                                                 "TPCdesign")) |>
            do.call(what = rbind)
    })


    return(out)
}


dtemps_list <- goo(100)

colMeans(dtemps_list)
# [1]  9.1889 15.9090 22.8040 29.6880 36.5703  << on cluster
apply(dtemps_list, 2, min)
# [1]  2.91 11.09 17.88 24.95 31.80  << on cluster
apply(dtemps_list, 2, max)
# [1] 12.22 21.18 28.21 34.26 40.84  << on cluster



temps <- colMeans(dtemps_list)

set.seed(1)
one_test_fit(1L, temps, n_reps = 3L, obs_cv = 0.1, ctmin = 5, ctmax = 40,
             a = 1, b = 1)



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Function ace_design_temps gives different outputs on cluster vs locally.
# Testing where the difference starts to occur.

n_temps = 5L
ctmin = 5
ctmax = 40
a = 1
b = 1
ctmin_err = 2.50
ctmax_err = 1.50
logb_err = 0.26
min_sep = 1
n_filler = 1L
n_draws = 250L
n_starts = 7L
digits = 2L
n_threads = 1L



n_optimal <- n_temps - n_filler

# genus-level relatives was 2.50°C for Tmin, 1.50°C for Tmax,
# and 0.26 for log(b)
prior_ctmin <- ctmin + c(-1,1) * ctmin_err
prior_ctmax <- ctmax + c(-1,1) * ctmax_err
prior_lb <- log(b) + c(-1,1) * logb_err


set.seed(456789)
theta_draws <- cbind(ctmin = runif(n_draws, prior_ctmin[1], prior_ctmin[2]),
                     ctmax = runif(n_draws, prior_ctmax[1], prior_ctmax[2]),
                     b = exp(runif(n_draws, prior_lb[1], prior_lb[2])),
                     a = rep(a, n_draws))
colMeans(theta_draws)




min_temp <- min(prior_ctmin)
max_temp <- max(prior_ctmax)

if ((n_optimal - 1L) * min_sep > (max_temp - min_temp)) {
    stop("min_sep is too large for given min and max temperatures and the ",
         "number of desired optimized temps")
}

set.seed(3456789)
start_temp_list <- lapply(1:n_starts, function(i) {
    lhs_min <- ctmin + min_sep
    lhs_max <- ctmax - min_sep
    d <- lhs::improvedLHS(n = n_optimal, k = 1) * (lhs_max - lhs_min) + lhs_min
    d <- d[order(d),,drop=FALSE]
    colnames(d) <- "temp"
    return(d)
})

start_temp_list |> do.call(what = cbind) |> rowMeans()


# Limit temps to have a minimum separation:
limits_minsep <- function(d, i, j) {
    grid <- seq(from = min_temp, to = max_temp, length.out = 2000)
    other_points <- as.vector(d)[-i]
    for (s in other_points) {
        grid <- grid[(grid < (s - min_sep)) | (grid > (s + min_sep))]
    }
    return(grid)
}



set.seed(23467890)
design_robust <- acebayes::pace(utility = TPCdesign:::utility_briere2D,
                                start.d = start_temp_list,
                                B = rep(list(theta = theta_draws), 2),
                                lower = min_temp,
                                upper = max_temp,
                                limits = limits_minsep,
                                N2 = 0L,
                                deterministic = TRUE,
                                mc.cores = n_threads)
opt_temps <- sort(round(design_robust$d, digits))

final_temps <- TPCdesign:::gap_filler(n_filler, opt_temps, min_temp, max_temp, digits)



curve(briere2_tpc(x, 5, 40, 1, 1), 0, 45)
abline(v = opt_temps, col ="red", lty = "22")
abline(v = final_temps[!final_temps %in% opt_temps], col ="blue", lty = "22")


TPCdesign:::utility_briere2D(start_temp_list[[1]], theta_draws)
utility_briere2D(start_temp_list[[1]], theta_draws)


one_test_fit(1L, final_temps, n_reps = 3L, obs_cv = 0.1, ctmin = 5, ctmax = 40,
             a = 1, b = 1)

one_test_fit(1L, c(9.37, 15.40, 21.43, 28.70, 32.44), n_reps = 3L,
             obs_cv = 0.1, ctmin = 5, ctmax = 40,
             a = 1, b = 1)



derivs <- list(
    # dy / dctmin:
    \(temp) -a * temp * (ctmax - temp)^b,
    # dy / dctmax:
    \(temp) a * temp * (temp - ctmin) * b * (ctmax - temp)^(b - 1.0),
    # dy / db:
    \(temp) a * temp * (temp - ctmin) * (ctmax - temp)^b * log((ctmax - temp)),
    # dy / da:
    \(temp) temp * (temp - ctmin) * (ctmax - temp)^b)

temps <- c(9.37, 15.40, 21.43, 28.70, 32.44)


# Gradient
grad_at_temps <- function(temps) {
    G <- matrix(0, length(temps), 4)
    inside <- which(temps > ctmin & temps < ctmax)
    for (i in 1:4) {
        d <- derivs[[i]]
        G[inside,i] <- d(temps[inside])
    }
    return(G)
}
info_at_temps <- function(temps) {
    G <- grad_at_temps(temps)
    # TPCdesign:::logdet(t(G) %*% G)

    lapply(seq_len(nrow(G)), function(i) outer(G[i, ], G[i, ]))
}





greedy_d_optimal <- function(T_candidates, k, n_start = 4) {

    I_candidates <- info_at_temps(T_candidates)

    # seed with a few points spread across the range to avoid singular start
    chosen_idx <- round(seq(1, length(T_candidates), length.out = n_start))
    I_total <- Reduce(`+`, I_candidates[chosen_idx])

    while (length(chosen_idx) < k) {
        best_det <- -Inf
        best_idx <- NA
        for (i in seq_along(T_candidates)) {
            if (i %in% chosen_idx) next
            I_try <- I_total + I_candidates[[i]]
            d <- det(I_try)
            if (d > best_det) {
                best_det <- d
                best_idx <- i
            }
        }
        chosen_idx <- c(chosen_idx, best_idx)
        I_total <- I_total + I_candidates[[best_idx]]
    }

    sort(T_candidates[chosen_idx])
}

greedy_d_optimal(seq(5, 40, 0.1), k = 4, n_start = 2)


with(list(G = grad_at_temps(c(5.0, 27.9, 39.1, 40.0))), {TPCdesign:::logdet(t(G) %*% G)})

with(list(G = grad_at_temps(final_temps)), {TPCdesign:::logdet(t(G) %*% G)})

with(list(G = grad_at_temps(c(9.37, 15.40, 21.43, 28.70, 32.44))), {TPCdesign:::logdet(t(G) %*% G)})








score <- function(temp) {
    warm_gap = ctmax - temp
    warm_gap_b = warm_gap^b
    cold_gap = temp - ctmin
    out <- numeric(4)
    # dy / dctmin:
    out[1] = -a * temp * warm_gap_b
    # dy / dctmax:
    out[2] = a * temp * cold_gap * b * warm_gap^(b - 1.0)
    # dy / db:
    out[3] = a * temp * cold_gap * warm_gap_b * log(warm_gap)
    # dy / da:
    out[4] = temp * cold_gap * warm_gap_b
    return(out)
}


map(opt_temps, score) |> do.call(what = rbind)
map(c(9.37, 21.43, 28.70, 32.44), score) |> do.call(what = rbind)

map(opt_temps, score) |> do.call(what = rbind) |> rowSums()
map(c(9.37, 21.43, 28.70, 32.44), score) |> do.call(what = rbind) |> rowSums()





# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




# Takes ~2 min
x0 <- one_combo_fits(1L, tibble(n_temps = 5L, n_reps = 3L, b = 1, obs_cv = 0.1,
                                ctmin = 5, ctmax = 40, a = 1,
                                n_filler = 1L, n_threads = 1L))
x0 |>
    filter(method != "real") |>
    group_by(method) |>
    summarize(rmse = mean(rmse)) |>
    (\(xx) {print(xx); return(xx)})() |>
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "uniform"]))

# # A tibble: 2 × 2
#   method   rmse
#   <chr>   <dbl>
# 1 design   354.
# 2 uniform  382.
# # A tibble: 1 × 1
#      rel
#    <dbl>
# 1 -0.113



 # Takes ~1 min
x2 <- one_combo_fits(1L, tibble(n_temps = 8L, n_reps = 3L, b = 1, obs_cv = 0.1,
                                 ctmin = 5, ctmax = 40, a = 1),
                    .opt_temp_p = 1, min_sep = 2, n_threads = 5L)
x2 |>
    filter(method != "real") |>
    group_by(method) |>
    summarize(rmse = mean(rmse)) |>
    (\(xx) {print(xx); return(xx)})() |>
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "uniform"]))
# # A tibble: 2 × 2
#   method  rmse
#   <chr>  <dbl>
# 1 design  181.
# 2 uniform    228.
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
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "uniform"]))
# # A tibble: 2 × 2
#   method  rmse
#   <chr>  <dbl>
# 1 design  98.5
# 2 uniform   233.
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
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "uniform"]))
# # A tibble: 2 × 2
#   method  rmse
#   <chr>  <dbl>
# 1 design  143.
# 2 uniform    180.
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
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "uniform"]))


one_combo_fits(1L, crossing(n_temps = 10L, n_reps = 3L, b = 0.2, obs_cv = 0.4,
                            ctmin = 5, ctmax = 40, a = 1),
               .opt_temp_p = 1, n_threads = 5L) |>
    filter(method != "real") |>
    group_by(method) |>
    summarize(rmse = mean(rmse)) |>
    (\(xx) {print(xx); return(xx)})() |>
    summarize(rel = log2(rmse[method == "design"] / rmse[method == "uniform"]))
