# designing informative TPCs

# issues remaining: 
#' 1) is the uncertainty in the priors incorporated correctly? Since this is simulated data, we can pretend we know the 1 true curve, yet we are sampling a different set from the priors each time we propose a design? It seems like it's overcomplicating our demonstration of the idea. But when someone puts our process into practice, they will not know the true curve and will need to specify priors.

# Details about this code:
# Breire 2 functional curve
# equal replicates across temperatures
# randomly generated designs
# full curve RMSE
# only generate positive observed values of the trait (so not normally distributed), keep variance relatively small

# the process is: 1) Define a TPC model and define priors 2) generate designs, simulate data, fit TPC, and evaluate error metric, 3) check convergence and pick top 10% of designs 4) compare RMSE to a uniform design (does it do better?).

# set functions for: 
# briere2_tpc with parameters (T, a, CTmin, CTmax, b), keeping a constant
#   performance <- a * T * (T - CTmin) * (pmax(CTmax - T, 0))^b
# draw_priors
# random_design_generator
# simulate_and_eval_rmse

# Create Figure 1 
# 3 Panels (Symmetric, Mild Skew, Severe Skew)
# Density plots of top 10% designs + overlay best design as points + overlay curve shape

# Temperature–Trait Experiment Design Simulation

# Required packages
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(nls.multstart)

set.seed(111)

# 1) Define the Brière-2 thermal performance curve (TPC) ####
briere2_tpc <- function(T, a, CTmin, CTmax, b) {
  a * T * (T - CTmin) * pmax(CTmax - T, 0)^b
}

# 2) Draw parameter sets from informative priors (excluding 'a', since we'll scale curves)
draw_priors <- function(n, b_val) {
  tibble(
    CTmin = rnorm(n, mean = 5, sd = 2),
    CTmax = rnorm(n, mean = 40, sd = 2),
    b     = b_val
  )
}

# 3) Generate a random design of K temperatures
generate_design <- function(K, temp_range = c(5, 40)) {
  # Define all possible integers
  all_ints <- seq(temp_range[1], temp_range[2], by = 1)
  # Sample K distinct integers (without replacement)
  sort(sample(all_ints, size = K, replace = FALSE))
}


# 4) Simulate data + compute RMSE using Gamma noise, with peak scaling
simulate_and_eval_rmse <- function(design_temps, priors, n_reps , obs_cv = 0.2) {
  # 1) Sample one parameter set from priors
  params <- priors %>% dplyr::sample_n(1)
  
  # 2) Compute "true" unscaled performance at design points (with a = 1)
  true_unscaled <- briere2_tpc(
    T     = design_temps,
    a     = 1,
    CTmin = params$CTmin,
    CTmax = params$CTmax,
    b     = params$b
  )
  scale_factor <- 1 / max(true_unscaled)  # scale peak to 1
  true_vals    <- true_unscaled * scale_factor
  
  # 3) Simulate gamma data at these design points
  obs <- purrr::map2_dfr(design_temps, true_vals, ~ {
    mu    <- .y
    shape <- 1 / (obs_cv^2)
    scale <- mu * (obs_cv^2)
    tibble::tibble(
      T = rep(.x, each = n_reps),
      y = rgamma(n_reps, shape = shape, scale = scale)
    )
  })
  
  # 4) Fit model to these K points 
  formula_str <- y ~ a * T * (T - CTmin) * pmax(CTmax - T, 0)^b
  
  fit <- tryCatch({
    nls_multstart(
      formula     = as.formula(formula_str),
      data        = obs,
      start_lower = c(a = 0,  CTmin = 0,  CTmax = 30, b = 0.01),
      start_upper = c(a = 2,  CTmin = 15, CTmax = 50, b = 3),
      
      iter        = 50,
      supp_errors = "Y"
    )
  }, error = function(e) NULL)
  
  if (is.null(fit) || !fit$convInfo$isConv) {
    # Return infinite RMSE if fit fails
    return(Inf)
  }
  
  # 5) Evaluate RMSE across a dense temperature grid =====#
  # Construct a dense grid
  Tgrid <- seq(5, 40, by = 0.1)  #probably doesn't need to be this fine
  
  # Calculate "True" curve across that grid, using the same param set & scaling
  #    We'll recompute unscaled with a=1, then apply the same scale_factor
  true_unscaled_grid <- briere2_tpc(
    T     = Tgrid,
    a     = 1,
    CTmin = params$CTmin,
    CTmax = params$CTmax,
    b     = params$b
  )
  # Ensure it doesn't blow up if max=0:
  scale_factor2 <- ifelse(max(true_unscaled_grid) > 0,
                          1 / max(true_unscaled_grid),
                          1)
  true_scaled_grid <- true_unscaled_grid * scale_factor2
  
  # 6) Fitted curve across Tgrid
  #    We must pass each T value to predict(), but we need a new data frame
  #    with columns that match the formula's variables (esp. T).
  #    Also note we need a column for 'a', 'CTmin', 'CTmax','b' if they're in the formula.
  #    But nls won't require them if they're coefficients. We'll do partial approach:
  
  # Extract fitted coefficients
  coefs <- coef(fit)
  # a, CTmin, CTmax
  # We'll compute predicted curve by replicating the same Brière expression:
  # y = a * T * (T - CTmin) * pmax(CTmax - T,0)^b
  pred_scaled_grid <- coefs["a"] * Tgrid * (Tgrid - coefs["CTmin"]) *
    pmax(coefs["CTmax"] - Tgrid, 0)^coefs["b"]
  
  #===== Compare entire grid =====#
  sqrt(mean((true_scaled_grid - pred_scaled_grid)^2))
}

# 8) Experimental settings
K <- 6             # number of temperature treatments
N_search <- 1000    # random designs to sample
n_reps   <- 5      # replicates per temperature

scenario_bs <- c(
  severe_skew = 0.2,
  mild_skew   = 0.5,
  near_sym    = 2
)

# 9) Evaluate designs
all_results <- map_df(names(scenario_bs), function(scn) {
  priors <- draw_priors(100, b_val = scenario_bs[[scn]])
  map_dfr(seq_len(N_search), function(i) {
    temps <- generate_design(K)
    rmse  <- simulate_and_eval_rmse(temps, priors, n_reps)
    tibble(
      scenario  = scn,
      design_id = i,
      temps     = list(temps),
      rmse      = rmse
    )
  })
})

# Replace any NaN rmse with Inf just in case for when nls fails
all_results <- all_results %>% mutate(rmse = ifelse(is.nan(rmse), Inf, rmse))

# 10) Select top 5% designs per scenario
top_designs <- all_results %>%
  group_by(scenario) %>%
  filter(rmse <= quantile(rmse, 0.05, na.rm = TRUE)) %>%
  ungroup()

# 11) Check convergence of N_search and plot
conv_df <- all_results %>%
  group_by(scenario) %>%
  arrange(scenario, design_id) %>%
  mutate(best_rmse = cummin(rmse)) %>%
  ungroup()

# Drop initial NaN/Inf rows for plotting convergence
conv_plot_df <- conv_df %>% filter(is.finite(best_rmse))

ggplot(conv_plot_df, aes(x = design_id, y = best_rmse, color = scenario)) +
  geom_line() +
  labs(
    title = "Convergence of Best RMSE over Search Iterations",
    x     = "Iteration", y = "Best RMSE so far"
  ) +
  theme_minimal()


# 12) Figure 1: density of top designs + best design points + equalized TPC curves
# Prepare density data
dens_df <- top_designs %>%
  unnest_longer(temps, values_to = "T")

# Best design points
best_df <- top_designs %>%
  group_by(scenario) %>%
  slice_min(rmse) %>%
  unnest_longer(temps, values_to = "T")

# Prepare equalized TPC curves (peak = 1)
tpc_curves <- map_df(names(scenario_bs), function(scn) {
  b_val <- scenario_bs[[scn]]
  Tgrid <- seq(5, 40, length.out = 200)
  perf   <- briere2_tpc(Tgrid, a = 1, CTmin = 5, CTmax = 40, b = b_val)
  perf   <- perf / max(perf)
  tibble(scenario = scn, T = Tgrid, perf = perf)
})

# Build and save Figure 1
plot_list <- map(names(scenario_bs), function(scn) {
  df_d <- dens_df %>% filter(scenario == scn)
  df_b <- best_df %>% filter(scenario == scn)
  df_c <- tpc_curves %>% filter(scenario == scn)

  density_max <- max(density(df_d$T)$y)
  df_c <- df_c %>% mutate(scaled = perf * density_max)

  ggplot(df_d, aes(x = T)) +
    geom_density(color = "grey60") +
    geom_point(data = df_b, aes(x = T, y = 0), shape = 4, size = 3) +
    geom_line(data = df_c, aes(x = T, y = scaled), linetype = "dashed") +
    labs(title = scn, x = "Temperature", y = "Density / Scaled Performance") +
    theme_minimal()
})

figure1 <- plot_grid(plotlist = plot_list, ncol = 3)
figure1



