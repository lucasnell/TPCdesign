test_that("sim_gamma_data returns correctly structured output", {
    temp <- seq(10, 35, 5)
    n_reps <- 4L
    out <- sim_gamma_data(temp, n_reps, obs_cv = 0.1,
                          ctmin = 8, ctmax = 38, a = 1, b = 0.2)

    expect_s3_class(out, "data.frame")
    expect_named(out, c("temp", "y"))
    expect_equal(nrow(out), length(temp) * n_reps)
    expect_equal(sort(unique(out$temp)), sort(temp))
    expect_equal(as.vector(table(out$temp)), rep(n_reps, length(temp)))
})


test_that("sim_gamma_data validates its arguments", {
    args <- list(temp = seq(10, 35, 5), n_reps = 2L, obs_cv = 0.1,
                 ctmin = 8, ctmax = 38, a = 1, b = 0.2)

    expect_error(do.call(sim_gamma_data, modifyList(args, list(n_reps = 0L))))
    expect_error(do.call(sim_gamma_data, modifyList(args, list(obs_cv = 0))))
    expect_error(do.call(sim_gamma_data, modifyList(args, list(obs_cv = -0.1))))
    expect_error(do.call(sim_gamma_data, modifyList(args, list(a = -1))))
    expect_error(do.call(sim_gamma_data, modifyList(args, list(b = -0.2))))
})


test_that("sim_gamma_data is reproducible when the RNG seed is set", {
    args <- list(temp = seq(10, 35, 5), n_reps = 5L, obs_cv = 0.2,
                 ctmin = 8, ctmax = 38, a = 1, b = 0.2)

    set.seed(42)
    out1 <- do.call(sim_gamma_data, args)
    set.seed(42)
    out2 <- do.call(sim_gamma_data, args)

    expect_equal(out1, out2)
})


test_that("sim_gamma_data's simulated means approximate the underlying TPC", {
    temp <- seq(10, 35, 5)
    ctmin <- 8; ctmax <- 38; a <- 1; b <- 0.2
    true_y <- briere2_tpc(temp, ctmin, ctmax, a, b)

    set.seed(123)
    out <- sim_gamma_data(temp, n_reps = 20000L, obs_cv = 0.2,
                          ctmin = ctmin, ctmax = ctmax, a = a, b = b)

    obs_means <- tapply(out$y, out$temp, mean)
    obs_means <- as.numeric(obs_means[as.character(temp)])

    expect_equal(obs_means, true_y, tolerance = 0.05)
})


test_that("sim_gamma_data respects scale_tpc", {
    temp <- seq(10, 35, 5)
    ctmin <- 8; ctmax <- 38; a <- 1; b <- 0.2
    true_y_scaled <- briere2_tpc(temp, ctmin, ctmax, a, b, scale = TRUE)

    set.seed(99)
    out <- sim_gamma_data(temp, n_reps = 20000L, obs_cv = 0.2,
                          ctmin = ctmin, ctmax = ctmax, a = a, b = b,
                          scale_tpc = TRUE)

    obs_means <- tapply(out$y, out$temp, mean)
    obs_means <- as.numeric(obs_means[as.character(temp)])

    expect_equal(obs_means, true_y_scaled, tolerance = 0.05)
})


test_that("sim_gamma_data returns 0 where expected performance is zero", {
    # at or outside ctmin and ctmax
    temp <- c(5, 8, 20, 38, 40)
    set.seed(7)
    out <- sim_gamma_data(temp, n_reps = 3L, obs_cv = 0.1,
                          ctmin = 8, ctmax = 38, a = 1, b = 0.2)
    expect_true(all(out$y[out$temp != 20] == 0))
    expect_true(all(out$y[out$temp == 20] > 0))

    # a = 0 makes the whole curve zero
    out_a0 <- sim_gamma_data(c(10, 20, 30), n_reps = 3L, obs_cv = 0.1,
                             ctmin = 8, ctmax = 38, a = 0, b = 0.2)
    expect_true(all(out_a0$y == 0))
})


test_that("trunc_rnorm draws from a normal truncated below at zero", {
    mu <- 1; sigma <- 2
    set.seed(11)
    x <- trunc_rnorm(1e5, mu, sigma)

    # mean and SD of a normal truncated below at 0
    alpha <- -mu / sigma
    lambda <- stats::dnorm(alpha) / (1 - stats::pnorm(alpha))
    true_mean <- mu + sigma * lambda
    true_sd <- sigma * sqrt(1 + alpha * lambda - lambda^2)

    expect_length(x, 1e5)
    expect_true(all(x >= 0))
    expect_equal(mean(x), true_mean, tolerance = 0.01)
    expect_equal(stats::sd(x), true_sd, tolerance = 0.01)
})


test_that("trunc_rnorm handles degenerate parameters", {
    # mu more than 5 SDs below zero returns all zeros
    expect_equal(as.vector(trunc_rnorm(10, -6, 1)), rep(0, 10))
    # sigma = 0 returns mu
    expect_equal(as.vector(trunc_rnorm(10, 3, 0)), rep(3, 10))
})


test_that("trunc_rnorm is reproducible when the RNG seed is set", {
    set.seed(3)
    x1 <- trunc_rnorm(20, 1, 1)
    set.seed(3)
    x2 <- trunc_rnorm(20, 1, 1)
    expect_identical(x1, x2)
})


test_that("trunc_rnorm_range returns the minimum and maximum draws", {
    set.seed(5)
    r <- as.vector(trunc_rnorm_range(1e4, 1, 2))
    expect_length(r, 2)
    expect_gte(r[1], 0)
    expect_gt(r[2], r[1])

    set.seed(5)
    x <- trunc_rnorm(1e4, 1, 2)
    expect_equal(r, range(x))

    expect_equal(as.vector(trunc_rnorm_range(10, -6, 1)), c(0, 0))
    expect_equal(as.vector(trunc_rnorm_range(10, 200, 0)), c(200, 200))
})
