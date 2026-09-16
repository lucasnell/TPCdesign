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
