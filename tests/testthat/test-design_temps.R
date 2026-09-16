test_that("design_temps works", {

    # design_temps relies on R's RNG (via lhs::randomLHS and stats::runif),
    # so fix the seed to keep this stochastic optimization test reproducible
    set.seed(1)

    test_dt <- function(n_temps, min_sep = 1, n_filler = 1) {
        design_temps(n_temps = n_temps, min_sep = min_sep,
                     ctmin = 5, ctmax = 40, a = 1, b = 0.5,
                     n_draws = 50L, n_filler = n_filler,
                     n_starts = 1L,
                     digits = 2L)
    }

    for (L in 2:10) expect_length(test_dt(L), L)

    for (ms in c(0.1, 0.5, 1, 2)) {
        for (L in 2:10) {
            dt <- test_dt(L, min_sep = ms, n_filler = 0)
            expect_gte(min(diff(dt)), ms)
        }
    }

})


test_that("design_temps validates its arguments", {
    base_args <- list(n_temps = 5, ctmin = 5, ctmax = 40, a = 1, b = 0.5,
                      n_draws = 10L, n_starts = 1L)

    expect_error(do.call(design_temps, modifyList(base_args, list(n_temps = 1))))
    expect_error(do.call(design_temps, modifyList(base_args, list(ctmax = 5))))
    expect_error(do.call(design_temps, modifyList(base_args, list(a = -1))))
    expect_error(do.call(design_temps, modifyList(base_args, list(b = 0))))
    expect_error(do.call(design_temps, modifyList(base_args, list(n_filler = 10))))
    expect_error(do.call(design_temps, modifyList(base_args, list(ctmin_err = -1))))
})


test_that("design_temps errors when min_sep is incompatible with n_temps", {
    expect_error(
        design_temps(n_temps = 20, min_sep = 5, ctmin = 5, ctmax = 40, a = 1, b = 0.5,
                     n_draws = 10L, n_filler = 0, n_starts = 1L),
        regexp = "min_sep is too large"
    )
})
