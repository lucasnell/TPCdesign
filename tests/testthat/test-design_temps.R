test_that("design_temps works", {

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
