test_that("briere2_tpc_Topt works", {
    briere2_tpc_Topt_R <- function(ctmin, ctmax, b) {
      (2 * ctmax + ctmin + b * ctmin + sqrt(-4 * (2 + b) * ctmax * ctmin +
                                                (-2 * ctmax - ctmin - b * ctmin)^2)) /
          (2 * (2 + b))
    }
    topt_tests <- lapply(1:1000, \(i) {
        ctmin <- runif(1, 0, 10)
        ctmax <- runif(1, 30, 50)
        b <- runif(1, 0.01, 4)
        x <- briere2_tpc_Topt_R(ctmin, ctmax, b)
        y <- briere2_tpc_Topt(ctmin, ctmax, b)
        return(all.equal(x, y))
    }) |> do.call(what = c)
    expect_all_true(topt_tests)
})


test_that("briere2_tpc_deriv works", {
    briere2_tpc_deriv_R <- function(temp, ctmin, ctmax, a, b) {
        a * (ctmax - temp)^b * temp + a * (ctmax - temp)^b * (temp - ctmin) -
            a * b * (ctmax - temp)^(b - 1) * temp * (temp - ctmin)
    }
    deriv_tests <- lapply(1:1000, \(i) {
        ctmin <- runif(1, 0, 10)
        ctmax <- runif(1, 30, 50)
        a <- runif(1, 0.01, 1)
        b <- runif(1, 0.01, 4)
        temps <- seq(ctmin+1, ctmax-1, length.out = 101)
        x <- briere2_tpc_deriv_R(temps, ctmin, ctmax, a, b)
        y <- briere2_tpc_deriv(temps, ctmin, ctmax, a, b)
        return(all.equal(x, y))
    }) |> do.call(what = c)
    expect_all_true(deriv_tests)
})


test_that("briere2_tpc works", {
    briere2_tpc_R <- function(temp, ctmin, ctmax, a, b) {
        a * temp * pmax(temp - ctmin, 0.0) * pmax(ctmax - temp, 0.0)^b
    }
    tpc_tests <- lapply(1:1000, \(i) {
        ctmin <- runif(1, 0, 10)
        ctmax <- runif(1, 30, 50)
        a <- runif(1, 0.01, 1)
        b <- runif(1, 0.01, 4)
        temps <- seq(ctmin+1, ctmax-1, length.out = 101)
        x <- briere2_tpc_R(temps, ctmin, ctmax, a, b)
        y <- briere2_tpc(temps, ctmin, ctmax, a, b)
        return(all.equal(x, y))
    }) |> do.call(what = c)
    expect_all_true(tpc_tests)
})



