test_that("briere2_tpc_Topt works", {
    briere2_tpc_Topt_R <- function(ctmin, ctmax, b) {
      (2 * ctmax + ctmin + b * ctmin + sqrt(-4 * (2 + b) * ctmax * ctmin +
                                                (-2 * ctmax - ctmin - b * ctmin)^2)) /
          (2 * (2 + b))
    }
    topt_tests <- do.call(c, lapply(1:1000, \(i) {
        ctmin <- runif(1, 0, 10)
        ctmax <- runif(1, 30, 50)
        b <- runif(1, 0.01, 4)
        x <- briere2_tpc_Topt_R(ctmin, ctmax, b)
        y <- briere2_tpc_Topt(ctmin, ctmax, b)
        return(all.equal(x, y))
    }))
    expect_all_true(topt_tests)
})


test_that("briere2_tpc_deriv works", {
    briere2_tpc_deriv_R <- function(temp, ctmin, ctmax, a, b) {
        a * (ctmax - temp)^b * temp + a * (ctmax - temp)^b * (temp - ctmin) -
            a * b * (ctmax - temp)^(b - 1) * temp * (temp - ctmin)
    }
    deriv_tests <- do.call(c, lapply(1:1000, \(i) {
        ctmin <- runif(1, 0, 10)
        ctmax <- runif(1, 30, 50)
        a <- runif(1, 0.01, 1)
        b <- runif(1, 0.01, 4)
        temps <- seq(ctmin+1, ctmax-1, length.out = 101)
        x <- briere2_tpc_deriv_R(temps, ctmin, ctmax, a, b)
        y <- briere2_tpc_deriv(temps, ctmin, ctmax, a, b)
        return(all.equal(x, y))
    }))
    expect_all_true(deriv_tests)
})


test_that("briere2_tpc_Topt recycles length-1 arguments", {
    out1 <- briere2_tpc_Topt(ctmin = c(5, 6, 7), ctmax = 40, b = 0.5)
    out2 <- sapply(c(5, 6, 7), function(cm) briere2_tpc_Topt(cm, 40, 0.5))
    expect_equal(as.vector(out1), as.vector(out2))
})


test_that("briere2_tpc_Topt errors on incompatible argument lengths", {
    expect_error(briere2_tpc_Topt(ctmin = c(5, 6), ctmax = c(30, 35, 40), b = 0.5))
})


test_that("briere2_tpc works", {
    briere2_tpc_R <- function(temp, ctmin, ctmax, a, b) {
        a * temp * pmax(temp - ctmin, 0.0) * pmax(ctmax - temp, 0.0)^b
    }
    tpc_tests <- do.call(c, lapply(1:1000, \(i) {
        ctmin <- runif(1, 0, 10)
        ctmax <- runif(1, 30, 50)
        a <- runif(1, 0.01, 1)
        b <- runif(1, 0.01, 4)
        temps <- seq(ctmin+1, ctmax-1, length.out = 101)
        x <- briere2_tpc_R(temps, ctmin, ctmax, a, b)
        y <- briere2_tpc(temps, ctmin, ctmax, a, b)
        return(all.equal(x, y))
    }))
    expect_all_true(tpc_tests)
})


test_that("briere2_tpc with scale = TRUE rescales the curve to a maximum of 1", {
    temp <- seq(8.01, 37.99, length.out = 500)
    y <- briere2_tpc(temp, ctmin = 8, ctmax = 38, a = 1, b = 0.5, scale = TRUE)
    expect_equal(max(y), 1)
    expect_true(all(y <= 1))
})



