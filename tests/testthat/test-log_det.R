test_that("log_det matches determinant() for nonsingular matrices", {
    set.seed(1)
    for (n in 2:6) {
        A <- crossprod(matrix(stats::rnorm(n * n), n))
        expect_equal(log_det(A), as.numeric(determinant(A)$modulus))
    }
    # negative determinant: log of absolute value
    B <- matrix(c(0, 1, 1, 0), 2)
    expect_equal(log_det(B), 0)
    expect_equal(log_det(diag(c(2, 3, 0.5))), log(3))
})


test_that("log_det returns -Inf for singular matrices", {
    # col3 = 2 * col2 - col1; determinant() can return a finite value here
    expect_equal(log_det(matrix(1:9 / 10, 3)), -Inf)
    expect_equal(log_det(diag(c(1, 0, 2))), -Inf)
})


test_that("log_det errors on non-square matrices", {
    expect_error(log_det(matrix(1:6, 2)), regexp = "Matrix must be square")
})
