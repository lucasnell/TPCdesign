test_that("single_integer validates correctly", {
    expect_null(single_integer(5L, "x"))
    expect_null(single_integer(5, "x"))

    expect_error(single_integer(5.5, "x"))
    expect_error(single_integer(c(1L, 2L), "x"))
    expect_error(single_integer(NA_integer_, "x"))
    expect_error(single_integer("5", "x"))
    expect_error(single_integer(NULL, "x"))

    expect_null(single_integer(5L, "x", .min = 5L))
    expect_error(single_integer(4L, "x", .min = 5L))
    expect_null(single_integer(5L, "x", .max = 5L))
    expect_error(single_integer(6L, "x", .max = 5L))

    expect_error(single_integer(5.5, "n_temps"), regexp = "n_temps.*must be a single integer")
})


test_that("single_number validates correctly", {
    expect_null(single_number(5.5, "x"))

    expect_error(single_number(c(1, 2), "x"))
    expect_error(single_number(NA_real_, "x"))
    expect_error(single_number("5", "x"))
    expect_error(single_number(NULL, "x"))

    expect_null(single_number(5, "x", .min = 5))
    expect_error(single_number(4.9, "x", .min = 5))
    expect_null(single_number(5, "x", .max = 5))
    expect_error(single_number(5.1, "x", .max = 5))
})


test_that("single_string validates correctly", {
    expect_null(single_string("a", "x"))

    expect_error(single_string(1, "x"))
    expect_error(single_string(c("a", "b"), "x"))
    expect_error(single_string(NA_character_, "x"))
    expect_error(single_string(NULL, "x"))
})


test_that("single_logical validates correctly", {
    expect_null(single_logical(TRUE, "x"))
    expect_null(single_logical(FALSE, "x"))

    expect_error(single_logical(1, "x"))
    expect_error(single_logical(c(TRUE, FALSE), "x"))
    expect_error(single_logical(NA, "x"))
    expect_error(single_logical(NULL, "x"))
})


test_that("is_type validates correctly", {
    expect_null(is_type(1:5, "x", type = "integer",
                        len_min = 5, len_max = 5, .min = 1, .max = 5))

    expect_error(is_type(1:5, "x", type = "integer", len_min = 6))
    expect_error(is_type(1:5, "x", type = "integer", len_max = 4))
    expect_error(is_type(1:5, "x", type = "integer", .min = 2))
    expect_error(is_type(1:5, "x", type = "integer", .max = 4))
    expect_error(is_type(1:5, "x", type = "numeric"))
    expect_error(is_type(NULL, "x", type = "integer"))

    expect_null(is_type(1:5, "x", type = is.integer))
    expect_error(is_type("a", "x", type = is.numeric))
})
