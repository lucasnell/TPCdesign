test_that("gap_filler returns opt_temps unchanged when n_filler is 0", {
    pts <- c(10, 20, 30)
    expect_equal(gap_filler(0, pts, min_temp = 5, max_temp = 35, digits = 2), pts)
})


test_that("gap_filler fills the largest gap first", {
    pts <- c(10, 30)
    out <- gap_filler(1, pts, min_temp = 5, max_temp = 35, digits = 2)

    # boundaries 5 and 35 create gaps of 5, 20, 5; the 20-wide gap between
    # 10 and 30 should be split first
    expect_length(out, 3)
    expect_true(20 %in% out)
})


test_that("gap_filler adds filler points in order of decreasing gap size", {
    pts <- c(10, 30)
    out <- gap_filler(2, pts, min_temp = 5, max_temp = 35, digits = 2)

    # after adding 20, remaining gaps are 5, 10, 10, 5; one of the two
    # 10-wide gaps (5-10 or 30-35) gets split next
    expect_length(out, 4)
    expect_true(20 %in% out)
})


test_that("gap_filler drops domain boundaries that weren't already optimal", {
    pts <- c(10, 20)
    out <- gap_filler(1, pts, min_temp = 5, max_temp = 35, digits = 0)

    expect_false(5 %in% out)
    expect_false(35 %in% out)
})


test_that("gap_filler retains a boundary that was already an optimal point", {
    pts <- c(5, 20)
    out <- gap_filler(1, pts, min_temp = 5, max_temp = 35, digits = 2)

    expect_true(5 %in% out)
    expect_false(35 %in% out)
})
